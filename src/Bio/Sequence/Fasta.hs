{- |
   Module: Bio.Sequence.Fasta

   This module incorporates functionality for reading and writing
   sequence data in the Fasta format.
   Each sequence consists of a header (with a '>' prefix)
   and a set of lines containing the sequence data..
-}

module Bio.Sequence.Fasta
    ( Sequence(..)
    -- * Reading and writing plain FASTA files
    , readFasta, writeFasta, hReadFasta, hWriteFasta
    -- * Reading and writing quality files
    , readQual, writeQual, hWriteQual
    -- * Combining FASTA and quality files
    , readFastaQual, hWriteFastaQual, writeFastaQual
    -- * Counting sequences in a FASTA file
    , countSeqs
    -- * Helper function for reading your own sequences
    , mkSeqs
    -- * Other
    , toStr, seqid, seqheader, seqdata, seqlength
) where


import Data.Char (isSpace)
import Data.List (groupBy, intersperse)
import System.IO

import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.ByteString.Lazy as BB
import Data.ByteString.Lazy.Char8 (ByteString)

import Bio.Core.Sequence

splitsAt :: Offset -> ByteString -> [ByteString]
splitsAt n s = let (s1,s2) = B.splitAt (unOff n) s
               in if B.null s2 then [s1] else s1 : splitsAt n s2

data Sequence = Seq SeqLabel SeqData (Maybe QualData)
                deriving (Show, Eq)

instance BioSeq Sequence where
  seqid     (Seq lab _sq _mqual) = SeqLabel {unSL = B.takeWhile (/= ' ') $ unSL lab}
  seqheader (Seq lab _sq _mqual) = lab
  seqdata   (Seq _lab sq _mqual) = sq
  seqlength (Seq _lab sq _mqual) = Offset {unOff = B.length $ unSD sq}

toStr :: SeqData -> String
toStr  = B.unpack . unSD

-- | Lazily read sequences from a FASTA-formatted file
readFasta :: FilePath -> IO [Sequence]
readFasta f = (mkSeqs . B.lines) `fmap` B.readFile f

-- | Read quality data for sequences to a file.
readQual :: FilePath -> IO [Sequence]
readQual f = (mkQual . B.lines) `fmap` B.readFile f

-- | Write quality data for sequences to a file.
writeQual :: FilePath -> [Sequence] -> IO ()
writeQual f ss = do
  h <- openFile f WriteMode
  hWriteQual h ss
  hClose h

-- | Read sequence and associated quality.  Will error if
--   the sequences and qualites do not match one-to-one in sequence.
readFastaQual :: FilePath -> FilePath -> IO [Sequence]
readFastaQual  s q = do
  ss <- readFasta s
  qs <- readQual q
  -- warning: assumes correct qual file!
  let mkseq s1@(Seq x y _) s2@(Seq _ _ (Just z))
          | seqid s1 /= seqid s2 = error ("mismatching sequence and quality: "
                                                ++B.unpack (unSL $ seqid s1)++","++B.unpack (unSL $ seqid s2))
          | B.length (unSD y) /= B.length (unQD z)   = error ("mismatching sequence and quality lengths:"
                                                ++ B.unpack (unSL $ seqid s1)++","++show (B.length (unSD y),B.length (unQD z)))
          | otherwise   = Seq x y (Just z)
      mkseq _ _ = error "readFastaQual: could not combine Fasta and Qual information"
      zipWith' f (x:xs) (y:ys) = f x y : zipWith' f xs ys
      zipWith' _ [] [] = []
      zipWith' _ _ _ = error "readFastaQual: Unbalanced sequence and quality"
  return $ zipWith' mkseq ss qs

-- | Write sequence and quality data simulatnously
--   This may be more laziness-friendly.
writeFastaQual :: FilePath -> FilePath -> [Sequence] -> IO ()
writeFastaQual f q ss = do
  hf <- openFile f WriteMode
  hq <- openFile q WriteMode
  hWriteFastaQual hf hq ss
  hClose hq
  hClose hf

hWriteFastaQual :: Handle -> Handle -> [Sequence] -> IO ()
hWriteFastaQual hf hq = mapM_ wFQ
    where wFQ s = wFasta hf s >> wQual hq s

-- | Write sequences to a FASTA-formatted file.
--   Line length is 60.
writeFasta :: FilePath -> [Sequence] -> IO ()
writeFasta f ss = do
  h <- openFile f WriteMode
  hWriteFasta h ss
  hClose h

-- | Lazily read sequence from handle
hReadFasta :: Handle -> IO [Sequence]
hReadFasta h = (mkSeqs . B.lines) `fmap` B.hGetContents h

-- | Write sequences in FASTA format to a handle.
hWriteFasta :: Handle -> [Sequence] -> IO ()
hWriteFasta h = mapM_ (wFasta h)

wHead :: Handle -> SeqLabel -> IO ()
wHead h l = do
   B.hPut h $ B.pack ">"
   B.hPut h (unSL l)
   B.hPut h $ B.pack "\n"

wFasta :: Handle -> Sequence -> IO ()
wFasta h (Seq l d _) = do
  wHead h l
  let ls = splitsAt 60 (unSD d)
  mapM_ (B.hPut h) $ intersperse (B.pack "\n") ls
  B.hPut h $ B.pack "\n"

mkQual :: [ByteString] -> [Sequence]
mkQual = map f . blocks
    where f [] = error "mkQual on empty input - this is impossible"
          f (l:ls) = Seq (SeqLabel (B.drop 1 l)) (SeqData B.empty)
                     (Just $ QualData $ BB.pack $ go 0 ls)
          isDigit x = x <= 58 && x >= 48 
          go i (s:ss) = case BB.uncons s of Just (c,rs) -> if isDigit c then go (c - 48 + 10*i) (rs:ss)
                                                           else let rs' = BB.dropWhile (not . isDigit) rs
                                                                in if BB.null rs' then i : go 0 ss else i : go 0 (rs':ss)
                                            Nothing -> i : go 0 ss
          go _ [] = []

hWriteQual :: Handle -> [Sequence] -> IO ()
hWriteQual h = mapM_ (wQual h)

wQual :: Handle -> Sequence -> IO ()
wQual h (Seq l _ (Just (QualData q))) = do
  wHead h l
  let qls = splitsAt 20 q
      qs = B.pack . unwords . map show . BB.unpack
  mapM_ ((\l' -> B.hPut h l' >> B.hPut h (B.pack "\n")) . qs) qls
wQual _ (Seq _ _ Nothing) = return ()

-- | Convert a list of FASTA-formatted lines into a list of sequences.
--   Blank lines are ignored.
--   Comment lines start with "#" are allowed between sequences (and ignored).
--   Lines starting with ">" initiate a new sequence.
mkSeqs :: [ByteString] -> [Sequence]
mkSeqs = map mkSeq . blocks

mkSeq :: [ByteString] -> Sequence
mkSeq (l:ls) = Seq (SeqLabel (B.drop 1 l))
               (SeqData (B.filter (not . isSpace) $ B.concat $ takeWhile isSeq ls))
               Nothing
               where isSeq s = (not . B.null) s &&
                               (flip elem (['A'..'Z']++['a'..'z']) . B.head) s
mkSeq []     = error "empty input to mkSeq"

-- | Split lines into blocks starting with '>' characters
--   Filter out # comments (but not semicolons?)
blocks :: [ByteString] -> [[ByteString]]
blocks = groupBy (const (('>' /=) . B.head))
         . filter ((/='#') . B.head)
         . dropWhile (('>' /=) . B.head)
         . filter (not . B.null)

countSeqs :: FilePath -> IO Int
countSeqs f = do
  ss <- B.readFile f
  let hdrs = filter (('>'==).B.head) $ filter (not . B.null) $ B.lines ss
  return (length hdrs)
