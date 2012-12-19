import Test.QuickCheck
import Text.Printf (printf)

main = mapM_ (\(s,a) -> printf "%-25s: " s >> a) tests

prop_test_00 s = s == s
  where _ = s :: Int

tests = [("test_00", quickCheck prop_test_00)]
