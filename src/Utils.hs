module Utils(
    fst3,
    snd3,
    þrd3,
    fst4,
    frþ4,
    mergeBy,
    uncurry3,
    map2,
    map3,
    zipWithLong,
    cmdHasOption,
    labHead
) where

import System.Environment

fst3 :: (a, b, c) -> a
fst3 (x, _, _) = x

snd3 :: (a, b, c) -> b
snd3 (_, x, _) = x

þrd3 :: (a, b, c) -> c
þrd3 (_, _, x) = x

fst4 :: (a, b, c, d) -> a
fst4 (x, _, _, _) = x

frþ4 :: (a, b, c, d) -> d
frþ4 (_, _, _, x) = x

mergeBy :: (a -> a -> Maybe a) -> [a] -> [a]
mergeBy f = mergeBy' []
    where mergeBy' ys [] = ys
          mergeBy' ys (x:xs) = mergeBy' (mergeInto ys x) xs
          mergeInto [] x = [x]
          mergeInto (y:ys) x = case f y x of
              Just y' -> y':ys
              Nothing -> y: mergeInto ys x

uncurry3 :: (a -> b -> c -> d) -> (a,b,c) -> d
uncurry3 f (a,b,c) = f a b c

map2 :: (a -> b) -> [[a]] -> [[b]]
map2 = map . map
map3 :: (a -> b) -> [[[a]]] -> [[[b]]]
map3 = map . map2

-- Extend the shorter list to the length of the longer.
zipWithLong :: a -> b -> (a -> b -> c) -> [a] -> [b] -> [c]
zipWithLong xr yr f [] ys = map (f xr) ys
zipWithLong xr yr f xs [] = map (flip f yr) xs
zipWithLong xr yr f (x:xs) (y:ys) = f x y : zipWithLong xr yr f xs ys

cmdHasOption :: Char -> IO Bool
cmdHasOption c = (elem c . concat . map tail . filter ((=='-') . head)) <$> getArgs

labHead :: String -> [a] -> a
labHead s [] = error ("labHead: empty list at " ++ s)
labHead _ (x:_) = x
