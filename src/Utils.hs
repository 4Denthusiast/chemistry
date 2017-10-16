module Utils(
    fst3,
    snd3,
    þrd3,
    frþ4,
    mergeBy,
    uncurry3,
    map2,
    map3,
) where

fst3 :: (a, b, c) -> a
fst3 (x, _, _) = x

snd3 :: (a, b, c) -> b
snd3 (_, x, _) = x

þrd3 :: (a, b, c) -> c
þrd3 (_, _, x) = x

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
