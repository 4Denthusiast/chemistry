module Cache(
    makeAtomUsingCache,
    cacheAtomEA,
    CacheMode(..),
    getCacheMode
) where

import System.Environment
import System.IO.Error
import System.IO
import Data.Bifunctor (first)
import Data.List
import Data.Maybe

import Utils
import Orbital
import Atom

data CacheMode = NoCache | UseCache | TrustCache

loadCache :: IO [[[String]]]
loadCache = map (map (splitBy ',') . splitBy ';') <$> catchIOError (withFile "cache.csv" ReadMode hGetLines) (\e -> if isDoesNotExistError e then return [] else ioError e)
    where splitBy c = map tail . groupBy (const (/= c)) . (c:)
          hGetLines h = catchIOError ((:) <$> hGetLine h <*> hGetLines h) (\e -> if isEOFError e then return [] else ioError e)

readEA :: String -> Maybe [(N,L,Int)]
readEA "" = Nothing
readEA s = Just $ map readPart $ words s
    where readPart (n:l:'^':o) = let l' = fromMaybe (error ("Invalid angular momentum label in cache file: "++[l])) $ elemIndex l angularMomentumLabels in (read [n] - l', l', read o)
          readPart s = error ("Invalid electron arrangement sytax in cache file: " ++ s)

-- The bool output is whether the returned e.a. is for exactly the requested z, a & charge.
getCachedEA :: Int -> Int -> Int -> IO (Bool, [(N,L,Int)])
getCachedEA z a c = findEA <$> loadCache
    where findEA ls = fromMaybe (False, []) $ fmap (\(v0,(v1,(v2,c))) -> (v0&&v1&&v2, c)) $ search z $ map (search (a-z+1) . map (search (z-c+1) . map readEA)) ls
          search :: Int -> [Maybe a] -> Maybe (Bool, a)
          search n xs = find True $ uncurry interleave $ first reverse $ splitAt n (xs ++ [Nothing])
          interleave [] xs = xs
          interleave xs [] = xs
          interleave (x:xs) ys = x: interleave ys xs
          find :: Bool -> [Maybe a] -> Maybe (Bool, a)
          find b [] = Nothing
          find b (Nothing:xs) = find False xs
          find b (Just x:_) = Just (b,x)

-- I don't understand why strictly evaluating s first helps, but otherwise the cache file keeps getting emptied.
putCachedEA :: Int -> Int -> Int -> String -> IO ()
putCachedEA z a c s = seq s $ writeFile "cache.csv" =<< (unlines <$> map (intercalate ";" . map (intercalate ",")) <$> update (z-1) [] (update (a-z) [] (update (z-c) "" (const s))) <$> loadCache)

update :: Int -> a -> (a -> a) -> [a] -> [a]
update n r f [] = replicate n r ++ [f r]
update 0 r f (x:xs) = f x : xs
update n r f (x:xs) = x : update (n-1) r f xs

cacheAtomEA :: Atom -> IO ()
cacheAtomEA a = putCachedEA (atomicNumber a) (massNumber a) (charge a) (prettyElectronArrangement' a)

prettyElectronArrangement' :: Atom -> String
prettyElectronArrangement' = (\s -> if null s then " " else s) . unwords . reverse . filter ((/= '!') . last) . reverse . words . prettyElectronArrangement

getCacheMode :: IO CacheMode
getCacheMode = do
    t <- cmdHasOption 't'
    i <- cmdHasOption 'i'
    return $ if t
        then TrustCache
    else if i
        then NoCache
    else
        UseCache

makeAtomUsingCache :: Int -> Int -> Int -> IO Atom
makeAtomUsingCache z a c = do
    cm <- getCacheMode
    (t,ea) <- getCachedEA z a c
    let cm' = case (t,cm) of {(False, TrustCache) -> UseCache; (_,x) -> x}
    (\at -> cacheAtomEA at >> return at) $ case cm' of
        NoCache    -> makeAtom           z a c
        UseCache   -> makeAtomWithEA     z a c ea
        TrustCache -> makeAtomTrustingEA z a c ea
