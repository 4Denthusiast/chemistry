module Main where

import System.Environment
import Lib

main :: IO ()
main = do
    args <- getArgs
    case filter ((/= '-') . head) args of
        "graph":args1 -> case args1 of
            [] -> putStrLn "Argument required: atomic number"
            [z]       -> processAtom (read z) (2 * read z) 0
            [z, a]    -> processAtom (read z) (read a) 0
            [z, a, c] -> processAtom (read z) (read a) (parseSignedNumber c)
            _         -> putStrLn "too many arguments"
        "ea":args1 -> case args1 of
            []     -> printEaTable Nothing
            [z]    -> printEaTable $ Just $ read z
            _      -> putStrLn "too many arguments"
        [] -> putStrLn "Argument required: action (options: \"graph\", \"ea\")"
        x:_ -> putStrLn ("Argument not recognised: " ++ x)


-- Starting commend-line aarguments with "-" interferes with Stack, and putting the sign at the end is standard chemical notation anyway.
parseSignedNumber :: String -> Int
parseSignedNumber s
    | last s == '-' = negate $ read $ reverse $ tail $ reverse $ s
    | otherwise     = read s
