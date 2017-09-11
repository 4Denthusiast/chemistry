module Main where

import System.Environment
import Lib

main :: IO ()
main = do
    args <- getArgs
    case args of
        "graph":args1 -> case args1 of
            [] -> putStrLn "Argument required: atomic number"
            [z]    -> processAtom (read z) 0
            [z, c] -> processAtom (read z) (parseSignedNumber c)
            _      -> putStrLn "too many arguments"
        "ea":args1 -> case args1 of
            []     -> putStrLn "Argument required: max atomic number"
            [z]    -> printEaTable (read z)
            _      -> putStrLn "too many arguments"
        [] -> putStrLn "Argument required: action"
        x:_ -> putStrLn ("Argument not recognised: " ++ x)


parseSignedNumber :: String -> Int
parseSignedNumber s
    | last s == '-' = negate $ read $ reverse $ tail $ reverse $ s
    | otherwise     = read s
