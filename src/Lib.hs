module Lib(
    processAtom,
    drawGraph,
    printEaTable
) where

import Atom
import AtomEnergy
import Graphing

processAtom :: Int -> Int -> IO ()
processAtom z c = let atom = makeAtom z c in do
    putStrLn $ prettyElectronArrangement atom
    let energy = totalEnergy atom
        ions   = reverse (cationsOf atom) ++ anionsOf atom
    putStrLn $ show $ energy
    putStrLn $ show $ energies' atom
    graphAtom atom
    putStrLn $ unlines $ map (showIon energy) ions

anionsOf :: Atom -> [Atom]
anionsOf = takeWhile atomFull . tail . iterate anionise

cationsOf :: Atom -> [Atom]
cationsOf atom = takeWhile (\ion -> charge ion <= atomicNumber ion && totalEnergy ion < totalEnergy atom + 0.07) $ iterate cationise atom

showIon :: Double -> Atom -> String
showIon e0 ion = (if charge ion >= 0 then "+" else "") ++ show (charge ion) ++ ": " ++ show (e0 - totalEnergy ion)

printEaTable :: Int -> IO ()
printEaTable = mapM_ processElement . flip take aperiodicTable
    where processElement a = do
            putStr $ show $ atomicNumber a
            putStr ":  "
            putStrLn $ prettyElectronArrangement a
            graphAtom a

--eaTable = map (\a -> show (atomicNumber a) ++ ": " ++ prettyElectronArrangement a) aperiodicTable
