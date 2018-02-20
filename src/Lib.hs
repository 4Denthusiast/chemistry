module Lib(
    processAtom,
    printEaTable
) where

import Atom
import AtomEnergy
import Graphing
import Cache

processAtom :: Int -> Int -> IO ()
processAtom z c = do
    atom <- makeAtomUsingCache z c
    putStrLn $ prettyElectronArrangement atom
    let energy = totalEnergy atom
        ions   = reverse $ reverse (anionsOf atom) ++ cationsOf atom
    putStrLn $ show $ energy
    putStrLn $ show $ energies atom
    graphAtom atom
    putStrLn $ unlines $ map (showIon energy) ions

anionsOf :: Atom -> [Atom]
anionsOf = takeWhile ((==0) . incorrectCharge) . iterate anionise

cationsOf :: Atom -> [Atom]
cationsOf atom = takeWhile (\ion -> charge ion <= atomicNumber ion && totalEnergy ion < totalEnergy atom + 0.07) $ tail $ iterate cationise atom

showIon :: Double -> Atom -> String
showIon e0 ion = (if charge ion >= 0 then "+" else "") ++ show (charge ion) ++ ": " ++ show (e0 - totalEnergy ion)

printEaTable :: Maybe Int -> IO ()
printEaTable n = mapM_ processElement =<< (maybe id take n <$> atoms)
    where processElement a' = do
              a <- a'
              putStr $ show $ atomicNumber a
              putStr ":  "
              putStrLn $ prettyElectronArrangement a
              graphAtom a
          atoms = do
              cm <- getCacheMode
              return $ case cm of
                  NoCache -> map return aperiodicTable
                  _       -> map (flip makeAtomUsingCache 0) [1..]
