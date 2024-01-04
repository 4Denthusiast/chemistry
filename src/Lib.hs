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
    let energy = totalEnergy atom
    putStrLn $ prettyElectronArrangement atom
    putStrLn $ show $ energy
    putStrLn $ show $ energies atom
    graphAtom atom
    anions <- anionsOf atom
    cations <- cationsOf atom
    let ions = reverse $ reverse anions ++ [atom] ++ cations
    putStrLn $ unlines $ map (showIon energy) ions

anionsOf :: Atom -> IO [Atom]
anionsOf atom = do
    cm <- getCacheMode
    case cm of
        NoCache -> return $ takeWhile ((==0) . incorrectCharge) $ tail $ iterate anionise atom
        _ -> takeWhileM ((==0) . incorrectCharge) $ map (makeAtomUsingCache (atomicNumber atom) . (charge atom -)) [1..]

cationsOf :: Atom -> IO [Atom]
cationsOf atom = do
        cm <- getCacheMode
        case cm of
            NoCache -> return $ takeWhile valid $ tail $ iterate cationise atom
            _ -> takeWhileM valid $ map (makeAtomUsingCache (atomicNumber atom)) [charge atom+1..atomicNumber atom]
    where valid ion = charge ion < atomicNumber ion && totalEnergy ion < totalEnergy atom + 0.07

takeWhileM :: Monad m => (a -> Bool) -> [m a] -> m [a]
takeWhileM f [] = return []
takeWhileM f (x:xs) = x >>= (\x' -> if f x' then (x':) <$> takeWhileM f xs else return [])

showIon :: Double -> Atom -> String
showIon e0 ion = (if charge ion >= 0 then "+" else "") ++ show (charge ion) ++ ": " ++ show (e0 - totalEnergy ion)

printEaTable :: Maybe Int -> IO ()
printEaTable n = mapM_ processElement =<< (maybe id take n <$> atoms)
    where processElement a' = do
              a <- a'
              putStr $ show $ atomicNumber a
              putStr ":  "
              putStrLn $ prettyElectronArrangement a
              putStrLn $ show $ energies a
              putStrLn $ show $ orbitalRadii a
              graphAtom a
          atoms = do
              cm <- getCacheMode
              return $ case cm of
                  NoCache -> map return aperiodicTable
                  _       -> map (flip makeAtomUsingCache 0) [1..]
