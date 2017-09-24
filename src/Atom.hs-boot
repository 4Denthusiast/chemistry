module Atom(
    Atom(..),
    energies',
    electronArrangement,
    prevElectronArrangement,
    indices
) where

import Orbital

data Atom = Atom
    {
        atomicNumber :: Int,
        charge :: Int,
        energies :: [[Energy]],
        orbitals :: [[Orbital]],
        shieldingPotentials :: [[Potential]],
        atomGrid :: Grid,
        prevOccs :: [[Int]],
        forcedEA :: Maybe [(N, L, Int)]
    }

energies' :: Atom -> [[Energy]]

electronArrangement :: Atom -> [(N, L, Int)]
prevElectronArrangement :: Atom -> [(N, L, Int)]

indices :: [[(N,L)]]
