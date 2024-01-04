module Atom(
    Atom(..),
    PerOrbital,
    getPO,
    indices,
    electronArrangement
) where

import Orbital

data PerOrbital a = PerOrbital a [[a]]
instance Functor PerOrbital
instance Applicative PerOrbital

getPO :: PerOrbital a -> N -> L -> a

indices :: PerOrbital (Maybe (N,L))

data Atom = Atom
    {
        atomicNumber :: Int,
        massNumber :: Int,
        charge :: Int,
        atomGrid :: Grid,
        forcedOccs :: PerOrbital Int,
        energies :: PerOrbital Energy,
        prevEnergies :: PerOrbital Energy,
        smoothingFactors :: PerOrbital Double,
        orbitals :: PerOrbital Orbital,
        occupations :: PerOrbital Double
    }


electronArrangement :: Atom -> [(N, L, Int)]
