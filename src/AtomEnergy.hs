module AtomEnergy(
    totalEnergy,
    orbitalEnergy,
    approxOrbitalEnergy,
    splitSpin,
    carefulEnergies
) where

import {-# SOURCE #-} Atom
import Orbital
import Utils

import Data.List


totalEnergy :: Atom -> Energy
totalEnergy atom = sum nuclearEnergies + sum eeEnergies
    where ea              = concatMap splitSpin $ electronArrangement atom
          nuclearEnergies = map (nuclearEnergy atom) ea
          eeEnergies      = concat $ zipWith (zipWith (eeEnergy atom) . repeat) ea (tails ea)

splitSpin :: (N, L, Int) -> [(N, L, Spin, Int)]
splitSpin (n, l, o)
    | o <= ll   = [(n, l, Up, o)]
    | otherwise = [(n, l, Up, ll), (n, l, Down, o-ll)]
        where ll = (l+1)^2

nuclearEnergy :: Atom -> (N, L, Spin, Int) -> Energy
nuclearEnergy atom (n, l, _, o) = fromIntegral o * sum (zipWith4 (\r r' hψ ψ -> r*r*r*(r'-r)*hψ*ψ) rs (tail rs) (oneElectronHamiltonian rs z l ψs) ψs)
    where rs = atomGrid atom
          z  = atomicNumber atom
          ψs = getPO (orbitals atom) n l

oneElectronHamiltonian :: Grid -> Int -> L -> Orbital -> Orbital
oneElectronHamiltonian rs z l ψs = zipWith (+) potentialTerm laplacianTerm
    where vs            = fst $ basePotential rs (fromIntegral l) (fromIntegral z)
          potentialTerm = zipWith (*) vs ψs
          rds           = zipWith (-) (tail rs) rs
          dup (x:xs)    = x:x:xs
          dψs           = dup $ zipWith (/) (zipWith (-) (tail  ψs)  ψs) rds -- The first term would be 1, but these are normalised.
          ddψs          =       zipWith (/) (zipWith (-) (tail dψs) dψs) rds
          laplacianTerm = zipWith3 (\ddψ dψ r -> (-ddψ - 3/r * dψ)/2) ddψs dψs rs

eeEnergy :: Atom -> (N, L, Spin, Int) -> (N, L, Spin, Int) -> Energy
eeEnergy atom (n0, l0, s0, o0) (n1, l1, s1, o1) = coulumb'sConstant * sum (zipWith4 (\r r' c x -> r*r*r*(r'-r)*(c-x)) rs (tail rs) coulumbTerm exchangeTerm)
    where rs   = atomGrid atom
          ifDiag a b = if n0 == n1 && l0 == l1 && s0 == s1 then a else b
          (o0', o1') = (fromIntegral o0, fromIntegral o1)
          ψs0  = getPO (orbitals atom) n0 l0
          ψs1  = getPO (orbitals atom) n1 l1
          ψs00 = zipWith (*) ψs0 ψs0
          ψs01 = zipWith (*) ψs0 ψs1
          ψs11 = zipWith (*) ψs1 ψs1
          coulumb01 = zipWith (*) ψs00 $ zipWith (\r q -> q/(r*r)) rs $ scanl (+) 0 $ zipWith3 (\r r' d -> r*r*r*(r'-r)*d) rs (tail rs) ψs11
          coulumb10 = zipWith (*) ψs11 $ zipWith (\r q -> q/(r*r)) rs $ scanl (+) 0 $ zipWith3 (\r r' d -> r*r*r*(r'-r)*d) rs (tail rs) ψs00
          coulumbTerm = map (* ifDiag (o0'*(o1'-1)/2) (o0'*o1')) $ zipWith (+) coulumb01 coulumb10
          ls    = [l0 + l1, l0 + l1 - 2 .. abs (l0 - l1)]
          exchange l = zipWith (*) ψs01 $ zipWith (\r q -> q/(r*r*r^l)) rs $ scanl (+) 0 $ zipWith3 (\r r' d -> r*r*r*r^l*(r'-r)*d) rs (tail rs) ψs01
          exchangeTerm = if s0 /= s1 then repeat 0 else map (* ((o0'*o1' - ifDiag o0' 0)*ifDiag 1 2/(fromIntegral ((l0+1)*(l1+1))))) $ map sum $ transpose $ map exchange ls

orbitalEnergy :: Atom -> (N, L, Spin) -> Energy
orbitalEnergy atom (n, l, s) = nuclearEnergy atom (n, l, s, 1)
    + sum (map (eeEnergy atom (n,l,s,1)) $ concatMap splitSpin $ electronArrangement atom)

carefulEnergies :: Atom -> PerOrbital [Energy]
carefulEnergies atom = map (orbitalEnergy atom . \(n,l,s,o)->(n,l,s)) <$> splitSpin <$> (maybe (const (0,0,0)) (uncurry (,,)) <$> indices <*> (round <$> occupations atom))

approxOrbitalEnergy :: Atom -> (N, L, Spin) -> Energy
approxOrbitalEnergy atom (n, l, s) = getPO (energies atom) n l
