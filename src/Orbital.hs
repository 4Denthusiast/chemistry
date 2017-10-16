module Orbital(
    logGrid,
    basePotential,
    trimmedOrbital,
    realOrbital,
    normalizedOrbital,
    principalNumber,
    findEnergy,
    findEnergyHinted,
    coulumb'sConstant,
    Orbital,
    Potential,
    Grid,
    Energy,
    L,N,Spin(..),
    dropWhileEndExcept,
) where

import Dual
import Data.List
import Data.Bifunctor (first)
import Data.Bits.Floating.Ulp

import Debug.Trace
traceShow1 :: Show a => a -> b -> b
traceShow1 = const id

type Potential = ([Double], [Double]) -- The first list is the actual potential. The second if for the exchange energy, which can't be described simply as a potential.
type Grid = [Double] -- radius sampling points
type DOrbital = [(DDouble, DDouble)]
type Orbital = [Double]
type Energy = Double
type L = Int
type N = Int
data Spin = Up | Down deriving (Eq, Ord, Enum)
instance Show Spin where{show Up = "↑"; show Down = "↓"}

logGrid :: Double -> Double -> Grid
logGrid spacing z = iterate (* (1 + spacing)) (log (slimConstant/coulumb'sConstant) - 3/(sqrt z +1))

-- The effective potential purely from nuclear charge and angular momentum
basePotential :: Grid -> Double -> Double -> Potential
basePotential grid l charge = (map (\r -> l*(l+2) * r^^(-2) * 0.5 + coulumbForce r) grid, repeat 0)
    where coulumbForce r = charge * (-coulumb'sConstant + slimConstant / exp r) / (r*r)

coulumb'sConstant, slimConstant :: Double
coulumb'sConstant = 0.5
slimConstant = 4

-- Calculates the orbital with the fixed energy value.
singleOrbital :: Grid -> Potential -> Energy -> DOrbital
singleOrbital rs (vs,xs) e = scanl progress (0, 1) (zip4 rs (zipWith (-) (tail rs) rs) vs xs)
    where progress (ψ, dψ) (r, dr, v, x) = (ψ', dψ')
              where ddψ = ((dual v) - (Dual e 1))*ψ*2 - (dual x)*2 - 3/(dual r) * dψ
                    dψ' = (dual dr) * ddψ  + dψ
                    ψ'  = (dual dr) *  dψ' +  ψ

trimmedOrbital :: Grid -> Potential -> Energy -> DOrbital
trimmedOrbital rs vs e = trimOrbital rs vs $ singleOrbital rs vs e

trimOrbital :: Grid -> Potential -> DOrbital -> DOrbital
trimOrbital rs (vs,_) = map snd . trim3 . trim2 . trim1 . zip rs
    where trim1 rψs = takeWhile (\(r, (ψ, _)) -> abs (std ψ) * r < 3*threshold rψs) rψs
          thresholdPoint = (2*) $ snd $ head $ filter (\(v, r) -> v < 0 || r > 8) $ zip vs rs
          threshold = maximum . map (\(r, (ψ, _)) -> abs (std ψ) * r) . takeWhile ((< thresholdPoint) . fst)
          trim2 = takeWhile (\(r, (ψ, dψ)) -> r < 2 || r^2 * (3*(std ψ)^2 + (std ψ + std dψ * r)^2) > 1e-4)
          trim3 = dropWhileEndExcept (\(r, (ψ, dψ)) -> ψ * dψ > 0 && r > 3) 150

dropWhileEndExcept :: (a -> Bool) -> Int -> [a] -> [a]
dropWhileEndExcept f n = (\(t, h) -> reverse h ++ take n (reverse t)) . span f . reverse

principalNumber :: Grid -> DOrbital -> N
principalNumber rs ψs = fst $ foldl
    (\(n, f) ψ ->
        if f && ψ < (-0.01) then
            (n+1, False)
        else if not f && ψ > 0.01 then
            (n+1, True)
        else
            (n, f)
    )
    (0, False) (zipWith ((*) . std . fst) ψs rs)

orbitalExists :: Grid -> Potential -> N -> Bool
orbitalExists rs vs n = principalNumber rs orb >= n && (not $ null $ filter (<0) $ map fst $ orb)
    where orb = trimmedOrbital rs vs 0

findEnergy :: Grid -> Potential -> N -> Energy
findEnergy = findEnergyHinted (-1)

findEnergyHinted :: Energy -> Grid -> Potential -> N -> Energy
findEnergyHinted e0 rs vs n = let
        e0' = min e0 (-1e-16)
        approxEq a b = abs (a-b) < 1e-12 * min (abs a) (abs b)
        iter :: Energy -> Energy -> Energy -> Energy
        iter e ll hh = traceShow1 e $ let
                orbital      = trimmedOrbital rs vs e
                ((ψ, dψ), r) = last $ zip orbital rs
                roundError   = doubleUlp e * ((\(Dual _ x) -> x) ψ)
                roundWarning = if roundError * r > 0.01 then trace "Rounding error." else id
                err          = ψ + (dual r*0.2)*dψ
                e'           = let (Dual era ere) = err in max (e*1.6) $ min (e*0.6) $ e - era/ere
                n1           = principalNumber rs orbital
                e''
                    | e > -1e-20             = 0
                    | approxEq e' e          = e
                    | n1 < n                 = traceShow1 n1 $ traceShow1 n $ iter (0.6 * max e ll) (0.8 * max e ll) hh
                    | n1 > n                 = traceShow1 n1 $ traceShow1 n $ iter (1.5 * min e hh) ll (1.2 * min e hh)
                    | otherwise              = iter e' ll hh
            in roundWarning e''
    in
        if orbitalExists rs vs n then iter e0' (-1/0) 0 else 0

realOrbital :: Grid -> Potential -> Energy -> Orbital
realOrbital rs vs e = map (std . fst) $ trimmedOrbital rs vs e

normalizedOrbital :: Grid -> Potential -> Energy -> Orbital
normalizedOrbital rs vs e = map (/ scale) o
    where o     = realOrbital rs vs e
          o'    = zipWith3 (\ψ r r' -> ψ*ψ*r*r*r*(r'-r)) o rs (tail rs)
          scale = sqrt $ sum o'
