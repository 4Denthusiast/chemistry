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
import Utils
import Data.List
import Data.Bifunctor (first)
import Data.Bits.Floating.Ulp

import Debug.Trace
import Data.IORef
import System.IO.Unsafe
traceShow1 :: Show a => a -> b -> b
traceShow1 = const id

type Potential = ([Double], [Double], Double, Double) -- The first list is the actual potential. The second is for the exchange energy, which can't be described simply as a potential. The last two elements are the values A and B such that V = A/rr + B log r + O(1) near the origin.
type Grid = [Double] -- radius sampling points
type DOrbital = [(DDouble, DDouble)]
type Orbital = [Double]
type Energy = Double
type L = Int
type N = Int
data Spin = Up | Down deriving (Eq, Ord, Enum)
instance Show Spin where{show Up = "↑"; show Down = "↓"}

logGrid :: Double -> Double -> Grid
logGrid spacing z = iterate (* (1 + spacing)) ((a*c*spacing^3)**c)
    where a = z*coulumb'sConstant
          c = 1/(2*sqrt(1+2*a)+1)
-- The energy error due to the r<r0 approximate analytic solution is about B r0^(2f+3), in the notation of the comment below. The discretization error is about A^2 spacing^3 / 12. These are optimally balanced when (2f+3)B r0^(2f+3) ~= A^2 spacing^3. This is further approximated with A=B=kZ.
-- TODO: This indicates that the discretization error increases quite rapidly with Z. Consider using a smaller spacing for high Z.

-- The effective potential purely from nuclear charge and angular momentum
basePotential :: Grid -> Double -> Double -> Double -> Potential
basePotential grid l charge sCharge = (map (\r -> l*(l+2) * r^^(-2) * 0.5 + coulumbForce r) grid, repeat 0, sCharge*slimConstant - charge*coulumb'sConstant + l*(l+2)/2, sCharge*slimConstant/2)
    where coulumbForce r = -charge * coulumb'sConstant / (r*r) + sCharge * slimConstant * besselK1 r / r

coulumb'sConstant, slimConstant :: Double
coulumb'sConstant = 0.9
slimConstant = coulumb'sConstant

{-
Let V = A/rr + B log r + C + O(r log r) for small r,
Let f = sqrt(1 + 2A) - 1 (so that A = ff/2 + f), and h = B/(2f+4)
Let Ψ = r^(f + hrr) = r^f e^(hrr log r)
Ψ' = (f/r + h (2r log r + r)) Ψ
Ψ'' = (-f/rr + h (2log r + 3) + ff/rr + 2fh(2log r + 1) + hhrr(4(log r)^2 + 4log r + 1)) Ψ
∇^2 Ψ/2 = 1/2 (Ψ'' + 3/r Ψ') = VΨ - E
energy error = 1/2 (3h(2log r + 1) + h(2log r + 3) + 2fh(2log r + 1)) - B log r - C + E + O(r log r)
= h/2 (8+4f) log r - B log r + h (3 + f) - C + E + O(r log r)
~= B/2 - C + E
The B and C terms are similar in magnitude, adding up to about 0.9 k_s * mass number. This neglects screening, but that should also be much smaller.
-}

-- Calculates the orbital with the fixed energy value.
singleOrbital :: Grid -> Potential -> Energy -> DOrbital
singleOrbital rs (vs,xs,a,b) e = scanl progress (1, dual $ f/r0 + h*(2*r0*log r0 + r0)) (zip5 rs drs (head drs:drs) vs xs)
    where f = sqrt (1 + 2*a) - 1
          h = b/(2*f+4)
          r0 = head rs
          drs = zipWith (-) (tail rs) rs
          progress (ψ, dψ) (r, dr, prevDr, v, x) = (ψ', dψ')
              where ddψ = (((dual v) - (Dual e 1))*ψ*2 - (dual x)*2 - 3/(dual r) * dψ)/(1 + 3/2*dual(dr/r)) --The division at the end implies the estimate of the derivative at r is about (dψ+dψ')/2.
                    dψ' = (dual (dr+prevDr)/2) * ddψ  + dψ
                    ψ'  = (dual dr) *  dψ' +  ψ

trimmedOrbital :: Grid -> Potential -> Energy -> DOrbital
trimmedOrbital rs vs e = trimOrbital rs vs $ singleOrbital rs vs e

trimOrbital :: Grid -> Potential -> DOrbital -> DOrbital
trimOrbital rs (vs,_,_,_) = map snd . trim3 . trim2 . trim1 . zip rs
    where trim1 rψs = takeWhile (\(r, (ψ, _)) -> abs (std ψ) * r < 3*threshold rψs) rψs
          thresholdPoint = (2*) $ þrd3 $ labHead "trimOrbital" $ filter (\(v, v', r) -> r > 1 && v < v' || r > 8) $ zip3 vs (tail vs) rs
          threshold = maximum . map (\(r, (ψ, _)) -> abs (std ψ) * r) . takeWhile ((< thresholdPoint) . fst)
          trim2 = takeWhile (\(r, (ψ, dψ)) -> r < 2 || r^2 * (3*(std ψ)^2 + (std ψ + std dψ * r)^2) > 1e-4)
          trim3 = dropWhileEndExcept (\(r, (ψ, dψ)) -> ψ * dψ > 0 && r > 3) 150

dropWhileEndExcept :: (a -> Bool) -> Int -> [a] -> [a]
dropWhileEndExcept f n = (\(t, h) -> reverse h ++ take n (reverse t)) . span f . reverse

principalNumber :: Grid -> Energy -> DOrbital -> N
principalNumber rs e = length . group . map (\(r,(ψ, dψ)) -> inf ψ + (1 / (sqrt (-2*e)+1.5/r))* inf dψ > 0) . zip rs

orbitalExists :: Grid -> Potential -> N -> Bool
orbitalExists rs vs n = principalNumber rs 0 orb >= n && (not $ null $ filter (<0) $ map fst $ orb)
    where orb = trimmedOrbital rs vs 0

findEnergy :: Grid -> Potential -> N -> Energy
findEnergy = findEnergyHinted (-1)

findEnergyHinted :: Energy -> Grid -> Potential -> N -> Energy
findEnergyHinted e0 rs vs n = traceShow1 (e0, n) $ seq (stashPotential vs) $ let
        e0' = min e0 (-1e-16)
        approxEq a b = abs (a-b) < 1e-12 * min (abs a) (abs b)
        iter :: Energy -> Energy -> Energy -> Energy
        iter e ll hh = traceShow1 [e,ll,hh] $ let
                orbital      = trimmedOrbital rs vs e
                ((ψ, dψ), r, v) = last $ zip3 orbital rs (fst4 vs)
                roundError   = doubleUlp e * (inf ψ)
                roundWarning = if roundError * r > 0.01 then trace "Rounding error." else id
                llhhError    = if ll > hh then error "Energy interval of negative size: ll > hh" else ()
                err          = ψ + dual (1 / (sqrt (-2*e)+1.5/r))*dψ
                e'           = let (Dual era ere) = err in max (e*1.6) $ min (e*0.6) $ e - era/ere
                (ll', hh')   = case mconcat [ --Calculate by cases whether the energy is too high or too low.
                                   markOrd "n" $ compare n1 n,
                                   --markOrd "s" $ max EQ $ compare ( ψ * ((-1)^n)) 0,
                                   --markOrd "s" $ min EQ $ compare (dψ * ((-1)^n)) 0,
                                   markOrd "e" $ compare e e'
                               ] of
                                   LT -> (max ll e,hh)
                                   GT -> (ll,min hh e)
                mid          = if isInfinite ll then 2*hh' else (hh'+ll')/2
                n1           = principalNumber rs e orbital
                e''
                    | e > -1e-20     = 0
                    | approxEq e' e  = e
                    | approxEq ll hh = ll
                    | max e e' >= hh = traceShow1 'h' $ iter mid ll' hh' -- If it's no better than it was before, default to binary search.
                    | min e e' <= ll = traceShow1 'l' $ iter mid ll' hh'
                    | n1 < n         = traceShow1 n1 $ traceShow1 n $ iter (0.3^(n-n1) * max e ll) ll' hh'
                    | n1 > n         = traceShow1 n1 $ traceShow1 n $ iter (2.5^(n1-n) * min e hh) ll' hh'
                    | otherwise      = iter e' ll' hh'
            in seq llhhError $ roundWarning e''
    in
        if orbitalExists rs vs n then iter e0' (-1/0) 0 else 0

markOrd :: String -> Ordering -> Ordering
markOrd s o  = case o of {LT -> traceShow1 (s++"<") LT; GT -> traceShow1 (s++">") GT; EQ -> EQ}

potentialBox :: IORef Potential
potentialBox = unsafePerformIO $ newIORef ([],[],0,0)

stashPotential :: Potential -> ()
stashPotential = unsafePerformIO . writeIORef potentialBox

getStashPotential :: () -> Potential
getStashPotential () = unsafePerformIO $ readIORef potentialBox

realOrbital :: Grid -> Potential -> Energy -> Orbital
realOrbital rs vs e = map (std . fst) $ trimmedOrbital rs vs e

normalizedOrbital :: Grid -> Potential -> Energy -> Orbital
normalizedOrbital rs vs e = map (/ scale) o
    where o     = realOrbital rs vs e
          o'    = zipWith3 (\ψ r r' -> ψ*ψ*r*r*r*(r'-r)) o rs (tail rs)
          scale = sqrt $ sum o'

aimLengths :: DOrbital -> [Double]
aimLengths = map (std . uncurry (/))
