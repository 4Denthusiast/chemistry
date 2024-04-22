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

type Potential = ([Double], [Double], Double, Double) -- The first list is the actual potential. The second if for the exchange energy, which can't be described simply as a potential. The last two elements are the values A and B such that V = (A+rB+O(r^2))r^-2 near the origin.
type Grid = [Double] -- radius sampling points
type DOrbital = [(DDouble, DDouble)]
type Orbital = [Double]
type Energy = Double
type L = Int
type N = Int
data Spin = Up | Down deriving (Eq, Ord, Enum)
instance Show Spin where{show Up = "↑"; show Down = "↓"}

logGrid :: Double -> Double -> Grid
logGrid spacing z = iterate (* (1 + spacing)) ((4*c*spacing^3)**(1/c))
    where c = 2*sqrt(1+2*z*coulumb'sConstant)+2
-- Ideally, r0 ~= (spacing^3*2(2f+4)A/B)^(1/(2f+4)) (in the notation of the comment below), for an optimal combination of step size and start point. This is further approximated with A=B/2=kZ.

-- The effective potential purely from nuclear charge and angular momentum
basePotential :: Grid -> Double -> Double -> Double -> Potential
basePotential grid l charge sCharge = (map (\r -> l*(l+2) * r^^(-2) * 0.5 + coulumbForce r) grid, repeat 0, sCharge*slimConstant - charge*coulumb'sConstant + l*(l+2)/2, -sCharge*slimConstant)
    where coulumbForce r = (-charge * coulumb'sConstant + sCharge * slimConstant / exp r) / (r*r)

coulumb'sConstant, slimConstant :: Double
coulumb'sConstant = 0.9
slimConstant = coulumb'sConstant

{-Assume that, in the vicinity of 0, Ψ can be approximated as
Ψ = r^f(1 + g r + O(r^2))
and V can be approximated as
V = (A+rB)r^-2
(i.e. A = (s charge-charge)*k+l(l+2)/2, and B = -s charge*k)
then

Ψ' = fr^(f-1)(1+gr+...) + r^f(g+...) = r^(f-1)(f+(f+1)gr+...)
Ψ'' = (f-1)r^(f-2)(f+(f+1)gr...) + r^(f-1)((f+1)g+...) = r^(f-2)(f(f-1)+f(f+1)gr+...)

∇^2 Ψ/2 = VΨ (the energy term is negligible here)
1/2 (Ψ'' + 3/r Ψ') = (A+rB)r^-2 Ψ
1/2 (f(f-1) + f(f+1)gr + 3(f + (f+1)gr)) = (A+rB)(1+gr)
ff + 2f + r(f+3)(f+1)g = 2A + r(2B+2Ag) + O(r^2)

ff+2f-2A = 0, f = (-2 +- sqrt(4+8A))/2 = -1+-sqrt(1+2A)
The positive root is the correct one: f = sqrt(1 + 2A) - 1
(ff+4f+3)g-2Ag = 2B
g(2f+3) = 2B
g = 2B/(2f+3)

therefore near the origin, Ψ'/Ψ = f/r + f+g + O(r)

Alternatively, if we use Ψ = r^f e^(hr+O(r^2))
Ψ' = (f+hr)r^(f-1)e^hr
Ψ'' = h r^(f-1)e^hr + (f+hr)(f-1+hr)r^(f-2)e^hr
  = (ff-f+2fhr+3(f+hr)+O(rr))r^(f-2)e^hr
1/2 (ff+2f+2fhr+3hr) = A+Br
2fh+3h=2B
h=2B/(2f+3)

For Ψ = r^f e^hr exactly, the potential implied (for 0 energy) is
∇^2 Ψ/2 = VΨ
1/2 (ff+2f+2fhr+3hr+hhrr)r^(f-2)e^hr = V r^f e^hr
V = (A+Br+hhrr/2)r^-2
= (A+Br)r^-2 + 2BB/(2f+3)^2
= (A+Br)r^-2 + 2BB/(4ff+12f+9)
= (A+Br)r^-2 + 2BB/(8A+4sqrt(1+2A)+5)
whereas the true V-E is
(A+Br)r^-2 - B + screening - E + O(r)
The implied error in energy is therefore
2BB/(8A+4sqrt(1+2A)+5) + B - screening + E + O(r)
The B term is dominant (or at least similar) except in the case of high atomic number (> about 3) with very few neutrons, which this program deals with poorly for other reasons anyway.
-}

-- Calculates the orbital with the fixed energy value.
singleOrbital :: Grid -> Potential -> Energy -> DOrbital
singleOrbital rs (vs,xs,a,b) e = scanl progress (1, dual $ f/head rs + g) (zip5 rs drs (head drs:drs) vs xs)
    where f = sqrt (1 + 2*a) - 1
          g = 2*b/(2*f+3)
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
