module OrbitalShape(
    asymptoticChart
) where

import Polynomial
import Orbital
import Atom
import Data.Complex
import Data.List

normalise :: CPoly4 -> CPoly4
normalise p = (p *) $ constant $ recip $ sqrt $ evaluate [1] $ sphericalIntegral $ p * conjPoly4 p

orbitalShape :: L -> L -> L -> CPoly4
orbitalShape l m1 m2 = normalise $ (*angle) $ sum $ zipWith3 ((.) (*) . (*)) xys wzs $ unfoldr nextCoeff (1, 0)
    where [m1', m2'] = map abs [m1, m2]
          l'         = if mod (l-m1'-m2') 2 == 0 then div (l-m1'-m2') 2 else error ("Invalid magnetic quantum numbers: l="++show l++", m1="++show m1++", m2="++show m2)
          nextCoeff (x, i) = if i <= l'
              then Just $ (constant x, (-x/fromIntegral ((i+1)*(i+m1'+1))*fromIntegral ((l'-i)*(l'+m2'-i)), i+1))
              else Nothing
          angle = (variable (if m1 > 0 then 1 else 0))^m1' * (variable (if m2 > 0 then 3 else 2))^m2'
          xys = iterate (* (variable 0 * variable 1)) 1
          wzs = reverse $ take (l'+1) $ iterate (* (variable 2 * variable 3)) 1

exponential :: Double -> Double -> CPoly4
exponential a r = toHopfBasis $ sum $ take n terms
    where n     = min 20 $ floor $ 2 + exp 1 * r * a
          terms = scanl (\x i -> x*(variable 0)* constant ((a:+0)/fromIntegral i)) 1 [1..]

-- if (asymptoticDecline atom n l == (a, k)), the n, l orbital in the atom has asymptotic form ke^(-ar)r^(-1.5).
asymptoticDecline :: Atom -> N -> L -> (Double, Double)
asymptoticDecline atom n l = (a, k)
    where a      = sqrt $ (-2) * energies atom !! l !! (n-1)
          orb    = orbitals atom !! l !! (n-1)
          (ψ, r) = (!!30) $ reverse $ zip orb (atomGrid atom)
          k      = ψ / (r**(-1.5) * exp (-a*r))

-- Uses the expression for the asymptotic decline to approximate what the orbital would do. (Tweaking the energy at regular intervals would be a much better way of doing this, but also far more complicated.
extendedOrbital :: Atom -> N -> L -> Orbital
extendedOrbital atom n l = replace (dropEnd 30 orb) orbTail
    where (a, k)      = asymptoticDecline atom n l
          orbTail     = map (\r -> k * r**(-1.5) * exp (-a*r)) $ atomGrid atom
          orb         = orbitals atom !! l !! (n-1)
          replace x y = zipWith id (map const x ++ repeat id) y
          dropEnd n x = zipWith const x (drop n x)

asymptoticChart :: Atom -> N -> L -> [(Double, Double)]
asymptoticChart atom n l = map (\r -> (r, k*exp (-a*r) * r**(-0.5))) $ atomGrid atom
    where (a, k) = asymptoticDecline atom n l

exponentialOverlap :: Grid -> Orbital -> CPoly4 -> (Double, Double) -> [Double]
exponentialOverlap rs ψs ψl (a, k) = zipWith valueAtRadius rs ψs
    where ψl' = exponential a (snd $ last $ zip ψs $ takeWhile (<400) rs) * constant (k :+ 0)
          ψψl = ψl * ψl'
          l   = degree ψl
          shell = sphericalIntegral ψψl
          valueAtRadius r ψ = realPart (evaluate [r :+ 0] shell) * ψ * r^^(-l-3)

asymptoticOverlap :: Atom -> N -> L -> Atom -> N -> L -> L -> L -> (Double, [Double])
asymptoticOverlap atm0 n0 l0 atm1 n1 l1 m1 m2 = (a, exponentialOverlap rs ψs ψl (a, k))
    where rs = atomGrid atm1
          ψs = extendedOrbital atm1 n1 l1
          ψl = orbitalShape l1 m1 m2
          (a, k) = asymptoticDecline atm0 n0 l0

-- TODO: trim the integral over the extended orbital down to finite size.
asymOverlap :: Atom -> N -> L -> Atom -> N -> L -> L -> L -> (Double, Double, Double)
asymOverlap atm0 n0 l0 atm1 n1 l1 m1 m2 = (a, integrate ov, integrate (zipWith (*) ov vs))
    where (a, ov)   = asymptoticOverlap atm0 n0 l0 atm1 n1 l1 m1 m2
          rs        = atomGrid atm1
          vs        = map (min 0) $ getPotential atm1 n1 l1 Nothing -- Exclude the negative part as (to a crude approximation), the other orbital will be repelled and won't penetrate. It will have additional kinetic energy from this, but that's ignored.
          integrate = trimSum 0 0 0 . zipWith3 (\r r' x -> r^3*x*(r'-r)) rs (tail rs)
          trimSum s m a [] = s
          trimSum s m a (x:xs) = let a' = a + 0.03*(abs x - a)
                                     m' = max m (abs x)
                                     s' = s + x
                                 in if a < 0.01*m then s' else trimSum s' m' a' xs

bondEnergy :: Bool -> Atom -> N -> L -> L -> L -> Atom -> N -> L -> L -> L -> (Double, Double, Double)
bondEnergy sgnB atm0 n0 l0 m10 m20 atm1 n1 l1 m11 m21 = ((a0+a1)/2, sgn*hov/2, hov * ov)
    where (a0, ov0, hov0) = asymOverlap atm1 n1 l1 atm0 n0 l0 m10 m20
          (a1, ov1, hov1) = asymOverlap atm0 n0 l0 atm1 n1 l1 m11 m21
          ov  = ov0 + ov1
          hov = hov0 + hov1
          sgn = if sgnB then 1 else -1
