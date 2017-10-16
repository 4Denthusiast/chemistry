{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE KindSignatures      #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE FlexibleInstances   #-}
module Polynomial(
    Polynomial,
    constant,
    variable,
    degree,
    evaluate,
    laplacian,
    Poly4,
    CPoly4,
    CPoly1,
    toHopfBasis, fromHopfBasis,
    sphericalIntegral,
    conjPoly4,
    poly4Norm,
    assumePolyReal
) where

import Utils
import Data.List
import Data.Complex
import GHC.TypeLits --(natVal, Nat, KnownNat)

-- Polynomials in n variables.
data Polynomial (n :: Nat) a = Polynomial [Monomial n a] deriving (Eq)
data Monomial   (n :: Nat) a = Monomial{
    monCoefficient :: a,
    monExponents   :: [Int]
} deriving (Eq)

monomialTimes :: Num a => Monomial n a -> Monomial n a -> Monomial n a
monomialTimes (Monomial x ex) (Monomial y ey) = Monomial (x*y) (zipWith (+) ex ey)

instance Functor (Monomial n) where
    fmap f (Monomial a e) = Monomial (f a) e

simplifyPoly :: (Eq a, Num a) => [Monomial n a] -> Polynomial n a
simplifyPoly = Polynomial . filter ((/=0) . monCoefficient) . mergeBy mergeMon . sortBy compareExponents
    where compareExponents (Monomial _ e0) (Monomial _ e1) = mappend (compare (sum e1) (sum e0)) (compare e1 e0)
          mergeMon (Monomial x ex) (Monomial y ey)
              | ex == ey  = Just $ Monomial (x+y) ex
              | otherwise = Nothing

data NatProxy (n::Nat)
constant :: forall a n. (KnownNat n, Eq a, Num a) => a -> Polynomial n a
constant x = simplifyPoly [Monomial x $ genericReplicate (natVal (undefined :: NatProxy n)) 0]

variable :: forall a n. (KnownNat n, Eq a, Num a) => Int -> Polynomial n a
variable i = Polynomial [Monomial 1 es]
    where n  = fromInteger $ natVal (undefined :: NatProxy n)
          es = replicate i 0 ++ 1 : replicate (n-i-1) 0

degree :: Polynomial n a -> Int
degree (Polynomial []) = -1 -- negative infinity is the usual convention, but making this an Int seems more appropriate.
degree (Polynomial (Monomial _ es : _)) = sum es

instance (Eq a, Num a, KnownNat n) => Num (Polynomial n a) where
    (Polynomial xs) + (Polynomial ys) = simplifyPoly (xs ++ ys)
    negate (Polynomial xs)            = Polynomial (map (fmap negate) xs)
    (Polynomial xs) * (Polynomial ys) = simplifyPoly (monomialTimes <$> xs <*> ys)
    abs _ = error "The typeclass Num doesn't contain a method: conjugate"
    signum _ = error "Signum doesn't make much sense for polynomials."
    fromInteger = constant . fromInteger

evaluate :: Num a => [a] -> Polynomial n a -> a
evaluate xs (Polynomial ms) = sum (map evaluateMonomial ms)
    where evaluateMonomial (Monomial a es) = a*product (zipWith (^) xs es)

substitute :: (Eq a, Num a, KnownNat n, KnownNat m) => [Polynomial n a] -> Polynomial m a -> Polynomial n a
substitute xs (Polynomial ms) = evaluate xs (Polynomial ms')
    where ms' = map (fmap constant) ms

laplacian :: (KnownNat n, Eq a, Num a) => Polynomial n a -> Polynomial n a
laplacian (Polynomial ms) = sum $ map monLap ms
    where monLap (Monomial a es) = Polynomial $ zipWith3 (\h t e -> Monomial (a*fromIntegral (e*(e-1))) (h++(e-2):t)) (inits es) (tail $ tails es) es


type Poly4 = Polynomial 4
type CPoly4 = Poly4 (Complex Double)
type CPoly1 = Polynomial 1 (Complex Double)

toHopfBasis, fromHopfBasis :: CPoly4 -> CPoly4
[toHopfBasis, fromHopfBasis] = map substitute [map (*constant 0.5) [x+y, i*(x-y), z+w, i*(z-w)], [x-i*y, x+i*y, z-i*w, z+i*w]]
    where [x,y,z,w] = map variable [0..3]
          i         = constant (0 :+ 1)

-- The integral over a sphere of radius r of the polynomial, which should be in the hopf basis.
sphericalIntegral :: CPoly4 -> CPoly1
sphericalIntegral (Polynomial ms) = constant ((2*pi)^2) * sum (map monIntegral ms)
    where monIntegral :: Monomial 4 (Complex Double) -> CPoly1
          monIntegral (Monomial x [a,a',b,b'])
              | a == a' && b == b' = (variable 0)^(2*a+2*b+3) * constant (x / fromIntegral (binomial a b) / 2 / fromIntegral (1+a+b))
              | otherwise          = 0

binomial :: Integral a => a -> a -> a
binomial a b
    | a > b     = binomial b a
    | a == 0    = 1
    | otherwise = div (binomial (a-1) (b+1) * (b+1)) a

conjPoly4 :: CPoly4 -> CPoly4
conjPoly4 (Polynomial ms) = Polynomial $ map conj ms
    where conj (Monomial a [x, y, z, w]) = Monomial (conjugate a) [y, x, w, z]

poly4Norm :: CPoly4 -> Polynomial 1 Double
poly4Norm p = assumePolyReal $ sphericalIntegral $ p * conjPoly4 p

assumePolyReal :: (Eq a, Num a) => Polynomial n (Complex a) -> Polynomial n a
assumePolyReal (Polynomial ms) = Polynomial $ map (fmap assumeReal) ms
    where assumeReal (x :+ 0) = x
          assumeReal _        = error "Real number required."


class NamedDimensions (n :: Nat) where
    dimName :: forall proxy. proxy n -> Int -> String

instance NamedDimensions 4 where
    dimName _ 0 = "x"
    dimName _ 1 = "y"
    dimName _ 2 = "z"
    dimName _ 3 = "w"

instance NamedDimensions 1 where
    dimName _ 0 = "r" --"r" is more likely to get used, but "x" is more general.

instance (Show a, Num a, Eq a, NamedDimensions n, KnownNat n) => Show (Polynomial n a) where
    show 0 = "0"
    show (Polynomial ms) = intercalate " + " $ map show ms

instance (Show a, Num a, Eq a, NamedDimensions n, KnownNat n) => Show (Monomial n a) where
    show (Monomial a es) = showCoeff ++ concatMap (\(i, e) -> dimName (undefined :: NatProxy n) i ++ showPower e) (filter ((>0) . snd) $ zip [0..] es)
        where showCoeff   = if a == 1 && any (/=0) es then "" else show a
              showPower e = if e == 1 then "" else '^' : show e
