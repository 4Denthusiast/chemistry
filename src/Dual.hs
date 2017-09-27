module Dual(
    Dual(..),
    std,
    dual,
    DDouble,
) where

data Dual n = Dual n n deriving (Eq)
std :: Dual n -> n
std (Dual a e) = a

type DDouble = Dual Double

dual :: (Num n) => n -> Dual n
dual n = Dual n 0

instance (Eq n, Show n, Num n) => Show (Dual n) where
    show (Dual a e)
        | e == 0    = show a
        | a == 0    = show e ++ "ε"
        | otherwise = show a ++ " + " ++ show e ++ "ε"

instance (Ord a) => Ord (Dual a) where
    compare (Dual a e) (Dual a' e') = mappend (compare a a') (compare e e')

instance (Eq n, Num n) => Num (Dual n) where
    (Dual a e) + (Dual a' e') = Dual (a + a') (e + e')
    negate (Dual a e) = Dual (- a) (- e)
    (Dual a e) * (Dual a' e') = Dual (a * a') (a * e' + e * a')
    signum (Dual a e)
        | a == 0    = dual $ signum e
        | otherwise = dual $ signum a
    abs (Dual a e)
        | a == 0    = Dual 0 (abs e)
        | otherwise = let sa = signum a in Dual (a*sa) (e*sa) -- This assumes that sa^2 = 1.
    fromInteger = dual . fromInteger

instance (Real n) => Real (Dual n) where
    toRational = toRational . std

instance (Eq n, Fractional n) => Fractional (Dual n) where
    recip (Dual a e) = Dual (1/a) (-e/(a*a))
    fromRational = dual . fromRational
