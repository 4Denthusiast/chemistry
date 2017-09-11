module Atom(
    Atom(..),
    emptyAtom,
    makeAtom,
    updateAtom,
    energies',
    getPotential,
    electronArrangement,
    prevElectronArrangement,
    prettyElectronArrangement,
    angularMomentumLabels,
    occupations,
    aperiodicTable,
    anionise,
    cationise,
    atomFull,
    indices
) where

import Orbital
import AtomEnergy
import Utils
import Data.List
import Data.Function
import Data.Maybe
import Debug.Trace
{-traceShow = const id
trace = const id
traceShowId = id-}


angularMomentumLabels = "spdfghk" ++ repeat 'R'

data Atom = Atom
    {
        atomicNumber :: Int,
        charge :: Int,
        energies :: [[Energy]],
        orbitals :: [[Orbital]],
        shieldingPotentials :: [[Potential]],
        atomGrid :: Grid,
        prevOccs :: [[Int]]
    }

fnPower :: Int -> (a -> a) -> a -> a
fnPower 0 f x = x
fnPower n f x = f (fnPower (n-1) f x)

energies' :: Atom -> [[Energy]]
energies' = takeWhile (not . null) . energies

energies'' :: Atom -> [[Energy]]
energies'' atom = zipWith take (map (+ 1) occ) (energies atom)
    where occ = occupations atom ++ [0]

emptyAtom :: Int -> Int -> Atom
emptyAtom z c = Atom z c empty empty empty (logGrid 0.01 (fromIntegral z)) []
    where empty = repeat []

-- If the atomic numbers are similar, this should converge faster than starting from empty.
copyAtom :: Atom -> Int -> Int -> Atom
copyAtom atom z c = atom
    {
        atomicNumber = z,
        charge       = c,
        atomGrid     = logGrid 0.01 (fromIntegral z)
    }

makeAtom :: Int -> Int -> Atom
--makeAtom z c = fnPower iterations updateAtom $ emptyAtom z c
--    where iterations = 8 + floor (logBase (1/0.7) $ fromIntegral z)
makeAtom z c = updateAtomUntilConvergence $ emptyAtom z c

fillAtom :: Atom -> Atom
fillAtom atom = 
    let
        occ = occupations atom
        done = (all (\(es, n) -> not $ null $ drop n es) $ zip (energies atom) occ)
            && (length occ < length (energies' atom))
    in
        if done then atom else fillAtom $ addOrbital atom

{-electronArrangement :: Atom -> [(N, L, Int)]
electronArrangement atom = takeWhile ((>0) . þrd3) $ snd $
    mapAccumL (\a (n,l,o) -> (a-o, (n, l, min a o))) (atomicNumber atom - charge atom) $
    map (\(l, n) -> (n, l, 2*(l+1)^2)) $
    map snd $ sort $ filter ((<0) . fst) $ concat $
    zipWith zip (energies' atom) (map (map (\(a,b) -> (b,a))) indices)-}

{-smoothedElectronArrangement atom = map limitOrbital $ electronArrangement atom
    where limitOrbital x@(n, l, o) = (n,l,max (prev x - 1) $ min (prev x + 1) o)
          prev (n, l, _)           = (((prevOccs atom ++ repeat [])!! l) ++ repeat 0) !! (n-1)-}

electronArrangement :: Atom -> [(N, L, Int)]
electronArrangement atom = xA
    where x1 = concat $ zipWith3 (zipWith3 (const (,))) (energies' atom) indices (map (++ repeat 0) (prevOccs atom) ++ repeat (repeat 0))
          x2 = map (\((n, l), po) -> (n, l, po+1)) x1
          x3 = concatMap splitSpin x2
          x4 = map (\(n, l, s, po) -> (approxOrbitalEnergy atom (n, l, s), n, l, s, po-1)) x3
          x5 = filter (\(e, n, l, s, po) -> e<0 || po > 0) x4
          x6 = map (\(e, n, l, s, po) -> (e, n, l, po-1)) $ sort x5
          spareElectrons = atomicNumber atom - charge atom - (sum $ map (max 0 . frþ4) x6)
          change a ll po e = minimum [a, ll-po, if po < 0 then 1 else 2, if e == 0 then 0 else 2]
          x7 = map (\(e, n, l, po) -> (e, n, l, (l+1)^2, po)) x6
          x8 = snd $ mapAccumL (\a (e, n,l,ll,po) -> (a-change a ll po e, (n, l, max 0 po + change a ll po e))) spareElectrons x7
          x9 = takeWhile ((>0) . þrd3) x8
          xA = mergeBy (\(n0, l0, o0) (n1, l1, o1) -> if n0 == n1 && l0 == l1 then Just (n0,l0,o0+o1) else Nothing) x9

prevElectronArrangement :: Atom -> [(N, L, Int)]
prevElectronArrangement = concat . zipWith (zipWith (\(n, l) o -> (n, l, o))) indices . prevOccs

indices :: [[(N, L)]]
indices = [[(n, l) | n <- [1..]] | l <- [0..]]

occupations :: Atom -> [Int]
occupations atom = map length $ group $ sort $ map snd3 $ electronArrangement atom

addOrbital :: Atom -> Atom
addOrbital atom = addOrbitalAt atom (getNextOrbital $ energies' atom)

getNextOrbital :: [[Energy]] -> L
getNextOrbital [] = 0
getNextOrbital es
    | length (last es) > 1 = length es
    | otherwise            = snd $ minimum $ zip (map last es) [0..]

addOrbitalAt :: Atom -> L -> Atom
addOrbitalAt atom l = atom
        {
            energies = appendAt l e (energies atom),
            orbitals = appendAt l o (orbitals atom),
            shieldingPotentials = appendAt l p (shieldingPotentials atom)
        }
    where n  = 1 + length (energies atom !! l)
          (e, o, p) = genOrbitalAt atom n l (-1)

genOrbitalAt :: Atom -> N -> L -> Energy -> (Energy, Orbital, Potential)
genOrbitalAt atom n l e0 = (e, o, p)
    where rs = atomGrid atom
          vs = getPotential atom n l Nothing
          e  = {-trace ("n" ++ show n ++ "l" ++ show l) $-} findEnergyHinted e0 rs vs n
          o  = normalizedOrbital rs vs e
          pn = genShieldingPotential rs o
          po = ((map Just $ shieldingPotentials atom !! l) ++ repeat Nothing) !! (n-1)
          p  = (zipWith (\a b -> a + 0.7*(b-a)) pn) $ maybe (repeat 0) id po

updateAtom :: Atom -> Atom
updateAtom atom = fillAtom $ let
        (es', ψs', vs') = unzip3 $ map unzip3 $ zipWith (zipWith $ uncurry $ genOrbitalAt atom) indices (energies'' atom)
    in atom
        {
            energies = es' ++ repeat [],
            orbitals = ψs' ++ repeat [],
            shieldingPotentials = vs' ++ repeat [],
            prevOccs = map (map þrd3 . sort) $ groupBy (on (==) snd3) $ sortBy (on compare snd3) $ electronArrangement atom
        }

updateAtomUntilConvergence :: Atom -> Atom
updateAtomUntilConvergence atom = if atomsSimilar atom next then next else updateAtomUntilConvergence next
    where next = {-traceShow (energies' atom) $ {-traceShow (carefulEnergies atom) $-} trace (prettyElectronArrangement atom) $-} updateAtom atom

atomsSimilar :: Atom -> Atom -> Bool
atomsSimilar a0 a1 = electronArrangement a0 == electronArrangement a1 &&
                     (all (all id) $ zipWith3 (zipWith3 similarEnergies) (energies' a0) (energies' a1) boolOccs)
    where similarEnergies e0 e1 occupied = (not occupied && e0<1e-8 && e1==0) || abs (e0-e1) <= min 0.001 (0.03*abs (e0+e1))
          boolOccs = map ((++ repeat False) . flip replicate True) (occupations a0) ++ repeat (repeat False)

genShieldingPotential :: Grid -> Orbital -> Potential
genShieldingPotential rs ψs = map (* coulumb'sConstant) $ zipWith (+) vOut vIn
    where d     = zipWith3 (\ψ r r' -> ψ^2 * r^3 * (r' - r)) ψs rs (tail rs)
          invSq = zipWith (\r x -> x/(r^2)) rs
          vOut  = (scanr (+) 0 $ invSq d) ++ repeat 0
          vIn   = invSq $ scanl (+) 0 d ++ repeat 1

getPotential :: Atom -> N -> L -> Maybe Spin -> Potential
getPotential atom n0 l0 s0 = zipWith (+) baseV (getShieldingPotential atom n0 l0 s0)
    where rs    = atomGrid atom
          z     = atomicNumber atom
          baseV = basePotential rs (fromIntegral l0) (fromIntegral z)

getShieldingPotential :: Atom -> N -> L -> Maybe Spin -> Potential
getShieldingPotential atom n0 l0 s0 = if length (orbitals atom !! l0) < n0 then repeat 0 else foldr (zipWith (+)) (repeat 0) scaledVs
    where arr       = electronArrangement atom
          occ0'     = fromIntegral <$> þrd3 <$> find (\(n,l,_) -> n == n0 && l == l0) arr
          occ0      = maybe 1 id occ0'
          errQ      = atomicNumber atom - charge atom - sum (map þrd3 arr)
          arr'      = if isJust occ0' || errQ /= 0 then arr else init arr ++ [let (n, l, o) = last arr in (n,l,o-1), (n0, l0, 1)] --For unoccupied orbitals, the calculated energy should be what the energy would be if they were occupied.
          zeroEnergy= (((energies atom !! l0) ++ repeat 0) !! (n0-1)) == 0
          upOcc0    = case s0 of
                          Just Up   -> occ0
                          Just Down -> 0
                          Nothing   -> min occ0 $ fromIntegral (l0 + 1)^2
          ψs0       = orbitals atom !! l0 !! (n0-1)
          rs        = atomGrid atom
          orbOverlap n l = if zeroEnergy then 0 else
              (/fromIntegral (l+l0+1)) $ sum $ zipWith4 (\r r' ψ0 ψ1 -> abs (r*r*r*(r'-r)*ψ0*ψ1)) rs (tail rs) ψs0 $ orbitals atom !! l !! (n-1)
          occ n l o = {-trace (show n0 ++ "," ++ show l0 ++ ", " ++ show occ0 ++ ", " ++ show n ++ "," ++ show l ++ ", " ++ show o) $ traceShow (orbOverlap n l) $ traceShowId $-} if n == n0 && l == l0 
              then upOcc0*(occ0-upOcc0)/occ0 + max 0 (occ0 - 1)/2
              else o - orbOverlap n l * (o - o/occ0 * upOcc0 + (2*upOcc0/occ0 - 1) * min o (fromIntegral (l+1)^2))
          scaledVs  = map (\(n,l,o) -> map (* occ n l (fromIntegral o)) (shieldingPotentials atom !! l !! (n-1))) arr
-- Ideally, the exchange energy should be the integral as x1 and x2 range over the 4D space of ψ1ψ2(x1) ψ1ψ2(x2) (x1-x2)^-2, while the plain repulsion term is  ψ1ψ1(x1) ψ2ψ2(x2) (x1-x2)^-2. I've very crudely approximated the ratio of these terms in orbOverlap, as I currently can't be bothered doing it better.


appendAt :: Int -> a -> [[a]] -> [[a]]
appendAt n x xss = xssh ++ (xs ++ [x]) : xsst
    where (xssh, (xs:xsst)) = splitAt n xss

prettyElectronArrangement :: Atom -> String
--prettySmoothedElectronArrangement :: Atom -> String
prettyElectronArrangement = prettify electronArrangement
    where prettify genEA atom = (unwords $ map (\(n, l, o) -> show (n+l) ++ (angularMomentumLabels !! l) : '^' : show o) $ ea) ++ chargeLabel
              where ea = genEA atom
                    errorCharge = atomicNumber atom - charge atom - sum (map þrd3 ea)
                    chargeLabel = if errorCharge == 0 then "" else "  " ++ show errorCharge ++ "!!"

atomFull :: Atom -> Bool
atomFull atom = atomicNumber atom - charge atom == sum (map þrd3 $ electronArrangement atom)

incrementAtom :: Atom -> Atom
incrementAtom atom = updateAtomUntilConvergence $ fillAtom $ copyAtom atom (atomicNumber atom + 1) (charge atom)

cationise :: Atom -> Atom
cationise atom = {-trace ("Ion: +" ++ show (charge atom + 1)) $-} updateAtomUntilConvergence atom{charge = charge atom + 1}
anionise  :: Atom -> Atom
--anionise  atom = traceShow (charge atom - 1) $ updateAtomUntilConvergence atom{charge = charge atom - 1}
anionise atom = makeAtom (atomicNumber atom) (charge atom - 1)

aperiodicTable :: [Atom]
--aperiodicTable = map (flip makeAtom 0) [1..]
aperiodicTable = iterate incrementAtom (makeAtom 1 0)
