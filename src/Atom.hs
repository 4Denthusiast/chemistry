{-# LANGUAGE TupleSections #-}

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
import qualified Data.Set as S
import Data.Function
import Data.Maybe
import Debug.Trace


angularMomentumLabels = "spdfghi" ++ ['k'..'z'] ++ repeat 'R' -- "R" for "Rydberg".

data Atom = Atom
    {
        atomicNumber :: Int,
        charge :: Int,
        energies :: [[Energy]],
        orbitals :: [[Orbital]],
        shieldingPotentials :: [[Potential]],
        atomGrid :: Grid,
        prevOccs :: [[Double]],
        forcedEA :: Maybe [(N, L, Int)] --a value of Nothing represents that the electron arrangement should be determined from energies.
    }


energies' :: Atom -> [[Energy]]
energies' = takeWhile (not . null) . energies

energies'' :: Atom -> [[Energy]]
energies'' atom = zipWith take (map (+ 1) occ) (energies atom)
    where occ = occupations atom ++ [0]

prevOccs' :: Atom -> [[Double]]
prevOccs' = (++ repeat (repeat 0)) . map (++ repeat 0) . prevOccs

electronsRequired :: Atom -> Int
electronsRequired atom = atomicNumber atom - charge atom

incorrectCharge :: Atom -> Int
incorrectCharge atom = electronsRequired atom - sum (map þrd3 $ electronArrangement atom)

emptyAtom :: Int -> Int -> Atom
emptyAtom z c = Atom z c empty empty empty (logGrid 0.01 (fromIntegral z)) [] Nothing
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
makeAtom z c = relaxAtom $ fillAtom $ emptyAtom z c

fillAtom :: Atom -> Atom
fillAtom atom = 
    let
        occ = occupations atom
        done = (all (\(es, n) -> not $ null $ drop n es) $ zip (energies atom) occ)
            && (length occ < length (energies' atom))
    in
        if done then atom else fillAtom $ addOrbital occ atom

electronArrangement :: Atom -> [(N, L, Int)]
electronArrangement atom = maybe ea trimEA $ forcedEA atom
    where spareElectrons = electronsRequired atom
          ea = mergeSpins $
               takeWhile ((>0) . þrd3) $
               snd $ mapAccumL (\a (n,l,ll) -> (a-min a ll, (n, l, min a ll))) spareElectrons $
               map (\(e, n, l) -> (n, l, (l+1)^2)) $
               sort $ filter ((<0) . fst3) $
               (\s (n, l) -> (approxOrbitalEnergy atom (n, l, s), n, l)) <$> [Up, Down] <*>
               (concat $ zipWith (zipWith (const id)) (energies' atom) indices)
          mergeSpins = mergeBy (\(n0, l0, o0) (n1, l1, o1) -> if n0 == n1 && l0 == l1 then Just (n0,l0,o0+o1) else Nothing)
          trimEA = mergeSpins . map (\(n,l,s,o) -> (n,l,o)) . filter ((<0) . approxOrbitalEnergy atom . (\(n,l,s,o) -> (n,l,s))) . concatMap splitSpin

prevElectronArrangement :: Atom -> [(N, L, Int)]
prevElectronArrangement = concat . zipWith (zipWith (\(n, l) o -> (n, l, floor o))) indices . prevOccs

indices :: [[(N, L)]]
indices = [[(n, l) | n <- [1..]] | l <- [0..]]


occupations' :: [(N,L,Int)] -> [Int]
occupations' = map length . group . sort . map snd3
occupations :: Atom -> [Int]
occupations  = occupations' . electronArrangement

addOrbital :: [Int] -> Atom -> Atom
addOrbital occs atom = addOrbitalAt atom (getNextOrbital occs $ energies' atom)

getNextOrbital :: [Int] -> [[Energy]] -> L
getNextOrbital occs [] = 0
getNextOrbital occs es
    | length (last es) > 1 = length es
    | otherwise            = snd3 $ minimum $ filter (\(e, l, occ) -> null $ drop occ e) $ zip3 es [0..] (occs++repeat 0)

addOrbitalAt :: Atom -> L -> Atom
addOrbitalAt atom l = atom
        {
            energies = appendAt l e (energies atom),
            orbitals = appendAt l o (orbitals atom),
            shieldingPotentials = appendAt l p (shieldingPotentials atom)
        }
    where n  = 1 + length (energies atom !! l)
          (e, o, p) = genOrbitalAt atom 0.5 n l (-1)

genOrbitalAt :: Atom -> Double -> N -> L -> Energy -> (Energy, Orbital, Potential)
genOrbitalAt atom smooth n l e0 = (e, o, p)
    where rs = atomGrid atom
          vs = getPotential atom n l Nothing
          e  = findEnergyHinted e0 rs vs n
          o  = normalizedOrbital rs vs e
          pn = genShieldingPotential rs o
          po = ((map Just $ shieldingPotentials atom !! l) ++ repeat Nothing) !! (n-1)
          p  = (zipWith (\a b -> a + smooth*(b-a)) pn) $ maybe (repeat 0) id po

updateAtom :: Bool -> Double -> Atom -> Atom
updateAtom small smooth atom = (if small then id else fillAtom) $ let
        es = maybe id (zipWith take . occupations') (if small then forcedEA atom else Nothing) (energies'' atom)
        (es', ψs', vs') = unzip3 $ map unzip3 $ zipWith (zipWith $ uncurry $ genOrbitalAt atom smooth) indices es
        lerp x a b = if abs (a-b) < 0.05 then b else b*x + a*(1-x)
    in atom
        {
            energies = es' ++ repeat [],
            orbitals = ψs' ++ repeat [],
            shieldingPotentials = vs' ++ repeat [],
            prevOccs = zipWith (zipWith (lerp smooth)) (prevOccs' atom) $ map (map (fromIntegral . þrd3) . sort) $ groupBy (on (==) snd3) $ sortBy (on compare snd3) $ electronArrangement atom
        }

updateAtomUntilConvergence :: Atom -> Atom
updateAtomUntilConvergence = uauc' 0
    where uauc' n atom
              | atomsSimilar atom next           = next
              | on (>) incorrectCharge next atom = if n >= 3 then next else uauc' (n+1) next
              | otherwise                        = uauc' n next
                  where next = traceShow (energies' atom) $ traceShow (prevOccs atom) $ {-traceShow (carefulEnergies atom) $-} traceShow (prettifyElectronArrangement atom <$> forcedEA atom) $ trace (prettyElectronArrangement atom) $ updateAtom True (1/ fromIntegral (n+1)) atom

atomsSimilar :: Atom -> Atom -> Bool
atomsSimilar a0 a1 = prevOccs a0 == prevOccs a1 &&
                     (all (all id) $ zipWith3 (zipWith3 similarEnergies) (energies' a0) (energies' a1) boolOccs)
    where similarEnergies e0 e1 occupied = (not occupied && e0<1e-8 && e1==0) || abs (e0-e1) <= min 0.001 (0.01*abs (e0+e1))
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
    where arr       = concat $ zipWith (zipWith (uncurry (,,))) indices $ prevOccs atom
          occ0'     = þrd3 <$> find (\(n,l,_) -> n == n0 && l == l0) arr
          occ0      = maybe 1 id occ0'
          errQ      = fromIntegral (electronsRequired atom) - sum (map þrd3 arr)
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
          occ n l o = if n == n0 && l == l0 
              then upOcc0*(occ0-upOcc0)/occ0 + max 0 (occ0 - 1)/2
              else o - orbOverlap n l * (o - o/occ0 * upOcc0 + (2*upOcc0/occ0 - 1) * min o (fromIntegral (l+1)^2))
          scaledVs  = map (\(n,l,o) -> map (* occ n l o) (shieldingPotentials atom !! l !! (n-1))) arr
-- Ideally, the exchange energy should be the integral as x1 and x2 range over the 4D space of ψ1ψ2(x1) ψ1ψ2(x2) (x1-x2)^-2, while the plain repulsion term is  ψ1ψ1(x1) ψ2ψ2(x2) (x1-x2)^-2. I've very crudely approximated the ratio of these terms in orbOverlap, as I currently can't be bothered doing it better.


appendAt :: Int -> a -> [[a]] -> [[a]]
appendAt n x xss = xssh ++ (xs ++ [x]) : xsst
    where (xssh, (xs:xsst)) = splitAt n xss

prettyElectronArrangement :: Atom -> String
prettyElectronArrangement atom = prettifyElectronArrangement atom (electronArrangement atom)

prettifyElectronArrangement :: Atom -> [(N, L, Int)] -> String
prettifyElectronArrangement atom ea = (unwords $ map (\(n, l, o) -> show (n+l) ++ (angularMomentumLabels !! l) : '^' : show o) $ ea) ++ chargeLabel
    where errorCharge = electronsRequired atom - sum (map þrd3 ea)
          chargeLabel = if errorCharge == 0 then "" else "  " ++ show errorCharge ++ "!!"

atomFull :: Atom -> Bool
atomFull atom = incorrectCharge atom == 0

incrementAtom :: Atom -> Atom
incrementAtom atom = relaxAtom $ copyAtom atom (atomicNumber atom + 1) (charge atom)

cationise :: Atom -> Atom
cationise atom = relaxAtom atom{charge = charge atom + 1}
anionise  :: Atom -> Atom
anionise  atom = relaxAtom atom{charge = charge atom - 1}

aperiodicTable :: [Atom]
aperiodicTable = iterate incrementAtom (makeAtom 1 0)


atomAdjEAs :: Atom -> [[(N, L, Int)]]
atomAdjEAs atom = adjacentElectronArrangements (electronArrangement atom) (electronsRequired atom)

adjacentElectronArrangements :: [(N, L, Int)] -> Int -> [[(N, L, Int)]]
adjacentElectronArrangements ea eReq = traceShowId $ map (sort . filter ((>0) . þrd3)) $ case signum (eReq - sum oByL) of
        -1 -> positiveOthers ++ others
        0  -> others
        1  -> negativeOthers ++ others
    where eaByL = groupBy (on (==) snd3) $ sortOn snd3 ea
          maxL  = snd3 $ head $ last $ [(1,-1,0)] : eaByL
          nextL [] _ = ([], [])
          nextL (xs@((_,l0,_):_):xss) l = if l0 == l then (xss, xs) else (xs:xss, [])
          eaByL' = snd $ mapAccumL nextL eaByL [0..maxL+1]
          oByL   = map (sum . map þrd3) eaByL'
          newOrbs = flip (zipWith (,,0)) [0..] $ map ((+1) . fst3 . last . ((0,0,0):)) eaByL'
          additionFrontier = unionBy (on (==) snd3) (filter (\(n,l,o) -> o < 2*(l+1)^2) ea) newOrbs
          removalFrontier  = mergeBy (\x0@(_,l0,_) x1@(_,l1,_) -> if l0 == l1 then Just (max x0 x1) else Nothing) ea
          modifiedEAs (na,la,oa) (nr,lr,or) = valid >> (flip (++) ea' <$> newXs)
              where steps l = (*(l+1)^2) <$> [0,1,2]
                    oa' = filter (>0) $ map (+(-oa)) $ steps la
                    or' = filter (>0) $ map (or-   ) $ steps lr
                    o'  = filter (<= min (maximum oa') (maximum or')) $ union oa' or'
                    ea' = ea \\ [(na,la,oa),(nr,lr,or)]
                    newXs = (\o -> [(na,la,oa+o), (nr,lr,or-o)]) <$> o'
                    valid = if la >= lr && na >= nr then [] else [()]
          others = concat $ modifiedEAs <$> additionFrontier <*> removalFrontier
          chargedEA o' (na,la,oa) = (na,la,oa+o') : (delete (na,la,oa) ea)
          positiveOthers = map (chargedEA (-1))  removalFrontier
          negativeOthers = map (chargedEA   1 ) additionFrontier

forceEA :: Atom -> [(N, L, Int)] -> Atom
forceEA atom ea = atom{forcedEA = Just ea}

correctEA :: Atom -> Bool
correctEA atom = case forcedEA atom of
    Just ea -> electronArrangement atom == ea
    Nothing -> True

relaxAtom :: Atom -> Atom
relaxAtom atom = relax' atom' S.empty
    where atom' = updateAtomUntilConvergence $ forceEA atom (electronArrangement atom)
          relax' bestA triedEAs = case as of {[] -> bestA; (bestA':_) -> relax' bestA' triedEAs'}
              where nextEAs = filter (flip S.notMember triedEAs) $ atomAdjEAs bestA
                    nextAtomInits = sortOn totalEnergy $ map (forceEA $ updateAtom False 0.5 bestA) nextEAs
                    nextAtoms = map updateAtomUntilConvergence nextAtomInits
                    (failedAtoms, as) = span (not . better) nextAtoms
                    better a = correctEA a && on (<) (\a' -> (abs (incorrectCharge a'), totalEnergy a')) a bestA
                    triedEAs' = S.union triedEAs $ S.fromList $ map (fromJust . forcedEA) $ bestA : failedAtoms
