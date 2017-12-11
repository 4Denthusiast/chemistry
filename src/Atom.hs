{-# LANGUAGE TupleSections #-}

module Atom(
    PerOrbital(..),
    Atom(..),
    AtomCache(..),
    getPO,
    trimPO,
    poToList,
    indices,
    emptyAtom,
    makeAtom,
    electronsRequired,
    incorrectCharge,
    updateAtom,
    getPotential,
    electronArrangement,
    prettyElectronArrangement,
    angularMomentumLabels,
    aperiodicTable,
    anionise,
    cationise
) where

import Orbital
import AtomEnergy
import Utils

import Data.List
import Control.Applicative
import qualified Data.Set as S
import Data.Function
import Data.Maybe
import Data.Bifunctor (first, second)
import Debug.Trace


angularMomentumLabels = "spdfghi" ++ (['k'..'z'] \\ "sp") ++ repeat 'R' -- "R" for "Rydberg".

data PerOrbital a = PerOrbital a [[a]]

instance Functor PerOrbital where
    fmap f (PerOrbital r xss) = PerOrbital (f r) (map (map f) xss)
instance Foldable PerOrbital where -- This instance is only appropriate where p == (PerOrbital z _).
    foldr f z p = foldr f z (concat $ poToList p)
instance Applicative PerOrbital where
    pure = emptyPO
    (PerOrbital r xss) <*> (PerOrbital r' xss') = PerOrbital (r r') (zipWithDefault (zipWithDefault id r r') [] [] xss xss')
        where zipWithDefault g d d' [] xs = map (     g d ) xs
              zipWithDefault g d d' xs [] = map (flip g d') xs
              zipWithDefault g d d' (x:xs) (x':xs') = g x x' : zipWithDefault g d d' xs xs'
liftPO2Short :: (a -> b -> c) -> PerOrbital a -> PerOrbital b -> PerOrbital c
liftPO2Short f (PerOrbital r xss) (PerOrbital r' xss') = PerOrbital (f r r') $ zipWith (zipWith f) xss xss'
instance (Show a) => Show (PerOrbital a) where
    show (PerOrbital _ xss) = (++"]") $ unlines $ ("[":) $ zipWith (\l xs -> l : ": " ++ show xs) angularMomentumLabels xss
instance (Eq a) => Eq (PerOrbital a) where
    p == p' = let (PerOrbital r x, PerOrbital r' x') = (trimPO p, trimPO p') in r == r' && x == x'

getPO :: PerOrbital a -> N -> L -> a
getPO (PerOrbital r xss) n l = case (drop l xss) of
    []     -> r
    (xs:_) -> case (drop (n-1) xs) of
        []    -> r
        (x:_) -> x

emptyPO :: a -> PerOrbital a
emptyPO x = PerOrbital x []

poToList :: PerOrbital a -> [[a]]
poToList (PerOrbital _ xss) = xss

trimPO :: Eq a => PerOrbital a -> PerOrbital a
trimPO (PerOrbital r xss) = PerOrbital r $ dropWhileEnd null $ map (dropWhileEnd (r ==)) xss

untrimPO :: PerOrbital a -> PerOrbital a
untrimPO (PerOrbital r xss) = PerOrbital r $ map (++[r]) xss ++ [[r]]

indices :: PerOrbital (Maybe (N, L))
indices = PerOrbital Nothing [[Just (n, l) | n <- [1..]] | l <- [0..]]



data Atom = Atom
    {
        atomicNumber :: Int,
        charge :: Int,
        atomGrid :: Grid,
        forcedOccs :: PerOrbital Int,
        energies :: PerOrbital Energy,
        prevEnergies :: PerOrbital Energy, -- For detecting rate of convergence.
        orbitals :: PerOrbital Orbital,
        occupations :: PerOrbital Double
    }

data AtomCache = AtomCache
    {
        atomInCache :: Atom,
        totalPotential :: [Double], --doesn't include the nuclear potential.
        shieldingPotentials :: PerOrbital [Double] --used to exclude each electron's self-shielding. TODO: try just using the exchange terms for this: it is equivalent after convergence (not before), and doing so would simplify some things considerably.
    }

electronsRequired :: Atom -> Int
electronsRequired atom = atomicNumber atom - charge atom

incorrectCharge :: Atom -> Double
incorrectCharge atom = fromIntegral (electronsRequired atom) - sum (occupations atom)

emptyAtom :: Int -> Int -> Atom
emptyAtom z c = Atom z c (logGrid 0.01 (fromIntegral z)) (emptyPO 0) (emptyPO 0) (emptyPO 0) (emptyPO []) (emptyPO 0)

-- If the atomic numbers are similar, this should converge faster than starting from empty.
copyAtom :: Atom -> Int -> Int -> Atom
copyAtom atom z c = atom {
        atomicNumber = z,
        charge = c,
        atomGrid = logGrid 0.01 (fromIntegral z)
    }

copyAtom' :: Atom -> Int -> Int -> Atom
copyAtom' atom z' c' = copyAtom atom (z' + atomicNumber atom) (c' + charge atom)

incrementAtom :: Atom -> Atom
incrementAtom atom = relaxAtom $ copyAtom' atom 1 0

cationise :: Atom -> Atom
cationise atom = relaxAtom $ copyAtom' atom 0 1
anionise  :: Atom -> Atom
anionise  atom = relaxAtom $ copyAtom' atom 0 (-1)

aperiodicTable :: [Atom]
aperiodicTable = iterate incrementAtom (makeAtom 1 0)



lOccupations :: AtomCache -> [Int]
lOccupations atom = map (length . dropWhileEnd (==0)) $ poToList $ occupations $ atomInCache atom

electronArrangement :: Atom -> [(N, L, Int)]
electronArrangement atom = occsToEA $ round <$> occupations atom
doubleElectronArrangement :: Atom -> [(N, L, Double)]
doubleElectronArrangement atom = occsToEA $ occupations atom

occsToEA :: Eq a => PerOrbital a -> [(N, L, a)]
occsToEA = concat . poToList . liftPO2Short (uncurry (,,) . fromJust) indices . trimPO
eaToOccs :: [(N, L, Int)] -> PerOrbital Int
eaToOccs = PerOrbital 0 . map (map (sum . map þrd3) . compileOn ((+(-1)) . fst3)) . compileOn snd3
    where compileOn :: (a -> Int) -> [a] -> [[a]]
          compileOn = compileOn' 0
          compileOn' n f [] = []
          compileOn' n f xs = uncurry (:) $ second (compileOn' (n+1) f) $ partition ((n==) . f) xs

forceEA :: Atom -> [(N, L, Int)] -> Atom
forceEA atom ea = atom{forcedOccs = eaToOccs ea}

forcedEA atom = occsToEA $ forcedOccs atom

correctEA :: Atom -> Bool
correctEA atom = occupations atom == (fromIntegral <$> forcedOccs atom)

prettyElectronArrangement :: Atom -> String
prettyElectronArrangement atom = prettifyElectronArrangement atom (electronArrangement atom)

prettifyElectronArrangement :: Atom -> [(N, L, Int)] -> String
prettifyElectronArrangement atom ea = (unwords $ map (\(n, l, o) -> show (n+l) ++ (angularMomentumLabels !! l) : '^' : show o) $ ea) ++ chargeLabel
    where errorCharge = electronsRequired atom - sum (map þrd3 ea)
          chargeLabel = if errorCharge == 0 then "" else "  " ++ show errorCharge ++ "!!"



genShieldingPotential :: Grid -> Orbital -> [Double]
genShieldingPotential rs ψs = map (* coulumb'sConstant) $ zipWith (+) vOut vIn
    where d     = zipWith3 (\ψ r r' -> ψ^2 * r^3 * (r' - r)) ψs rs (tail rs)
          invSq = zipWith (\r x -> x/(r^2)) rs
          vOut  = (scanr (+) 0 $ invSq d) ++ repeat 0
          vIn   = invSq $ scanl (+) 0 d ++ repeat 1

--Calculate further details of the same atom.
genAtomCache :: Atom -> AtomCache
genAtomCache atom = AtomCache
        atom
        (foldr (zipWith (+)) (repeat 0) ((map . (*)) <$> (trimPO $ occupations atom) <*> vs))
        vs
    where vs = genShieldingPotential (atomGrid atom) <$> orbitals atom

getShieldingPotential :: AtomCache -> N -> L -> Maybe Spin -> Potential
getShieldingPotential atomCache n0 l0 s0 = (vs, xs)
    where atom      = atomInCache atomCache
          arr       = doubleElectronArrangement atom
          occ0      = getPO (occupations atom) n0 l0
          zeroEnergy= getPO (energies atom) n0 l0 == 0
          upOcc0    = case s0 of
                          Just Up   -> occ0
                          Just Down -> 0
                          Nothing   -> min occ0 $ fromIntegral (l0 + 1)^2
          ψs0       = getPO (orbitals atom) n0 l0
          rs        = atomGrid atom
          xch n l o = let
                  ψs = getPO (orbitals atom) n l
                  orbOverlap = zipWith4 (\r r' ψ0 ψ -> (r'-r)*ψ*ψ0*r^3) rs (tail rs) ψs0 ψs
                  dl = abs (l-l0) -- The lowest term with non-zero coefficient
                  bias i = zipWith (*) (map (^^i) rs)
                  inner = bias   dl    $ scanr (+) 0 $ bias (-dl-2) orbOverlap
                  outer = bias (-dl-2) $ scanl (+) 0 $ bias   dl    orbOverlap
                  o' = xchOcc n l o
              in zipWith3 (\ψ inn out -> o'*ψ*(inn+out)/fromIntegral ((l+1)*(l0+1))) ψs inner outer ++ repeat 0
          xchOcc n l o
              | occ0 == 0          = min o (fromIntegral (l+1)^2)
              | n == n0 && l == l0 = 2*upOcc0^2 / occ0 + occ0 - 2*upOcc0 - 1
              | otherwise          = ((occ0 - upOcc0)*o + (2*upOcc0 - occ0) * min o (fromIntegral (l+1)^2))/occ0
          vs        = (if zeroEnergy then id else zipWith (flip (-)) (getPO (shieldingPotentials atomCache) n0 l0)) $ totalPotential atomCache
          xs        = foldr (zipWith (+)) (repeat 0) $ if zeroEnergy then [] else map (uncurry3 xch) arr

getPotential :: AtomCache -> N -> L -> Maybe Spin -> Potential
getPotential atomC n0 l0 s0 = first (zipWith (+) $ fst baseV) (getShieldingPotential atomC n0 l0 s0)
    where atom  = atomInCache atomC
          rs    = atomGrid atom
          z     = atomicNumber atom
          baseV = basePotential rs (fromIntegral l0) (fromIntegral z)



genOrbitalAt :: AtomCache -> N -> L -> Energy -> (Energy, Orbital)
genOrbitalAt atom n l e0 = (e, o)
    where rs = atomGrid $ atomInCache $ atom
          vs = getPotential atom n l Nothing
          e  = findEnergyHinted e0 rs vs n
          o  = normalizedOrbital rs vs e

updateAtom :: AtomCache -> Double -> Bool -> Atom
updateAtom prevAC smooth expand = let
        atom = atomInCache prevAC
        es = liftPO2Short const (untrimPO $ energies atom) ((if expand then untrimPO else id) $ trimPO $ forcedOccs atom)
        unzipPO p = (fst <$> p, snd <$> p)
        (es', ψs') = unzipPO $ liftPO2Short (maybe (const (0,[])) $ uncurry $ genOrbitalAt prevAC) indices es
        smoothOcc o' o = if abs (o'-o)*(1-smooth) < 0.1 then o else smooth*o + (1-smooth)*o'
        os' = trimPO $ liftA3 (\e o' o -> if e==0 then 0 else smoothOcc o' (fromIntegral o)) es' (occupations atom) (forcedOccs atom)
    in atom
        {
            energies = es',
            prevEnergies = energies atom,
            orbitals = ψs',
            occupations = os'
        }
-- Given that f(xi)=yi, estimate an x such that (lerp x y1 y2) is a fixed point of f, assuming it's linear.
linearFix x1 y1 x2 y2 = max 0.05 $ min 1 $ (x1-y1)/(x1-y1-x2+y2)

mixAtoms :: AtomCache -> AtomCache -> AtomCache
mixAtoms ac0 ac1 = let
        [atom0, atom1] = map atomInCache [ac0, ac1]
        xs = linearFix <$> prevEnergies atom0 <*> energies atom0 <*> prevEnergies atom1 <*> energies atom1
        lerp x a b = case x of {0 -> a; 1 -> b; _ -> b*x + a*(1-x)}
        lerpxs f = lerp <$> xs <*> f atom0 <*> f atom1
        lerpsxs f = (zipWithLong 0 0 . lerp) <$> xs <*> f ac0 <*> f ac1
        occs = lerpxs occupations
        shields = lerpsxs shieldingPotentials
    in --trace ("mixing: "++show xs) $
    AtomCache
        atom1{
            energies = lerpxs energies,
            prevEnergies = lerpxs prevEnergies,
            occupations = occs,
            orbitals = lerpsxs (orbitals . atomInCache) -- This isn't normalised, but that shouldn't be a problem.
        }
        (foldr (zipWith (+)) (repeat 0) $ (map . (*)) <$> occs <*> shields)
        shields

startAtom :: AtomCache -- an atom such that (mixAtoms startAtom ~= id)
startAtom = AtomCache
    (Atom undefined undefined undefined undefined (emptyPO 0) (emptyPO (1e20)) (emptyPO []) (emptyPO 0))
    [] (emptyPO [])

traceAtom :: Atom -> a -> a
traceAtom atom = traceShow (energies atom) . traceShow (occupations atom) . {-traceShow (carefulEnergies atom) .-} trace (prettifyElectronArrangement atom $ occsToEA $ forcedOccs atom) . trace (prettyElectronArrangement atom ++ "\n")

atomsSimilar :: Atom -> Atom -> Bool
atomsSimilar a0 a1 = occupations a0 == occupations a1 &&
                     (all id (similarEnergies <$> (energies a0) <*> (energies a1) <*> boolOccs))
    where similarEnergies e0 e1 occupied = (not occupied && e0<1e-8 && e1==0) || abs (e0-e1) <= min 0.001 (0.01*abs (e0+e1))
          boolOccs = fmap (>0) (occupations a0)

updateAtomUntilConvergence :: Atom -> Atom
updateAtomUntilConvergence = (\a -> traceAtom a a) . uauc' 0 startAtom . removePrevEnergies
    where removePrevEnergies atom = atom{prevEnergies = emptyPO 0}
          uauc' n prev atom
              | atomsSimilar atom next           = next
              | on (>) incorrectCharge next atom = if n >= 3 then next else uauc' (n+1) ac next
              | otherwise                        = uauc' n ac next
                  where next = traceAtom atom $ updateAtom (mixAtoms prev ac) (1/ fromIntegral (n+1)) False
                        ac   = genAtomCache atom



adjacentElectronArrangements :: [(N, L, Int)] -> Int -> [[(N, L, Int)]]
adjacentElectronArrangements ea eReq = traceShowId $ map (sort . filter ((>0) . þrd3)) $ (if errCharge > 0 then veryNegativeOthers ++ negativeOthers else []) ++ others ++ positiveOthers
    where eaByL = groupBy (on (==) snd3) $ sortOn snd3 ea
          maxL  = snd3 $ head $ last $ [(1,-1,0)] : eaByL
          nextL [] _ = ([], [])
          nextL (xs@((_,l0,_):_):xss) l = if l0 == l then (xss, xs) else (xs:xss, [])
          eaByL' = snd $ mapAccumL nextL eaByL [0..maxL+1]
          oByL   = map (sum . map þrd3) eaByL'
          errCharge = eReq - sum oByL
          newOrbs = flip (zipWith (,,0)) [0..] $ map ((+1) . fst3 . last . ((0,0,0):)) eaByL'
          additionFrontier = unionBy (on (==) snd3) (filter (\(n,l,o) -> o < 2*(l+1)^2) ea) newOrbs
          removalFrontier  = mergeBy (\x0@(_,l0,_) x1@(_,l1,_) -> if l0 == l1 then Just (max x0 x1) else Nothing) ea
          modifiedEAs (na,la,oa) (nr,lr,or) = valid >> (flip (++) ea' <$> (newXs ++ discarding))
              where steps l = (*(l+1)^2) <$> [0,1,2]
                    oa' = filter (>0) $ map (+(-oa)) $ steps la
                    or' = filter (>0) $ map (or-   ) $ steps lr
                    o'  = filter (<= min (maximum oa') (maximum or')) $ union oa' or'
                    ea' = ea \\ [(na,la,oa),(nr,lr,or)]
                    newXs = (\o -> [(na,la,oa+o), (nr,lr,or-o)]) <$> o'
                    lla = 2*(la+1)^2
                    discarding = if lla < oa + or then [[(na,la,lla)]] else []
                    valid = if la >= lr && na >= nr then [] else [()]
          others = concat $ modifiedEAs <$> additionFrontier <*> removalFrontier
          chargedEA o' (na,la,oa) = (na,la,oa+o') : (delete (na,la,oa) ea)
          positiveOthers = map (chargedEA (-1))  removalFrontier
          negativeOthers = map (chargedEA   1 ) additionFrontier
          veryNegativeOthers = map (\(n,l,o,o') -> chargedEA o' (n,l,o)) $ filter ((>1) . frþ4) $ (\(n,l,o) -> (n,l,o,min errCharge $ 2*(l+1)^2-o)) <$> additionFrontier
          
atomAdjEAs :: Atom -> [[(N, L, Int)]]
atomAdjEAs atom = adjacentElectronArrangements (electronArrangement atom) (electronsRequired atom)

relaxAtom :: Atom -> Atom
relaxAtom atom = relax' (updateAtomUntilConvergence atom) S.empty
    where relax' bestA triedEAs = case as of {[] -> bestA; (bestA':_) -> relax' bestA' triedEAs'}
              where nextEAs = filter (flip S.notMember triedEAs) $ atomAdjEAs bestA
                    nextAtomInits = sortOn (\a -> (abs $ forcedIncorrectCharge a, totalEnergy a)) $ map (forceEA $ updateAtom (genAtomCache bestA) 0.3 True) nextEAs
                    nextAtoms = map updateAtomUntilConvergence nextAtomInits
                    (failedAtoms, as) = span (not . better) nextAtoms
                    better a = correctEA a && on (<) (\a' -> (max 0 (-incorrectCharge a'), totalEnergy a')) a bestA
                    triedEAs' = S.union triedEAs $ S.fromList $ map forcedEA $ bestA : failedAtoms
                    forcedIncorrectCharge a = electronsRequired a - sum (forcedOccs a)

makeAtom :: Int -> Int -> Atom
makeAtom z c = relaxAtom $ emptyAtom z c
