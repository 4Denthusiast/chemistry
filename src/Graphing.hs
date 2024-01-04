module Graphing(
    graphAtom,
    graphValues,
) where

import Orbital
import Dual
import Atom
import OrbitalShape

import Control.Monad
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Diagrams(toFile)
import Graphics.Rendering.Chart.Axis.Floating

graphAtom :: Atom -> IO ()
graphAtom atom = toFile def ("charts/Z" ++ show (atomicNumber atom) ++ "A" ++ show (massNumber atom) ++ ".svg") $ do
    layout_title .= "Element #" ++ show (atomicNumber atom) ++"("++ show (massNumber atom) ++") charge: "++ show (charge atom) ++",  "++prettyElectronArrangement atom
    layout_x_axis . laxis_generate .= autoScaledLogAxis (LogAxisParams (map show))
    mapM_
        (\(l, ψs) -> plot (line [l] $ map (denormalize . zip (atomGrid atom) . take 1500) ψs))
        (zip angularMomentumLabels $ poToList $ trimPO $ (\occ orb -> if occ>0 then orb else []) <$> (occupations atom) <*> (orbitals atom))
        --(zip angularMomentumLabels $ takeWhile (not . null) (orbitals atom))
        --(zip angularMomentumLabels $ takeWhile (not . null) $ take 1 $ (orbitals atom))

denormalize :: [(Double, Double)] -> [(Double, Double)]
denormalize or@((r0,ψ0):(r1,ψ1):_) = map (\(r, ψ) -> (r, ψ*s*r^1)) or
    where s = (r1-r0)/(ψ1-ψ0)

graphValues :: Double -> Grid -> [[Double]] -> IO ()
graphValues lim rs xs = toFile def "chart.svg" $ do
        layout_title .= "Thing"
        layout_x_axis . laxis_generate .= autoScaledLogAxis (LogAxisParams (map show))
        forM_ xs (plot . line "Values" . (:[]) . zip rs . map (min lim . max (-lim)))
