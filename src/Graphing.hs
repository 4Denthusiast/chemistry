module Graphing(
    graphAtom,
    drawGraph,
    graphValues,
) where

import Orbital
import Dual
import Atom
import OrbitalShape

import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Diagrams(toFile)
import Graphics.Rendering.Chart.Axis.Floating

graphAtom :: Atom -> IO ()
graphAtom atom = toFile def ("charts/" ++ show (atomicNumber atom) ++ ".svg") $ do
    layout_title .= "Element #" ++ show (atomicNumber atom) ++", charge: "++ show (charge atom) ++",  "++prettyElectronArrangement atom
    layout_x_axis . laxis_generate .= autoScaledLogAxis (LogAxisParams (map show))
    mapM_
        (\(l, ψs) -> plot (line [l] $ map (denormalize . zip (atomGrid atom) . take 800) ψs))
        (zip angularMomentumLabels $ zipWith take (occupations atom) (orbitals atom))
        --(zip angularMomentumLabels $ takeWhile (not . null) (orbitals atom))
        --(zip angularMomentumLabels $ takeWhile (not . null) $ take 1 $ (orbitals atom))

rs = logGrid 0.005 1
rst = (takeWhile (<10000) rs)
vs = basePotential rs 0 1
ψs = trimmedOrbital rs vs 0
ψs' = map (std . fst) ψs

denormalize :: [(Double, Double)] -> [(Double, Double)]
denormalize or@((r0,ψ0):(r1,ψ1):_) = map (\(r, ψ) -> (r, ψ*s*r^1)) or
    where s = (r1-r0)/(ψ1-ψ0)

drawGraph :: IO ()
drawGraph = putStrLn (show (principalNumber rs ψs)) >>
    (toFile def "chart.svg" $ do
        layout_title .= "Orbital"
        layout_x_axis . laxis_generate .= autoScaledLogAxis (LogAxisParams (map show))
        plot (line "H ∞s" [zipWith (\r ψ -> (r, r*ψ)) rst ψs'])
        plot (line "Decay" [zipWith (\r (ψ, dψ) -> (r, r^2 * (3*(std ψ)^2 + (std ψ + std dψ * r)^2))) rst ψs])
    )

graphValues :: Grid -> [Double] -> IO ()
graphValues rs xs = toFile def "chart.svg" $ do
        layout_title .= "Thing"
        layout_x_axis . laxis_generate .= autoScaledLogAxis (LogAxisParams (map show))
        plot (line "Values" [zip rs xs])
