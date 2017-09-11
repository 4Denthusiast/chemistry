# Ab initio 4D quantum chemistry simulation

A highly messy program for calculating properties of 4D atoms with an additional short-range repulsive force, using an approximation to the Hartree-Fock method.

## Building & Running

1. clone repository
2. Install stack: https://docs.haskellstack.org/en/stable/README/
3. Initialise stack
  $ stack setup
4. build
  $ stack build
5. run
  * To print properties of a specific atom (with charge in the format "n" or "n-"):
    $ stack run chemistry graph <atomic number> [<charge>]
  * To calculate the electron arrangements of all elements up to a specified atomic number
    $ stack run chemistry ea <atomic number>

## Features

Calculates electron arrangements, electron affinities/ionisation energies, asymptotic bond strengths (per orbital, as internuclear distance tends to infinity) and orbital shapes.

## Physics model

It's based on a realistic but very approximated and non-relativistic quantum model. Electrons are bound to the nucleus with an r^-2 potential, with a short range repulsive force of the form r^-2 * e^-r. The following physical constants are normalised to 1:

  * electron mass
  * electron charge
  * reduced Planck's constant
  * scale-factor for nuclear repulsive force
