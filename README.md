# Fundamentals of DLVO theory

Copyright (c) 2022 by Miguel A. Caro

These are materials for the teaching demo at Aalto PHYS about DLVO theory
of colloidal aggregation. Check a very brief introduction to the theory
under `slides/`. 

There is a code to perform simple molecular dynamics with a system made
of colloidal particles and solvent particles (positively and negatively
charged). The model is quite simple and based on a screened electrostatic
interaction plus Lennard-Jones potential in periodic boundary conditions
under the minimum image convention. To build the code just execute
`./build.sh` and add the `src/` directory to your path. You can then
use the `aggregation` binary which requires a local `input` file (an
example file can be found under `example/`). The
code is rather primitive and made for educational purposes. If you want
to use it for research you should probably spend some time figuring out
whether the different approximations are sound enough for your purposes.
