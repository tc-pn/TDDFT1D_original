# TDDFT1D_original
Here is conserved the original version of a Fortran TDDFT code of a 1D nuclear system. It was made to test its stability using a t_0, t_3 functional of the one-body density. I am not its original creator.

The EOS files are here to compute the equation of state of the fermionic system under the 1D infinite matter approximation.
The rest of the ils are the DFT code itself.

COMMON initializes all needed variables. For now, the mesh is composed of 250 points and max 40 single particle waves (around 10 000 values stored values at all times in the one body density matrix).

init_DFT1DMF.for computes the initial single-particle states in a iterative way, in a Wood-Saxon potential completed by a harmonic oscillator and an optical potential at the borders of the spatial lattice. The harmonic oscillator is there to prevent an initial state with particles in the continuum. There is a possibility to initialize the code using only an Harmonic Oscillator.

TDDFT1D.for contains the code for the time-propagation.

The other files are either modules containig necessary rountines to perform the calculation (routines_gan.for) or other specific initializations of variables.
