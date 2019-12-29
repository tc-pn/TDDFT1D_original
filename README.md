# TDDFT1D_original
Here is conserved the original version of a Fortran TDDFT code of a 1D nuclear system. It was made to test its stability using a t_0, t_3 functional of the one-body density. I am not its original creator.

The EOS files are here to compute the equation of state of the fermionic system under the 1D infinite matter approximation.
The rest of the ils are the DFT code itself.

COMMON initializes all needed variables.

init_DFT1DMF.for computes the initial single-particle states in a iterative way, in a Wood-Saxon potential completed by a harmonic oscillator and an optical potential at the borders of the spatial lattice.

TDDFT1D.for contains the code for the time-propagation.

The other files are either modules containig necessary rountines to perform the calculation (routines_gan.for) or other specific initializations of variables.
