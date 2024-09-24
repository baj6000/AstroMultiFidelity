# AstroMultiFidelity
Python software for multi-fidelity orbit uncertainty propagation.

This is the core software used for multi-fidelity uncertainty propagation presented in:

Jones and Weisman, "Multi-fidelity orbit uncertainty propagation", _Acta Astronautica_, Vol. 155, 2019, pp. 406-417, [https://doi.org/10.1016/j.actaastro.2018.10.023](https://doi.org/10.1016/j.actaastro.2018.10.023)

**NOTE:  You will need to use your own orbit propagation software.  An orbit propagtor is not provided in this package.**

This software only provides a minimal example based on a harmonic oscillator.  A multi-fidelity solution is not really needed for that test case, but it is used as a simply illustration of the software.  You will need to provide a wrapper for an orbit propagator of your choice to duplicate the results of the paper.  

## Dependencies:

This software assumes that you have packages that are standard in Conda (for example).  Specifically, you will need:

- NumPy
- MatPlotLib

## Acknolwedgement

If you use this software as part of a publication, we ask that you cite our paper listed above.  It is also prudent to cite the papers upon which we developed this tool:

A. Narayan, C. Gittelson, and D. Xiu, "A stochastic collocation algorithm with multi-fidelity models", _SIAM Journal of Scientific Computation_, Vol. 36, No. 22, 2014, pp. A495–A521, https://doi.org/10.1137/130929461.

X. Zhu, A. Narayan, and D. Xiu, "Computational aspects of stochastic collocation with multifidelity models", _SIAM/ASA Journal of Uncertainty Quantification_ Vol. 2, No. 1, 2014, pp. 444–463, https://doi.org/10.1137/130949154.
