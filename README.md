# simple-2d-ns

Hello World to incompressible Navier--Stokes. The basic features of this solver are:
- Two-dimensional flow
- Second-order finite difference, staggered grid
- Fractional-step method
- Pressure-Poisson equation solved via FFTs
- OpenMP acceleration
- Fully explicit Euler time-integration
- Uniform grid spacing

Various branches have been created in what is (generally) increasing order of implementation complexity:
- `main`: Doubly-periodic domain, Pressure-Poisson equation solved via 2D FFTs
- `walls`: Implementation of wall-boundary conditions in one direction. The Poisson equation is solved with 1D FFTs in the periodic x-direction and a tridiagonal inversion in the remaining wall-normal direction. A scalar field is also added so that 2D Rayleigh-Benard flow can be simulated.
- `implicit`: Implicit time-stepping treatment of the diffusive terms using an Alternating Direct Implicit (ADI) approach or Helmholtz solver. Currently a WIP.
- `RK3`: Implementation of RK3 sub-stepping. TODO.
- `non-uniform-grid`: Grid-stretching in the wall-normal direction. TODO
## Dependencies

- FFTW3
- Fortran compiler