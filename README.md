# Poisson 2D
Solvers for Poisson equation in 1D and 2D

## To do list

### Laplace operator

- [ ] Extend to metric w/ non-diagonal terms (see Christoffel symbols w/ Jibril)
- [ ] Extend to pencil=5 (i.e. the cross 9-point stencil)...
- [ ] ... or to the square 9-point stencil (see (3.17) in 3.5 in LeVeque)
- [ ] Consider writing a GR version to solve Einstein equation

### BCs

- [ ] For NGC=1, to see if it helps, implement a type of BCs where value in GCs is set such as Laplacian in the neighboring inside cell = source term

### Other

- [ ] Modify get_residual such as it validates the solution when eps>eps0 only on the (inside) sides (because it is often the case if BCs are zero or cont, although the solver has converged toward a fairly accurate solution)
- [ ] Clean BiCGSTAB subroutine such as there is no fold/unfold inside
- [ ] Speed up the code by replacing matrix product w/ the actual non-zero products w/ closest neighbors
- [ ] get_avg should probably account for the determinant of the metric
- [ ] Write a proper read_par subroutine like in AMRVAC
- [ ] Write a proper makefile like in AMRVAC

### Solver

- [ ] Put back Jacobi to compare performance
- [ ] Consider using a pre-conditioner
- [ ] Consider implementing a multigrid solver

### Grid

- [ ] Enable stretching in direction 1 only
- [ ] Enable cos stretching for direction 2 (in spherical)
- [ ] For strong stretching (i.e. 40 cells from 1 to 100), it does not always converge => can it be solved?

### MPI

- [ ] Reactivate the calls to MPI library
- [ ] Replace crash subroutine w/ mpistop
- [ ] Parallelize the code
