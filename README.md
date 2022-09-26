# SphereDiffProblem

Model for solving the linear diffusion equation over sphere surface

### Features:

- Second order aproximation in time and space    
- Calculation parallelization using OpenMP
- Mass conservative
- Numerical squeme:
    - Crank Nicolson squeme (Unconditionally stable)
    - Implicit (no  iterations required)
    - L2 norm conservative
