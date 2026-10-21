<!-- # Exploration of Numerical Methods for PDEs
This repository contains Python files in which numerical methods for various PDEs have been implemented.

See `example_sims` for the simulations that were generated using the Python code (mp4s and pngs).  
See `FDM` for the Python scripts use to generate the simulations. These are all based on the finite difference method.  
See `Notes` for detailed notes I took while learning about FDM. These include derivations of some finite difference equations 
with some basic error analysis, as well as an analysis of the stability of various numerical schemes using Fourier analysis.

The notes are based off of
- the 6th chapter of Haberman’s Applied PDEs with Fourier Series and Boundary Value Problems (4th ed.)
- lectures from MIT’s graduate-level course on Numerical Methods for Partial Differential Equations (16.920J)

The following are the dependencies of the scripts:
```
numpy
matplotlib
scipy
```

Here is an example simulation:
![example sim](/example_sims/2D_poisson.png) -->


# Exploration of Numerical Methods for PDEs

This repository contains Python implementations and detailed notes exploring numerical methods for solving Partial Differential Equations (PDEs) using Finite Difference Methods (FDM).

## Repository Structure
```
├── notes/
│   └── Numerical_Methods_for_PDEs.pdf   # Detailed notes on FDM, Taylor series, and stability
├── pdes/
│   ├── __init__.py
│   ├── heat_1D.py                       # 1D heat equation solver
│   ├── heat_2D.py                       # 2D heat equation solver
│   └── poisson_2D.py                    # 2D Poisson equation solver
├── simulations/
│   ├── heat_1D_sim.mp4                  # Generated 1D heat simulation video
│   ├── heat_2D_sim.mp4                  # Generated 2D heat simulation video
│   └── poisson_2D.png                   # Generated 2D Poisson solution surface plot
├── main.py                              # Entry point script
└── README.md                            # Repository documentation
```

## Theoretical Background

The derivations, error analyses, and stability proofs (including von Neumann stability analysis) documented in notes/Numerical_Methods_for_PDEs.pdf are based on:

- Chapter 6 of Haberman’s Applied PDEs with Fourier Series and Boundary Value Problems (4th ed.)
- Lectures from MIT’s graduate-level course Numerical Methods for Partial Differential Equations (16.920J)

## Dependencies

The code relies on the following Python packages:
```
numpy
matplotlib
scipy
numba
```

Here is an example simulation for the 2D heat equation:
<video src="simulations/heat_2D_sim.mp4" controls width="640" height="640">
</video>
<!-- Visualizations & Simulations

Simulation Files

1D Heat Equation Simulation: simulations/heat_1D_sim.mp4

2D Heat Equation Simulation: simulations/heat_2D_sim.mp4

2D Poisson Equation Solution

Here is an example steady-state numerical solution generated for the 2D Poisson equation: -->

<!-- Future Additions

The following enhancements and extensions are planned for future development:

JAX JIT for GPU Acceleration: Transitioning core numerical kernels and grid updates from CPU-bound Numba implementations to JAX. This will enable automatic differentiation, seamless multi-GPU/accelerator scaling, and significantly faster execution times for high-resolution spatiotemporal meshes.

Wave Equation & Nonhomogeneous Heat Cases: Expanding the PDE solver suite to include the hyperbolic wave equation (modeling vibrating strings/membranes) alongside fully generalized nonhomogeneous heat equation formulations (incorporating time- and space-dependent internal source terms).

Memory-Efficient Grid Management for Simulations: Optimizing video generation and grid history storage. Instead of keeping full 3D/4D arrays ([time, x, y]) in RAM across thousands of iterations, future iterations will implement streaming writes, memory-mapped files, or on-the-fly frame rendering to drastically lower memory footprints for large-scale simulations. -->
