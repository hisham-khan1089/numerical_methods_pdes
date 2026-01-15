# Exploration of Numerical Methods for PDEs
This repository contains Python files in which various numerical methods for PDEs have been implemented.

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
![example sim](/example_sims/nonhomogeneous_2D_heat_sim.mp4)
