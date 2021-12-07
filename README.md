# ACSE-4-SPH

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

### Installation Guide

Move into the directory into which you want to clone the repository. Then execute 
```
git clone https://github.com/acse-2019/acse-4-sph-garten.git
```
in the terminal. Then change directory into the root of the repo
```
cd acse-4-sph-garten
```

### User instructions

To run the simulation to generate files for post-processing, 
```
cd src
```
and build and run 
```
SPH_Snippet.cpp
```
Post-processing can be performed using the 
```
garten_plot_shp.ipynb
```
file, located inside the tests directory.

### Documentation

The code documentation can be found inside the 
```
Documentation.txt
```
file.

### Testing

The tool includes tests, which you can use to check its operation on your system. With the code compiled, these can be run
with

```
python run_tests.py
```
### Results

Analysis of the results of the simulation can be found inside the 
```
data_imgs
```
folder inside the 
```
pj4_data_git
```
folder.

