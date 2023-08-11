# colltree
[![DOI](https://zenodo.org/badge/504667062.svg)](https://zenodo.org/badge/latestdoi/504667062)

colltree is a package that helps you better visualize the formation histories of planets in N-body simulations with imperfect collisions. It takes output files from the simulation and generates a collision history tree for each planet, then creates an interactive dashboard to display both these individual collision histories and the overall collision statistics from all of the simulated runs. colltree is built to run with the 9 collision types used in SyMBA, but can be altered for codes with different numbers of collision types.

## Features
- an interactive dashboard that plots the collision statistics of the final planets in all simulations and will plot out the collision history tree (see [Scora et al. (2022)](https://doi.org/10.3847/1538-4357/ac9cda) for a full description of the collision history plot) and collision histogram of any planet once clicked on
- creates a .csv table with collision history information for all planets
- creates a .csv table with some basic statistics on the collisions in each simulation (i.e. numbers of each type of collision)
- generates a table for each specific planet's collision tree on demand

## Installation

You can download the source code from github and run `python setup.py` in the directory. 


## Usage

The code is set up so that everything can be run from the Jupyter notebook `interactive-collisions.ipynb` inside of Jupyter lab. 


### Important note: Input files

This code reads in Fortran-generated text input files with spaces as dividers. The input files required are based off of the output files generated from the version of SyMBA [(Duncan et al 1998)](https://iopscience.iop.org/article/10.1086/300541/pdf) with imperfect collisions [(Scora et al 2020)](https://academic.oup.com/mnras/advance-article/doi/10.1093/mnras/staa568/5762783?guestAccessKey=10e0ee18-dd24-4987-bfcd-a76aad9fec12). 

The functions to read in the input files are all in util.py, and can be modified to read in different input files as needed. 

The input files needed are:
- impact information for each simulation (impact velocity/escape velocity and impact angle at minimum)
- collision information for each simulation (time of collision, collision type, and masses of both initial and final bodies)
- output dumps for each simulation (time of dump, mass, semi-major axis, eccentrity and inclination of embryos and debris, if included in simulation, that exist at that time)
- input parameter file with minimum embryo mass (if no specific minimum embryo mass is input)

Example files for each of the above-mentioned files are included in the `examples` directory. Core-mass fractions (CMFs) are included throughout the code, and required for these output files as well. If this information is not required for your analysis, then you can include columns of zeros for the relevant CMF columns and remove CMF from the plotting options. 
