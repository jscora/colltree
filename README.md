# colltree

colltree is a code that takes collision information from an N-body code with imperfect collisions and generates a collision history tree for each planet and creates an interactive dashboard to display these collision histories and overall collision statistics from all of the simulated runs.

## Features
- generates a .csv table with collision history information for all planets, one with some basic statistics on the collisions in each simulation (i.e. numbers of each type of collision), and a table for a specific planet's collision tree
- generates an interactive dashboard that plots the final planets in all simulations and will plot out the collision history tree and histogram of any planet once clicked on

## Installation



## Usage

The code is set up so that everything can be run from the Jupyter notebook `interactive-collisions.ipynb` inside of Jupyter lab. 


### Important note: Input files

This code reads in Fortran-generated text input files with spaces as dividers. The input files required are based off of the output files generated from the version of SyMBA (Duncan et al 1998) with imperfect collisions (Scora et al 2020). 

The functions to read in the input files are all in util.py, and can be modified to read in different input files as needed. 

The input files needed are:
- impact information for each simulation (impact velocity/escape velocity and impact angle at minimum)
- collision information for each simulation (time of collision, collision type, and masses of both initial and final bodies)
- output dumps for each simulation (time of dump, mass, semi-major axis, eccentrity and inclination of embryos and debris (if included in simulation) that exist at that time)
- input parameter file with minimum embryo mass (if no specific minimum embryo mass is input)

Example files for each of the above-mentioned files are included in the `examples` directory. Core-mass fractions (CMFs) are included throughout the code, and required for these output files as well. If this information is not required for your analysis, then you can include columns of zeros for the relevant CMF columns and remove CMF from the plotting options. 
