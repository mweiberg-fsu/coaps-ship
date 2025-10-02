# Flux Uncertainty Calculations and Visualizations
This directory contains the files to calculate flux uncertainties using the Monte Carlo Method. Input values are retrieved from 37 vessels through Erddap and stored by year in csv files.

## Project Overview
Generate uncertainty values for calculated latent and sensible heat fluxes, stress magnitude, and dimensionless obukhov as a potential future parameter.

## Instructions
1. `python run.py` will execute the code.
2. individual functions (e.g., run_parallel_downloads, process_data) can be commented out

## File Overview
- **run.py** - This is the parent file to run *data_get*, *data_prc*, and *data_plt* functions.
- **imports.py** - Stores necessary Python modules and libraries to run the program
- **config.py** - Config options for running the program (e.g., ships, years, processor cores, and directories)
- **constants.py** - Contains constants necessary to run calculations
- **data_get.py** - Functions to download input files for each vessel (interated over by each year)
- **data_prc.py** - Functions to calculate output variables for uncertainty
- **data_plt.py** - Functions to plot one-to-one scatter plots to check/compare the differences between runs
- **utils.py** - Contains generalized functions for creating directories, logging, and other multi-purpose utility
- **MFT23.py** - Functions for flux calculation by Dr. Bourassa

## Directory Overview
- **logs/data_get.log** - Log file for debugging and tracking functions and output in *data_get.py*
- **logs/data_prc.log** - Log file for debugging and tracking functions and output in *data_prc.py*
- **data/input/SHIP_ID/*.csv** - Contains individual folders for each vessel with csv files for each active year ship was in service

## Workflow
1. Program is executed by 'python run.py'.
1. 'data_get.py' downloads the necessary files for given config options.
1. 'data_prc.py' then processes each individual csv file.
1. 'MFT23.py' calculates fluxes for permutated values.
1. STD is calculated for each set of permutated values.
1. 'data_plt.py' generates one-to-one scatter plots.

