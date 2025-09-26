# Flux Uncertainty Calculations and Visualizations

This directory contains the files to calculate flux uncertainties using the Monte Carlo Method using input values retrieved from Erddap, as well as the files to generate quality-checking plots.

## Project Overview

Generate uncertainty values for calculated latent and sensible heat fluxes, stress magnitude, and dimensionless obukhov as a potential future parameter to be available in Erddap.

--- 

## Directory Overview

- **/home/lfox/rdr24/csv_trial/erddap_analysis/Uncertainty_Calculations/compare.py** – Create one to one plots scatter to check and compare the differences between two runs
- **/home/lfox/rdr24/csv_trial/erddap_analysis/Uncertainty_Calculations/download.py** – Download large data input files and delete files with no data 
- **/home/lfox/rdr24/csv_trial/erddap_analysis/Uncertainty_Calculations/extremeTimeSeries.py** – Time series plots for isolated periods with an extreme value
- **/home/lfox/rdr24/csv_trial/erddap_analysis/Uncertainty_Calculations/MFT23.py** – Functions for flux calculation by Mark A. Bourassa
- **/home/lfox/rdr24/csv_trial/erddap_analysis/Uncertainty_Calculations/README.md** – README file
- **/home/lfox/rdr24/csv_trial/erddap_analysis/Uncertainty_Calculations/uncertainty_calc.py** – Main Code that performs all uncertainty calculations

---

## Workflow 

1. CSV files are downloaded in 'download.py'.
2. 'uncertainty_calc.py' performs permutations bounded by the values in the moving window.
3. 'MFT23.py' calculates a flux for each set of permutated values.
4. A Standard deviation is calculted for each set of permutated values.
5. 'compare.py' generates one to one scatter plots based on multiple runs of the code.

---

## Outputs 

### **CSV Files**

#### **Input CSV Files**
- **Location:** '/Net/mdc/people/lfox/ship_csvs_qc'

#### **Output CSV Files**
- **Location:** '/Net/mdc/people/lfox/csv_ship_outputs_qc'

### **Plots**

#### **One-to-one Scatter Plots**

#### **Time Series Plots**

---

## Instructions

### Generate Input CSVs
1. Run 'download.py' with user's specifications in lines 26-42.
2. Alternatively, use dowloaded quality controlled input CSVs from: '/Net/mdc/people/lfox/ship_csvs_qc'.

## Generate Uncertainty CSVs
1. Run 'uncertainty_calc.py' to generate CSVs with uncertainty values, specify input/output file directories in lines 117-121.
'MFT23.py' must be in the same directory for code to work.

## Generate one to one plots for multiple runs
1. Generate 2 output uncertainty CSVs for two separate runs with 'uncertainty_calc.py'.
2. Run compare.py with these two output files, specify their locations in lines 14-18.

## Generate time series plots for extreme case
1. Download a sample CSV file of data in which a known extreme value occurs
2. Insert this file in line 9 of 'extremeTimeSeries.py'
3. Run 'extremeTimeSeries.py'

---

**To Do:** Create a filter for extreme flux standard deviation values that prompts a recalculation.

---

**Estimated Runtime:** ~1 hour per year of a single ship's data

