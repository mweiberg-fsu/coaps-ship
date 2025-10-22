import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from config import *
from utils import *

create_directory(imgs_dir)

def percent_difference(a, b):
    return ((b - a) / ((a + b) / 2)) * 100

######################################################################################################################################################################################

# Input csv file for one run of uncertainty calculations
file_run1 = 'data/output/KAOU/KAOU2025_processed_v1.csv'

# Input csv file for a second run of uncertainty calculations
file_run2 = 'data/output/KAOU/KAOU2025_processed_v2.csv'

# Insert number of permutations done on the uncertainty calculations
n = 100

######################################################################################################################################################################################

dataRun1 = pd.read_csv(file_run1)
dataRun2 = pd.read_csv(file_run2)

time1 = dataRun1['time'].to_numpy() 
call1 = dataRun1['platform_call_sign'].to_numpy()
lat1 = dataRun1['latitude'].to_numpy()
lon1 = dataRun1['longitude'].to_numpy()
shfStd1 = dataRun1['hfss stdv'].to_numpy()
shfMean1 = dataRun1['hfss mean'].to_numpy()
lhfStd1 = dataRun1['hfls stdv'].to_numpy()
lhfMean1 = dataRun1['hfls mean'].to_numpy()
tauStd1 = dataRun1['tau stdv'].to_numpy()
tauMean1 = dataRun1['tau mean'].to_numpy()
dmoMean1 = dataRun1['mean dmo'].to_numpy()
Tstdv1 = dataRun1['T stdv'].to_numpy()
TSstdv1 = dataRun1['TS stdv'].to_numpy()
Pstdv1 = dataRun1['P stdv'].to_numpy()
RHstdv1 = dataRun1['RH stdv'].to_numpy()
SPDstdv1 = dataRun1['SPD stdv'].to_numpy()

time2 = dataRun2['time'].to_numpy() 
call2 = dataRun2['platform_call_sign'].to_numpy()
lat2 = dataRun2['latitude'].to_numpy()
lon2 = dataRun2['longitude'].to_numpy()
shfStd2 = dataRun2['hfss stdv'].to_numpy()
shfMean2 = dataRun2['hfss mean'].to_numpy()
lhfStd2 = dataRun2['hfls stdv'].to_numpy()
lhfMean2 = dataRun2['hfls mean'].to_numpy()
tauStd2 = dataRun2['tau stdv'].to_numpy()
tauMean2 = dataRun2['tau mean'].to_numpy()
dmoMean2 = dataRun2['mean dmo'].to_numpy()
Tstdv2 = dataRun2['T stdv'].to_numpy()
TSstdv2 = dataRun2['TS stdv'].to_numpy()
Pstdv2 = dataRun2['P stdv'].to_numpy()
RHstdv2 = dataRun2['RH stdv'].to_numpy()
SPDstdv2 = dataRun2['SPD stdv'].to_numpy()

row1 = []
row2 = []
for i in range(len(time2)):
    for j in range(len(time1)):
        if time2[i] == time1[j]:
            if lat2[i] == lat1[j]:
                if lon2[i] == lon1[j]:
                    row2.append(
                        [time2[i], call2[i], lat2[i], lon2[i], shfStd2[i], shfMean2[i], lhfStd2[i], lhfMean2[i],
                        tauStd2[i], tauMean2[i], dmoMean2[i], Tstdv2[i], TSstdv2[i], Pstdv2[i], RHstdv2[i], SPDstdv2[i]]
                     )
                    row1.append(
                        [time1[j], call1[j], lat1[j], lon1[j], shfStd1[j], shfMean1[j], lhfStd1[j], lhfMean1[j],
                        tauStd1[j], tauMean1[j], dmoMean1[j], Tstdv1[j], TSstdv1[j], Pstdv1[j], RHstdv1[j], SPDstdv1[j]]
                     )

shfStdv, shfMean, lhfStdv, lhfMean, tauStdv, tauMean, meanDMO = [],[],[],[],[],[],[]
Tstdv, TSstdv, Pstdv, RHstdv, SPDstdv = [],[],[],[],[]

shfSa, shfMa, lhfSa, lhfMa, tauSa, tauMa, da = [],[],[],[],[],[],[]
Ta, TSa, Pa, RHa, SPDa = [],[],[],[],[]

shfSb, shfMb, lhfSb, lhfMb, tauSb, tauMb, db = [],[],[],[],[],[],[]
Tb, TSb, Pb, RHb, SPDb = [],[],[],[],[]

for i in range(len(row2)):
    shfSa.append(row1[i][4])
    shfSb.append(row2[i][4])
    shfMa.append(row1[i][5])
    shfMb.append(row2[i][5])
    lhfSa.append(row1[i][6])
    lhfSb.append(row2[i][6])
    lhfMa.append(row1[i][7])
    lhfMb.append(row2[i][7])
    tauSa.append(row1[i][8])
    tauSb.append(row2[i][8])
    tauMa.append(row1[i][9])
    tauMb.append(row2[i][9])
    da.append(row1[i][10])
    db.append(row2[i][10])
    Ta.append(row1[i][11])
    Tb.append(row2[i][11])
    TSa.append(row1[i][12])
    TSb.append(row2[i][12])
    Pa.append(row1[i][13])
    Pb.append(row2[i][13])
    RHa.append(row1[i][14])
    RHb.append(row2[i][14])
    SPDa.append(row1[i][15])
    SPDb.append(row2[i][15])

# Debug: Check range of shfSa to verify threshold
print(f"SHF Standard Deviation (Run 1) - Min: {min(shfSa):.2f}, Max: {max(shfSa):.2f}, Mean: {np.mean(shfSa):.2f}")

extreme_speed = []
extreme_speedT = []
extreme_speed_time = []
extreme_indices = []  # To store indices for cross-referencing
for i in range(len(shfSa)):
    if shfSa[i] > 15:
        extreme_speed.append(shfSa[i])
        extreme_speedT.append(row1[i][0][11:16])
        extreme_speed_time.append(row1[i][0])
        extreme_indices.append(i)
        print(f"Extreme SHF Std Found at index {i}: {shfSa[i]:.2f} W/m², Time: {row1[i][0]}")

# For checking outliers: Print cross-references with other variables for extreme cases
print("\nExtreme SHF Std Dev Cases (Run 1 > 15):")
if extreme_indices:
    for idx in extreme_indices:
        print(f"Time: {row1[idx][0]}")
        print(f"SHF Std (Run 1): {shfSa[idx]:.2f} W/m², Wind Speed Std (Run 1): {SPDa[idx]:.2f} m/s, Temp Std (Run 1): {Ta[idx]:.2f} °C")
        print(f"SHF Std (Run 2): {shfSb[idx]:.2f} W/m², Wind Speed Std (Run 2): {SPDb[idx]:.2f} m/s, Temp Std (Run 2): {Tb[idx]:.2f} °C")
        print("---")
else:
    print("No SHF Standard Deviation values > 15 found in Run 1.")

#####################################################################################
########################## Create One-to-One Scatter Plots ##########################
#####################################################################################

def create_scatter_plot(x, y, xlabel, ylabel, title, filename):
    plt.figure()
    plt.scatter(x, y, alpha=0.6)
    min_val = min(min(x), min(y))
    max_val = max(max(x), max(y))
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', label='y=x')  # Add diagonal line
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.savefig(filename)
    plt.close()

create_scatter_plot(shfSa, shfSb, 'Run 1 (W/m²)', 'Run 2 (W/m²)', f'SHF Standard Deviation n={n}', f'{imgs_dir}/shfStdv{n}.png')
create_scatter_plot(shfMa, shfMb, 'Run 1 (W/m²)', 'Run 2 (W/m²)', f'SHF Mean n={n}', f'{imgs_dir}/shfMean{n}.png')
create_scatter_plot(lhfSa, lhfSb, 'Run 1 (W/m²)', 'Run 2 (W/m²)', f'LHF Standard Deviation n={n}', f'{imgs_dir}/lhfStdv{n}.png')
create_scatter_plot(lhfMa, lhfMb, 'Run 1 (W/m²)', 'Run 2 (W/m²)', f'LHF Mean n={n}', f'{imgs_dir}/lhfMean{n}.png')
create_scatter_plot(tauSa, tauSb, 'Run 1 (N/m²)', 'Run 2 (N/m²)', f'Tau Standard Deviation n={n}', f'{imgs_dir}/tauStdev{n}.png')
create_scatter_plot(tauMa, tauMb, 'Run 1 (N/m²)', 'Run 2 (N/m²)', f'Tau Mean n={n}', f'{imgs_dir}/tauMean{n}.png')
create_scatter_plot(da, db, 'Run 1', 'Run 2', f'Mean DMO n={n}', f'dmoMean{n}.png')
create_scatter_plot(Ta, Tb, 'Run 1 (°C)', 'Run 2 (°C)', f'Temperature Standard Deviation n={n}', f'{imgs_dir}/Tstdv{n}.png')
create_scatter_plot(TSa, TSb, 'Run 1 (°C)', 'Run 2 (°C)', f'Sea Surface Temperature Standard Deviation n={n}', f'{imgs_dir}/TSstdv{n}.png')
create_scatter_plot(Pa, Pb, 'Run 1 (hPa)', 'Run 2 (hPa)', f'Pressure Standard Deviation n={n}', f'{imgs_dir}/Pstdv{n}.png')
create_scatter_plot(RHa, RHb, 'Run 1 (%)', 'Run 2 (%)', f'Relative Humidity Standard Deviation n={n}', f'{imgs_dir}/RHstdv{n}.png')
create_scatter_plot(SPDa, SPDb, 'Run 1 (m/s)', 'Run 2 (m/s)', f'Wind Speed Standard Deviation n={n}', f'{imgs_dir}/SPDstdv{n}.png')

#####################################################################################
##################### Calculate Percent Difference Between Runs #####################
#####################################################################################

for i in range(len(row2)):
    shfStdv.append(percent_difference(row2[i][4], row1[i][4]))
    shfMean.append(percent_difference(row2[i][5], row1[i][5]))
    lhfStdv.append(percent_difference(row2[i][6], row1[i][6]))
    lhfMean.append(percent_difference(row2[i][7], row1[i][7]))
    tauStdv.append(percent_difference(row2[i][8], row1[i][8]))
    tauMean.append(percent_difference(row2[i][9], row1[i][9]))
    meanDMO.append(percent_difference(row2[i][10], row1[i][10]))
    Tstdv.append(percent_difference(row2[i][11], row1[i][11]))
    TSstdv.append(percent_difference(row2[i][12], row1[i][12]))
    Pstdv.append(percent_difference(row2[i][13], row1[i][13]))
    RHstdv.append(percent_difference(row2[i][14], row1[i][14]))
    SPDstdv.append(percent_difference(row2[i][15], row1[i][15]))

# Plot histograms of percent differences
def create_histogram(data, title, filename):
    plt.figure()
    plt.hist(data, bins=20, edgecolor='black')
    plt.title(title)
    plt.xlabel('Percent Difference (%)')
    plt.ylabel('Frequency')
    plt.savefig(filename)
    plt.close()

create_histogram(shfStdv, f'Percent Difference: SHF Standard Deviation n={n}', f'{imgs_dir}/shfStdv_diff_hist{n}.png')
create_histogram(shfMean, f'Percent Difference: SHF Mean n={n}', f'{imgs_dir}/shfMean_diff_hist{n}.png')
create_histogram(lhfStdv, f'Percent Difference: LHF Standard Deviation n={n}', f'{imgs_dir}/lhfStdv_diff_hist{n}.png')
create_histogram(lhfMean, f'Percent Difference: LHF Mean n={n}', f'{imgs_dir}/lhfMean_diff_hist{n}.png')
create_histogram(tauStdv, f'Percent Difference: Tau Standard Deviation n={n}', f'{imgs_dir}/tauStdv_diff_hist{n}.png')
create_histogram(tauMean, f'Percent Difference: Tau Mean n={n}', f'{imgs_dir}/tauMean_diff_hist{n}.png')
create_histogram(meanDMO, f'Percent Difference: Mean DMO n={n}', f'{imgs_dir}/dmoMean_diff_hist{n}.png')
create_histogram(Tstdv, f'Percent Difference: Temperature Standard Deviation n={n}', f'{imgs_dir}/Tstdv_diff_hist{n}.png')
create_histogram(TSstdv, f'Percent Difference: Sea Surface Temperature Standard Deviation n={n}', f'{imgs_dir}/TSstdv_diff_hist{n}.png')
create_histogram(Pstdv, f'Percent Difference: Pressure Standard Deviation n={n}', f'{imgs_dir}/Pstdv_diff_hist{n}.png')
create_histogram(RHstdv, f'Percent Difference: Relative Humidity Standard Deviation n={n}', f'{imgs_dir}/RHstdv_diff_hist{n}.png')
create_histogram(SPDstdv, f'Percent Difference: Wind Speed Standard Deviation n={n}', f'{imgs_dir}/SPDstdv_diff_hist{n}.png')

# Calculate and print Pearson correlation coefficients
def print_correlation(x, y, var_name):
    corr, p_value = pearsonr(x, y)
    print(f"Pearson Correlation for {var_name}: {corr:.4f} (p-value: {p_value:.4f})")

print("\nCorrelation Coefficients Between Runs:")
print_correlation(shfSa, shfSb, "SHF Standard Deviation")
print_correlation(shfMa, shfMb, "SHF Mean")
print_correlation(lhfSa, lhfSb, "LHF Standard Deviation")
print_correlation(lhfMa, lhfMb, "LHF Mean")
print_correlation(tauSa, tauSb, "Tau Standard Deviation")
print_correlation(tauMa, tauMb, "Tau Mean")
print_correlation(da, db, "Mean DMO")
print_correlation(Ta, Tb, "Temperature Standard Deviation")
print_correlation(TSa, TSb, "Sea Surface Temperature Standard Deviation")
print_correlation(Pa, Pb, "Pressure Standard Deviation")
print_correlation(RHa, RHb, "Relative Humidity Standard Deviation")
print_correlation(SPDa, SPDb, "Wind Speed Standard Deviation")