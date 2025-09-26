import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time


def percent_difference(a,b):
    return((b-a) / ((a+b) / 2)) * 100
    

######################################################################################################################################################################################

# Input csv file for one run of uncertainty calculations
file_run1 = '.csv'

# Input csv file for a second run of uncertainty calculations
file_run2 = '.csv'

#Insert number of permutations done on the uncertainty calculations
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

extreme_speed = []
extreme_speedT = []
extreme_speed_time = []
for i in range(len(SPDa)):
    if shfSa[i] > 15:
        extreme_speed.append(shfSa[i])
        # # print(row1[i][0][11:16])
        extreme_speedT.append(row1[i][0][11:16])
        extreme_speed_time.append(row1[i][0])
        print(row1[i][0])



#####################################################################################
########################## Create One-to-One Scatter Plots ##########################
#####################################################################################

plt.figure()
plt.scatter(shfSa, shfSb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'SHF Standard Deviation n={n}')
plt.savefig(f'shfStdv{n}.png')

plt.figure()
plt.scatter(shfMa, shfMb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'SHF Mean n={n}')
plt.savefig(f'shfMean{n}.png')

plt.figure()
plt.scatter(lhfSa, lhfSb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'LHF Standard Deviation n={n}')
plt.savefig(f'lhfStdv{n}.png')

plt.figure()
plt.scatter(lhfMa, lhfMb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'LHF Mean n={n}')
plt.savefig(f'lhfMean{n}.png')

plt.figure()
plt.scatter(tauSa, tauSb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'tau Standard Deviation n={n}')
plt.savefig(f'tauStdev{n}.png')

plt.figure()
plt.scatter(tauMa, tauMb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title('tau Mean n=150')
plt.savefig('tauMean150.png')

plt.figure()
plt.scatter(da, db)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'Mean dmo n={n}')
plt.savefig(f'dmoMean{n}.png')

plt.figure()
plt.scatter(Ta, Tb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'T Standard Deviation n={n}')
plt.savefig(f'Tstdv{n}.png')

plt.figure()
plt.scatter(TSa, TSb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'TS Standard Deviation n={n}')
plt.savefig(f'TSstdev{n}.png')

plt.figure()
plt.scatter(Pa, Pb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'P Standard Deviation n={n}')
plt.savefig(f'Pstdv{n}.png')

plt.figure()
plt.scatter(RHa, RHb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'RH Standard Deviation n={n}')
plt.savefig(f'RHstdv{n}.png')

plt.figure()
plt.scatter(SPDa, SPDb)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title(f'SPD Standard Deviation n={n}')
plt.savefig(f'SPDstdv{n}.png')

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
    meanDMO.append(percent_difference(row2[i][10],row1[i][10]))
    Tstdv.append(percent_difference(row2[i][11], row1[i][11]))
    TSstdv.append(percent_difference(row2[i][12], row1[i][12]))
    Pstdv.append(percent_difference(row2[i][13], row1[i][13]))
    RHstdv.append(percent_difference(row2[i][14], row1[i][14]))
    SPDstdv.append(percent_difference(row2[i][15], row1[i][15]))

