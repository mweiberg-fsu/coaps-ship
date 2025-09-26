import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

###############################################################################################################################################################################################################
# Insert file name for csv with input data for a time period of outlier data
extreme_file = '.csv'

#Insert figure title
figure_title = 'Variable, Time frame'

#Variable to plot
variable = 'in_T'
###############################################################################################################################################################################################################

dataExtreme = pd.read_csv(extreme_file)
t = dataExtreme['time'].to_numpy() 
lat = dataExtreme['latitude'].to_numpy()
lon = dataExtreme['longitude'].to_numpy()
var = dataExtreme[variable].to_numpy()


for i in range(1, len(SPD)):
    var[i] = float(var[i])


matchingVar = []

for i in range(len(extreme_speed_time)):
    for j in range(len(t)):
        if extreme_speed_time[i] == t[j]:
            matchingVar.append(round(var[j],4))


x = extreme_speedT
y = matchingVar
plt.figure()
plt.plot(x,y,color='blue',linestyle='-',marker='o')
plt.xlabel('hour:minute')
plt.ylabel(varibale)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title(figure_title)
plt.tight_layout()
plt.savefig(f'Extreme_{variable}.png')