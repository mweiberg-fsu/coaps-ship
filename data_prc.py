## File: data_prc.py
## Description: This file handles data processing and flux calculations for ship data, including uncertainty estimations and output generation.

from imports import *
from config import *
from utils import *

import numpy as np
from MFT23 import *
from datetime import datetime, timedelta, timezone
import gc
import warnings

logger = setup_logger('data_prc', f'{logs_dir}/data_prc.log', level=logging.DEBUG)

# Calculate mean of data lists for usable data samples
def Mean(values):
    values = np.array(values)
    valid = values > -1111.0
    if np.sum(valid) > 0.5 * len(values):
        return np.mean(values[valid])
    else:
        return -1111.0

# Calculate the variance and standard deviation of data lists, uses Mean function defined above
def VarStd(values):
    values = np.array(values)
    valid = values > -1111.0
    if np.any(valid):
        return np.var(values[valid]), np.std(values[valid])
    else:
        return -1111.0, -1111.0

# Finds the uncertainty of data samples using random gaussian permutations
def Uncertainty(values):
    meanVal = Mean(values)
    _, stdv = VarStd(values)
    if meanVal < -1000.0 or stdv < 0 or stdv == -1111.0:
        return -9999.0
    return np.random.normal(meanVal, stdv)

# Removes matching indices for all elements in the lists for each variable
def remove_matching_indices(list1, list2, list3, list4, list5, target_value):
    indices_to_remove = []
    for i in range(len(list1)):
        if list1[i] == target_value or list2[i] == target_value or list3[i] == target_value or list4[i] == target_value or list5[i] == target_value:
            indices_to_remove.append(i)
    
    for index in sorted(indices_to_remove, reverse=True):
        del list1[index]
        del list2[index]
        del list3[index]
        del list4[index]
        del list5[index]
    
    return list1, list2, list3, list4, list5

# Checks that each element in the data array is within reasonable bounds for the variable
def check_array(arr, min_val, max_val, replacement_val):
    modified_arr = np.where(arr < min_val, replacement_val, arr)
    modified_arr = np.where(modified_arr > max_val, replacement_val, modified_arr)
    return modified_arr

def is_valid_input(*args):
    return all(x != -9999.0 for x in args)

############################################################################################################
#############################################  Constants  ##################################################
############################################################################################################

CONV_CRIT = 0.005	
convect = 0.0 		
warn = 1   		#0
eqv_neut = 0		#Height adjustment is not modified, normal winds. 
z_wanted = 10.0
ref_ht_wind = 17.9	#Height given by SAMOS
ref_ht_tq = 17.9    #Height given by SAMOS     
net_heat_flux = 	5.0		
sst_prm = 1
z0_mom_prm = 6		
z0_theta_q_prm = 0
stable_prm = 0				
wave_ang = 0		
dyn_in_prm = 0		#Wind speed
dyn_in_val2 = 0.0
astab = 1		
ss_prm = 2 			#Wind-wave Stability Param.
ss_val = 48.0       
air_moist_prm = 1	#Relative Humidity in Air
sfc_moist_prm = 1	#Relative Humidity at Surface
sfc_moist_val = 0.98
salinity = 0.0349	#Average salinity
sfc_current1 = 0.0
sfc_current2 = 0.0
oil_fract_area = 0.0
missing = -9999.0
dimensionless_m_o_length = 0.0
zo_m = 0.0 

##############################################################################################################################################################################################
# Select number of permutations per calculation
n = 100 

# Select chunk size for reading csv file
chunk_size = 1000

# Select number of records per rolling window
window_size = 11

# Input directory from config
input_csvs = directory_destination

# Output directory (define here as it's new)
output_csvs = 'data/output/'
create_directory(output_csvs)

#Select proper number for flux model
flux_model = -1		#Blair 23, use 1 for S88
##############################################################################################################################################################################################

# Headers printed to output file
headers = ['time', 'platform_call_sign', 'latitude', 'longitude', 'hfss stdv', 'hfss mean', 'hfls stdv', 'hfls mean', 
           'tau stdv', 'tau mean', 'mean dmo', 'T stdv', 'TS stdv', 'RH stdv', 'P stdv', 'SPD stdv']

def process_data():
    logger.info("Starting the data processing and flux calculation")
    
    start_time = time.time()  # Start timer
    
    file_count = 0
    for ship in ships:
        ship_start = time.time()
        directory_path = f'{input_csvs}/{ship}'
        # Iterate through each file in each ship directory
        for filename in os.listdir(directory_path):
            file_count += 1
            loop_start = time.time()
            file_path = os.path.join(directory_path, filename) 
            data1 = pd.read_csv(file_path, skiprows=[2])
            airTemp1 = data1['in_T'].to_numpy() 
            number_of_chunks = len(airTemp1) // chunk_size + 1
            logger.debug(f"Processing file {file_path} with {number_of_chunks} chunks")
            csv_filename = f'{output_csvs}/{filename[:-4]}_processed.csv'
            with open(csv_filename, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(headers)
            with open(csv_filename, "a", newline="") as file:
                writer = csv.writer(file)
                
                count = 0
                chunk_index = 0
                chunk_count = 0
                column_names = pd.read_csv(file_path, nrows=0).columns.tolist()
                overlap = window_size // 2
                prev_tail = None
                for i, chunk in enumerate(pd.read_csv(file_path, skiprows=2, chunksize=chunk_size - 2*overlap, names=column_names)):
                    if i > 0 and prev_tail is not None:
                        chunk = pd.concat([prev_tail, chunk], ignore_index=True)
                    try:
                        next_chunk = pd.read_csv(file_path, skiprows=2 + (i+1)*(chunk_size - 2*overlap), nrows=overlap, names=column_names)
                    except ValueError:
                        next_chunk = pd.DataFrame()
                    chunk = pd.concat([chunk, next_chunk], ignore_index=True)
                    if len(chunk) >= overlap:
                        prev_tail = chunk.iloc[-overlap:]
                    else:
                        prev_tail = None
                    chunk_count += 1
                    chunk_start = time.time()
                    logger.debug(f'{chunk_count}/{number_of_chunks}')

                    airTemp = chunk['in_T'].to_numpy().astype(float)
                    hum = chunk['in_RH'].to_numpy().astype(float) / 100
                    pressure = chunk['in_P'].to_numpy().astype(float) * 100
                    waterTemp = chunk['in_TS'].to_numpy().astype(float)
                    truew = chunk['in_SPD'].to_numpy().astype(float)
                    time_list = chunk['time'].to_numpy()
                    latitude = chunk['latitude'].to_numpy()
                    longitude = chunk['longitude'].to_numpy()
                    ship_id = chunk['platform_call_sign'].to_numpy()

                    airTemp = check_array(airTemp, -60, 100, -9999.0)
                    waterTemp = check_array(waterTemp, -4, 100, -9999.0)
                    hum = check_array(hum, 0, 1, -9999.0)
                    pressure = check_array(pressure, 0, 120000, -9999.0)
                    truew = check_array(truew, 0, 100, -9999.0)

                    row_buffer = []
                    for k in range(5, len(truew)-5):
                        g = k - 5
                        h = k + 6
                        lhf, shf, m_o, tau = [], [], [], []
                        a = truew[g:h]
                        b = pressure[g:h]
                        c = hum[g:h]
                        d = airTemp[g:h]
                        e = waterTemp[g:h]
                        truew_new = list(a)
                        press_new = list(b)
                        hum_new = list(c)
                        airT_new = list(d)
                        waterT_new = list(e)
                        truew_new, press_new, hum_new, airT_new, waterT_new = remove_matching_indices(truew_new, press_new, hum_new, airT_new, waterT_new, -9999.0)
                        z = len(truew_new)  
                        if z >= 3:
                            j = 0
                            while j < n:
                                tspd_gauss = Uncertainty(truew_new)
                                press_gauss = Uncertainty(press_new)
                                hum_gauss = Uncertainty(hum_new)
                                airT_gauss = Uncertainty(airT_new)
                                waterT_gauss = Uncertainty(waterT_new)
                                
                                if is_valid_input(tspd_gauss, press_gauss, hum_gauss, airT_gauss, waterT_gauss):
                                    try:
                                        flux = mft_fluxes(dyn_in_prm, tspd_gauss, dyn_in_val2, sfc_current1, sfc_current2, convect, press_gauss, air_moist_prm, hum_gauss, sfc_moist_prm, sfc_moist_val, salinity, ss_prm, ss_val, airT_gauss, sst_prm, waterT_gauss, ref_ht_wind, ref_ht_tq, z_wanted, astab, eqv_neut, net_heat_flux, warn, flux_model, z0_mom_prm, z0_theta_q_prm, stable_prm, oil_fract_area, dimensionless_m_o_length, zo_m, missing)
                                        m_o.append(flux[7])
                                        lhf.append(flux[2])
                                        shf.append(flux[1])
                                        tau.append(flux[3][0])
                                        j += 1
                                    except Exception as e:
                                        logger.error(f'error: {e} occurred at index {k}')
                                else:
                                    logger.warning(f'Skipped iteration due to invalid inputs at index {k}')
                                    # Do not increment j, but since always invalid if bad window, break to avoid infinite
                                    if j == 0:  # If first try invalid, assume all will be, skip
                                        break
                        if len(shf) > 0:
                            row = [
                                time_list[k], ship_id[k], latitude[k], longitude[k],
                                np.std(shf), np.mean(shf),
                                np.std(lhf), np.mean(lhf),
                                np.std(tau), np.mean(tau),
                                np.mean(m_o),
                                np.std(airT_new), np.std(waterT_new),
                                np.std(hum_new), np.std(press_new), np.std(truew_new)
                            ]
                        else:
                            row = [
                                time_list[k], ship_id[k], latitude[k], longitude[k],
                                -9999, -9999, -9999, -9999,
                                -9999, -9999, -9999,
                                -9999, -9999, -9999, -9999, -9999
                            ]
                        row_buffer.append(row)
                    
                    if row_buffer:
                        writer.writerows(row_buffer)
                    count += 1
                    chunk_index += 1000
                    chunk_end = time.time()
                    chunk_total = chunk_end - chunk_start
                    logger.info(f'Time to complete chunk: {chunk_total}')
            
            loop_end = time.time()
            loop_total = loop_end - loop_start
            logger.info(f'Time to complete file: {loop_total}')
            logger.debug(len(airTemp1))
        
        ship_end = time.time()
        ship_total = ship_end - ship_start
        logger.info(f'Time to complete ship {ship}: {ship_total}')
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logger.info(f"Data processing completed in {int(minutes)} min {seconds:.2f} sec")