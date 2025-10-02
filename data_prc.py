## File: data_prc.py
## Description: This file handles data processing and flux calculations for ship data, including uncertainty estimations and output generation.

from imports import *
from config import *
from constants import *
from utils import *
from MFT23 import *

logger = setup_logger('data_prc', f'{logs_dir}/data_prc.log', level=logging.DEBUG) # Set to DEBUG for detailed trace

# Calculate mean of data lists for usable data samples
def Mean(values):
    values = np.array(values)
    valid = values > -1111.0
    if np.sum(valid) > 0.5 * len(values):
        return np.mean(values[valid])
    else:
        return -1111.0

# Calculate the variance and standard deviation of data lists, uses mean function defined above
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
    
    # Remove from end to start to avoid index shifting
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

# Validates that all inputs are not NaN
def is_valid_input(*args):
    return all(x != -9999.0 for x in args)

# Set input and output directories
input_csvs = directory_destination
create_directory(output_csvs)


# Main function to process data and calculate fluxes
def process_data():
    logger.info("Starting the data processing and flux calculation")
    
    start_time = time.time()  # Start timer
    file_count = 0  # Initialize file counter

    for ship in ships:  # Iterate through each ship
        ship_start = time.time()  # Start timer for ship
        directory_path = f'{input_csvs}/{ship}'  # Directory for current ship
       
        # Iterate through each file in ship directory
        for filename in os.listdir(directory_path):
            file_count += 1
            loop_start = time.time()
            file_path = os.path.join(directory_path, filename)
            
            # Read the CSV file once, skipping metadata rows
            try:
                data = pd.read_csv(file_path, skiprows=[1, 2])  # Adjust skiprows based on file structure
                logger.info(f"Loaded file {file_path} with {len(data)} rows")
            except Exception as e:
                logger.error(f"Failed to read {file_path}: {e}")
                continue

            # Define overlap as in original code
            overlap = window_size // 2  # Define overlap size
            # Validate chunk_size and overlap
            if chunk_size <= 2 * overlap:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}. "
                             f"Ensure chunk_size > 2 * overlap.")
                continue
            if chunk_size <= 0 or overlap < 0:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}. "
                             f"chunk_size must be positive and overlap non-negative.")
                continue

            # Calculate number of chunks
            effective_chunk_size = chunk_size - 2 * overlap
            number_of_chunks = max(1, (len(data) // effective_chunk_size) + (1 if len(data) % effective_chunk_size > 0 else 0))
            logger.info(f"Processing file {file_path} with {number_of_chunks} chunks")

            # Initialize output CSV with headers
            csv_filename = f'{output_csvs}/{filename[:-4]}_processed.csv'
            with open(csv_filename, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(headers)

            # Process file in chunks with overlap
            chunk_count = 0
            row_buffer = []  # Buffer to hold rows before writing to CSV

            # Manually slice the DataFrame into chunks
            for start_idx in range(0, len(data), effective_chunk_size):
                chunk_count += 1
                chunk_start = time.time()
                logger.info(f'Processing chunk {chunk_count}/{number_of_chunks}')

                # Define chunk boundaries with overlap
                chunk_start_idx = max(0, start_idx - overlap)  # Include previous overlap
                chunk_end_idx = min(len(data), start_idx + chunk_size + overlap)  # Include next overlap
                chunk = data.iloc[chunk_start_idx:chunk_end_idx].copy()  # Use .copy() to avoid SettingWithCopyWarning

                # Extract variables from chunk
                airTemp = chunk['in_T'].to_numpy().astype(float)
                hum = chunk['in_RH'].to_numpy().astype(float) / 100
                pressure = chunk['in_P'].to_numpy().astype(float) * 100
                waterTemp = chunk['in_TS'].to_numpy().astype(float)
                truew = chunk['in_SPD'].to_numpy().astype(float)
                time_list = chunk['time'].to_numpy()
                latitude = chunk['latitude'].to_numpy()
                longitude = chunk['longitude'].to_numpy()
                ship_id = chunk['platform_call_sign'].to_numpy()

                # Validate and clean data arrays
                airTemp = check_array(airTemp, -60, 100, -9999.0)
                waterTemp = check_array(waterTemp, -4, 100, -9999.0)
                hum = check_array(hum, 0, 1, -9999.0)
                pressure = check_array(pressure, 0, 120000, -9999.0)
                truew = check_array(truew, 0, 100, -9999.0)

                # Process each data point in the chunk
                for k in range(5, len(truew) - 5):
                    g = k - 5  # Define window for moving calculations
                    h = k + 6  # Define window for moving calculations
                    lhf, shf, m_o, tau = [], [], [], []  # Initialize lists for flux calculations
                    a = truew[g:h]
                    b = pressure[g:h]
                    c = hum[g:h]
                    d = airTemp[g:h]
                    e = waterTemp[g:h]

                    # Vectorized removal of invalid indices
                    valid_mask = (a != -9999.0) & (b != -9999.0) & (c != -9999.0) & (d != -9999.0) & (e != -9999.0)
                    truew_new = a[valid_mask]
                    press_new = b[valid_mask]
                    hum_new = c[valid_mask]
                    airT_new = d[valid_mask]
                    waterT_new = e[valid_mask]
                    z = len(truew_new)  # Length of valid data in window

                    # Perform flux calculations if enough valid data points
                    if z >= 3:  # Require at least 3 valid points in window
                        flux_start = time.time()
                        # Precompute Gaussian samples for all iterations
                        tspd_gauss = np.random.normal(Mean(truew_new), VarStd(truew_new)[1], n)
                        press_gauss = np.random.normal(Mean(press_new), VarStd(press_new)[1], n)
                        hum_gauss = np.random.normal(Mean(hum_new), VarStd(hum_new)[1], n)
                        airT_gauss = np.random.normal(Mean(airT_new), VarStd(airT_new)[1], n)
                        waterT_gauss = np.random.normal(Mean(waterT_new), VarStd(waterT_new)[1], n)

                        for j in range(n):
                            if is_valid_input(tspd_gauss[j], press_gauss[j], hum_gauss[j], airT_gauss[j], waterT_gauss[j]):
                                try:
                                    flux = mft_fluxes(
                                        dyn_in_prm, tspd_gauss[j], dyn_in_val2, sfc_current1, sfc_current2,
                                        convect, press_gauss[j], air_moist_prm, hum_gauss[j], sfc_moist_prm,
                                        sfc_moist_val, salinity, ss_prm, ss_val, airT_gauss[j], sst_prm,
                                        waterT_gauss[j], ref_ht_wind, ref_ht_tq, z_wanted, astab, eqv_neut,
                                        net_heat_flux, warn, flux_model, z0_mom_prm, z0_theta_q_prm, stable_prm,
                                        oil_fract_area, dimensionless_m_o_length, zo_m, missing
                                    )
                                    m_o.append(flux[7])
                                    lhf.append(flux[2])
                                    shf.append(flux[1])
                                    tau.append(flux[3][0])
                                except Exception as e:
                                    logger.error(f'Error in flux calculation at index {k}: {e}')
                            else:
                                logger.warning(f'Skipped iteration due to invalid inputs at index {k}')
                                if j == 0:  # If first try invalid, assume all will be, skip
                                    break

                        flux_end = time.time()
                        flux_total = flux_end - flux_start
                        logger.debug(f'Time for flux calculations at index {k}: {flux_total}')
                    
                    # Prepare row for output
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

                # Write buffered rows to CSV after processing chunk
                if row_buffer:
                    with open(csv_filename, "a", newline="") as file:
                        writer = csv.writer(file)
                        writer.writerows(row_buffer)
                    row_buffer = []  # Clear buffer after writing
                chunk_end = time.time()
                chunk_total = chunk_end - chunk_start
                logger.info(f'Time to complete chunk: {chunk_total}')

            # End of file processing
            loop_end = time.time()
            loop_total = loop_end - loop_start
            logger.info(f'Time to complete file: {loop_total}')

        # End of ship processing
        ship_end = time.time()
        ship_total = ship_end - ship_start
        logger.info(f'Time to complete ship {ship}: {ship_total}')
    
    # End of all processing
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logger.info(f"Data processing completed in {int(minutes)} min {seconds:.2f} sec")