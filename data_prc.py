from config import *
from constants import *
from utils import *
from MFT23 import *
from multiprocessing import Pool

logger = setup_logger('data_prc', f'{logs_dir}/data_prc.log', level=logging.INFO)

# Step 1 : Checks that each element in the data array is within reasonable bounds, validate it
def check_array(arr, min_val, max_val, replacement_val):
    modified_arr = np.where(arr < min_val, replacement_val, arr)
    modified_arr = np.where(modified_arr > max_val, replacement_val, modified_arr)
    return modified_arr

# Validates that all inputs are not NaN
def is_valid_input(*args):
    return all(not np.isnan(x) for x in args)

# Set input and output directories
input_csvs = directory_destination
create_directory(output_csvs)

# Step 2. Function to process a single chunk (to be run in parallel)
def process_chunk(args):
    chunk, chunk_idx, number_of_chunks, window_radius, headers = args
    logger.info(f'Processing chunk {chunk_idx+1}/{number_of_chunks} in process {os.getpid()}')
    chunk_start = time.time()
    
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
    airTemp = check_array(airTemp, -60, 100, np.nan)
    waterTemp = check_array(waterTemp, -4, 100, np.nan)
    hum = check_array(hum, 0, 1, np.nan)
    pressure = check_array(pressure, 0, 120000, np.nan)
    truew = check_array(truew, 0, 100, np.nan)

    # Step 3: Compute rolling means and stds
    rolling_window = 2 * window_radius + 1
    min_periods = 3

    airT_series = pd.Series(airTemp)
    airT_mean = airT_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    airT_std = airT_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    waterT_series = pd.Series(waterTemp)
    waterT_mean = waterT_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    waterT_std = waterT_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    hum_series = pd.Series(hum)
    hum_mean = hum_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    hum_std = hum_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    press_series = pd.Series(pressure)
    press_mean = press_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    press_std = press_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    truew_series = pd.Series(truew)
    truew_mean = truew_series.rolling(rolling_window, center=True, min_periods=min_periods).mean().to_numpy()
    truew_std = truew_series.rolling(rolling_window, center=True, min_periods=min_periods).std(ddof=0).to_numpy()

    # Step 4: Process each data point in the chunk
    row_buffer = []
    for k in range(window_radius, len(truew) - window_radius):
        lhf, shf, m_o, tau = [], [], [], []

        mean_tspd = truew_mean[k]
        std_tspd = truew_std[k]

        # Skip if mean_tspd is NaN or std_tspd is negative
        if np.isnan(mean_tspd) or std_tspd < 0:
            row = [
                time_list[k], ship_id[k], latitude[k], longitude[k],
                -9999, -9999, -9999, -9999,
                -9999, -9999, -9999,
                -9999, -9999, -9999, -9999, -9999
            ]
        else: # Only compute fluxes if mean_tspd is valid
            tspd_gauss = np.random.normal(mean_tspd, std_tspd, n)
            press_gauss = np.random.normal(press_mean[k], press_std[k], n)
            hum_gauss = np.random.normal(hum_mean[k], hum_std[k], n)
            airT_gauss = np.random.normal(airT_mean[k], airT_std[k], n)
            waterT_gauss = np.random.normal(waterT_mean[k], waterT_std[k], n)

            for j in range(n): # Loop over n Gaussian samples
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
                    if j == 0:
                        break

        # Prepare output row    
        if len(shf) > 0:
            row = [
                time_list[k], ship_id[k], latitude[k], longitude[k],
                np.std(shf), np.mean(shf),
                np.std(lhf), np.mean(lhf),
                np.std(tau), np.mean(tau),
                np.mean(m_o),
                airT_std[k], waterT_std[k],
                hum_std[k], press_std[k], truew_std[k]
            ]
        else:
            row = [
                time_list[k], ship_id[k], latitude[k], longitude[k],
                -9999, -9999, -9999, -9999,
                -9999, -9999, -9999,
                -9999, -9999, -9999, -9999, -9999
            ]
        row_buffer.append(row)

    chunk_end = time.time()
    chunk_total = chunk_end - chunk_start
    logger.info(f'Time to complete chunk {chunk_idx+1}: {chunk_total}')
    return row_buffer

# Step 5: Main function to process data and calculate fluxes
def process_data():
    logger.info("Starting the data processing and flux calculation")
    
    start_time = time.time()
    file_count = 0

    # Ensure the base output directory exists
    create_directory(output_csvs)

    for ship in ships:
        ship_start = time.time()
        # Input directory for the ship
        directory_path = f'{input_csvs}/{ship}'
        # Output directory for the ship
        ship_output_dir = f'{output_csvs}/{ship}'
        # Create ship-specific output directory
        create_directory(ship_output_dir)
       
        for filename in os.listdir(directory_path):
            file_count += 1
            loop_start = time.time()
            file_path = os.path.join(directory_path, filename)
            
            try:
                data = pd.read_csv(file_path, skiprows=[1, 2])
                logger.info(f"Loaded file {file_path} with {len(data)} rows")
            except Exception as e:
                logger.error(f"Failed to read {file_path}: {e}")
                continue

            overlap = window_size // 2
            if chunk_size <= 2 * overlap:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}.")
                continue
            if chunk_size <= 0 or overlap < 0:
                logger.error(f"Invalid chunk_size ({chunk_size}) or overlap ({overlap}) for file {file_path}.")
                continue

            effective_chunk_size = chunk_size - 2 * overlap
            number_of_chunks = max(1, (len(data) // effective_chunk_size) + (1 if len(data) % effective_chunk_size > 0 else 0))
            logger.info(f"Processing file {file_path} with {number_of_chunks} chunks")

            # Construct output filename: SHIP_ID[year]_processed.csv
            output_filename = f"{filename[:-4]}_processed.csv"  # Append _processed before .csv
            csv_filename = os.path.join(ship_output_dir, output_filename)

            # Write header to CSV
            with open(csv_filename, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(headers)

            # Define window_radius before chunk preparation
            window_radius = 5  # Hardcoded as in original code; adjust if window_size changes

            # Prepare chunks for parallel processing
            chunks = []
            for start_idx in range(0, len(data), effective_chunk_size):
                chunk_start_idx = max(0, start_idx - overlap)
                chunk_end_idx = min(len(data), start_idx + chunk_size + overlap)
                chunk = data.iloc[chunk_start_idx:chunk_end_idx].copy()
                chunks.append((chunk, len(chunks), number_of_chunks, window_radius, headers))

            # Process chunks in parallel using Pool
            with Pool(processes=proc_num) as pool:
                results = pool.map(process_chunk, chunks)

            # Write results to CSV in order
            with open(csv_filename, "a", newline="") as file:
                writer = csv.writer(file)
                for row_buffer in results:
                    writer.writerows(row_buffer)

            loop_end = time.time()
            loop_total = loop_end - loop_start
            logger.info(f'Time to complete file: {loop_total}')

        ship_end = time.time()
        ship_total = ship_end - ship_start
        logger.info(f'Time to complete ship {ship}: {ship_total}')
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logger.info(f"Data processing completed in {int(minutes)} min {seconds:.2f} sec")