## Filename: data_get.py
## Description: This file handles downloading and validating ship data files.

from imports import *
from config import *
from utils import *

logger = setup_logger('data_get', f'{logs_dir}/data_get.log', level=logging.DEBUG)

# Years for which existing files will be forcibly overwritten
force_download_years = [2025]

def download_ship_year_data(ship, year):
    try:
        ship_directory = f'{directory_destination}/{ship}'
        create_directory(ship_directory)

        file_name = f'{ship_directory}/{ship}{year}.csv'

        # Skip if file exists and year is not in force_download_years
        if os.path.exists(file_name) and year not in force_download_years:
            logger.info(f"Skipping {ship}-{year}, file already exists")
            return (ship, year, "skipped")

        logger.debug(f"Downloading data for ship {ship} in year {year}")

        # ERDDAP data URL
        csv_url = (
            f'https://erddap-samos.coaps.fsu.edu/erddap/tabledap/SAMOS_Fluxes_B23_v301.csv?'
            f'time%2Cplatform_call_sign%2Clatitude%2Clongitude%2Cin_T%2Cin_RH%2Cin_TS%2Cin_SPD%2Cin_P'
            f'&time%3E={year}-01-01'
            f'&time%3C={year}-12-31T23%3A59%3A00Z'
            f'&platform_call_sign=%22{ship}%22'
            f'&latitude%3E=-78.64944&latitude%3C=89.99979'
            f'&longitude%3E=0&longitude%3C=359.9999'
            f'&in_T_qc=1&in_RH_qc=1&in_TS_qc=1&in_SPD_qc=1&in_P_qc=1'
        )

        # Download the CSV file
        req = requests.get(csv_url)
        req.raise_for_status()  # Raises error for bad status codes
        url_content = req.content

        with open(file_name, 'wb') as csv_file:
            csv_file.write(url_content)

        logger.debug(f"Saved file to {file_name}")
        return (ship, year, "downloaded")

    # Throw exceptions for logging
    except Exception as e:
        logger.error(f"Failed to download {ship}-{year}: {e}")
        return (ship, year, "failed")

# Function to validate downloaded files
def validate_downloaded_files():
    logger.info("Starting file validation...")
    for ship in ships: # Iterate over each ship
        for year in years: # Iterate over each year
            file_path = f'{directory_destination}/{ship}/{ship}{year}.csv' # Construct file path
            logger.debug(f'Checking file {file_path} for variable {variable_name}') # Debug log

            # Check if file exists
            if check_variable(file_path, variable_name):
                logger.info(f'{file_path} is good')
            else:
                logger.warning(f'Removing empty or invalid file {file_path}')
                try:
                    os.remove(file_path)
                except Exception as e:
                    logger.error(f"Failed to remove {file_path}: {e}")

# Function to run downloads in parallel
def run_parallel_downloads(max_threads=proc_num):
    tasks = [(ship, year) for ship in ships for year in years]
    total = len(tasks)

    logger.info(f"Starting parallel downloads with {proc_num} threads")
    completed = {"downloaded": 0, "skipped": 0, "failed": 0}

    # Use ThreadPoolExecutor for parallel downloads
    with ThreadPoolExecutor(max_workers=proc_num) as executor:
        future_to_task = {
            executor.submit(download_ship_year_data, ship, year): (ship, year)
            for ship, year in tasks
        }

        # Process completed futures as they finish
        for count, future in enumerate(as_completed(future_to_task), start=1):
            ship, year = future_to_task[future]
            try:
                _, _, status = future.result()
                completed[status] += 1

                # Log the result based on status
                if status == "downloaded":
                    logger.info(f'{count}/{total} - Downloaded {ship}-{year}')
                elif status == "skipped":
                    logger.info(f'{count}/{total} - Skipped {ship}-{year} (already exists)')
                else:
                    logger.warning(f'{count}/{total} - Failed {ship}-{year}')
            except Exception as e:
                logger.error(f'{count}/{total} - Exception downloading {ship}-{year}: {e}')
                completed["failed"] += 1

    logger.info("All downloads completed.")
    logger.info(f"Summary: {completed['downloaded']} downloaded, {completed['skipped']} skipped, {completed['failed']} failed")

    validate_downloaded_files() # Validate files after all downloads
