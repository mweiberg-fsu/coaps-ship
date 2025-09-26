# data_get.py
from imports import *
from config import *
from utils import *

logger = setup_logger('data_get', f'{logs_dir}/data_get.log', level=logging.DEBUG)

def download_ship_year_data(ship, year):
    try:
        ship_directory = f'{directory_destination}/{ship}'
        create_directory(ship_directory)
        logger.debug(f'Downloading data for ship {ship} in year {year}')
        
        csv_url = f'https://erddap-samos.coaps.fsu.edu/erddap/tabledap/SAMOS_Fluxes_B23_v301.csv?time%2Cplatform_call_sign%2Clatitude%2Clongitude%2Cin_T%2Cin_RH%2Cin_TS%2Cin_SPD%2Cin_P&time%3E={year}-01-01&time%3C={year}-12-31T23%3A59%3A00Z&platform_call_sign=%22{ship}%22&latitude%3E=-78.64944&latitude%3C=89.99979&longitude%3E=0&longitude%3C=359.9999&in_T_qc=1&in_RH_qc=1&in_TS_qc=1&in_SPD_qc=1&in_P_qc=1'

        file_name = f'{ship_directory}/{ship}{year}.csv'
        req = requests.get(csv_url)
        url_content = req.content

        with open(file_name, 'wb') as csv_file:
            csv_file.write(url_content)

        return (ship, year, True)
    except Exception as e:
        logger.error(f"Failed to download {ship}-{year}: {e}")
        return (ship, year, False)

def validate_downloaded_files():
    for ship in ships:
        for year in years:
            file_path = f'{directory_destination}/{ship}/{ship}{year}.csv'
            logger.debug(f'Checking file {file_path} for variable {variable_name}')
            if check_variable(file_path, variable_name):
                logger.info(f'{file_path} is good')
            else:
                logger.warning(f'Removing empty or invalid file {file_path}')
                os.remove(file_path)

def run_parallel_downloads(max_threads=8):
    tasks = [(ship, year) for ship in ships for year in years]
    total = len(tasks)

    logger.info(f"Starting parallel downloads with {max_threads} threads")

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        future_to_task = {executor.submit(download_ship_year_data, ship, year): (ship, year) for ship, year in tasks}
        
        for count, future in enumerate(as_completed(future_to_task), start=1):
            ship, year = future_to_task[future]
            try:
                _, _, success = future.result()
                logger.info(f'{count}/{total} - Downloaded {ship}-{year}' if success else f'{count}/{total} - Failed {ship}-{year}')
            except Exception as e:
                logger.error(f'{count}/{total} - Error downloading {ship}-{year}: {e}')

    logger.info("All downloads completed. Starting file validation.")
    validate_downloaded_files()
