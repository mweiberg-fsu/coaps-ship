## File: run.py
## Description: Main script to run data retrieval and processing

from imports import *
from config import *
from data_get import *
from data_plt import *



## File: run.py
## Description: Main script to run data retrieval and processing

from imports import *
from config import *
from data_get import *
from data_prc import process_data  # Import process_data from data_prc.py

def main():
    logger = setup_logger('run', f'{logs_dir}/run.log', level=logging.INFO)
    logger.info("Starting the data retrieval and processing pipeline")
    print("Starting the data retrieval and processing pipeline.")
    
    start_time = time.time()  # Start timer for entire process
    print(f"Processing ships: {ships}")

    # Step 1: Download data
    logger.info("Initiating data download")
    # run_parallel_downloads(max_threads=proc_num)  # Download data for ships in config.ships
    print("Data download complete. Proceeding to data processing.")
    
    # Step 2: Process downloaded data
    logger.info("Initiating data processing")
    process_data()  # Process data for ships in config.ships
    print("Data processing complete.")
    
    # Log total execution time
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logger.info(f"Data retrieval and processing completed in {int(minutes)} min {seconds:.2f} sec")

if __name__ == "__main__":
    main()
