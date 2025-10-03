## File: run.py
## Description: Main script to run data retrieval and processing

from imports import *
from config import *
from data_get import *
from data_plt import *
from data_prc import process_data_prl
from data_prc_nonprl import process_data  # Import process_data from data_prc.py

def main():
    # Step 1: Assign logging
    logger = setup_logger('run', f'{logs_dir}/run.log', level=logging.INFO)
    logger.info("Starting the data retrieval and processing pipeline")
    print("Starting the data retrieval and processing pipeline.")
    
    start_time = time.time()  # Start timer for entire process
    print(f"Processing ships: {ships}")

    # Step 2: Download data
    logger.info("Initiating data download")
    run_parallel_downloads(max_threads=proc_num)  # Main function for downloading data
    print("Data download complete. Proceeding to data processing.")
    
    # Step 3: Process downloaded data
    logger.info("Initiating data processing")
    process_data_prl() # Process ship data for uncertainty calculations
    print("Data processing complete.")

    # Step 4: Plot data (next step)
    ## plot_data()  # temp


    # Log total execution time
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logger.info(f"Data retrieval and processing completed in {int(minutes)} min {seconds:.2f} sec")

if __name__ == "__main__":
    main()
