from imports import *
from data_get import *


def main():
    logger.info("Starting the data retrieval process")
    
    start_time = time.time()  # ⏱️ Start timer
    
    run_parallel_downloads(max_threads=4)
    
    end_time = time.time()    # ⏱️ End timer
    elapsed_time = end_time - start_time
    minutes, seconds = divmod(elapsed_time, 60)

    logger.info(f"Data retrieval and cleaning completed in {int(minutes)} min {seconds:.2f} sec")
    logger.debug("No data processing implemented yet in data_prc.py")

if __name__ == "__main__":
    main()
