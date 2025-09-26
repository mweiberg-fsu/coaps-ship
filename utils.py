# utils.py
from imports import *

## Create directory if it doesn't exist
def create_directory(path):
    os.makedirs(path, exist_ok=True)

## Initialize and configure logger
def setup_logger(name, log_file, level=logging.INFO):
    create_directory(os.path.dirname(log_file)) # Ensure log directory exists
    handler = logging.FileHandler(log_file) # Append mode
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s') # Standard format
    handler.setFormatter(formatter) # Set formatter for handler
    logger = logging.getLogger(name) # Create logger
    logger.setLevel(level) # Set logging level
    logger.addHandler(handler) # Add handler to logger
    return logger # Return the configured logger

# Function to clear log file contents on startup
def clear_log_file(log_file):
    try:
        if os.path.exists(log_file): # Check if log file exists
            os.remove(log_file) # Remove the existing log file on run
            print(f"Cleared log file: {log_file}")
        else:
            print(f"No log file found at: {log_file}")
    except Exception as e:
        print(f"Error clearing log file {log_file}: {e}")

# Function to check if a specific variable exists in the CSV file header
def check_variable(file_path, variable_name):
    try:
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            header = next(reader)  # Get the header row
            return variable_name in header
    except FileNotFoundError:
        print(f"Error: File not found at path: {file_path}")
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False