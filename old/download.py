from urllib.request import urlretrieve
import requests
import time
import csv
import pandas as pd
import io 
import shutil
import os


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



##################################################################################################################################################################################################################

# Select ships from this list
# ships = ["WTEA", "KAOU", "KAQP", "KCEJ", "KNBD", "KTDQ", "NEPP", "NRUO", "VLHJ", "VLMJ", "VMIC", "VNAA", "WARL", "WBP3210", "WCX7445", "WDA7827", "WDC9417",
#           "WDD6114", "WDG7520", "WDN7246", "WECB", "WKWB", "WSAF", "WSQ2674", "WTDF", "WTDH", "WTDK", "WTDL", "WTDM", "WTDO", "WTEB", "WTEC", "WTED", "WTEE",
#            "WTEF", "WTEG", "WTEJ", "WTEK", "WTEO", "WTEP", "WTER", "WTEU", "WTEY", "WXAQ", "ZCYL5", "ZGOJ7", "ZMFR"]


ships = ["WTEA"]
#Select years from this range
years = [2025]
for i in range(2005,2025):
    years.append(i)

# Select csv file destination (Written now to have separate directories for each ship)
# Delete the first {ship} directory in fileName and file_path to not download to individual ship directories
directory_destination = 'data/input/'

#############################################################################################################################################################################################




################################################################################################################
################################################ Download CSV Files ############################################
################################################################################################################

total = len(ships) * len(years)
count = 1
for s in range(len(ships)):
    ship = ships[s]
    
    # Create ship-specific directory if it doesn't exist
    ship_directory = f'{directory_destination}/{ship}'
    os.makedirs(ship_directory, exist_ok=True)  # Creates directory if it doesn't exist
    
    for year in years:
        print(f'{count}/{total}')
        csv_url = f'https://erddap-samos.coaps.fsu.edu/erddap/tabledap/SAMOS_Fluxes_B23_v300.csv?time%2Cplatform_call_sign%2Clatitude%2Clongitude%2Cin_T%2Cin_RH%2Cin_TS%2Cin_SPD%2Cin_P&time%3E={year}-01-01&time%3C={year}-12-31T23%3A59%3A00Z&platform_call_sign=%22{ship}%22&latitude%3E=-78.64944&latitude%3C=89.99979&longitude%3E=0&longitude%3C=359.9999&in_T_qc=1&in_RH_qc=1&in_TS_qc=1&in_SPD_qc=1&in_P_qc=1'
        
        fileName = f'{ship_directory}/{ship}{year}.csv'
        req = requests.get(csv_url)
        url_content = req.content
        csv_file = open(fileName, 'wb')
        csv_file.write(url_content)
        count += 1
        csv_file.close()


################################################################################################################
################################################ Delete Empty Files ############################################
################################################################################################################

variable_name = 'in_T'
for ship in ships:
    for year in years:
        file_path = f'{directory_destination}/{ship}/{ship}{year}.csv'
        if check_variable(file_path,variable_name):
            print(f'{file_path} is good')
        else:
            os.remove(file_path)
