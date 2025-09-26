# config.py

# List of available ships
# WTEA", "KAOU", "KAQP", "KCEJ", "KNBD", "KTDQ", "NEPP", "NRUO", "VLHJ", "VLMJ", "VMIC", "VNAA", "WARL", "WBP3210", "WCX7445", "WDA7827", "WDC9417",
# "WDD6114", "WDG7520", "WDN7246", "WECB", "WKWB", "WSAF", "WSQ2674", "WTDF", "WTDH", "WTDK", "WTDL", "WTDM", "WTDO", "WTEB", "WTEC", "WTED", "WTEE",
# "WTEF", "WTEG", "WTEJ", "WTEK", "WTEO", "WTEP", "WTER", "WTEU", "WTEY", "WXAQ", "ZCYL5", "ZGOJ7", "ZMFR"


# Define ships to process
# ships = ["WTEA", "KAOU", "KAQP", "KCEJ", "KNBD", "KTDQ", "NEPP", "NRUO", "VLHJ", "VLMJ", "VMIC", "VNAA", "WARL", "WBP3210", "WCX7445", "WDA7827", "WDC9417",
# "WDD6114", "WDG7520", "WDN7246", "WECB", "WKWB", "WSAF", "WSQ2674", "WTDF", "WTDH", "WTDK", "WTDL", "WTDM", "WTDO", "WTEB", "WTEC", "WTED", "WTEE",
# "WTEF", "WTEG", "WTEJ", "WTEK", "WTEO", "WTEP", "WTER", "WTEU", "WTEY", "WXAQ", "ZCYL5", "ZGOJ7", "ZMFR"]

ships = ["WTEA"]  

# Define years from 2005 to 2025
years = []
for i in range(2005, 2026):
    years.append(i)

proc_num = 4 # Number of parallel threads for downloading     

# Define directory paths
directory_destination = 'data/input/'
variable_name = 'in_T'
logs_dir = 'logs'
     