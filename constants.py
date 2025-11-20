# constants.py
# Description: This file contains constants necessary to run calculations

# General constants for flux calculations
CONV_CRIT = 0.005               # Critical value for convection	
convect = 0.0                   # Convection value		
warn = 1   		                # Warning flag
eqv_neut = 0		            # Height adjustment is not modified, normal winds. 
z_wanted = 10.0                 # Height to which wind is adjusted
ref_ht_wind = 17.9	            # Height given by SAMOS for wind
ref_ht_tq = 17.9                # Height given by SAMOS for air temp and humidity     
net_heat_flux = 5.0	            # Net heat flux (W/m^2)	
sst_prm = 1                     # SST parameter)
z0_mom_prm = 6		            # Momentum roughness parameter
z0_theta_q_prm = 2              # Temperature and humidity roughness parameter
stable_prm = 0                  # Stability parameter				
wave_ang = 0                    # Wave angle		
dyn_in_prm = 0		            # Wind speed
dyn_in_val2 = 0.0               # Wind direction
astab = 1                       # Atmospheric stability		
ss_prm = 2 			            # Wind-wave Stability Param
ss_val = 48.0                   # Wind-wave Stability Value     
air_moist_prm = 1	            # Relative Humidity in Air
sfc_moist_prm = 1	            # Relative Humidity at Surface
sfc_moist_val = 0.98            # Relative Humidity at Surface Value
salinity = 0.0349	            # Average salinity
sfc_current1 = 0.0              # Surface current u-component
sfc_current2 = 0.0              # Surface current v-component
oil_fract_area = 0.0            # Oil fraction of area
missing = -9999.0               # NaN value
dimensionless_m_o_length = 0.0  # Dimensionless momentum roughness length
zo_m = 0.0                      # Momentum roughness length

# Program specific constants
n = 500                         # Select number of permutations to calculate
chunk_size = 1000               # Chunk size for processing data in pieces
window_size = 11                # Window size for rolling calculations
flux_model = -1                 # Flux model selection (Blair 23, use 1 for S88)