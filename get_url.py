csv_url = (
    f'https://erddap-samos.coaps.fsu.edu/erddap/tabledap/SAMOS_Fluxes_B23_v301.csv?'
    f'time%2Cplatform_call_sign%2Clatitude%2Clongitude%2Cin_T%2Cin_RH%2Cin_TS%2Cin_SPD%2Cin_P'
    f'&time%3E=2011-01-01'
    f'&time%3C=2011-12-31T23%3A59%3A00Z'
    f'&platform_call_sign=%22KAOU%22'
    f'&latitude%3E=-78.64944&latitude%3C=89.99979'
    f'&longitude%3E=0&longitude%3C=359.9999'
    f'&in_T_qc=1&in_RH_qc=1&in_TS_qc=1&in_SPD_qc=1&in_P_qc=1'
)

print(csv_url)
