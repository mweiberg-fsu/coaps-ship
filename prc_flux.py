#!/usr/bin/env python3
"""
mft_sensible_flux.py
====================
Sensible heat flux (SHF) calculator based on BVW/MFT23 model.

Input file : input.csv  (hard-coded)
Output file: output.csv (hard-coded)

CSV format:
    Row 1: column names
    Row 2: units (ignored)
    Row 3+: data

Author: Grok (based on MFT23.py by Mark Bourassa)
"""

import pandas as pd
import numpy as np
import math
import os

# ===================================================================
# 1. HARD-CODED FILE NAMES
# ===================================================================
INPUT_FILE = "KAOU2011.csv"
OUTPUT_FILE = "output.csv"

# ===================================================================
# 2. PHYSICAL CONSTANTS
# ===================================================================
G = 9.81          # gravity                     [m/s^2]
KV = 0.40          # von Karman constant
R = 287.05        # gas constant dry air        [J/kg*K]
Prt = 1.0           # turbulent Prandtl number
cp0 = 1004.0        # specific heat dry air       [J/kg*K]

# ===================================================================
# 3. DEFAULT PARAMETERS
# ===================================================================
Z_REF_WIND = 10.0      # wind measurement height               [m]
Z_REF_TQ = 10.0      # temperature/humidity height           [m]
CONV_CRIT = 5e-5      # convergence criterion

# ===================================================================
# 4. HELPER FUNCTIONS
# ===================================================================


def sat_vap_press(Tc):
    """Saturation vapor pressure over water (Pa) – augmented Tetens."""
    return 611.2 * math.exp(17.67 * Tc / (Tc + 243.5))


def spec_hum_from_rh(RH_frac, Tc, P_pa):
    """Specific humidity [kg/kg] from RH (0–1), T (°C), P (Pa)."""
    es = sat_vap_press(Tc)
    e = RH_frac * es
    return 0.622 * e / (P_pa - 0.378 * e)


def air_density(P_pa, Tc, q):
    Tk = Tc + 273.15
    return P_pa / (R * Tk * (1.0 + 0.61 * q))


def moist_specific_heat(q_sfc):
    return cp0 * (1.0 + 0.9 * q_sfc)


def kinematic_viscosity_air(Tc):
    """Kinematic viscosity of air (m²/s)."""
    Tk = Tc + 273.15
    mu = 1.327e-5 * (Tk / 273.15)**1.81
    rho = air_density(101325, Tc, 0.0)
    return mu / rho


# ===================================================================
# 5. ROUGHNESS & STABILITY
# ===================================================================
B = 0.019
sfcten = 0.074     # surface tension [N/m]
denwat = 1025.0    # seawater density [kg/m³]


def z0_momentum_bvw(ustar, nu, wave_age=28.0):
    term1 = 0.11 * nu / ustar
    term2 = B * sfcten / (ustar**2 * denwat)
    term3 = 0.016 * ustar**2 / G
    return term1 + math.sqrt(term2**2 + term3**2)


def renewal_time_shear(ustar, z0):
    return 53.32 * math.sqrt(1.4e-5 * z0 / ustar**3)


def stanton_from_renewal(ustar, renewal):
    alpha_t = 2.1e-5
    return math.sqrt(alpha_t / (ustar**2 * renewal))


def z0_thermal_cfc(z0m, ustar, Ch):
    return z0m * math.exp(KV * (5.0 - 1.0 / (Prt * Ch)))


def psi_theta(z, zeta):
    if zeta < 0:  # unstable
        x = (1.0 - 9.0 * zeta)**0.5
        return 2.0 * math.log((1.0 + x) / 2.0)
    else:  # stable
        return -5.0 * zeta

# ===================================================================
# 6. CORE: Compute SHF for one row
# ===================================================================


def compute_shf_row(row):
    T_air = float(row["in_T"])
    RH = float(row["in_RH"]) / 100.0
    T_skin = float(row["in_TS"])
    U10 = float(row["in_SPD"])
    P_mb = float(row["in_P"])
    P_pa = P_mb * 100.0

    q_air = spec_hum_from_rh(RH, T_air, P_pa)
    q_sfc = spec_hum_from_rh(1.0, T_skin, P_pa)
    rho = air_density(P_pa, T_air, q_air)
    cp = moist_specific_heat(q_sfc)

    ustar = max(0.001 * U10, 1e-6)
    L_monin = 1e6
    wave_age = 28.0

    for _ in range(40):
        nu = kinematic_viscosity_air(T_air)
        z0m = z0_momentum_bvw(ustar, nu, wave_age)

        # Wind profile with stability
        zeta_u = Z_REF_WIND / L_monin
        psi_u = psi_theta(Z_REF_WIND, zeta_u)
        ln_u = math.log(Z_REF_WIND / z0m)
        U_eff = U10 + (ustar / KV) * psi_u
        if U_eff < 0.1:
            U_eff = 0.1

        ustar_new = KV * U_eff / ln_u
        ustar_new = max(ustar_new, 1e-6)

        # Thermal roughness
        renewal = renewal_time_shear(ustar_new, z0m)
        Ch = stanton_from_renewal(ustar_new, renewal)
        z0t = z0_thermal_cfc(z0m, ustar_new, Ch)
        zzt = Z_REF_TQ / z0t

        # T*
        psi_t = psi_theta(Z_REF_TQ, Z_REF_TQ / L_monin)
        T_star = (T_skin - T_air) / (math.log(zzt) - psi_t) * (KV / Prt)

        # Monin-Obukhov length
        Tv = (T_air + 273.15) * (1.0 + 0.61 * q_air)
        bflux = -G * T_star / Tv
        if abs(bflux) < 1e-12:
            L_new = 1e9 * np.sign(L_monin)
        else:
            L_new = -ustar_new**3 / (KV * bflux)
        L_new = np.clip(L_new, -1e8, 1e8)

        # Convergence
        if (abs(ustar_new - ustar) / ustar < CONV_CRIT and
                abs(L_new - L_monin) / abs(L_monin) < CONV_CRIT):
            ustar = ustar_new
            L_monin = L_new
            break

        ustar = 0.7 * ustar + 0.3 * ustar_new
        L_monin = 0.7 * L_monin + 0.3 * L_new

    shf = -rho * cp * ustar * T_star
    return round(shf, 3)

# ===================================================================
# 7. MAIN: Read → Compute → Write
# ===================================================================


def main():
    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: Input file '{INPUT_FILE}' not found.")
        return

    # FIX: skip first 2 rows (header + units), use row 3 as header
    df = pd.read_csv(INPUT_FILE, skiprows=2, header=0)

    required = {"in_T", "in_RH", "in_TS", "in_SPD", "in_P"}
    missing = required - set(df.columns)
    if missing:
        print(f"ERROR: Missing columns in CSV: {missing}")
        return

    print(f"Loaded {len(df)} data rows from '{INPUT_FILE}'")
    df["SHF_Wm2"] = df.apply(compute_shf_row, axis=1)

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Results saved to '{OUTPUT_FILE}'")
    print("\nFirst 5 rows:")
    print(df[["time", "in_T", "in_RH", "in_TS", "in_SPD", "in_P", "SHF_Wm2"]].head())


if __name__ == "__main__":
    main()
