# MFT23_numba.py
import numpy as np
from numba import jit, prange
import math as m

# ------------------------------------------------------------------
# 1. Helper functions – copy the *pure-math* part from MFT23.py
# ------------------------------------------------------------------


@jit(nopython=True)
def bstar_f(neutral, qmixa, t_air, tstar, qstar, min_val):
    # original formula – replace with the real one from MFT23.py
    return (9.81 / t_air) * (tstar + 0.61 * t_air * qstar)


@jit(nopython=True)
def psi_u_f(zref_minus_dspl, zzu, monin_inv, stab_prm):
    # very short version of the stability function
    if monin_inv == 0.0:
        return 0.0
    x = (1.0 - 16.0 * zref_minus_dspl / monin_inv) ** 0.25
    if monin_inv > 0:
        return -2.0 * np.log((1.0 + x) / 2.0)
    else:
        return 2.0 * np.log((1.0 + x) / 2.0) + np.log((1.0 + x * x) / 2.0) - 2.0 * np.arctan(x) + np.pi / 2.0


@jit(nopython=True)
def ustar_f(fixed_uen, neutral, bstar, ustar, convect_adj, psi_u, monin_inv, wind_vect_i, zzu, ten):
    # simple COARE-like iteration step
    return (wind_vect_i * 0.4) / (np.log(ten / zzu) - psi_u)


@jit(nopython=True)
def z0_and_waves(...):   # <-- **copy the whole body from MFT23.py** (replace math → np)
    # For brevity I keep it as a stub – you only need the math, no logging.
    # Return the 12-tuple exactly as the original function.
    return trouble, ustar, zzu, dspl, hsig, wave_len, orbital_vel, wave_age, dom_wave_phs_spd, betac, betag, betas

# ------------------------------------------------------------------
# 2. Vectorised find_ustar (the only function you call from outside)
# ------------------------------------------------------------------


@jit(nopython=True, parallel=True, cache=True)
def find_ustar_numba(
    wind_vect,          # (N,2)  – wind speed (u) and 0 (v)
    qmixa, t_air, tstar, qstar, ss_val, oil_fract_area,
    fixed_hsig, neutral, convect_adj, zref, ss_prm, z0_prm, stab_prm,
    betac_prime, betag_prime, denwat, nu, sfcten, ww, ww_eql,
    fixed_dom_wave_phs_spd, fixedwa, no_capw, no_sfcten, use_dh, use_orb_vel
):
    N = wind_vect.shape[0]
    dsplcmnt_ht = np.zeros(N)
    monin_inv = np.zeros(N)
    ustar_out = np.zeros((N, 2))
    zzu_out = np.zeros(N)
    betac_out = np.zeros(N)
    betag_out = np.zeros(N)
    betas_out = np.zeros(N)

    for k in prange(N):
        # ---- initialise per-sample variables -----------------------
        ustar = np.array([0.0, 0.0])
        ustar_old = np.array([0.0, 0.0])
        ucount = 0
        done = False
        betac = betag = betas = dspl = monin = 0.0
        hsig = period = wave_len = orbital_vel = wave_age = dom_wave_phs_spd = zzu = 0.0

        while not done and ucount < 30:
            ucount += 1
            ustar_old[:] = ustar[:]

            # ----- buoyancy flux ----------------------------------------
            bstar = bstar_f(neutral, qmixa[k],
                            t_air[k], tstar[k], qstar[k], 1e-6)

            for i in range(2):
                if m.fabs(wind_vect[k, i]) < 0.001:
                    ustar[i] = 0.0
                    continue

                trouble, ustar, zzu, dspl, hsig, wave_len, orbital_vel, wave_age, \
                    dom_wave_phs_spd, betac, betag, betas = z0_and_waves(
                        fixed_dom_wave_phs_spd, fixedwa, neutral,
                        fixed_hsig and i == 0, stab_prm, no_capw, no_sfcten,
                        use_dh, use_orb_vel, betac_prime, betag_prime,
                        betac, betag, betas, denwat, nu, sfcten, zref,
                        dspl, monin, wave_age, dom_wave_phs_spd, hsig,
                        period, wave_len, orbital_vel, ww, ww_eql,
                        ss_prm, ss_val[k], i, z0_prm, zzu, ustar, oil_fract_area[k]
                    )

                psi_u = psi_u_f(zref - dspl, zzu, monin, stab_prm)
                ustar[i] = ustar_f(False, neutral, bstar, ustar, convect_adj,
                                   psi_u, monin, wind_vect[k, i], zzu, 10.0)

            # ----- convergence test ------------------------------------
            if ustar[0] <= 1e-6:
                done = m.fabs(ustar[0] - ustar_old[0]) < 1e-4
            else:
                done = m.fabs((ustar[0] - ustar_old[0]) / ustar[0]) < 0.01

            if ustar[1] <= 1e-6:
                done = done and m.fabs(ustar[1] - ustar_old[1]) < 1e-4
            else:
                done = done and m.fabs(
                    (ustar[1] - ustar_old[1]) / ustar[1]) < 0.01

        dsplcmnt_ht[k] = dspl
        monin_inv[k] = monin
        ustar_out[k] = ustar
        zzu_out[k] = zzu
        betac_out[k] = betac
        betag_out[k] = betag
        betas_out[k] = betas

    return dsplcmnt_ht, monin_inv, ustar_out, zzu_out, betac_out, betag_out, betas_out
