"""
 'Python' library for combined sea state and atmospheric boundary layer model

  Developed by Mark A. Bourassa

  Based on Original Bourassa-Vincent-Wood combined sea state and atmospheric
  boundary layer model developed at Purdue University.  Several upgrades and
  options have been added while at the Center for Ocean-Atmospheric
  Prediction Studies (COAPS) at Florida State University.

  Dr. Bourassa is currently (May 2023) a Professor at the Florida
  State University's Department of Earth, Ocean and Atmospheric Science, and is
  the Associate Director as the Center for Ocean-Atmospheric Prediction Studies.
  Please send comments, suggestions, and requests for added features to
      Dr. Mark Bourassa
      Center for Ocean-Atmospheric Prediction Studies
      Florida State University
      Tallahassee, FL 32306-2840

      Email: bourassa@coaps.fsu.edu
      Phone: (850) 644-6923

  The bulk of 'in code' documentation is in the subroutine pmix_. Documentation
  is also on the COAPS website: https://www.coaps.fsu.edu/~bourassa/MFT_code/MFT.py

"""
import math as m
import logging

# import numpy as np

# TRUE and FALSE
TRUE = 1
FALSE = 0

# # # Constants # # #
B = 0.0190000  # Wu's equilibrium const for capillary waves  []
KV = 0.40000  # von Karman's constant []
G = 9.81000  # gravitational acc. (at sea level)  [m/s^2]
H = 1000.0  # approximate height of boundary layer (at low wind speed)  [m]
PI2 = 6.28318  # two PI  []
R = 287.05000  # specific gas constant for dry air  [J/kg]
RW = 441.0000  # specific gas constant for water vapour  [J/kg]
DTR = 3.14159 / 180.0  # degrees to radians conversion factor

stefanBoltzmann = 5.67e-08
Prt = 1.00  # turbulent Prandtl number []
Prair = 0.7
PrWater = 7.0
Sc = 1.00  # turbulent Schmidt number []
thermalDiff = 2e-05  # or air from Gill p. 71 [m^2/s]
molecularDiff = 2.4e-05  # for H2O in air from Gill p. 68 [m^2/s]

# Next value determines water type for solar absorption
# 1 = Type I, 1.1 = Type IA,  1.2 = Type IB, 2.0 = Type II,
# 3.0 = Type III, 0 = no solar absorption
waterType = 1

logger = logging.getLogger("MFT_log")


class MFTError(Exception):
    ...


class InvalidFluxModelParameterError(MFTError):
    ...


class UnrealisticGeophysicalValueError(MFTError):
    ...


def pmix_(
    dyn_in_prm,
    dyn_in_val,
    dyn_in_val2,
    sfc_current,
    convect,
    conv_crit,
    pressure,
    air_moist_prm,
    air_moist_val_dum,
    sfc_moist_prm,
    sfc_moist_val_dum,
    salinity,
    ss_prm,
    ss_val_dum,
    t_air,
    sst_prm,
    t_skin,
    ref_ht_wind,
    ref_ht_tq,
    astab,
    net_long_wave_flux,
    flux_model,
    z0_mom_prm,
    z0_theta_q_prm,
    stab_prm,
    zo_m,
    dimensionless_m_o_length,
    oil_fract_area,
):
    """
    'main' subroutine for combined sea state and atmospheric flux model

    the original BVW model Included the following, but there are now many more options.
       BVW sea state parameterization,
       BVW atmospheric boundary layer (roughness length) model (optional),
       Toba's 2/3 wave relation,
       Monin-Obhukov atmospheric stability parameterization,
       Godfrey and Beljaar's convective parameterization.

    Options for dynamic input parameters
       Most users will want to input wind speeds.  However, there is a
       growing need to convert scatterometer observations of stress or
       equivalent neutral wind speed to winds and/or stresses.
       dyn_in_prm    Parameter treated as known
       0             Wind speed,
       1             Friction velocity,
       2             Surface stress,
       3             Equivalent neutral wind speed. Note 'equivalent neutral'
                     is NOT equivalent to 'neutral'.

    dyn_in_val and dyn_in_val2 are the input values matching the physical parameter
      selected via input through the dyn_in_prm. This can be input as a speed in
      dyn_in_val (with 0.0 input for dyn_in_val2) or as vector components.

    Options for seastate parameterizations (or parameters):
       There are three possible seastate assumptions: any one of the
       following can be treated as known: wind-wave stability parameter,
       phase speed, or wave age.
       ss_prm     Parameter treated as known
       0          Wind-wave stability parameter (set to 1.0 for local
                  equilibrium),
       1          Phase speed of the dominant waves  [m/s],
       2          Wave age the dominant waves (dom_wave_phs_spd/ustar)  [],
       3          Significant Wave Height  [m],
       4          Significant slope  [],
       5          Period of the dominant waves [s],
       6          Orbital velocity of dominant waves [m/s].

    Options for atmospheric moisture input:
       Choose the moisture parameter that is easiest for you to deal with:
       air_moist_prm  Parameter for moisture of air
       0          Absolute humidity at the reference height of the thermometer
                  and humidity sensor [g vapor / g air],
       1          Relative Humidity [fraction],
       2          Dew point temperature [C],
       3          Wet bulb temperature [C].

    Options for surface moisture input:
       Choose the moisture parameter that is easiest for you to deal with:
       sfc_moist_prm  Parameter for moisture of air
       0          Absolute humidity 'at' (near) the surface [g vapor / g air],
       1          Relative Humidity [fraction],
       2          Dew point temperature [C],
       3          Wet bulb temperature [C].


    Definitions:
    Input parameters  (passed through call):
        CON_P       Convective parameter (unitless).  Recommended value between
                    0.7 and 1.25.  For details see TOGA NOTES #4.
        crit        Convergence criterion (unitless fraction).
        air_moist_prm   Moisture parameter index.
        air_moist_val   Value of the parameter corresponding to the above index.
        press       Atmospheric surface pressure (Pa).
        s           Salinity (unitless).
        ss_prm      Seastate parameter index.
        ss_val      Value of the parameter corresponding to the above index.
        t_air         Air temperature at the reference height of the thermometer
                    and humidity sensor (degrees C).
        u           Mean wind speed (m/s) at the height (zref) of the anemometer.
        warn        Normally this should be set to zero (0). If it is set to one (1),
                    then in cases where the code fails to converge for very small fluxes,
                    the output will be set to missing. This might be useful for cal/val
                    of flux algorithms, as this would remove output that is not a suitably
                    close match to the parameterizations. However, normally assuming a
                    flux of zero is a small error.
        zref        Height (metres) of the wind observations.
        zrefq       Height (metres) of the humidity observations.  See zreft.
        zreft       Height (metres) of the temperature observations.  Note: in
                    the current version of the code this must equal to height
                    of the humidity observations.
        astab       Option for atmospheric stability condition:
                        0) atmospheric stability is assumed to be neutral,
                        1) assumption of a neutral atmospheric stability is
                           not made.
        oil_fract_area       fraction of surface covered by oil. Normally the input value for
                    this parameter is zero.  Larger value will suppress capillary waves
                    and mimic the suppression of shorter gravity waves.

    Other parameters:
        bstar       Component of buoyancy flux (-bstar * ustar).
        cmin        Minimum phase speed of water waves - determined through
                    the phase relation   [m/s].
        air_specific_heat          Specific heat of air [J/kg].
        denair      density of air  [kg/m^3].
        denwat      density of water  [kg/m^3].
        latent_heat_vaporization          Latent heat of vaporization  [J/kg].
        sfcten      surface tension of air-water interface [N/m].
        ww          Wind-wave stability parameter (unitless).  Equal to one for
                    local equilibrium; greater than one for rising seas; less
                    than one (and greater than zero) for falling seas.
        ww_eql      Value of U(Hs)/dom_wave_phs_spd for local wind-wave equilibrium (0.83).
        zzu          wind reference height divided by momentum roughness length [],
        zzq         moisture reference height divided by moisture roughness
                    length  [],
        zzt         temperature reference height divided by temperature roughness
                    length  [],

    Output parameters:
        dom_wave_phs_spd          dominant phase speed of gravity waves   [m/s],
        hsig        significant wave height  [m],
        inv_monin   inverse Monin-Obhukov scale length  [m^-1],
        lhf         latent heat flux  [W/m^2],
        qstar       scaling parameter for moisture  [],
        shf         sensible heat flux  [W/m^2],[W/m^2],
        tau         stress vector   [N/m^2];  tau[0] is parallel the direction
                    of wave propagation, and tau[1] is perpendicular tau[0] in
                    a right-handed coordinate system (with the positive
                    vertical axis pointing upward,
        tstar       scaling term for potential temperature  [degrees C]
        ustar       friction velocity  [m/s],
        wave_age    wave age,  dom_wave_phs_spd/ustar  [],
        zo_m        momentum roughness length (two components: parallel and
                    perpendicular direction of wave propagation) [],
                    :rtype: object
    """


    # global ustar

    no_capw = FALSE
    use_dh = FALSE
    use_orb_vel = FALSE
    no_sfcten = FALSE

    wind_vect = [0.0, 0.0]
    orbital_vel = [0.0, 0.0]
    ustar = [0.0, 0.0]
    tau = [0.0, 0.0]

    if flux_model < 0:
        z0_prm = z0_mom_prm
        theta_q_zo_prm = z0_theta_q_prm
    else:
        theta_q_zo_prm = 0
        stab_prm = 0
        if flux_model == 0:  # BVW
            z0_prm = 0
        elif flux_model == 1:  # Smith 1988
            z0_prm = 0
            ss_prm = 2
            ss_val_dum = 43.64
            no_capw = TRUE
            stab_prm = 1
            sst_prm = 0
        elif flux_model == 2:
            # BVW without interaction between waves and atmospheric stab
            z0_prm = 0
        elif flux_model == 3:
            # BVW with Smith (1988) stability parameterization
            z0_prm = 0
            stab_prm = 1
        elif flux_model == 4:
            # BVW without capillary waves
            z0_prm = 0
            no_capw = TRUE
        elif flux_model == 5:
            # BVW without sfcten in phase relations
            z0_prm = 0
            no_sfcten = TRUE
        elif flux_model == 6:
            # BVW without capillary waves and without sfcten in the phase relation
            z0_prm = 0
            no_sfcten = TRUE
            no_capw = TRUE
        elif flux_model == 7:
            # Taylor and Yelland (2001) roughness length
            z0_prm = 2
            no_sfcten = TRUE
            no_capw = TRUE
        elif flux_model == 8:
            # Taylor and Yelland (2001) roughness length + capillary
            z0_prm = 3
        elif flux_model == 9:
            # Bourassa 2006 friction velocity with CFC zot and zoq
            z0_prm = 1
            theta_q_zo_prm = 1
        elif flux_model == 10:
            # case 9 with displacement height estimated within the code
            z0_prm = 1
            theta_q_zo_prm = 1
            use_dh = TRUE
        elif flux_model == 11:
            # Bourassa 2006 friction velocity with Zilitinkevich et al. zot and zoq
            z0_prm = 1
            theta_q_zo_prm = 2
        elif flux_model == 12:
            # Bourassa 2006 friction velocity with LKB zot and zoq
            z0_prm = 1
            theta_q_zo_prm = 3
        elif flux_model == 13:
            # Bourassa 2006 friction velocity with COARE3.0 zot and zoq
            z0_prm = 1
            theta_q_zo_prm = 4
        elif flux_model == 14:
            # Bourassa 2006 friction velocity with wall theory for moisture
            z0_prm = 1
            theta_q_zo_prm = 0
        elif flux_model == 15:
            # Bourassa 2006 friction velocity with CFC zot and zoq and oil spill z0m and zero LHF
            z0_prm = 1
            theta_q_zo_prm = 1
        else:
            raise InvalidFluxModelParameterError("Invalid choice of flux_model parameter: " + str(flux_model))

    if ref_ht_wind <= 0:
        raise UnrealisticGeophysicalValueError("Unrealistic value of reference height for wind: " + str(ref_ht_wind))

    if ref_ht_tq <= 0:
        raise UnrealisticGeophysicalValueError(
            "Unrealistic value of reference height for temperature and humidity: " + str(ref_ht_tq)
        )

    molec_mass_salt = 58.443000  # molecular mass of salt  [kg/kmol]
    molec_mass_water = 18.01600  # molecular mass of water [kg/kmol]
    betag_prime = 1.0
    betac_prime = 1.0
    if no_capw:
        betac_prime = 0.0
    betas_prime = 0.0 + float(no_capw)
    # betac = [betac_prime, betac_prime]
    betag = [betag_prime, betag_prime]
    # betas = [betas_prime, betas_prime]
    ww_eql = 0.345 / KV

    crit = float(conv_crit)
    convect_adj = float(convect)
    press = float(pressure)
    air_moist_val = float(air_moist_val_dum)
    sfc_moist_val = float(sfc_moist_val_dum)
    s = float(salinity)
    ss_val = float(ss_val_dum)
    t_air = float(t_air)
    tskin = float(t_skin)
    zref = float(ref_ht_wind)
    # /* the height of the temperature sensor must be equal to the height of
    # the humidity sensor.  */
    zrefq = float(ref_ht_tq)
    zreft = zrefq

    if z0_prm != 7:
        zzu = 10000.0
    else:
        zzu = zreft / zo_m  # this should be zreft minus the displacement height

    # set up the constants and initial guesses for the seastate.
    wave_age = [28.0, 28.0]  # initial guess
    hsig = [0.01, 0.01]  # initial guess assumes waves exist
    wave_len = 10.0  # initial guess assumes waves exist
    dom_wave_phs_spd = [4.0, 4.0]  # initial guess assuming wave exist.

    fixed_dom_wave_phs_spd = FALSE
    fixedwa = FALSE
    fixed_hsig = FALSE

    if ss_prm == 0:
        # the wind-wave stability parameter is specified
        pass
    elif ss_prm == 1:
        # phase speed is known
        fixed_dom_wave_phs_spd = TRUE
        dom_wave_phs_spd[0] = ss_val
    elif ss_prm == 2:
        # wave age (ustar/dom_wave_phs_spd) is known
        fixedwa = TRUE
        wave_age[0] = ss_val
    elif ss_prm == 3:
        # significant wave height is known
        fixed_hsig = TRUE
        hsig[0] = ss_val
    elif ss_prm == 4:
        # significant slope is known
        pass
    elif ss_prm == 5:
        # period of the dominant waves is known
        pass
    elif ss_prm == 6:
        # Orbital velocity is known
        orbital_vel[0] = ss_val
    else:
        raise InvalidFluxModelParameterError("Invalid choice of ss_prm: " + f"{ss_prm}")

    # set up the atmospheric stability parameters
    monin_inv = 0.0
    if astab == 0:  # treat the stratification as being neutral (even if that isn't realistic)
        neutral = TRUE
        monin_inv = 0.0
    elif astab == 1:  # treat the atmosphere as being non-neutral
        neutral = FALSE
    elif astab == 2:  # use an input value of z/l
        neutral = FALSE
        monin_inv = dimensionless_m_o_length / zref
    else:
        raise InvalidFluxModelParameterError(
            "Invalid atmospheric stability option: "
            + f"{astab}"
            + " allowable options are: 0 (neutral), 1 (non-neutral), and 2 (input z/L)"
        )

    betag_set = betag_prime
    betac_set = betac_prime
    betas_set = betas_prime

    # determine the surface tension
    sfcten = 7.6100000e-2 - 1.55000e-4 * tskin + 2.77000e-5 * s * molec_mass_water / molec_mass_water

    # determine the latent heat of vaporization
    latent_heat_vaporization = 4186.8 * (597.31 - 0.56525 * tskin)

    # use the atmospheric moisture parameter to determine the specific humidity [kg/kg].
    qmixa = find_q_(air_moist_prm, air_moist_val, press, t_air)
    if qmixa == -1.0:
        raise InvalidFluxModelParameterError("Invalid choice of air_moist_prm: " + f"{air_moist_prm}")

    # use the surface moisture parameter to determine the specific humidity [kg/kg].
    qmixw = find_q_(sfc_moist_prm, sfc_moist_val, press, tskin)
    if qmixw == -1.0:
        raise InvalidFluxModelParameterError("Invalid choice of sfc_moist_prm: " + f"{sfc_moist_prm}")
    if z0_prm == 4:
        qmixw = qmixa

    # determine the heat capacity at constant pressure
    air_specific_heat = 1004.0 * (1.0 + 0.9 * qmixw)

    # determine the density of the moist air
    denair = press / (R * (t_air + 273.14) * (1.00000000 + 0.6100000 * qmixa))

    # determine the density of pure water
    denwat = (999.8396 + 18.224944 * tskin - 0.007922210 * tskin * tskin) / (1.0000000 + 0.018159725 * tskin)

    # determine the density of saline water
    va = (
        (12.97000 + 0.23400000 * tskin - 4.2100000e-3 * tskin * tskin + 2.8570000e-5 * tskin * tskin * tskin)
        * 10.00000e-3
        + m.sqrt(s * denwat * 10.00000e-3) * 2.9820000e-3
        - 4.9700000e-5 * tskin
        + 6.0320000e-7 * tskin * tskin
    )
    denwat = denwat * (1.0000000 + s) / (1.000000 + va * denwat * s / molec_mass_salt)

    # determine the viscosity of air
    nu = (
        1.3260000e-5 * (1.00000 + 0.006542 * t_air + 8.301000e-6 * t_air * t_air + 4.840000e-9 * t_air * t_air * t_air)
    ) / denair

    # setup known inputs for the dynamic input (typically wind speed)
    fixed_ustar = FALSE
    fixed_uen = FALSE

    if dyn_in_prm == 0:  # wind speed is specified
        # split the wind into component parallel and perpendicular to the mean
        #      direction of wave propagation.
        wind_vect[0] = float(dyn_in_val)
        wind_vect[1] = float(dyn_in_val2)
    elif dyn_in_prm == 1:  # friction velocity is specified
        fixed_ustar = TRUE
        # split the wind into component parallel and perpendicular to the mean
        #      direction of wave propagation.
        ustar[0] = float(dyn_in_val)
        ustar[1] = float(dyn_in_val2)
        if m.fabs(ustar[0]) < 0.00001:
            ustar[0] = 0.00001
        if m.fabs(ustar[1]) < 0.00001:
            ustar[1] = 0.00001
    elif dyn_in_prm == 2:  # surface stress (magnitude) is known
        fixed_ustar = TRUE
        # split the wind into component parallel and perpendicular to the mean
        #     direction of wave propagation.
        ustar[0] = m.sqrt(float(dyn_in_val) / denair)
        ustar[1] = m.sqrt(float(dyn_in_val2) / denair)
        if m.fabs(ustar[0]) < 0.00001:
            ustar[0] = 0.00001
        if m.fabs(ustar[1]) < 0.00001:
            ustar[1] = 0.00001
    elif dyn_in_prm == 3:  # equivalent neutral wind speed is known
        fixed_uen = TRUE
        # split the wind into component parallel and perpendicular to the mean
        #     direction of wave propagation.
        wind_vect[0] = float(dyn_in_val)
        wind_vect[1] = float(dyn_in_val2)
    else:
        raise InvalidFluxModelParameterError("Invalid choice of dyn_in_prm: " + f"{dyn_in_prm}")

    # One the first pass assume than waves exist
    hsig[1] = 0.01

    if z0_prm != 7:
        zzu = 10000.0

    (
        count,
        zzu,
        zzt,
        zzq,
        ustar,
        tstar,
        qstar,
        dsplcmnt_ht,
        monin_inv,
        hsig,
        wave_len,
        orbital_vel,
        wave_age,
        dom_wave_phs_spd,
        orbital_vel,
        sst_bulk_adjustment,
    ) = solve(
        fixed_dom_wave_phs_spd,
        fixedwa,
        fixed_hsig,
        fixed_uen,
        fixed_ustar,
        neutral,
        use_dh,
        use_orb_vel,
        convect_adj,
        wave_len,
        crit,
        denwat,
        nu,
        qmixa,
        qmixw,
        sfcten,
        ss_prm,
        ss_val,
        wind_vect,
        sfc_current,
        t_air,
        sst_prm,
        tskin,
        net_long_wave_flux,
        zref,
        zreft,
        zrefq,
        ww_eql,
        betag_set,
        betac_set,
        betas_set,
        z0_prm,
        oil_fract_area,
        astab,
        denair,
        air_specific_heat,
        latent_heat_vaporization,
        theta_q_zo_prm,
        stab_prm,
        monin_inv,
        no_capw,
        no_sfcten,
        zzu,
        ustar,
    )

    if (
        count <= 1
        or ustar[0] == 10.0
        or ustar[1] == 10.0
        or (wave_age[0] == 1.08 and betag[0] > 0.5 and (flux_model == 7 or flux_model == 8))
        or (wave_age[0] == 250.0 and betag[0] > 0.5 and (flux_model == 7 or flux_model == 8))
        or m.fabs(monin_inv) == 10.0
    ):
        # count = -1
        logger.warning(
            f"non-conv diags: {count} {ustar[0]} {ustar[1]} {tstar} {qstar} {wave_age[0]} {monin_inv} {betag[0]}"
        )
        if ustar[0] == 10.0 or ustar[1] == 10.0:
            ustar[0] = 0.0  # usually a good approximation
            ustar[1] = 0.0
            logger.warning("friction velocity set to zero.")

    ustar_mag = m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1])
    tau[0] = denair * ustar_mag * ustar[0]
    tau[1] = denair * ustar_mag * ustar[1]
    shf = float(-denair * air_specific_heat * ustar_mag * tstar)
    lhf = (1.0 - oil_fract_area) * float(-denair * latent_heat_vaporization * ustar_mag * qstar)
    dom_wave_phs_spd = float(betag[0] * dom_wave_phs_spd[0])
    hsig = float(betag[0] * hsig[0])

    return (
        count,
        tau,
        lhf,
        shf,
        ustar,
        tstar,
        qstar,
        zzu,
        zzt,
        zzq,
        monin_inv,
        dsplcmnt_ht,
        wave_age,
        dom_wave_phs_spd,
        hsig,
        orbital_vel,
        sst_bulk_adjustment,
        stab_prm,
    )
    # end of routine pmix_()


def find_q_(moist_prm, moist_val, press, temperature):
    # use the atmospheric moisture parameter to determine the specific humidity [kg/kg].
    if moist_prm == 0:  # input is specific humidity
        spec_hum = moist_val
    elif moist_prm == 1:  # relative humidity (fraction) is used as input
        # determine saturation vapour pressure over a smooth surface of water, at the temperature
        satvp = (1.000700 + 3.460000e-8 * press) * 611.2100 * m.exp(17.50200 * temperature / (240.97000 + temperature))
        spec_hum = moist_val * satvp * 0.6220 / (press - 0.3780 * moist_val * satvp)
    elif moist_prm == 2:  # dew point temperature is used as input
        # determine the latent heat of vaporization
        satvp = (1.000700 + 3.460000e-8 * press) * 611.2100 * m.exp(17.50200 * moist_val / (240.97000 + moist_val))
        spec_hum = 0.6220 * satvp / (press - 0.3780 * satvp)
    elif moist_prm == 3:  # wet bulb temperature is used as input
        spec_hum = 0.0000066 * (1.0 + 0.00115 * moist_val)  # dummy use of spec_hum
        satvp = (1.000700 + 3.460000e-8 * press) * 611.2100 * m.exp(
            17.50200 * moist_val / (240.97000 + moist_val)
        ) - spec_hum * press * (temperature - moist_val)
        spec_hum = 0.6220 * satvp / (press - 0.3780 * satvp)
    else:
        spec_hum = -1.0
    return spec_hum


def mft_fluxes(
    dyn_in_prm,
    dyn_in_val,
    dyn_in_val2,
    sfc_current1,
    sfc_current2,
    convect,
    pressure,
    air_moist_prm,
    air_moist_val,
    sfc_moist_prm,
    sfc_moist_val,
    salinity,
    ss_prm,
    ss_val,
    t_air,
    sst_prm,
    t_skin,
    ref_ht_wind,
    ref_ht_tq,
    z_wanted,
    astab,
    eqv_neut_prm,
    net_long_wave_flux,
    warn,
    flux_model,
    z0_mom_prm,
    z0_theta_q_prm,
    stab_prm,
    oil_fract_area,
    dimensionless_m_o_length,
    zo_input,
    missing,
):
    dyn_in_prm = validate_parameter(dyn_in_prm, 0, 10, "dyn_in_prm is outside acceptable range")
    dyn_in_val = validate_parameter(dyn_in_val, -100.0, 100.0, "dy_in_val is outside acceptable range")
    dyn_in_val2 = validate_parameter(dyn_in_val2, -100.0, 100.0, "dyn_in_val2 is outside acceptable range")
    sfc_current1 = validate_parameter(sfc_current1, -10.0, 10.0, "sfc_current is outside acceptable range")
    sfc_current2 = validate_parameter(sfc_current2, -10.0, 10.0, "sfc_current2 is outside acceptable range")
    convect = validate_parameter(convect, 0.0, 10.0, "convect is outside acceptable range")
    pressure = validate_parameter(pressure, 0.0, 120000.0, "pressure is outside acceptable range")
    air_moist_prm = validate_parameter(air_moist_prm, 0, 5, "air_moist_prm is outside acceptable range")
    air_moist_val = validate_parameter(air_moist_val, -100.0, 100.0, "air_moist_val is outside acceptable range")
    sfc_moist_prm = validate_parameter(sfc_moist_prm, 0, 5, "sfc_moist_prm is outside acceptable range")
    sfc_moist_val = validate_parameter(sfc_moist_val, -100.0, 100.0, "sfc_moist_val is outside acceptable range")
    salinity = validate_parameter(salinity, 0.0, 100.0, "salinity is outside acceptable range")
    ss_prm = validate_parameter(ss_prm, -0, 10, "ss_prm is outside acceptable range")
    ss_val = validate_parameter(ss_val, -1000.0, 1000.0, "ss_val is outside acceptable range")
    t_air = validate_parameter(t_air, -60.0, 100.0, "t_air is outside acceptable range " + str(t_air))
    sst_prm = validate_parameter(sst_prm, 0, 5, "sst_prm is outside acceptable range" + str(sst_prm))
    t_skin = validate_parameter(t_skin, -4.0, 100.0, "t_skin is outside acceptable range")
    ref_ht_wind = validate_parameter(ref_ht_wind, 0.0, 100.0, "ref_ht_wind is outside acceptable range")
    ref_ht_tq = validate_parameter(ref_ht_tq, 0.0, 100.0, "ref_ht_tq is outside acceptable range")
    z_wanted = validate_parameter(z_wanted, 0.0, 500.0, "z_wanted is outside acceptable range")
    astab = validate_parameter(astab, 0, 5, "astab is outside acceptable range")
    eqv_neut_prm = validate_parameter(eqv_neut_prm, 0, 5, "eqv_neut_prm is outside acceptable range")
    net_long_wave_flux = validate_parameter(
        net_long_wave_flux, 0.0, 1000.0, "net_long_wave_flux is outside acceptable range"
    )
    # warn = validate_parameter(warn, TRUE, FALSE, 'warn is outside acceptable range')
    flux_model = validate_parameter(flux_model, -10, 100, "flux_model is outside acceptable range")
    z0_mom_prm = validate_parameter(z0_mom_prm, 0, 10, "z0_mom_prm is outside acceptable range")
    z0_theta_q_prm = validate_parameter(z0_theta_q_prm, 0, 10, "z0_theta_q_prm is outside acceptable range")
    stab_prm = validate_parameter(stab_prm, 0, 10, "stab_prm is outside acceptable range")
    oil_fract_area = validate_parameter(oil_fract_area, 0.0, 1.0, "oil_fract_area is outside acceptable range")
    dimensionless_m_o_length = validate_parameter(
        dimensionless_m_o_length, -4.0, 4.0, "dimensionless_m_o_length is outside acceptable range"
    )
    zo_input = validate_parameter(zo_input, -4.0, 4.0, "zo_input is outside acceptable range")
    missing = validate_parameter(missing, -10000.0, 10000.0, "missing is outside acceptable range")

    conv_crit = 0.00005  # convergence criterion (fractional change)  [unitless]

    sfc_current = [sfc_current1, sfc_current2]
    monin_inv = dimensionless_m_o_length / ref_ht_wind

    # The physics breaks down for zero wind speed.  If the wind speed is zero, make an adjustment
    low_lim = 0.00001
    zero_speed_flag = FALSE
    dyn_in_mag = (dyn_in_val**2 + dyn_in_val2**2) ** 0.5
    if dyn_in_mag <= low_lim:
        zero_speed_flag = TRUE
        dyn_in_val = low_lim

    if abs(dyn_in_val) <= low_lim:
        dyn_in_val = low_lim

    if abs(dyn_in_val2) <= low_lim:
        dyn_in_val2 = low_lim

    if flux_model == 9 or flux_model == 10:
        z0_theta_q_prm = 1

    (
        bvw_flag,
        tau,
        lhf,
        shf,
        ustar,
        tstar,
        qstar,
        zzu,
        zzt,
        zzq,
        monin_inv,
        dsplcmnt_ht,
        wave_age,
        dom_wave_phs_spd,
        hsig,
        orbital_vel,
        sst_bulk_adjustment,
        stab_prm,
    ) = pmix_(
        dyn_in_prm,
        dyn_in_val,
        dyn_in_val2,
        sfc_current,
        convect,
        conv_crit,
        pressure,
        air_moist_prm,
        air_moist_val,
        sfc_moist_prm,
        sfc_moist_val,
        salinity,
        ss_prm,
        ss_val,
        t_air,
        sst_prm,
        t_skin,
        ref_ht_wind,
        ref_ht_tq,
        astab,
        net_long_wave_flux,
        flux_model,
        z0_mom_prm,
        z0_theta_q_prm,
        stab_prm,
        zo_input,
        monin_inv,
        oil_fract_area,
    )

    # record neutral roughness length
    zo_m = (ref_ht_wind - dsplcmnt_ht) / zzu
    zo_neut_eqv = zo_m

    q_sfc = find_q_(sfc_moist_prm, sfc_moist_val, pressure, t_skin)
    # determine the heat capacity at constant pressure
    air_specific_heat = 1004.0 * (1.0 + 0.9 * q_sfc)

    # correct to the height of z_wanted
    ht_adj_ratio = (z_wanted - dsplcmnt_ht) / (ref_ht_wind - dsplcmnt_ht)
    psi_u = psi_u_f(z_wanted - dsplcmnt_ht, ht_adj_ratio * zzu, monin_inv, stab_prm)
    if eqv_neut_prm == 0:
        # make adjustment normally for wind
        zzt = zzt * (z_wanted - dsplcmnt_ht) / (ref_ht_tq - dsplcmnt_ht)
        zzq = zzq * (z_wanted - dsplcmnt_ht) / (ref_ht_tq - dsplcmnt_ht)
        psi_theta = psi_theta_f(z_wanted - dsplcmnt_ht, zzt, monin_inv, stab_prm)
        psi_q = psi_q_f(z_wanted - dsplcmnt_ht, zzq, monin_inv, stab_prm)
    elif eqv_neut_prm == 1:
        # Tang and Liu eqv neut winds (use the calculated friction velocity, but not stability adjustment)
        psi_u = 0.0
        psi_theta = 0.0
        psi_q = 0.0

    elif eqv_neut_prm == 2:
        # if option was selected, calculate equivalent friction velocities as
        #     defined by Geernaert and Katsaros (1986)
        psi_u = psi_u_f(z_wanted - dsplcmnt_ht, z_wanted / zo_m[0], monin_inv, stab_prm)
        ustar[0] = pow(dyn_in_val / ustar[0] - psi_u / KV - m.log(zo_neut_eqv[0] / zo_m[0]) / KV, -1.0) * dyn_in_val
        ustar[1] = pow(dyn_in_val2 / ustar[1] - psi_u / KV - m.log(zo_neut_eqv[1] / zo_m[1]) / KV, -1.0) * dyn_in_val2
        psi_theta = 0.0
        psi_q = 0.0
    else:
        raise InvalidFluxModelParameterError("Invalid choice of eqv_neut_prm: " + f"{eqv_neut_prm}")

    ht_adj_ratio = (z_wanted - dsplcmnt_ht) / (ref_ht_tq - dsplcmnt_ht)
    t_at_z = t_skin + Prt * tstar * (m.log(zzt * ht_adj_ratio) - psi_theta) / KV - G * z_wanted / air_specific_heat
    q_at_z = q_sfc + Sc * qstar * (m.log(zzq * ht_adj_ratio) - psi_q) / KV

    u_at_z = [0.0, 0.0]
    # change the output vectors to the coordinate system used to input the data
    dyn_in_mag = m.sqrt(dyn_in_val**2 + dyn_in_val2**2)
    if dyn_in_mag != 0.0:
        vect_cmp1 = dyn_in_val / dyn_in_mag
        vect_cmp2 = dyn_in_val2 / dyn_in_mag
    else:
        vect_cmp1 = 0.0
        vect_cmp2 = 0.0
    if zzu != missing:
        if dyn_in_val != 0.0:
            relative_wind_speed = (
                ustar[0] * (m.log(zzu * (z_wanted - dsplcmnt_ht) / (ref_ht_wind - dsplcmnt_ht)) - psi_u) / KV
            )
            u_at_z[0] = relative_wind_speed * vect_cmp1 + sfc_current1 + orbital_vel[0]
            u_at_z[1] = relative_wind_speed * vect_cmp2 + sfc_current2 + orbital_vel[1]
        else:
            u_at_z[0] = 0.0
            u_at_z[1] = 0.0
        ustar_mag = (ustar[0] ** 2 + ustar[1] ** 2) ** 0.5
        if ustar_mag != 0.0:
            ustar[0] = ustar_mag * vect_cmp1
            ustar[1] = ustar_mag * vect_cmp2
        else:
            tau[0] = 0.0
            tau[1] = 0.0
        tau_mag = (tau[0] ** 2 + tau[1] ** 2) ** 0.5
        if tau_mag != 0.0:
            tau[0] = tau_mag * vect_cmp1
            tau[1] = tau_mag * vect_cmp2
        else:
            tau[0] = 0.0
            tau[1] = 0.0

    modified_sst = t_skin + sst_bulk_adjustment
    # print(
    #     "In mft_flxues: log term = ",
    #     ustar[0] * m.log(zzu + 1.0) / KV,
    #     "psi_term = ",
    #     ustar[0]*psi_u/KV,
    # )

    dimensionless_m_o_length = monin_inv * ref_ht_wind
    if warn and bvw_flag < 0:
        # Either something went noticeably wrong with the calculation, or an error in the input data was identified
        tau = [missing, missing]
        shf = missing
        lhf = missing
        ustar = [missing, missing]
        tstar = missing
        qstar = missing
        dimensionless_m_o_length = missing
        wave_age = [missing, missing]
        dom_wave_phs_spd = missing
        hsig = missing
        zo_m = missing
        u_at_z = missing
        t_at_z = missing
        q_at_z = missing
        modified_sst = missing

    return (
        bvw_flag,
        shf,
        lhf,
        tau,
        ustar,
        tstar,
        qstar,
        dimensionless_m_o_length,
        wave_age,
        dom_wave_phs_spd,
        hsig,
        zo_m,
        u_at_z,
        t_at_z,
        q_at_z,
        modified_sst,
    )
    # end of subroutine mft_fluxes


def bstar_f(neutral, qmixa, t_air, tstar, qstar, minimum_value):
    bstar = (
        float(not neutral)
        * G
        * ((1.00 + 0.610 * qmixa) * tstar + 0.6100 * (t_air + 273.15) * qstar)
        / ((t_air + 273.15) * (1.0 + 0.610 * qmixa))
    )
    if m.fabs(bstar) < minimum_value:
        if bstar == 0.0:
            bstar = minimum_value
        else:
            bstar = minimum_value * bstar / m.fabs(bstar)
    return bstar


# # # # # solve # # # # #
def solve(
    fixed_dom_wave_phs_spd,
    fixedwa,
    fixed_hsig,
    fixed_uen,
    fixed_ustar,
    neutral,
    use_dh,
    use_orb_vel,
    convect_adj,
    wave_len,
    crit,
    denwat,
    nu,
    qmixa,
    qmixw,
    sfcten,
    ss_prm,
    ss_val,
    wind_vect_original,
    sfc_current,
    t_air,
    sst_prm,
    tskin,
    net_long_wave_flux,
    zref,
    zreft,
    zrefq,
    ww_eql,
    betag_set,
    betac_set,
    betas_set,
    z0_prm,
    oil_fract_area,
    astab,
    denair,
    air_specific_heat,
    latent_heat_vaporization,
    theta_q_zo_prm,
    stab_prm,
    monin_inv,
    no_capw,
    no_sfcten,
    zzu,
    ustar,
):
    #    convert the temperatures to equivalent potential temperature
    #    estimated through the dry static energy
    # logger.debug('in solve, astab =', astab, ' neutral = ', neutral, 'TRUE = ', TRUE, 'FALSE = ', FALSE )

    ustar_old = [0.0, 0.0]
    directional_vect_comp = [
        1.0,
        0.0,
    ]  # direction of orbital velocity in a coordinate system relative to the wind shear (first guess)
    tskin_ept = tskin  # technically should add G * displacement height / air_specific_heat
    ta_ept = t_air + G * zref / air_specific_heat
    dsplcmnt_ht = 0.0
    wind_vect = wind_vect_original
    sst_bulk_adjustment = 0

    ww = 1.0
    dom_wave_phs_spd = [4.0, 0.0]
    wave_age = [28.0, 28.0]
    hsig = [1.0, 0.0]
    period = [1.0, 1.0]
    orbital_vel = [0.0, 0.0]
    if ss_prm == 0:
        ww = ss_val
    elif ss_prm == 1:
        dom_wave_phs_spd = [ss_val, 0.0]
    elif ss_prm == 2:
        wave_age = [ss_val, 0.0]
    elif ss_prm == 3:
        hsig = [ss_val, 0.0]
    elif ss_prm == 4:
        ...
    elif ss_prm == 5:
        period = [ss_val, 1.0]
    elif ss_prm != 6:
        orbital_vel = [ss_val, 0.0]
    else:
        raise InvalidFluxModelParameterError("Invalid choice of ss_prm: " + f"{ss_prm}")
    adjusted_orbital_vel = orbital_vel  # seems like this should be 0.8 * orb_vel: check

    betag_prime = betag_set
    betag = [betag_prime, betag_prime]
    betag[1] = betag_prime
    betac_prime = betac_set
    betac = [betac_prime, betac_prime]
    betas_prime = betas_set
    betas = [betas_prime, betas_prime]

    #    intial guess for value of friction velocity (ustar)
    #    component parallel direction of wave motion
    ustar[0] = float(fixed_ustar) * ustar[0] + float(not fixed_ustar) * 0.003
    ustar[1] = 0.0  # component perpendicular direction of wave motion
    viscosity = 1.4e-5
    alpha = 1.0 / (t_air + 273.15)
    zo_mag = zref / zzu  # first guess (should use zref - displacement height)
    ustar_mag = m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1])
    renewal_time = get_renewal_time(zo_mag, ustar_mag, denair, air_specific_heat, alpha, net_long_wave_flux, viscosity)
    stanton = get_stanton(ustar_mag, renewal_time)
    dalton = get_dalton(ustar_mag, renewal_time)
    zzt = z_o_zt(zreft - dsplcmnt_ht, nu, ustar, zo_mag, stanton, theta_q_zo_prm)
    zzq = z_o_zq(zrefq - dsplcmnt_ht, nu, ustar, zo_mag, dalton, theta_q_zo_prm)

    tstar = 0.0003 * (t_air - tskin)
    if m.fabs(tstar) < 0.001:
        tstar = -0.001
    qstar = 0.0003 * (qmixa - qmixw)
    if m.fabs(qstar) < 0.0001:
        qstar = -0.0001

    ustar_old[0] = ustar[0] / 2.00
    ustar_old[1] = ustar[1] / 2.00
    cmin = m.sqrt(2.0 * m.sqrt(G * sfcten / denwat))

    count_u_loop = 0
    unreasonable = TRUE
    while unreasonable == TRUE:
        # iterate over atmospheric stability and friction velocity
        if astab != 2:
            monin_inv_old2 = -0.001  # initial guess
        else:
            monin_inv_old2 = monin_inv
        ustar_mag_old = 0.0
        ustar_mag = m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1])
        min_ustar = 0.00001
        if ustar_mag < min_ustar:
            ustar_mag = min_ustar
        count_u_loop = 0
        while (
            (
                m.fabs(monin_inv - monin_inv_old2) > crit * m.fabs(monin_inv_old2)
                or m.fabs(ustar_mag - ustar_mag_old) > crit * m.fabs(ustar_mag_old)
                or count_u_loop < 2
            )
            and count_u_loop < 30
            and not (fixed_ustar and count_u_loop >= 1)
        ):
            count_u_loop = count_u_loop + 1
            # determine roughness lengths for heat and moisture profiles
            viscosity = 1.4e-5
            alpha = 1.0 / (t_air + 273.15)
            zo_mag = (zref - dsplcmnt_ht) / zzu
            renewal_time = get_renewal_time(
                zo_mag, ustar_mag, denair, air_specific_heat, alpha, net_long_wave_flux, viscosity
            )
            stanton = get_stanton(ustar_mag, renewal_time)
            dalton = get_dalton(ustar_mag, renewal_time)
            zzt = z_o_zt(zreft - dsplcmnt_ht, nu, ustar, zo_mag, stanton, theta_q_zo_prm)
            zzq = z_o_zq(zrefq - dsplcmnt_ht, nu, ustar, zo_mag, dalton, theta_q_zo_prm)

            ustar_mag_old = ustar_mag
            monin_inv_old2 = monin_inv

            # iterate over solutions to tstar, qstar, and monin_inv
            tstar_old = 0.5 * tstar
            qstar_old = 0.5 * qstar
            monin_inv_old = monin_inv
            count_tq_loop = 0

            while (
                m.fabs(tstar - tstar_old) > crit * m.fabs(tstar_old)
                or m.fabs(qstar - qstar_old) > crit * m.fabs(qstar_old)
                or m.fabs(monin_inv - monin_inv_old) > 0.1 * crit * m.fabs(monin_inv_old)
                or count_tq_loop < 3
            ) and count_tq_loop < 40:
                count_tq_loop = count_tq_loop + 1

                tstar_old = tstar
                qstar_old = qstar
                monin_inv_old = monin_inv

                if ta_ept == tskin_ept:
                    tstar = 0.0
                else:
                    # determine the eqv. pot. temperature stability term
                    if neutral:
                        psi_theta = 0.0
                    else:
                        psi_theta = psi_theta_f(zreft, zzt, monin_inv, stab_prm)

                    # Determine the cool skin and warm layer adjustments to the sea minus air potential
                    # temperature difference
                    wind_mag = m.sqrt(wind_vect[0] ** 2 + wind_vect[1] ** 2)

                    if sst_prm == 0 or count_tq_loop < 2:
                        sst_bulk_adjustment = 0
                    elif sst_prm == 1:  # Sandra Castro's wind speed dependent exponential decay
                        sst_bulk_adjustment = -0.138 - 0.181 * m.exp(-0.135 * wind_mag)
                    elif sst_prm == 2:  # Sandra Castro's wind speed dependent polynomial fit
                        # Warning - this is flawed for wind speeds too much greater than 15 m/s
                        sst_bulk_adjustment = (
                            -0.389
                            + 0.075 * wind_mag
                            - 0.012 * wind_mag**2
                            + 0.0009 * wind_mag**3
                            - 0.000022 * wind_mag**4
                        )
                        if wind_mag > 15:
                            logger.warning("Cool skin, sst_prm = 2, wind speed outside valid range.")
                    elif sst_prm == 3:  # Sandra Castro's more complex version
                        sst_bulk_adjustment = 0
                    elif sst_prm == 4:
                        # CFC combo for warm layer and curl skin
                        shf = -denair * air_specific_heat * ustar_mag * tstar
                        lhf = -denair * latent_heat_vaporization * ustar_mag * qstar
                        rs = 30.0  # should be the net downward solar radiation
                        sst_bulk_adjustment = skin_temp_f(
                            zo_mag, ustar_mag, net_long_wave_flux, shf, lhf, rs, denair, air_specific_heat
                        )
                    else:
                        raise InvalidFluxModelParameterError("Invalid choice of sst_prm: " + f"{sst_prm}")

                    delta_theta = ta_ept - tskin_ept + sst_bulk_adjustment
                    if m.fabs(delta_theta) < 0.001:
                        delta_theta = -0.001
                    tstar = KV * delta_theta / (m.log(zzt + 1.0) - float(not neutral) * psi_theta) / Prt
                    # logger.debug("In solve, psi_theta = ', psi_theta, 'neutral = ', neutral, float(not neutral))
                if qmixa == qmixw:
                    qstar = 0
                else:
                    # determine the moisture stability term  */
                    if neutral:
                        psi_q = 0.0
                    else:
                        psi_q = psi_q_f(zreft - dsplcmnt_ht, zzq, monin_inv, stab_prm)

                    delta_q = qmixa - qmixw
                    if m.fabs(delta_q) < 0.0001:
                        delta_q = -0.0001
                    qstar = KV * delta_q / (m.log(zzq + 1.0) - float(not neutral) * psi_q) / Sc

                # determine the buoyancy flux
                bstar = bstar_f(neutral, qmixa, t_air, tstar, qstar, 0.000000)

                # !!!!This assumes temperature and humidity are measured at THE SAME HEIGHT

                # determine the inverse Monin-Obhukov scale length: 1/L
                ustar2_mag = ustar[0] * ustar[0] + ustar[1] * ustar[1]
                if astab != 2:
                    if ustar2_mag != 0.0:
                        monin_inv = KV * bstar / ustar2_mag
                    else:
                        monin_inv = 0.0
                # end of iteration on tstar, qstar, and monin_inv

            done = FALSE
            count = 0
            while (done == FALSE or count < 2) and count < 100:
                count = count + 1
                done = TRUE
                wave_age_old = wave_age
                zzu_old = zzu
                # Store the original vector for wind (or other input variable)
                # adjust to a coordinate system where the first component is parallel the wind shear
                if ss_prm == 6:  # orbital velocity is input in the original coordinate system
                    adjusted_orbital_vel = orbital_vel
                else:  # orbital velocity is transformed to the original wind coordinate system
                    adjusted_orbital_vel = [
                        orbital_vel[0] * directional_vect_comp[0] + orbital_vel[1] * directional_vect_comp[1],
                        -orbital_vel[0] * directional_vect_comp[1] + orbital_vel[1] * directional_vect_comp[0],
                    ]
                vect_cmp1 = wind_vect_original[0] - sfc_current[0] - float(use_orb_vel) * adjusted_orbital_vel[0]
                vect_cmp2 = wind_vect_original[1] - sfc_current[1] - float(use_orb_vel) * adjusted_orbital_vel[1]
                wind_shear_mag = (vect_cmp1**2 + vect_cmp2**2) ** 0.5
                if wind_shear_mag > 0.0:
                    directional_vect_comp[0] = vect_cmp1 / wind_shear_mag
                    directional_vect_comp[0] = vect_cmp2 / wind_shear_mag
                wind_vect = [wind_shear_mag, 0.0]
                # iteratively solve for u* for non-neutral conditions
                if not fixed_ustar:
                    dsplcmnt_ht, monin_inv, ustar, zzu, betac, betag, betas = find_ustar(
                        fixed_dom_wave_phs_spd,
                        fixedwa,
                        fixed_hsig,
                        fixed_uen,
                        neutral,
                        no_capw,
                        no_sfcten,
                        use_dh,
                        use_orb_vel,
                        betac_prime,
                        betag_prime,
                        betac,
                        betag,
                        betas,
                        convect_adj,
                        wave_age,
                        dom_wave_phs_spd,
                        hsig,
                        period,
                        wave_len,
                        orbital_vel,
                        denwat,
                        nu,
                        qmixa,
                        sfcten,
                        wind_vect,
                        t_air,
                        ww,
                        ww_eql,
                        ss_prm,
                        ss_val,
                        zref,
                        dsplcmnt_ht,
                        z0_prm,
                        zzu,
                        ustar,
                        tstar,
                        qstar,
                        stab_prm,
                        monin_inv,
                        oil_fract_area,
                    )
                else:
                    i = 0
                    if fixed_hsig and i == 0:
                        hsig_known = TRUE
                    else:
                        hsig_known = FALSE
                    (
                        trouble,
                        ustar,
                        zzu,
                        dsplcmnt_ht,
                        hsig,
                        wave_len,
                        orbital_vel,
                        wave_age,
                        dom_wave_phs_spd,
                        betac,
                        betag,
                        betas,
                    ) = z0_and_waves(
                        fixed_dom_wave_phs_spd,
                        fixedwa,
                        neutral,
                        hsig_known,
                        stab_prm,
                        no_capw,
                        no_sfcten,
                        use_dh,
                        use_orb_vel,
                        betac_prime,
                        betag_prime,
                        betac,
                        betag,
                        betas,
                        denwat,
                        nu,
                        sfcten,
                        zref,
                        dsplcmnt_ht,
                        monin_inv,
                        wave_age,
                        dom_wave_phs_spd,
                        hsig,
                        period,
                        wave_len,
                        orbital_vel,
                        ww,
                        ww_eql,
                        ss_prm,
                        ss_val,
                        i,
                        z0_prm,
                        zzu,
                        ustar,
                        oil_fract_area,
                    )

                    i = 1
                    hsig_known = FALSE
                    if ustar[1] != 0.0:
                        (
                            trouble,
                            ustar,
                            zzu,
                            dsplcmnt_ht,
                            hsig,
                            wave_len,
                            orbital_vel,
                            wave_age,
                            dom_wave_phs_spd,
                            betac,
                            betag,
                            betas,
                        ) = z0_and_waves(
                            fixed_dom_wave_phs_spd,
                            fixedwa,
                            neutral,
                            hsig_known,
                            stab_prm,
                            no_capw,
                            no_sfcten,
                            use_dh,
                            use_orb_vel,
                            betac_prime,
                            betag_prime,
                            betac,
                            betag,
                            betas,
                            denwat,
                            nu,
                            sfcten,
                            zref,
                            dsplcmnt_ht,
                            monin_inv,
                            wave_age,
                            dom_wave_phs_spd,
                            hsig,
                            period,
                            wave_len,
                            orbital_vel,
                            ww,
                            ww_eql,
                            ss_prm,
                            ss_val,
                            i,
                            z0_prm,
                            zzu,
                            ustar,
                            oil_fract_area,
                        )

                if (
                    m.fabs((wave_age[0] - wave_age_old[0])) > crit * wave_age[0]
                    or m.fabs((wave_age[1] - wave_age_old[1])) > crit * wave_age[1]
                    or m.fabs((zzu - zzu_old)) > crit * zzu
                ) and count < 60:
                    done = FALSE

            ustar_mag = m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1])
            if ustar_mag < min_ustar:
                ustar_mag = min_ustar
            # end of iteration over friction velocity and atmospheric stability

        # determine if the phase speed of the dominant waves is physically
        # reasonable. If the phase speed is less than the physical minimum,
        # then assume that no capillary waves are present. After the second
        # pass, if phase speed is still unacceptable then assume the surface
        # is smooth.

        if (
            fixedwa
            or fixed_dom_wave_phs_spd
            or (betag[0] != 0.0 and wave_age[0] * m.fabs(ustar[0]) > cmin)
            or (betag[1] != 0.0 and wave_age[1] * m.fabs(ustar[1]) > cmin)
            or betas_prime == 1.0
        ):
            unreasonable = FALSE
        else:
            if betag[0] != 0.0 and betag[1] != 0.0 and betac_prime != 0.0:
                if m.fabs(wind_vect[0]) > m.fabs(wind_vect[1]):
                    betag[1] = 0.0
                else:
                    betag[0] = 0.0 + float(no_capw)
            elif betac_prime > 0.0:
                betag[0] = 1.0
                betag[1] = 1.0
                betac_prime = 0.0
            else:
                if betag[0] > 0.0 and wave_age[0] * m.fabs(ustar[0]) < cmin:
                    betag[0] = 0.0 + float(no_capw)
                if betag[1] > 0.0 and wave_age[1] * m.fabs(ustar[1]) < cmin:
                    betag[1] = 0.0 + float(no_capw)
                if betag[0] == 0.0 and betag[1] == 0.0:
                    betag_prime = 0.0 + float(no_capw)
                    betas_prime = 1.0
        if no_capw:
            betac_prime = 0.0
            betag_prime = 1.0
            betas_prime = 1.0

    return (
        count_u_loop + fixed_ustar,
        zzu,
        zzt,
        zzq,
        ustar,
        tstar,
        qstar,
        dsplcmnt_ht,
        monin_inv,
        hsig,
        wave_len,
        orbital_vel,
        wave_age,
        dom_wave_phs_spd,
        adjusted_orbital_vel,
        sst_bulk_adjustment,
    )


# # # # # ustar_f # # # # #
def ustar_f(fixed_uen, neutral, bstar, ustar, convect_adj, psi_u, monin_inv, delta_u, zzu, max_value):
    low_lim = 0.000001
    wstar2 = pow(m.fabs(m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1]) * bstar) * H, 0.6666667)

    if monin_inv >= 0.0:
        wstar2 = 0.0
    zzu_dum = zzu
    if zzu_dum < 0.00001:
        zzu_dum = 0.00001
    if m.fabs(delta_u) < low_lim:
        ustar_test = 0.0
    elif m.log(zzu_dum + 1.0) - float(not fixed_uen) * float(not neutral) * psi_u >= low_lim:
        delta_u2 = delta_u * delta_u + convect_adj * wstar2
        if m.fabs(delta_u2) < low_lim:
            delta_u2 = low_lim
        delta_u = m.sqrt(delta_u2) * delta_u / m.fabs(delta_u)
        ustar_test = KV * delta_u / (m.log(zzu_dum + 1.0) - float(not fixed_uen) * float(not neutral) * psi_u)
        if ustar_test > max_value:
            ustar_test = 1.0  # max_value - 1.0
    else:
        ustar_test = max_value
    # print(
    #     "In ustar_f: delta_u = ",
    #     delta_u,
    #     "log term = ",
    #     m.log(zzu_dum + 1.0),
    #     "psi_term = ",
    #     float(not fixed_uen) * float(not neutral) * psi_u,
    # )

    return ustar_test


# # # # # psi_u_f # # #
def psi_u_f(zref, zzu_dum, monin_inv, stab_prm):
    # input height should be relative to the displacement height
    a_bh = 0.7
    b_bh = 0.75
    c_bh = 5.0
    d_bh = 0.35

    # unstable case
    if monin_inv < 0:
        if stab_prm == 0 or stab_prm == 2:  # Businger-Dyer parameterization
            # BVW
            # find zeta( (z-d) / L ) and zeta( zo / L )
            if zzu_dum < 0.0001:
                zzu_dum = 0.0001
            zeta = pow(1.0 - 15.0 * zref * monin_inv, 0.25)
            zeta0 = pow(1.0 - 15.0 * zref / zzu_dum * monin_inv, 0.25)
            psi_u = (
                2.0 * m.log((1.0 + zeta) / (1.0 + zeta0))
                - m.log((1.0 + zeta**2.0) / (1.0 + zeta0**2.0))
                + 2.0 * m.atan(zeta)
                - 2.0 * m.atan(zeta0)
            )
        elif stab_prm == 1:  # Smith 88 parameterization  */
            a_bh = m.sqrt(m.sqrt(1.0 - 16.0 * zref * monin_inv))
            psi_u = (
                2.0 * m.log(0.5 * (1.0 + a_bh)) + m.log(0.5 * (1.0 + a_bh * a_bh)) - 2.0 * m.atan(a_bh) + 3.14159 / 2.0
            )
        elif stab_prm == 3 and zzu_dum != 0.0:  # COAWST
            zeta = pow(1.000 - 16.000 * zref * monin_inv, 0.2500)
            zeta0 = pow(1.000 - 16.000 * zref / zzu_dum * monin_inv, 0.2500)
            psi_u = (
                2.0 * m.log((1.0 + zeta) / (1.0 + zeta0))
                - m.log((1.0 + zeta**2.0) / (1.0 + zeta0**2.0))
                + 2.0 * m.atan(zeta)
                - 2.0 * m.atan(zeta0)
            )
        elif stab_prm == 4:  # Paulson 1970
            zeta = pow(1.000 - 16.000 * zref * monin_inv, 0.2500)
            psi_u = (
                2.0 * m.log((1.0 + zeta) / 2.0) - m.log((1.0 + zeta**2.0) / 2.0) + 2.0 * m.atan(zeta) - 3.14159 / 2
            )
        else:
            raise InvalidFluxModelParameterError("Invalid choice of stab_prm: " + f"{stab_prm}")

    else:
        #  stable case  */
        if stab_prm == 0:
            # BVW
            psi_u = -(
                a_bh * zref * monin_inv
                + b_bh * (zref * monin_inv - c_bh / d_bh) * m.exp(-d_bh * zref * monin_inv)
                + b_bh * c_bh / d_bh
            )
        #    logger.warning('psi_u = ', psi_u, a_bh * zref * monin_inv,
        #     b_bh * (zref * monin_inv - c_bh / d_bh) * m.exp(-d_bh * zref * monin_inv),
        #     b_bh * c_bh / d_bh)
        elif stab_prm == 1 or stab_prm == 2 or stab_prm == 4:  # Hicks parameterization, also used in Smith 88
            psi_u = -5.0 * zref * monin_inv
        elif stab_prm == 3 and zzu_dum != 0.0:
            psi_u = -5.0 * (zref * monin_inv - zref / zzu_dum * monin_inv)
        else:
            raise InvalidFluxModelParameterError("Invalid choice of stab_prm: " + f"{stab_prm}")

    return psi_u


def psi_q_f(zrefq, zzq, monin_inv, stab_prm):
    # input height should be relative to the displacement height
    a_bh = 1.0
    b_bh = 0.667
    c_bh = 5.0
    d_bh = 0.35

    if monin_inv < 0:
        if stab_prm == 0 or stab_prm == 2:
            # BVW
            # find zeta( z / L ) and zeta( zo / L )
            qlamda = m.sqrt(1.000 - 9.000 * zrefq * monin_inv)
            q0lamda = m.sqrt(1.000 - 9.000 * zrefq / zzq * monin_inv)
            psi_q = 2.0000 * m.log((q0lamda + 1) / (qlamda + 1))
        elif stab_prm == 1:  # Smith 88 parameterization
            a_bh = m.sqrt(1.0 - 16.0 * zrefq * monin_inv)
            psi_q = 2.0 * m.log((1.0 + a_bh) / 2.0)
        elif stab_prm == 3:
            qlamda = pow(1.000 - 16.000 * zrefq * monin_inv, 0.500)
            q0lamda = pow(1.000 - 16.000 * zrefq / zzq * monin_inv, 0.500)
            psi_q = 2.0 * m.log((1.0 + qlamda) / (1.0 + q0lamda))
        else:
            raise InvalidFluxModelParameterError("Invalid choice of stab_prm: " + f"{stab_prm}")
    else:
        #  stable case
        if stab_prm == 0:  # BVW
            psi_q = -(
                pow(1.0 + 0.6667 * a_bh * zrefq * monin_inv, 1.5)
                + b_bh * (zrefq * monin_inv - c_bh / d_bh) * m.exp(-d_bh * zrefq * monin_inv)
                + b_bh * c_bh / d_bh
                - 1.0
            )
        elif stab_prm == 1 or stab_prm == 2:  # Smith 88 parameterization
            psi_q = -5.0 * zrefq * monin_inv
        elif stab_prm == 3:
            psi_q = -5.0 * (zrefq * monin_inv - zrefq / zzq * monin_inv)
        else:
            raise InvalidFluxModelParameterError("Invalid choice of stab_prm: " + f"{stab_prm}")

    return psi_q


# # # # # psi_theta_f # # # # #
def psi_theta_f(zreft, zzt, monin_inv, stab_prm):
    # input height should be relative to the displacement height
    a_bh = 1.0
    b_bh = 0.667
    c_bh = 5.0
    d_bh = 0.35

    if monin_inv < 0:
        if stab_prm == 0 or stab_prm == 2:  # BVW
            # find zeta( z / L ) and zeta( zo / L )  */
            tlamda2 = pow(1.000 - 9.000 * zreft * monin_inv, 0.500)
            t0lamda2 = pow(1.000 - 9.000 * zreft / zzt * monin_inv, 0.500)  # lamda squared
            psi_theta = 2.0 * m.log((1.0 + tlamda2) / 2.0) - 2.0 * m.log((1.0 + t0lamda2) / 2.0)
        elif stab_prm == 1:  # Smith 88 parameterization
            tlamda2 = pow(1.000 - 9.000 * zreft * monin_inv, 0.500)
            psi_theta = 2.0 * m.log((1.0 + tlamda2) / 2.0)
        elif stab_prm == 3:
            tlamda2 = pow(1.000 - 16.000 * zreft * monin_inv, 0.500)
            t0lamda2 = pow(1.000 - 16.000 * zreft / zzt * monin_inv, 0.500)
            psi_theta = 2.0000 * m.log((tlamda2 + 1.0) / (t0lamda2 + 1.0))
        else:
            raise InvalidFluxModelParameterError("Invalid choice of stab_prm: " + f"{stab_prm}")
    else:
        #  stable case
        if stab_prm == 0:  # BVW
            psi_theta = (
                -pow(1.0 + 0.6667 * a_bh * zreft * monin_inv, 1.5)
                - b_bh * (zreft * monin_inv - c_bh / d_bh) * m.exp(-d_bh * zreft * monin_inv)
                - b_bh * c_bh / d_bh
                + 1.0
            )
            # COARE fluxes - probably v3.0
            # psi_theta = ( ( 1 + b_bh * zreft * monin_inv )**1.5 + b_bh * ( zreft * monin_inv -
            #  14.28)/m.exp( zreft*monin_inv ) + 8.525 )
        elif stab_prm == 1 or stab_prm == 2:  # Smith 88 parameterization
            psi_theta = -5.0 * zreft * monin_inv
        elif stab_prm == 3:
            psi_theta = -5.0 * (zreft * monin_inv - zreft / zzt * monin_inv)
        else:
            raise InvalidFluxModelParameterError("Invalid choice of stab_prm: " + f"{stab_prm}")
    # print(
    #     "in psi_theta_f: zreft = ",
    #     zreft,
    #     " zzt = ",
    #     zzt,
    #     "1/L = ",
    #     monin_inv,
    #     "stab_prm = ",
    #     stab_prm,
    #     "psi_tehta = ",
    #     psi_theta,
    # )
    return psi_theta


# # # # # z_o_zq # # # # #
def z_o_zq(zrefq, nu, ustar, zo_mag, dalton, theta_q_zo_prm):
    #  zzq is  zref / z0Q
    ustar_mag = m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1])
    min_ustar = 0.00001
    if ustar_mag < min_ustar:
        ustar_mag = min_ustar
    if theta_q_zo_prm == 0:  # wall theory
        zzq = zrefq / (0.620 * nu / m.fabs(ustar_mag))
    elif theta_q_zo_prm == 1:  # CFC model
        zzq = zrefq / (zo_mag * m.exp(KV * (5.0 - 1.0 / (Sc * dalton))))
    elif theta_q_zo_prm == 2:  # Zilitinkevich et al. 2001
        zzq = zrefq / (zo_mag * m.exp(-4.0 * m.sqrt(zo_mag * ustar_mag / nu)))
    elif theta_q_zo_prm == 3:  # LKB 1979
        reynolds_rough_num = zo_mag * ustar_mag / nu
        if reynolds_rough_num < 0.11:
            const1 = 0.177
            const2 = 0
        elif reynolds_rough_num < 0.825:
            const1 = 1.376
            const2 = 0.929
        elif reynolds_rough_num < 3.0:
            const1 = 1.026
            const2 = -0.599
        elif reynolds_rough_num < 10.0:
            const1 = 1.625
            const2 = -1.018
        elif reynolds_rough_num < 30.0:
            const1 = 4.661
            const2 = -1.475
        elif reynolds_rough_num < 100.0:
            const1 = 34.904
            const2 = -2.067
        else:
            logger.warning("LKB zoq calculation failed due to large reynolds_rough_num: " + f"{reynolds_rough_num}")
            const1 = 34.904
            const2 = -2.067
        zzq = zrefq * ustar_mag * const1 * pow(reynolds_rough_num, const2) / nu
    elif theta_q_zo_prm == 4:  # COARE3.0
        zzq = 0.000055 * pow(zo_mag * ustar_mag / nu, -0.6)  # zo_q */
        if zzq < 0.00011:
            zzq = 0.00011
        zzq = zrefq / zzq
    elif theta_q_zo_prm == 5:  # Modified CFC model
        zzq = zrefq / (zo_mag * m.exp(KV * (0.5 - 1.0 / (Sc * dalton))))
    else:
        raise InvalidFluxModelParameterError("Invalid choice of theta_q_zo_prm: " + f"{theta_q_zo_prm}")
    return zzq


def z_o_zt(zreft, nu, ustar, zo_mag, stanton, theta_q_zo_prm):
    ustar_mag = m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1])
    min_ustar = 0.00001
    if ustar_mag < min_ustar:
        ustar_mag = min_ustar
    if theta_q_zo_prm == 0:
        # determine zreft / zT
        zzt = zreft / (0.40 * nu / m.fabs(ustar_mag))
    elif theta_q_zo_prm == 1:
        zzt = zreft / (zo_mag * m.exp(KV * (5.0 - 1.0 / (Prt * stanton))))
    elif theta_q_zo_prm == 2:
        # Zilitinkevich et al. 2001
        if zo_mag * m.exp(-4.0 * m.sqrt(zo_mag * ustar_mag / nu)) > 0.00001:
            zzt = zreft / (zo_mag * m.exp(-4.0 * m.sqrt(zo_mag * ustar_mag / nu)))
        else:
            zzt = 1000.0
    elif theta_q_zo_prm == 3:
        # LKB 1979
        reynolds_rough_num = zo_mag * ustar_mag / nu
        if reynolds_rough_num < 0.11:
            const1 = 0.292
            const2 = 0
        elif reynolds_rough_num < 0.825:
            const1 = 1.808
            const2 = 0.826
        elif reynolds_rough_num < 3.0:
            const1 = 1.393
            const2 = -0.528
        elif reynolds_rough_num < 10.0:
            const1 = 1.956
            const2 = -0.870
        elif reynolds_rough_num < 30.0:
            const1 = 4.994
            const2 = -1.297
        elif reynolds_rough_num < 100.0:
            const1 = 30.790
            const2 = -1.845
        else:
            logger.warning("LKB zot calculation failed due to large reynolds_rough_num: " + f"{reynolds_rough_num}")
            const1 = 30.790
            const2 = -1.845
        zzt = zreft * ustar_mag * const1 * pow(reynolds_rough_num, const2) / nu
    elif theta_q_zo_prm == 4:
        # COARE3.0
        zzt = 0.000055 * pow(zo_mag * ustar_mag / nu, -0.6)  # zo t
        if zzt < 0.00011:
            zzt = 0.00011
        zzt = zreft / zzt
    elif theta_q_zo_prm == 5:  # Modified CFC model
        zzt = zreft / (zo_mag * m.exp(KV * (0.5 - 1.0 / (Prt * stanton))))
    else:
        raise InvalidFluxModelParameterError("Invalid choice of theta_q_zo_prm: " + f"{theta_q_zo_prm}")
    return zzt


# # # # # z0_and_waves # # # # #
def z0_and_waves(
    fixed_dom_wave_phs_spd,
    fixedwa,
    neutral,
    hsig_known,
    stab_prm,
    no_capw,
    no_sfcten,
    use_dh,
    use_orb_vel,
    betac_prime,
    betag_prime,
    betac,
    betag,
    betas,
    denwat,
    nu,
    sfcten,
    zref,
    dsplcmnt_ht,
    monin_inv,
    wave_age,
    dom_wave_phs_spd,
    hsig,
    period,
    wave_len,
    orbital_vel,
    ww,
    ww_eql,
    ss_prm,
    ss_val,
    i,
    z0_prm,
    zzu,
    ustar,
    oil_fract_area,
):
    test = 0
    cmin = m.sqrt(2.0 * m.sqrt(G * sfcten / denwat))
    z0 = (zref - dsplcmnt_ht) / zzu
    # For the flux models (i.e., BVW) that calculate sea state
    if z0_prm == 0 or z0_prm == 2 or z0_prm == 3:
        b_toba = 0.0602
        # the sea state parameter is significant slope
        if ss_prm == 4 and i == 0:
            sig_slope = ss_val
            dom_wave_phs_spd[i] = 2.0 * 3.14159 * m.fabs(ustar[i]) * b_toba * b_toba / sig_slope / sig_slope
            wave_age[i] = dom_wave_phs_spd[i] / m.fabs(ustar[i])
            if wave_age[i] > 12000.0:
                wave_age[i] = 12000.0
            z0 = (
                betas[i] * 0.1100 * nu / m.fabs(ustar[i])
                + betac[i] * B * sfcten / (ustar[i] * ustar[i] * denwat)
                + betag[i] * 0.016 * ustar[i] * ustar[i] / G
            )
            psi_u = psi_u_f(betag[i] * hsig[i] - dsplcmnt_ht, hsig[i] / z0, monin_inv, stab_prm)
            ww = (m.log(hsig[i] / z0 + 1.0) - psi_u) / (KV * ww_eql * wave_age[i])
            # determine the wave characteristics
            period[i] = pow(pow(dom_wave_phs_spd[i] * sig_slope / b_toba, 2.0) / (G * ustar[i]), 0.333333)
            hsig[i] = b_toba * m.sqrt(G * m.fabs(ustar[i]) * pow(period[i], 3))
            wave_len = dom_wave_phs_spd[i] * period[i]
            orbital_vel[i] = PI2 * 0.5 * hsig[i] / period[i]
            # calculate displacement height, multiplied by zero or one depending on if it is used.
            # dsplcmnt_ht = float(use_dh) * hsig[0] * 0.8
        # the sea state parameter is significant wave height OR period
        if (hsig_known or ss_prm == 5) and i == 0:
            if ss_prm == 5:
                period_known = TRUE
            else:
                period_known = FALSE

            period[i] = (
                float(not period_known) * pow(hsig[i] * hsig[i] / (G * m.fabs(ustar[i]) * b_toba * b_toba), 0.333333)
                + float(period_known) * ss_val
            )
            hsig[i] = (
                float(period_known) * b_toba * m.sqrt(G * m.fabs(ustar[i]) * pow(period[i], 3))
                + float(not period_known) * ss_val
            )
            p1 = G * period[i] / (2.0 * 3.14159)
            p2 = pow(cmin, 4.0)
            p3 = 0.5 * p2 / p1 + pow(p1, 3.0) / 27.0
            p4 = 0.25 * p2 * p2 / (p1 * p1) + p1 * p1 * p2 / 27.0
            dom_wave_phs_spd[i] = pow(p3 + m.sqrt(p4), 0.333333) + pow(p3 - m.sqrt(p4), 0.333333) + p1 / 3.0
            wave_age[i] = dom_wave_phs_spd[i] / m.fabs(ustar[i])

            if wave_age[i] > 12000.0:
                wave_age[i] = 12000.0
            # determine the wave characteristics
            wave_len = dom_wave_phs_spd[i] * period[i]
            #  calculate displacement height, multiplied by zero or one depending on if it is used.
            dsplcmnt_ht = float(use_dh) * hsig[0] * 0.8
            # print("waves: i = ", i, " ustar = ", ustar, "dom_wave_phs_spd = ", dom_wave_phs_spd, " wa = ", wave_age)

            ustar, zzu, betac, betag, betas = z_over_z0(
                use_orb_vel,
                denwat,
                nu,
                sfcten,
                zref,
                dsplcmnt_ht,
                wave_age,
                hsig,
                wave_len,
                i,
                z0_prm,
                stab_prm,
                monin_inv,
                oil_fract_area,
                betac,
                betag,
                betas,
                zzu,
                ustar,
            )

            psi_u = psi_u_f(betag[i] * hsig[i] - dsplcmnt_ht, hsig[i] / z0, monin_inv, stab_prm)
            ww = (m.log(hsig[i] / z0 + 1.0) - psi_u) / (KV * ww_eql * wave_age[i])
            orbital_vel[i] = PI2 * 0.5 * hsig[i] / period[i]
        else:
            #  calculate displacement height, multiplied by zero or one depending on if it is used.
            dsplcmnt_ht = float(use_dh) * hsig[0] * 0.8

            ustar, zzu, betac, betag, betas = z_over_z0(
                use_orb_vel,
                denwat,
                nu,
                sfcten,
                zref,
                dsplcmnt_ht,
                wave_age,
                hsig,
                wave_len,
                i,
                z0_prm,
                stab_prm,
                monin_inv,
                oil_fract_area,
                betac,
                betag,
                betas,
                zzu,
                ustar,
            )

            # determine the wave characteristics
            if betag[i] > 0.1:
                dum = ustar[i] * ustar[i] * wave_age[i] * wave_age[i]
                if dum > cmin * cmin:
                    wave_len = 3.14159 * (dum + (m.sqrt(dum * dum - float(not no_sfcten) * pow(cmin, 4)))) / G
                    period[i] = pow(hsig[i] * hsig[i] / (G * m.fabs(ustar[i]) * b_toba * b_toba), 0.333333)
                else:
                    wave_len = PI2 * m.sqrt(sfcten / (G * denwat))

                # Toba's relation between wave height, friction velocity and period
                if ustar[i] > 0.0 and period[i] > 0:
                    hsig[i] = 0.062 * m.sqrt(G * m.fabs(ustar[i]) * pow(period[i], 3))
                else:
                    hsig[i] = 0.0
                # limit the height of breaking waves
                if hsig[i] > 0.137 * wave_len:
                    hsig[i] = 0.137 * wave_len
                if period[i] > 0:
                    orbital_vel[i] = PI2 * 0.5 * hsig[i] / period[i]
                else:
                    orbital_vel[i] = 0.0
            else:
                hsig[i] = 0.0
                period[i] = -9999.9
                dsplcmnt_ht = 0.0
                orbital_vel[i] = 0.0
                wave_len = -9999.0
        # end of block estimating seastate for BVW
    # determine beta-c
    betac[i] = betac_prime
    if z0_prm != 3:
        test_ht = dsplcmnt_ht + 10.0
        psi_u = psi_u_f(test_ht - dsplcmnt_ht, zzu, monin_inv, stab_prm)
        ueff = ustar[i] * (m.log(zzu + 1.0) - psi_u) / KV - float(use_orb_vel) * orbital_vel[i]
        # print("ueff = ", ueff, ustar, zzu, psi_u, dsplcmnt_ht, orbital_vel, i)
        if z0_prm == 0:
            if ueff < 1.8 and not no_capw:
                betac[i] = 0.0
                betag[i] = 0.0
            else:
                betac[i] = 1.0
                betag[i] = 1.0
        if z0_prm == 1:
            if ueff < 1.0:
                betac[i] = 0.0
                betag[i] = 0.0
            else:
                betac[i] = m.tanh(pow(0.4 * (ueff - 1.0), 3)) * betac_prime
                betag[i] = m.tanh(pow(0.3 * (ueff - 1.0), 3)) * betag_prime
        if z0_prm == 4:
            if ueff < 7.0:
                betac[i] = 0.0
            else:
                betac[i] = m.tanh(pow(0.4 * (ueff - 7.0), 3)) * betac_prime
                betag[i] = m.tanh(pow(0.3 * (ueff - 7.0), 3)) * betag_prime
        if betac[i] < 0.1 and (not no_capw):
            betag[i] = betac[i]
        if no_capw:
            betac[i] = 0.0
        betas[i] = 1.0 - betac[i]
        if no_capw:
            betas[i] = 1.0
        psi_u = psi_u_f(betag[i] * hsig[i] - dsplcmnt_ht, hsig[i] * zzu / zref, monin_inv, stab_prm)
        hzzg = hsig[i] * zzu / (zref - dsplcmnt_ht)  # can't be -ve
        if hzzg < 0.00001:
            hzzg = 0.00001
        if ss_prm == 0 or i == 1:
            # local equilibrium requires ww = 1
            wave_age[i] = (m.log(hzzg + 1.0) - float(not neutral) * psi_u) / (
                float((1.0 - i) * ww + float(i)) * KV * ww_eql
            )
            dom_wave_phs_spd[i] = wave_age[i] * m.fabs(ustar[i])

    if not fixedwa:
        # i must equal 0 to reach this point
        # wave age based on sea state observations (other than wave age)
        if m.fabs(ustar[i]) > 0.0:
            wave_age[i] = dom_wave_phs_spd[i] / m.fabs(ustar[i])
        else:
            wave_age[i] = 999.0
    else:
        dom_wave_phs_spd[i] = wave_age[i] * m.fabs(ustar[i])

    # ensure that wa is not horribly unreasonable
    if not (fixedwa or hsig_known or ss_prm == 4) and wave_age[i] > 250.0:
        wave_age[i] = 250.0
        test = -1
    else:
        test = test - test % 2
    if (not fixedwa) and wave_age[i] < 1.08:
        wave_age[i] = 1.08
        test = test - 2
    else:
        test = test - (test % 4 - test % 2)
    if not fixed_dom_wave_phs_spd:
        dom_wave_phs_spd[0] = wave_age[0] * m.fabs(ustar[0])

    # print("in z0_and_waves: betas =", betas, " betac = ", betac, "betag = ", betag)

    # determine zref / z0
    ustar, zzu, betac, betag, betas = z_over_z0(
        use_orb_vel,
        denwat,
        nu,
        sfcten,
        zref,
        dsplcmnt_ht,
        wave_age,
        hsig,
        wave_len,
        i,
        z0_prm,
        stab_prm,
        monin_inv,
        oil_fract_area,
        betac,
        betag,
        betas,
        zzu,
        ustar,
    )
    # ensure that the roughness length is not too absurd */
    if zref > zzu * 25:
        zzu = zref / 25.0
        test = test - 4
    else:
        test = test - (test % 8 - test % 4 - test % 2)

    return test, ustar, zzu, dsplcmnt_ht, hsig, wave_len, orbital_vel, wave_age, dom_wave_phs_spd, betac, betag, betas


def z_over_z0(
    use_orb_vel,
    denwat,
    nu,
    sfcten,
    zref,
    dsplcmnt_ht,
    wave_age,
    hsig,
    wave_len,
    i,
    z0_prm,
    stab_prm,
    monin_inv,
    oil_fract_area,
    betac,
    betag,
    betas,
    zzu,
    ustar,
):
    if not use_orb_vel:
        a_charnock = 0.018  # traditional value of Charnock's constant for wind-wave equilibrium
    else:
        a_charnock = 0.035  # used only when orbital velocity is vector subtracted from the wind shear

    oil_mod = 0.25  # crude guess for oil spills - not natural oil
    ustar_temp = m.sqrt(ustar[0] * ustar[0] + ustar[1] * ustar[1])

    if z0_prm == 0:
        # BVW
        zzu = zref / (
            (
                (betas[i] * 0.1100 * nu / m.fabs(ustar[i]))
                + m.sqrt(
                    pow(betac[i] * B * sfcten / (ustar_temp * m.fabs(ustar[i]) * denwat), 2.0)
                    + pow(betag[i] * 0.48 / wave_age[i] * ustar_temp * m.fabs(ustar[i]) / G, 2.0)
                )
            )
        )
        # * exp( -KV * Uorb / fabs( ustar[i] )
        # print("in z_over_z0: beta = ", betas, betac, betag, " ustar = ", ustar_temp, " zo = ", zref / zzu)

    elif z0_prm == 1:
        # Bourassa (2006)
        zzu = zref / (
            (
                (betas[i] * 0.1100 * nu / m.fabs(ustar[i]))
                + m.sqrt(
                    pow(betac[i] * B * sfcten / (ustar_temp * m.fabs(ustar[i]) * denwat), 2.0)
                    + pow(betag[i] * a_charnock * ustar_temp * m.fabs(ustar[i]) / G, 2.0)
                )
            )
        )
    elif z0_prm == 2:
        # BVW with Taylor & Yelland zog
        zzu = zref / (
            (betas[i] * 0.1100 * nu / m.fabs(ustar[i]))
            + m.sqrt(
                pow(betac[i] * B * sfcten / (ustar_temp * m.fabs(ustar[i]) * denwat), 2.0)
                + pow(betag[i] * hsig[i] * 1200.0 * pow(hsig[i] / wave_len, 4.5), 2.0)
            )
        )
    elif z0_prm == 3:
        # Taylor & Yelland
        zzu = zref / (betag[i] * hsig[i] * 1200.0 * pow(hsig[i] / wave_len, 4.5))
    elif z0_prm == 4:
        # Bourassa (2006) modified for oil slick
        zzu = zref / (
            (
                (betas[i] * 0.1100 * nu / m.fabs(ustar[i]))
                + oil_fract_area
                * oil_mod
                * m.sqrt(
                    pow(betac[i] * B * sfcten / (ustar_temp * m.fabs(ustar[i]) * denwat), 2.0)
                    + pow(betag[i] * a_charnock * ustar_temp * m.fabs(ustar[i]) / G, 2.0)
                )
            )
        )
    elif z0_prm == 5:
        # Smooth surface
        zzu = zref / (0.1100 * nu / m.fabs(ustar[i]))
    elif z0_prm == 6:
        # 2020 version of the oil slick code
        psi_u = psi_u_f(zref - dsplcmnt_ht, zzu, monin_inv, stab_prm)
        ueff = ustar[0] * (m.log(zzu + 1) - psi_u) / KV
        if ueff < 1.0:
            betac[i] = 0.0
            betag[i] = 0.0
        else:
            betac[i] = m.tanh((0.4 * (ueff - 1.0)) ** 3)
            betag[i] = m.tanh((0.2 * (ueff - 1.0)) ** 3)
        if ueff < 7.0:
            betac_oil = 0.0
            betag_oil = 0.0
        else:
            betac_oil = m.tanh((0.4 * (ueff - 7.0)) ** 3)
            betag_oil = m.tanh((0.3 * (ueff - 7.0)) ** 3)
        weighted_betac = (1.0 - oil_fract_area) * betac[i] + oil_fract_area * betac_oil * oil_mod
        weighted_betag = (1.0 - oil_fract_area) * betag[i] + oil_fract_area * betag_oil * oil_mod
        betas[i] = 1.0 - weighted_betag
        zzu = zref / (
            (betas[i] * 0.1100 * nu / m.fabs(ustar[i]))
            + m.sqrt(
                pow(weighted_betac * B * sfcten / (ustar[i] * ustar[i] * denwat), 2.0)
                + pow(weighted_betag * a_charnock * ustar_temp * m.fabs(ustar[i]) / G, 2.0)
            )
        )
    elif z0_prm == 7:
        # do not recalculate. Use the value that was initialized in the main program, based on
        # a value input to the program. The next line is not used.
        pass
    else:
        raise InvalidFluxModelParameterError("Unreasonable value of z0_prm: " + "f{z0_prm}")

    return ustar, zzu, betac, betag, betas


# # # # find_ustar # # # #
def find_ustar(
    fixed_dom_wave_phs_spd,
    fixedwa,
    fixed_hsig,
    fixed_uen,
    neutral,
    no_capw,
    no_sfcten,
    use_dh,
    use_orb_vel,
    betac_prime,
    betag_prime,
    betac,
    betag,
    betas,
    convect_adj,
    wave_age,
    dom_wave_phs_spd,
    hsig,
    period,
    wave_len,
    orbital_vel,
    denwat,
    nu,
    qmixa,
    sfcten,
    wind_vect,
    t_air,
    ww,
    ww_eql,
    ss_prm,
    ss_val,
    zref,
    dsplcmnt_ht,
    z0_prm,
    zzu,
    ustar,
    tstar,
    qstar,
    stab_prm,
    monin_inv,
    oil_fract_area,
):
    # global ustar
    ustar_old = [0.0, 0.0]
    ucount = 0
    # iteratively solve for u* for non-neutral conditions
    # determine the buoyancy flux
    done_outer = FALSE
    while not done_outer and ucount < 30:
        ucount = ucount + 1
        i = 0
        while i < 2:
            ustar_old[i] = ustar[i]
            if fixed_hsig and i == 0:
                hsig_known = TRUE
            else:
                hsig_known = FALSE
            bstar = bstar_f(neutral, qmixa, t_air, tstar, qstar, 0.000001)

            if m.fabs(wind_vect[i]) < 0.001:
                ustar[i] = 0.0
            else:
                (
                    trouble,
                    ustar,
                    zzu,
                    dsplcmnt_ht,
                    hsig,
                    wave_len,
                    orbital_vel,
                    wave_age,
                    dom_wave_phs_spd,
                    betac,
                    betag,
                    betas,
                ) = z0_and_waves(
                    fixed_dom_wave_phs_spd,
                    fixedwa,
                    neutral,
                    hsig_known,
                    stab_prm,
                    no_capw,
                    no_sfcten,
                    use_dh,
                    use_orb_vel,
                    betac_prime,
                    betag_prime,
                    betac,
                    betag,
                    betas,
                    denwat,
                    nu,
                    sfcten,
                    zref,
                    dsplcmnt_ht,
                    monin_inv,
                    wave_age,
                    dom_wave_phs_spd,
                    hsig,
                    period,
                    wave_len,
                    orbital_vel,
                    ww,
                    ww_eql,
                    ss_prm,
                    ss_val,
                    i,
                    z0_prm,
                    zzu,
                    ustar,
                    oil_fract_area,
                )

                psi_u = psi_u_f(zref - dsplcmnt_ht, zzu, monin_inv, stab_prm)
                ustar[i] = ustar_f(
                    fixed_uen,
                    neutral,
                    bstar,
                    ustar,
                    convect_adj,
                    psi_u,
                    monin_inv,
                    wind_vect[i],
                    zzu,
                    10.0,
                )
            i += 1
        if ustar[0] <= 0.000001:
            done_outer = m.fabs(ustar[0] - ustar_old[0]) < 0.0001
        else:
            done_outer = m.fabs((ustar[0] - ustar_old[0]) / ustar[0]) < 0.01
        if ustar[1] <= 0.000001:
            done_outer = done_outer and m.fabs(ustar[1] - ustar_old[1]) < 0.0001
        else:
            done_outer = done_outer and m.fabs((ustar[1] - ustar_old[1]) / ustar[1]) < 0.01

    return dsplcmnt_ht, monin_inv, ustar, zzu, betac, betag, betas


# # # get_renewal_time  # # #
def get_renewal_time(zo_mag, ustar_mag, density, air_specific_heat, alpha, net_long_wave_flux, viscosity):
    rfcr = -2.0e-04
    if density < 50:
        rfcr = -2.0e-01
    if density < 50:
        c_shear = 53.32
    else:
        c_shear = 209  # from Wick thesis
    if density < 50:
        c_conv = 0.80
    else:
        c_conv = 3.13  # from Wick thesis

    # get surface Richardson number

    r_fo = alpha * G * net_long_wave_flux * viscosity / (density * air_specific_heat * pow(ustar_mag, 4.0))
    if r_fo == m.fabs(r_fo):
        r_fo *= -1.0
    # now get time rate **/
    shear_contribution = c_shear * m.sqrt(viscosity * zo_mag / pow(ustar_mag, 3.0))
    conv_contribution = c_conv * m.sqrt(
        viscosity * density * air_specific_heat / (alpha * G * m.fabs(net_long_wave_flux))
    )
    renewal_time = (conv_contribution - shear_contribution) * m.exp(-rfcr / r_fo)
    renewal_time += shear_contribution
    return renewal_time


# # # # # get_solar # # # # #
def get_solar(solar):
    # this just determines the amount of solar radiation is
    #     deposited in the skin.

    if waterType == 0:
        return 0.0
    elif waterType == 1:
        r2 = 0.58
        gamma1 = 0.35
        gamma2 = 23
    elif waterType == 1.1:
        r2 = 0.62
        gamma1 = 0.60
        gamma2 = 20
    elif waterType == 1.2:
        r2 = 0.67
        gamma1 = 1.00
        gamma2 = 17
    elif waterType == 2:
        r2 = 0.77
        gamma1 = 1.50
        gamma2 = 14
    elif waterType == 3:
        r2 = 0.78
        gamma1 = 1.40
        gamma2 = 7.9
    else:
        r2 = 0.78
        gamma1 = 1.40
        gamma2 = 7.9

    z = 0.001
    solar_abs = solar - solar * (r2 * m.exp(-z / gamma1) + (1 - r2) * m.exp(-z / gamma2))
    return solar_abs


# # # # # get_stanton # # # # #
def get_stanton(ustar_mag, renewal_time):
    stanton = m.sqrt(thermalDiff / (ustar_mag * ustar_mag * renewal_time))
    return stanton


# # # # # get_dalton # # # # #
def get_dalton(ustar_mag, renewal_time):
    dalton = m.sqrt(molecularDiff / (ustar_mag * ustar_mag * renewal_time))
    return dalton


# # # # # skin_temp_f # # # # #
def skin_temp_f(zo_mag, ustar_mag, rl, shf, lhf, rs, denair, air_specific_heat):
    den_water = 1021.7
    heat_capacity_water = 4002.0
    kappawater = 1.4e-07  # from Gill page 71
    thermal_expansion = 3196e-07  # from Gill Table A3.1
    viscosity = 1e-06  # this is kinematic viscosity (Gill pg. 75)
    # get solar radiation absorbed over skin
    solar_radn_flux = get_solar(rs)
    # Get Net Heat Loss from ocean
    net_heat_flux = rl + shf + lhf + solar_radn_flux
    # get u* for water
    u_star_water = m.sqrt(denair / den_water) * ustar_mag
    # now get skin temperature
    renewal_time = get_renewal_time(
        zo_mag, u_star_water, denair, air_specific_heat, thermal_expansion, net_heat_flux, viscosity
    )

    skin_temperature_adjustment = m.sqrt(renewal_time)
    skin_temperature_adjustment *= net_heat_flux / (den_water * heat_capacity_water * m.sqrt(kappawater))

    return skin_temperature_adjustment


def validate_parameter(parameter, min_value, max_value, message=None):
    """Validates that a parameter value falls within a specified range.

    Args:
        parameter: The value to validate.
        min_value: The minimum allowed value.
        max_value: The maximum allowed value.
        message: An optional custom error message to raise if validation fails.

    Raises:
        ValueError: If the parameter value is not within the specified range,
            or if the types of the parameters are not compatible.

    Ensures that:
        - min_value and max_value are of the same type.
        - parameter is of the same type as min_value and max_value.
        - parameter's value falls within the range [min_value, max_value].
    """
    if type(min_value) is not type(max_value):
        raise ValueError("min_value and max_value must be of the same type")
    else:
        if type(parameter) is not type(min_value):
            # raise ValueError("Warning: parameter is not of correct type")  # should log a warning
            if message is None:
                message = f"{parameter} is the wrong variable type"
            if "%s" in message:
                logger.warning(message % parameter)
            logger.warning(message)
            if type(min_value) is int:
                parameter = int(float(parameter))
            if type(min_value) is float:
                parameter = float(parameter.item())
        if not min_value <= parameter <= max_value:
            # print(type(parameter), min_value, parameter, max_value)
            if message is None:
                message = f"{parameter} is not within the range [{min_value}, {max_value}]"
            if "%s" in message:
                raise ValueError(message % min_value % ' < ' % parameter % ' < ' % max_value)
            #raise ValueError(message)
    return parameter
