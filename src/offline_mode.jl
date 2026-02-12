export OfflineModeForcing

"""
    OfflineModeForcing(; name=:OfflineModeForcing)

Offline mode atmospheric forcing processor for the Community Land Model Urban (CLMU).
Computes derived atmospheric forcing variables from raw observational data for use in
uncoupled (offline) simulations, following Chapter 6 of Oleson et al. (2010).

This component takes total incident solar radiation, cosine of solar zenith angle,
atmospheric temperature, pressure, specific humidity, wind speed, and precipitation
rate as inputs and computes:

- Solar radiation partitioned into direct/diffuse and visible/near-infrared components (Eqs. 6.2-6.8)
- Atmospheric downwelling longwave radiation from the Idso (1981) formula (Eqs. 6.9-6.10)
- Rain and snow precipitation rates from temperature-dependent partitioning (Eqs. 6.11-6.13)
- Wind components, potential temperature, and reference height (p.149)

**Note**: Equation 6.1 describes temporal downscaling of 6-hourly solar radiation to
model time steps using cosine of solar zenith angle weighting. This discrete-time
algorithm is a data preprocessing step and is not implemented here; instead, `S_atm`
is taken as the already-interpolated total solar radiation at the model time step.

**Note**: Equations 6.14-6.15 describe alternative methods for computing specific humidity
from relative humidity or dew point temperature. These are optional data preprocessing
steps and are not implemented in this component.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp.
"""
@component function OfflineModeForcing(; name = :OfflineModeForcing)

    # Physical constants
    @constants begin
        σ_SB = 5.67e-8, [description = "Stefan-Boltzmann constant (Table 1.4)", unit = u"W/(m^2*K^4)"]
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
    end

    # Molecular weight ratio constants (same as CLMUAtmosphere, Table 1.4)
    @constants begin
        MW_da = 0.028966, [description = "Molecular weight of dry air (Table 1.4)", unit = u"kg/mol"]
        MW_wv = 0.018016, [description = "Molecular weight of water vapor (Table 1.4)", unit = u"kg/mol"]
        ε_ratio = MW_wv / MW_da, [description = "Ratio of molecular weights MW_wv/MW_da = 0.622 (dimensionless)"]
        one_minus_ε = 1.0 - MW_wv / MW_da, [description = "1 - MW_wv/MW_da = 0.378 (dimensionless)"]
    end

    # Solar radiation partitioning constants (Eq. 6.6)
    @constants begin
        α_vis = 0.5, [description = "Ratio of visible to total incident solar radiation (Eq. 6.6) (dimensionless)"]
    end

    # Polynomial coefficients for visible direct fraction (Eq. 6.7)
    @constants begin
        a_0 = 0.17639, [description = "Visible direct fraction polynomial coefficient a_0 (Eq. 6.7) (dimensionless)"]
        a_1 = 0.00380, [description = "Visible direct fraction polynomial coefficient a_1 (Eq. 6.7)", unit = u"m^2/W"]
        a_2 = -9.0039e-6, [description = "Visible direct fraction polynomial coefficient a_2 (Eq. 6.7)", unit = u"m^4/W^2"]
        a_3 = 8.1351e-9, [description = "Visible direct fraction polynomial coefficient a_3 (Eq. 6.7)", unit = u"m^6/W^3"]
    end

    # Polynomial coefficients for NIR direct fraction (Eq. 6.8)
    @constants begin
        b_0 = 0.29548, [description = "NIR direct fraction polynomial coefficient b_0 (Eq. 6.8) (dimensionless)"]
        b_1 = 0.00504, [description = "NIR direct fraction polynomial coefficient b_1 (Eq. 6.8)", unit = u"m^2/W"]
        b_2 = -1.4957e-5, [description = "NIR direct fraction polynomial coefficient b_2 (Eq. 6.8)", unit = u"m^4/W^2"]
        b_3 = 1.4881e-8, [description = "NIR direct fraction polynomial coefficient b_3 (Eq. 6.8)", unit = u"m^6/W^3"]
    end

    # Idso (1981) longwave radiation constants (Eq. 6.9)
    @constants begin
        idso_base = 0.70, [description = "Idso base emissivity factor (Eq. 6.9) (dimensionless)"]
        idso_coeff = 5.95e-7, [description = "Idso vapor pressure coefficient 5.95e-5 * 0.01 (Eq. 6.9)", unit = u"Pa^-1"]
        idso_T_ref = 1500.0, [description = "Idso temperature reference (Eq. 6.9)", unit = u"K"]
    end

    # Precipitation phase partitioning constant (Eq. 6.13)
    @constants begin
        f_P_slope = 0.5, [description = "Rain fraction slope (Eq. 6.13)", unit = u"K^-1"]
    end

    # Reference constants for clamping (dimensionless bounds)
    @constants begin
        R_min = 0.01, [description = "Minimum direct fraction (Eqs. 6.7-6.8) (dimensionless)"]
        R_max = 0.99, [description = "Maximum direct fraction (Eqs. 6.7-6.8) (dimensionless)"]
        zero_dimless = 0.0, [description = "Zero (dimensionless)"]
        one_dimless = 1.0, [description = "One (dimensionless)"]
    end

    # Reference constants for units
    @constants begin
        sqrt2_inv = 1.0 / sqrt(2.0), [description = "1/√2 for wind decomposition (dimensionless)"]
        z_prime_atm_val = 30.0, [description = "Atmospheric reference height (p.149)", unit = u"m"]
        zero_Wm2 = 0.0, [description = "Zero radiation", unit = u"W/m^2"]
        zero_precip = 0.0, [description = "Zero precipitation rate", unit = u"kg/(m^2*s)"]
    end

    # Raw atmospheric forcing inputs
    @parameters begin
        S_atm, [description = "Total incident solar radiation at model time step (Eq. 6.1 output)", unit = u"W/m^2"]
        T_atm, [description = "Atmospheric temperature", unit = u"K"]
        P_atm, [description = "Atmospheric pressure", unit = u"Pa"]
        q_atm, [description = "Atmospheric specific humidity", unit = u"kg/kg"]
        W_atm, [description = "Atmospheric wind speed", unit = u"m/s"]
        P_precip, [description = "Total precipitation rate", unit = u"kg/(m^2*s)"]
    end

    # Derived variables - solar radiation partitioning
    @variables begin
        R_vis_raw(t), [description = "Unclamped visible direct fraction (Eq. 6.7) (dimensionless)"]
        R_nir_raw(t), [description = "Unclamped NIR direct fraction (Eq. 6.8) (dimensionless)"]
        R_vis(t), [description = "Visible direct fraction, clamped to [0.01, 0.99] (Eq. 6.7) (dimensionless)"]
        R_nir(t), [description = "NIR direct fraction, clamped to [0.01, 0.99] (Eq. 6.8) (dimensionless)"]
        S_atm_dir_vis(t), [description = "Direct beam visible solar radiation (Eq. 6.2)", unit = u"W/m^2"]
        S_atm_dir_nir(t), [description = "Direct beam near-infrared solar radiation (Eq. 6.3)", unit = u"W/m^2"]
        S_atm_dif_vis(t), [description = "Diffuse visible solar radiation (Eq. 6.4)", unit = u"W/m^2"]
        S_atm_dif_nir(t), [description = "Diffuse near-infrared solar radiation (Eq. 6.5)", unit = u"W/m^2"]
    end

    # Derived variables - longwave radiation
    @variables begin
        e_atm(t), [description = "Atmospheric vapor pressure (Eq. 6.10)", unit = u"Pa"]
        L_atm_down(t), [description = "Atmospheric downwelling longwave radiation (Eq. 6.9)", unit = u"W/m^2"]
    end

    # Derived variables - precipitation phase
    @variables begin
        f_P_raw(t), [description = "Unclamped rain fraction (Eq. 6.13) (dimensionless)"]
        f_P(t), [description = "Rain fraction, clamped to [0, 1] (Eq. 6.13) (dimensionless)"]
        q_rain(t), [description = "Liquid precipitation rate (Eq. 6.11)", unit = u"kg/(m^2*s)"]
        q_snow(t), [description = "Solid precipitation rate (Eq. 6.12)", unit = u"kg/(m^2*s)"]
    end

    # Derived variables - wind and other
    @variables begin
        u_atm(t), [description = "Zonal wind component (p.149)", unit = u"m/s"]
        v_atm(t), [description = "Meridional wind component (p.149)", unit = u"m/s"]
        θ_atm(t), [description = "Potential temperature (p.149)", unit = u"K"]
        z_prime_atm(t), [description = "Atmospheric reference height (p.149)", unit = u"m"]
    end

    eqs = [
        # ---- Solar radiation partitioning (Eqs. 6.2-6.8) ----

        # Eq. 6.7 - Ratio of direct to total in visible (polynomial fit)
        # R_vis = a_0 + a_1*(α*S_atm) + a_2*(α*S_atm)^2 + a_3*(α*S_atm)^3
        R_vis_raw ~ a_0 + a_1 * (α_vis * S_atm) + a_2 * (α_vis * S_atm)^2 + a_3 * (α_vis * S_atm)^3,

        # Eq. 6.8 - Ratio of direct to total in NIR (polynomial fit)
        # R_nir = b_0 + b_1*(1-α)*S_atm + b_2*[(1-α)*S_atm]^2 + b_3*[(1-α)*S_atm]^3
        R_nir_raw ~ b_0 + b_1 * ((one_dimless - α_vis) * S_atm) +
                    b_2 * ((one_dimless - α_vis) * S_atm)^2 +
                    b_3 * ((one_dimless - α_vis) * S_atm)^3,

        # Clamp R_vis and R_nir to [0.01, 0.99]
        R_vis ~ max(R_min, min(R_max, R_vis_raw)),
        R_nir ~ max(R_min, min(R_max, R_nir_raw)),

        # Eq. 6.2 - Direct beam visible solar radiation
        # S_atm↓^μ_vis = R_vis * (α * S_atm)
        S_atm_dir_vis ~ R_vis * (α_vis * S_atm),

        # Eq. 6.3 - Direct beam NIR solar radiation
        # S_atm↓^μ_nir = R_nir * [(1-α) * S_atm]
        S_atm_dir_nir ~ R_nir * ((one_dimless - α_vis) * S_atm),

        # Eq. 6.4 - Diffuse visible solar radiation
        # S_atm↓_vis = (1 - R_vis) * (α * S_atm)
        S_atm_dif_vis ~ (one_dimless - R_vis) * (α_vis * S_atm),

        # Eq. 6.5 - Diffuse NIR solar radiation
        # S_atm↓_nir = (1 - R_nir) * [(1-α) * S_atm]
        S_atm_dif_nir ~ (one_dimless - R_nir) * ((one_dimless - α_vis) * S_atm),

        # ---- Longwave radiation (Eqs. 6.9-6.10) ----

        # Eq. 6.10 - Atmospheric vapor pressure from specific humidity
        # e_atm = P_atm * q_atm / (0.622 + 0.378 * q_atm)
        e_atm ~ P_atm * q_atm / (ε_ratio + one_minus_ε * q_atm),

        # Eq. 6.9 - Atmospheric downwelling longwave radiation (Idso 1981)
        # L_atm↓ = [0.70 + 5.95×10⁻⁵ × 0.01 × e_atm × exp(1500/T_atm)] × σ × T_atm⁴
        L_atm_down ~ (idso_base + idso_coeff * e_atm * exp(idso_T_ref / T_atm)) * σ_SB * T_atm^4,

        # ---- Precipitation phase partitioning (Eqs. 6.11-6.13) ----

        # Eq. 6.13 - Rain fraction: f_P = clamp(0.5 * (T_atm - T_f), 0, 1)
        f_P_raw ~ f_P_slope * (T_atm - T_f),
        f_P ~ max(zero_dimless, min(one_dimless, f_P_raw)),

        # Eq. 6.11 - Rain rate
        q_rain ~ P_precip * f_P,

        # Eq. 6.12 - Snow rate
        q_snow ~ P_precip * (one_dimless - f_P),

        # ---- Wind components and other derived quantities (p.149) ----

        # Wind decomposition: u_atm = v_atm = W_atm / √2
        u_atm ~ sqrt2_inv * W_atm,
        v_atm ~ sqrt2_inv * W_atm,

        # Potential temperature set equal to atmospheric temperature
        θ_atm ~ T_atm,

        # Atmospheric reference height set to 30 m
        z_prime_atm ~ z_prime_atm_val,
    ]

    return System(eqs, t; name)
end
