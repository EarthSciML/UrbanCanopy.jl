export CLMUAtmosphere

"""
    CLMUAtmosphere(; name=:CLMUAtmosphere)

Atmospheric interface for the Community Land Model Urban (CLMU) parameterization.
Computes atmospheric density, vapor pressure, and reference height from atmospheric
forcing variables, following Chapter 1 of Oleson et al. (2010).

The CLMU represents the urban surface using a canyon concept (Oke, 1987) consisting
of a canyon floor of width W bordered by two facing buildings of height H. The urban
canyon is divided into five columns: roof, sunlit wall, shaded wall, pervious road,
and impervious road. This component defines the atmospheric forcing interface and
the diagnostic equations that connect atmospheric model variables to the urban
canopy model.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp.
"""
@component function CLMUAtmosphere(; name = :CLMUAtmosphere)

    # Physical constants (Table 1.4)
    @constants begin
        κ_boltz = 1.38065e-23, [description = "Boltzmann constant", unit = u"J/K"]
        N_A = 6.02214e26, [description = "Avogadro's number", unit = u"mol^-1"]
        MW_da = 28.966, [description = "Molecular weight of dry air", unit = u"kg/mol"]
        MW_wv = 18.016, [description = "Molecular weight of water vapor", unit = u"kg/mol"]
    end

    # Derived constants
    @constants begin
        R_gas = 6.02214e26 * 1.38065e-23, [description = "Universal gas constant (N_A * κ)", unit = u"J/(K*mol)"]
        R_da = 6.02214e26 * 1.38065e-23 / 28.966, [description = "Dry air gas constant (R_gas / MW_da)", unit = u"J/(K*kg)"]
        R_wv = 6.02214e26 * 1.38065e-23 / 18.016, [description = "Water vapor gas constant (R_gas / MW_wv)", unit = u"J/(K*kg)"]
        ε_ratio = 18.016 / 28.966, [description = "Ratio of molecular weights MW_wv/MW_da (dimensionless)"]
        one_minus_ε = 1.0 - 18.016 / 28.966, [description = "1 - MW_wv/MW_da (dimensionless)"]
    end

    # Atmospheric forcing inputs (Table 1.1)
    @parameters begin
        z_prime_atm, [description = "Atmospheric reference height from atmospheric model", unit = u"m"]
        u_atm, [description = "Zonal wind at reference height", unit = u"m/s"]
        v_atm, [description = "Meridional wind at reference height", unit = u"m/s"]
        θ_atm, [description = "Atmospheric potential temperature", unit = u"K"]
        q_atm, [description = "Specific humidity at reference height", unit = u"kg/kg"]
        P_atm, [description = "Atmospheric pressure at reference height", unit = u"Pa"]
        T_atm, [description = "Atmospheric temperature at reference height", unit = u"K"]
        L_atm_down, [description = "Incident longwave radiation", unit = u"W/m^2"]
        q_rain, [description = "Liquid precipitation rate", unit = u"kg/(m^2*s)"]
        q_sno, [description = "Solid precipitation rate", unit = u"kg/(m^2*s)"]
        S_atm_dir_vis, [description = "Incident direct beam visible solar radiation", unit = u"W/m^2"]
        S_atm_dir_nir, [description = "Incident direct beam near-infrared solar radiation", unit = u"W/m^2"]
        S_atm_dif_vis, [description = "Incident diffuse visible solar radiation", unit = u"W/m^2"]
        S_atm_dif_nir, [description = "Incident diffuse near-infrared solar radiation", unit = u"W/m^2"]
    end

    # Surface parameters for reference height computation (Table 1.3, Footnote 1 p.16)
    @parameters begin
        z_0, [description = "Roughness length", unit = u"m"]
        z_d, [description = "Displacement height", unit = u"m"]
    end

    # Diagnostic variables
    @variables begin
        e_atm(t), [description = "Atmospheric vapor pressure", unit = u"Pa"]
        ρ_atm(t), [description = "Atmospheric air density", unit = u"kg/m^3"]
        z_atm(t), [description = "Reference height for flux computations", unit = u"m"]
    end

    eqs = [
        # Atmospheric vapor pressure (Chapter 1, p.17)
        # e_atm = q_atm * P_atm / (0.622 + 0.378 * q_atm)
        # where 0.622 = MW_wv/MW_da and 0.378 = 1 - MW_wv/MW_da
        e_atm ~ q_atm * P_atm / (ε_ratio + one_minus_ε * q_atm),

        # Air density (Chapter 1, p.17)
        # ρ_atm = (P_atm - 0.378 * e_atm) / (R_da * T_atm)
        ρ_atm ~ (P_atm - one_minus_ε * e_atm) / (R_da * T_atm),

        # Reference height for flux computations (Chapter 1, Footnote 1 p.16)
        # z_atm = z'_atm + z_0 + z_d
        z_atm ~ z_prime_atm + z_0 + z_d,
    ]

    return System(eqs, t; name)
end
