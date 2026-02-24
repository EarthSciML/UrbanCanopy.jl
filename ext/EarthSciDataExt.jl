module EarthSciDataExt

using EarthSciMLBase
using EarthSciMLBase: param_to_var, ConnectorSystem, CoupleType
using EarthSciData: GEOSFPCoupler
using UrbanCanopy: UrbanCanopyModelCoupler
using ModelingToolkit: t, ParentScope
using ModelingToolkit: @constants
using DynamicQuantities: @u_str

"""
    EarthSciMLBase.couple2(ucm::UrbanCanopyModelCoupler, g::GEOSFPCoupler)

Couple the [`UrbanCanopyModel`](@ref) to GEOS-FP meteorological reanalysis data.

Maps GEOS-FP surface fields to the urban canopy model's meteorological forcing
parameters. The Monin-Obukhov stability parameter ζ is computed from GEOS-FP's
friction velocity (USTAR) and sensible heat flux (HFLUX). Solar zenith angle
is computed from latitude, longitude, and time using astronomical formulas.
"""
function EarthSciMLBase.couple2(ucm::UrbanCanopyModelCoupler, g::GEOSFPCoupler)
    u, m = ucm.sys, g.sys

    # Convert top-level met forcing params to variables so coupling equations can drive them
    u = param_to_var(u, :S_atm, :T_atm, :P_atm, :q_atm, :W_atm, :P_precip, :ζ_in, :μ_zen)

    # ===== Constants for ζ (Monin-Obukhov stability) computation =====
    @constants begin
        C_p_cpl = 1004.64, [unit = u"J/(kg*K)", description = "Specific heat of dry air"]
        R_d_cpl = 287.058, [unit = u"J/(kg*K)", description = "Gas constant for dry air"]
        κ_cpl = 0.4, [description = "von Karman constant (dimensionless)"]
        g_cpl = 9.80616, [unit = u"m/s^2", description = "Gravitational acceleration"]
        α_cpl = 4.43, [description = "Displacement height coefficient (dimensionless)"]
        z_ref_cpl = 30.0, [unit = u"m", description = "Reference height"]
        ε_flux = 1.0e-4, [unit = u"W/m^2", description = "Min heat flux magnitude"]
        one_m_cpl = 1.0, [unit = u"m", description = "Unit length"]
        one_s_precip = 1.0, [unit = u"s", description = "Unit time for PRECTOT unit correction (GEOS-FP metadata labels precipitation rate as kg/m²/s² instead of kg/m²/s)"]
    end

    # Compute Monin-Obukhov length from GEOS-FP USTAR and HFLUX
    # L = -P * c_p * u*^3 / (R_d * κ * g * H)
    HFLUX_safe = ifelse(abs(m.A1₊HFLUX) > ε_flux, m.A1₊HFLUX, ε_flux)
    L_MO_expr = -m.I3₊PS * C_p_cpl * (m.A1₊USTAR / one_m_cpl * one_m_cpl)^3 /
        (R_d_cpl * κ_cpl * g_cpl * HFLUX_safe)

    # Compute d_canopy from geometry (Eq. 3.55, duplicated here since d is internal to flux)
    # d = H * (1 + α^(-W_roof) * (W_roof - 1))
    d_expr = u.H_canyon * (1 + (α_cpl^(-u.W_roof)) * (u.W_roof - 1))

    # Compute ζ = (z_ref - d) / L
    ζ_expr = (z_ref_cpl - d_expr) / L_MO_expr

    # ===== Solar zenith angle from lat/lon/time =====
    # cos(SZA) = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(hour_angle)
    τ = ParentScope(m.t_ref)  # absolute time reference [s]
    φ = ParentScope(m.lat)    # latitude [rad]
    λ = ParentScope(m.lon)    # longitude [rad]

    @constants begin
        sec_per_day = 86400.0, [unit = u"s", description = "Seconds per day"]
        days_per_year = 365.25, [description = "Days per year (dimensionless)"]
        obliquity = 0.4091, [unit = u"rad", description = "Earth obliquity 23.44 deg"]
        one_s = 1.0, [unit = u"s", description = "Unit time"]
        one_rad_cpl = 1.0, [unit = u"rad", description = "Unit angle"]
    end

    t_abs = (τ + t) / one_s  # dimensionless time in seconds
    # Day of year (approximate): fractional day within year
    day_frac = t_abs / (sec_per_day / one_s) # day number
    day_of_year = mod(day_frac, days_per_year)
    # Solar declination (simplified)
    dec = -obliquity * cos(2π * (day_of_year + 10) / days_per_year)
    # Hour angle from UTC time and longitude
    hour_angle = 2π * mod(day_frac, 1.0) - π + λ / one_rad_cpl
    # Cosine of solar zenith angle
    csza = sin(φ / one_rad_cpl) * sin(dec / one_rad_cpl) +
        cos(φ / one_rad_cpl) * cos(dec / one_rad_cpl) * cos(hour_angle)
    csza_clamped = max(csza, 0.01)  # avoid exactly zero or negative
    μ_expr = acos(csza_clamped) * one_rad_cpl

    return ConnectorSystem(
        [
            u.S_atm ~ m.A1₊SWGDN,
            u.T_atm ~ m.A1₊T2M,
            u.P_atm ~ m.I3₊PS,
            u.q_atm ~ m.A1₊QV2M,
            u.W_atm ~ sqrt(m.A1₊U10M^2 + m.A1₊V10M^2),
            u.P_precip ~ m.A1₊PRECTOT * one_s_precip,  # GEOS-FP PRECTOT metadata has kg/m²/s² instead of kg/m²/s
            u.ζ_in ~ ζ_expr,
            u.μ_zen ~ μ_expr,
        ],
        u,
        m
    )
end

end
