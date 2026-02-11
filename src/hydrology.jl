export SnowDensity, SnowIceContent, SnowWaterContent, SnowCompaction,
    SnowLayerCombination, SoilHydraulicProperties, SurfaceRunoffInfiltration,
    SoilWaterFlux, SoilWaterEquilibrium, RichardsEquation, GroundwaterDrainage,
    WaterTableDepth, AquiferWaterBalance, SnowCappingRunoff, SurfaceLayerUpdate,
    PerviousRoadWaterBalance, ImperviousWaterBalance

"""
    SnowDensity(; name=:SnowDensity)

Computes the bulk density of newly fallen snow as a function of atmospheric
temperature, following Anderson (1976) as parameterized in Eq. 5.10 of
Oleson et al. (2010).

The density depends on the atmospheric temperature relative to freezing:
- When T_atm > T_f + 2 K: ρ_sno = 50 + 1.7(17)^1.5
- When T_f - 15 < T_atm ≤ T_f + 2: ρ_sno = 50 + 1.7(T_atm - T_f + 15)^1.5
- When T_atm ≤ T_f - 15: ρ_sno = 50

**Reference**: Oleson et al. (2010), Chapter 5, Eq. 5.10, p. 115.
"""
@component function SnowDensity(; name = :SnowDensity)

    @constants begin
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
        one_K = 1.0, [description = "Unit temperature for nondimensionalization", unit = u"K"]
        ρ_base = 50.0, [description = "Base snow density", unit = u"kg/m^3"]
        ρ_coeff = 1.7, [description = "Snow density coefficient", unit = u"kg/m^3"]
        ρ_warm = 50.0 + 1.7 * 17.0^1.5, [description = "Snow density for warm conditions (T_atm > T_f+2)", unit = u"kg/m^3"]
        T_offset = 15.0, [description = "Temperature offset for density formula", unit = u"K"]
        T_thresh_high = 2.0, [description = "Upper temperature threshold", unit = u"K"]
        T_thresh_low = -15.0, [description = "Lower temperature threshold", unit = u"K"]
    end

    @parameters begin
        T_atm, [description = "Atmospheric temperature", unit = u"K"]
    end

    @variables begin
        ρ_sno(t), [description = "Bulk density of newly fallen snow (Eq. 5.10)", unit = u"kg/m^3"]
    end

    eqs = [
        # Eq. 5.10: Piecewise snow density as function of temperature
        # ρ_sno = 50 + 1.7*(T_diff/1K)^1.5 [kg/m³]
        ρ_sno ~ ifelse(T_atm > T_f + T_thresh_high,
            ρ_warm,
            ifelse(T_atm > T_f + T_thresh_low,
                ρ_base + ρ_coeff * ((T_atm - T_f + T_offset) / one_K)^1.5,
                ρ_base)),
    ]

    return System(eqs, t; name)
end


"""
    SnowIceContent(; name=:SnowIceContent)

Computes the conservation of ice mass in snow layers and new snowfall accumulation,
following Section 5.1.1 of Oleson et al. (2010).

This component handles:
- New snow depth rate from solid precipitation and snow density (derived from Eq. 5.9)
- Top layer ice content update from precipitation and frost/sublimation (Eqs. 5.12, 5.14)

The paper's Eq. 5.9 is Δz_sno = q_grnd,ice·Δt / ρ_sno, which is a per-timestep increment.
For continuous modeling we express this as a rate: dz_sno/dt = q_grnd,ice / ρ_sno [m/s].

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.1.1, Eqs. 5.6–5.14, pp. 114–115.
"""
@component function SnowIceContent(; name = :SnowIceContent)

    @parameters begin
        q_grnd_ice, [description = "Solid precipitation reaching the surface (Eq. 5.5)", unit = u"kg/(m^2*s)"]
        q_frost, [description = "Frost rate", unit = u"kg/(m^2*s)"]
        q_subl, [description = "Sublimation rate", unit = u"kg/(m^2*s)"]
        ρ_sno, [description = "Bulk density of newly fallen snow", unit = u"kg/m^3"]
    end

    @variables begin
        dz_sno_dt(t), [description = "Rate of new snow depth accumulation (from Eq. 5.9)", unit = u"m/s"]
        q_ice_top(t), [description = "Rate of ice accumulation for top snow layer (Eq. 5.7)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Derived from Eq. 5.9: dz_sno/dt = q_grnd_ice / ρ_sno [m/s]
        # (Eq. 5.9 is Δz_sno = q_grnd_ice·Δt / ρ_sno; dividing both sides by Δt gives the rate)
        dz_sno_dt ~ q_grnd_ice / ρ_sno,

        # Eq. 5.7: Ice accumulation rate for top layer
        q_ice_top ~ q_grnd_ice + (q_frost - q_subl),
    ]

    return System(eqs, t; name)
end


"""
    SnowWaterContent(; name=:SnowWaterContent)

Computes the flow of liquid water through snow layers, following Section 5.1.2
of Oleson et al. (2010).

The water flow between layers (Eq. 5.18) depends on the volumetric water and ice
contents, with an irreducible water saturation S_r = 0.033 representing capillary
retention. The flow is limited by the effective porosity of the underlying layer
(Eq. 5.21).

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.1.2, Eqs. 5.15–5.22, pp. 116–117.
"""
@component function SnowWaterContent(; name = :SnowWaterContent)

    @constants begin
        ρ_liq = 1000.0, [description = "Density of liquid water (Table 1.4)", unit = u"kg/m^3"]
        ρ_ice = 917.0, [description = "Density of ice (Table 1.4)", unit = u"kg/m^3"]
        S_r = 0.033, [description = "Irreducible water saturation (dimensionless)"]
        zero_kgpm2s = 0.0, [description = "Zero flux for comparison", unit = u"kg/(m^2*s)"]
        one_s = 1.0, [description = "Unit time for rate conversion", unit = u"s"]
    end

    @parameters begin
        w_liq, [description = "Mass of liquid water in snow layer", unit = u"kg/m^2"]
        w_ice, [description = "Mass of ice in snow layer", unit = u"kg/m^2"]
        Δz, [description = "Snow layer thickness", unit = u"m"]
    end

    @variables begin
        θ_ice(t), [description = "Volumetric ice content in snow layer (Eq. 5.19)", unit = u"m^3/m^3"]
        θ_liq(t), [description = "Volumetric liquid water content in snow layer (Eq. 5.20)", unit = u"m^3/m^3"]
        q_liq_out(t), [description = "Water flow out of snow layer (Eq. 5.18)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.19: Volumetric ice content
        θ_ice ~ min(w_ice / (Δz * ρ_ice), 1.0),

        # Eq. 5.20: Volumetric liquid water content (limited by 1 - θ_ice)
        θ_liq ~ min(w_liq / (Δz * ρ_liq), 1.0 - θ_ice),

        # Eq. 5.18: Water flow out of layer (≥ 0)
        # q_liq,i = ρ_liq[θ_liq,i - S_r(1-θ_ice,i)]Δz_i / Δt
        # Using Δt=1s since this gives the instantaneous rate in kg/(m²·s)
        q_liq_out ~ max(ρ_liq * (θ_liq - S_r * (1.0 - θ_ice)) * Δz / one_s, zero_kgpm2s),
    ]

    return System(eqs, t; name)
end


"""
    SnowCompaction(; name=:SnowCompaction)

Computes the total fractional compaction rate of snow layers, following Section 5.1.4
of Oleson et al. (2010).

Snow compaction includes three processes:
1. Destructive metamorphism (Eq. 5.26) - temperature dependent crystal breakdown
2. Overburden pressure (Eq. 5.28) - compression from weight of overlying snow
3. Melt-induced compaction (Eq. 5.31) - changes due to melt-freeze cycles

Compaction is not allowed if the layer is saturated (Eq. 5.25) or if ice content
is below a minimum value (w_ice,i ≤ 0.1 kg/m²).

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.1.4, Eqs. 5.24–5.33, pp. 118–120.
"""
@component function SnowCompaction(; name = :SnowCompaction)

    @constants begin
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
        c_3 = 2.777e-6, [description = "Fractional compaction rate at T_f (Eq. 5.26)", unit = u"s^-1"]
        c_4 = 0.04, [description = "Temperature coefficient for metamorphism (Eq. 5.26)", unit = u"K^-1"]
        η_0 = 9.0e5, [description = "Reference viscosity (Eq. 5.29)", unit = u"kg*s/m^2"]
        c_5 = 0.08, [description = "Viscosity temperature coefficient (Eq. 5.29)", unit = u"K^-1"]
        c_6 = 0.023, [description = "Viscosity density coefficient (Eq. 5.29)", unit = u"m^3/kg"]
        hundred_kgpm3 = 100.0, [description = "Reference density 100 kg/m^3", unit = u"kg/m^3"]
        thresh_liq_density = 0.01, [description = "Liquid water density threshold", unit = u"kg/m^3"]
        coeff_0046 = 0.046, [description = "Density coefficient (Eq. 5.27)", unit = u"m^3/kg"]
        one_per_s = 1.0, [description = "Unit rate for nondimensionalization", unit = u"s^-1"]
    end

    @parameters begin
        T_i, [description = "Snow layer temperature", unit = u"K"]
        w_ice, [description = "Mass of ice in layer", unit = u"kg/m^2"]
        w_liq, [description = "Mass of liquid water in layer", unit = u"kg/m^2"]
        Δz, [description = "Snow layer thickness", unit = u"m"]
        P_s, [description = "Snow load pressure on this layer (Eq. 5.30)", unit = u"kg/m^2"]
        f_ice_n, [description = "Ice fraction at time n (dimensionless)"]
        f_ice_n1, [description = "Ice fraction at time n+1 (dimensionless)"]
    end

    @variables begin
        C_R1(t), [description = "Destructive metamorphism compaction rate (Eq. 5.26)", unit = u"s^-1"]
        C_R2(t), [description = "Overburden compaction rate (Eq. 5.28)", unit = u"s^-1"]
        C_R3(t), [description = "Melt compaction rate (Eq. 5.31)", unit = u"s^-1"]
        C_R(t), [description = "Total fractional compaction rate (Eq. 5.24)", unit = u"s^-1"]
        η(t), [description = "Viscosity coefficient (Eq. 5.29)", unit = u"kg*s/m^2"]
        c_1(t), [description = "Density-dependent coefficient (Eq. 5.27) (dimensionless)"]
        c_2(t), [description = "Liquid-dependent coefficient (Eq. 5.27) (dimensionless)"]
    end

    eqs = [
        # Eq. 5.27: Density-dependent coefficient c_1
        c_1 ~ ifelse(w_ice / Δz > hundred_kgpm3,
            exp(-coeff_0046 * (w_ice / Δz - hundred_kgpm3)),
            1.0),

        # Eq. 5.27: Liquid-dependent coefficient c_2
        c_2 ~ ifelse(w_liq / Δz > thresh_liq_density, 2.0, 1.0),

        # Eq. 5.26: Destructive metamorphism compaction rate
        C_R1 ~ -c_3 * c_1 * c_2 * exp(-c_4 * (T_f - T_i)),

        # Eq. 5.29: Viscosity coefficient
        η ~ η_0 * exp(c_5 * (T_f - T_i) + c_6 * w_ice / Δz),

        # Eq. 5.28: Overburden compaction rate
        C_R2 ~ -P_s / η,

        # Eq. 5.31: Melt compaction rate
        # C_R3 = -(1/Δt) max(0, (f_ice^n - f_ice^{n+1}) / f_ice^n)
        # Using unit rate (1/s) since Δt factored into the f_ice difference
        C_R3 ~ -one_per_s * max(0.0, (f_ice_n - f_ice_n1) / max(f_ice_n, 1e-10)),

        # Eq. 5.24: Total compaction rate
        C_R ~ C_R1 + C_R2 + C_R3,
    ]

    return System(eqs, t; name)
end


"""
    SnowLayerCombination(; name=:SnowLayerCombination)

Computes the combined properties when two snow layers are merged, following
Section 5.1.5.1 of Oleson et al. (2010).

When two snow layers are combined, their thicknesses, water contents, and
temperatures are merged according to conservation principles (Eqs. 5.38–5.42).

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.1.5.1, Eqs. 5.38–5.42, pp. 122.
"""
@component function SnowLayerCombination(; name = :SnowLayerCombination)

    @constants begin
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
        L_f = 3.337e5, [description = "Latent heat of fusion (Table 1.4)", unit = u"J/kg"]
        C_ice = 2117.27, [description = "Specific heat of ice (Table 1.4)", unit = u"J/(kg*K)"]
        C_liq = 4188.0, [description = "Specific heat of liquid water (Table 1.4)", unit = u"J/(kg*K)"]
    end

    @parameters begin
        Δz_1, [description = "Thickness of layer 1", unit = u"m"]
        Δz_2, [description = "Thickness of layer 2", unit = u"m"]
        w_liq_1, [description = "Liquid water in layer 1", unit = u"kg/m^2"]
        w_liq_2, [description = "Liquid water in layer 2", unit = u"kg/m^2"]
        w_ice_1, [description = "Ice content in layer 1", unit = u"kg/m^2"]
        w_ice_2, [description = "Ice content in layer 2", unit = u"kg/m^2"]
        T_1, [description = "Temperature of layer 1", unit = u"K"]
        T_2, [description = "Temperature of layer 2", unit = u"K"]
    end

    @variables begin
        Δz_c(t), [description = "Combined thickness (Eq. 5.38)", unit = u"m"]
        w_liq_c(t), [description = "Combined liquid water (Eq. 5.39)", unit = u"kg/m^2"]
        w_ice_c(t), [description = "Combined ice content (Eq. 5.40)", unit = u"kg/m^2"]
        T_c(t), [description = "Combined temperature (Eq. 5.41)", unit = u"K"]
        h_1(t), [description = "Enthalpy of layer 1 (Eq. 5.42)", unit = u"J/m^2"]
        h_2(t), [description = "Enthalpy of layer 2 (Eq. 5.42)", unit = u"J/m^2"]
    end

    eqs = [
        # Eq. 5.38: Combined thickness
        Δz_c ~ Δz_1 + Δz_2,

        # Eq. 5.39: Combined liquid water
        w_liq_c ~ w_liq_1 + w_liq_2,

        # Eq. 5.40: Combined ice content
        w_ice_c ~ w_ice_1 + w_ice_2,

        # Eq. 5.42: Layer enthalpy
        h_1 ~ (C_ice * w_ice_1 + C_liq * w_liq_1) * (T_1 - T_f) + L_f * w_liq_1,
        h_2 ~ (C_ice * w_ice_2 + C_liq * w_liq_2) * (T_2 - T_f) + L_f * w_liq_2,

        # Eq. 5.41: Combined temperature
        T_c ~ T_f + (h_1 + h_2 - L_f * w_liq_c) / (C_ice * w_ice_c + C_liq * w_liq_c),
    ]

    return System(eqs, t; name)
end


"""
    SoilHydraulicProperties(; name=:SoilHydraulicProperties)

Computes the hydraulic properties of soil layers for the pervious road, following
Section 5.3.1 of Oleson et al. (2010).

Properties computed based on Clapp and Hornberger (1978) and Cosby et al. (1984):
- Saturated hydraulic conductivity k_sat (Eq. 5.70)
- Porosity θ_sat (Eq. 5.71)
- Clapp-Hornberger exponent B (Eq. 5.72)
- Soil matric potential ψ (Eq. 5.73)
- Saturated soil matric potential ψ_sat (Eq. 5.74)

Urban soils are assumed to have no organic matter, so simplified forms are used.
All values converted to SI units (m, m/s) from the paper's mm units.

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.3.1, Eqs. 5.69–5.74, pp. 131.
"""
@component function SoilHydraulicProperties(; name = :SoilHydraulicProperties)

    @constants begin
        # Paper gives k_coeff = 0.0070556 mm/s; multiply by 0.001 to get m/s
        k_coeff = 0.0070556e-3, [description = "Coefficient for saturated hydraulic conductivity (Eq. 5.70), converted from mm/s to m/s", unit = u"m/s"]
        # Paper gives ψ_coeff = -10.0 mm; multiply by 0.001 to get m
        ψ_coeff = -10.0e-3, [description = "Coefficient for saturated matric potential (Eq. 5.74), converted from mm to m", unit = u"m"]
        ψ_min = -1.0e5, [description = "Minimum soil matric potential (converted from -1e8 mm to m)", unit = u"m"]
        one_m = 1.0, [description = "Unit length for nondimensionalization", unit = u"m"]
    end

    @parameters begin
        pct_sand, [description = "Percent sand content (dimensionless)"]
        pct_clay, [description = "Percent clay content (dimensionless)"]
        θ_i, [description = "Total volumetric soil water content", unit = u"m^3/m^3"]
    end

    @variables begin
        k_sat(t), [description = "Saturated hydraulic conductivity (Eq. 5.70)", unit = u"m/s"]
        θ_sat(t), [description = "Volumetric water content at saturation / porosity (Eq. 5.71)", unit = u"m^3/m^3"]
        B(t), [description = "Clapp-Hornberger exponent (Eq. 5.72) (dimensionless)"]
        ψ_sat(t), [description = "Saturated soil matric potential (Eq. 5.74)", unit = u"m"]
        ψ(t), [description = "Soil matric potential (Eq. 5.73)", unit = u"m"]
    end

    eqs = [
        # Eq. 5.70: Saturated hydraulic conductivity
        k_sat ~ k_coeff * 10.0^(-0.884 + 0.0153 * pct_sand),

        # Eq. 5.71: Porosity
        θ_sat ~ 0.489 - 0.00126 * pct_sand,

        # Eq. 5.72: Clapp-Hornberger exponent
        B ~ 2.91 + 0.159 * pct_clay,

        # Eq. 5.74: Saturated matric potential
        ψ_sat ~ ψ_coeff * 10.0^(1.88 - 0.0131 * pct_sand),

        # Eq. 5.73: Soil matric potential
        # ψ_i = ψ_sat,i * (θ_i / θ_sat,i)^{-B_i}
        # Constrained: θ_i/θ_sat,i ∈ [0.01, 1] and ψ ≥ ψ_min
        ψ ~ max(ψ_sat * (max(θ_i, 0.01 * θ_sat) / θ_sat)^(-B) * one_m / one_m, ψ_min),
    ]

    return System(eqs, t; name)
end


"""
    SurfaceRunoffInfiltration(; name=:SurfaceRunoffInfiltration)

Computes surface runoff and infiltration for the pervious road, following
Section 5.2 of Oleson et al. (2010).

Implements the SIMTOP-based (Niu et al., 2005) runoff model (Eq. 5.49),
fractional saturated area (Eq. 5.50), impermeable fraction from frozen soil
(Eq. 5.51), and maximum infiltration capacity (Eq. 5.52).

All hydraulic properties are in SI units (m, m/s). The conversion between
kg/(m²·s) water flux and m/s uses ρ_liq: 1 kg/(m²·s) = 1/ρ_liq m/s.

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.2, Eqs. 5.47–5.58, pp. 124–128.
"""
@component function SurfaceRunoffInfiltration(; name = :SurfaceRunoffInfiltration)

    @constants begin
        α = 3.0, [description = "Scale-dependent parameter for frozen fraction (Eq. 5.51) (dimensionless)"]
        f_over = 0.5, [description = "Decay factor for saturated fraction (p. 126)", unit = u"m^-1"]
        θ_imp = 0.05, [description = "Minimum effective porosity (p. 127)", unit = u"m^3/m^3"]
        ρ_liq = 1000.0, [description = "Density of liquid water (Table 1.4)", unit = u"kg/m^3"]
        one_kgpm2 = 1.0, [description = "Unit mass per area", unit = u"kg/m^2"]
        zero_kgpm2s = 0.0, [description = "Zero flux", unit = u"kg/(m^2*s)"]
    end

    @parameters begin
        q_liq_0, [description = "Liquid water reaching surface from rain/snowmelt", unit = u"kg/(m^2*s)"]
        q_seva, [description = "Surface evaporation of liquid water", unit = u"kg/(m^2*s)"]
        f_max, [description = "Maximum saturated fraction (dimensionless)"]
        z_v, [description = "Water table depth", unit = u"m"]
        k_sat_1, [description = "Saturated hydraulic conductivity of top soil layer", unit = u"m/s"]
        B_1, [description = "Clapp-Hornberger exponent for top soil layer (dimensionless)"]
        ψ_sat_1, [description = "Saturated matric potential for top soil layer", unit = u"m"]
        θ_liq_1, [description = "Volumetric liquid water content of top soil layer", unit = u"m^3/m^3"]
        θ_sat_1, [description = "Porosity of top soil layer", unit = u"m^3/m^3"]
        θ_ice_1, [description = "Volumetric ice content of top soil layer", unit = u"m^3/m^3"]
        w_ice_1, [description = "Ice content of top soil layer", unit = u"kg/m^2"]
        w_liq_1, [description = "Liquid water content of top soil layer", unit = u"kg/m^2"]
        Δz_1, [description = "Thickness of top soil layer", unit = u"m"]
        has_snow, [description = "Whether snow layers are present (1=yes, 0=no) (dimensionless)"]
    end

    @variables begin
        f_frz_1(t), [description = "Impermeable fraction from frozen soil (Eq. 5.51) (dimensionless)"]
        f_sat(t), [description = "Fractional saturated area (Eq. 5.50) (dimensionless)"]
        s(t), [description = "Relative saturation (Eq. 5.53) (dimensionless)"]
        v(t), [description = "Infiltration variable v (Eq. 5.54) (dimensionless)"]
        q_infl_max(t), [description = "Maximum infiltration capacity (Eq. 5.52)", unit = u"kg/(m^2*s)"]
        q_over(t), [description = "Surface runoff (Eq. 5.49)", unit = u"kg/(m^2*s)"]
        q_infl(t), [description = "Infiltration rate (Eqs. 5.56–5.57)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.51: Impermeable fraction from frozen soil
        f_frz_1 ~ (exp(-α * (1.0 - w_ice_1 / max(w_ice_1 + w_liq_1, 1e-10 * one_kgpm2))) - exp(-α)) /
                   (1.0 - exp(-α)),

        # Eq. 5.50: Fractional saturated area
        f_sat ~ (1.0 - f_frz_1) * f_max * exp(-0.5 * f_over * z_v) + f_frz_1,

        # Eq. 5.53: Relative saturation of top soil layer
        s ~ max((θ_liq_1 / max(θ_imp, θ_sat_1 - θ_ice_1) - f_sat) / max(1.0 - f_sat, 0.01), 0.0),

        # Eq. 5.55: Matric potential derivative at saturation
        # (dψ/ds)|_{s=1} = -B_1 * ψ_sat_1
        # Eq. 5.54: Variable v (dimensionless)
        v ~ B_1 * ψ_sat_1 / (0.5 * Δz_1),

        # Eq. 5.52: Maximum infiltration capacity
        # k_sat_1 in m/s; convert to kg/(m²·s) by multiplying by ρ_liq
        q_infl_max ~ ρ_liq * k_sat_1 * (1.0 + v * (s - 1.0)),

        # Eq. 5.49: Surface runoff for pervious road (SIMTOP)
        q_over ~ f_sat * q_liq_0 + (1.0 - f_sat) * max(zero_kgpm2s, q_liq_0 - q_infl_max),

        # Eqs. 5.56-5.57: Infiltration
        q_infl ~ ifelse(has_snow > 0.5,
            q_liq_0 - q_over,                          # Eq. 5.57: with snow
            max(q_liq_0 - q_over - q_seva, zero_kgpm2s)),  # Eq. 5.56: no snow
    ]

    return System(eqs, t; name)
end


"""
    SoilWaterFlux(; name=:SoilWaterFlux)

Computes the soil water flux between two adjacent soil layers using the modified
Darcy's law following Zeng and Decker (2009), as described in Section 5.3.2 of
Oleson et al. (2010).

The flux is computed as Eq. 5.67:
  q = -k[z_h] * [(ψ_i - ψ_{i+1}) + (ψ_{E,i+1} - ψ_{E,i})] / (z_{i+1} - z_i)

All lengths and hydraulic properties in SI units (m, m/s).

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.3.2, Eqs. 5.87–5.88, p. 135.
"""
@component function SoilWaterFlux(; name = :SoilWaterFlux)

    @parameters begin
        k_interface, [description = "Hydraulic conductivity at layer interface", unit = u"m/s"]
        ψ_upper, [description = "Soil matric potential of upper layer", unit = u"m"]
        ψ_lower, [description = "Soil matric potential of lower layer", unit = u"m"]
        ψ_E_upper, [description = "Equilibrium matric potential of upper layer", unit = u"m"]
        ψ_E_lower, [description = "Equilibrium matric potential of lower layer", unit = u"m"]
        z_upper, [description = "Node depth of upper layer", unit = u"m"]
        z_lower, [description = "Node depth of lower layer", unit = u"m"]
    end

    @variables begin
        q(t), [description = "Soil water flux at interface (Eq. 5.88, positive upward)", unit = u"m/s"]
    end

    eqs = [
        # Eq. 5.88: Soil water flux (positive upward)
        q ~ -k_interface * ((ψ_upper - ψ_lower) + (ψ_E_lower - ψ_E_upper)) / (z_lower - z_upper),
    ]

    return System(eqs, t; name)
end


"""
    SoilWaterEquilibrium(; name=:SoilWaterEquilibrium)

Computes the equilibrium soil matric potential and volumetric water content
following Section 5.3.2.1 of Oleson et al. (2010).

The equilibrium state represents the hydrostatic distribution of soil moisture
above the water table (Eqs. 5.98–5.106).

All lengths in SI units (m).

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.3.2.1, Eqs. 5.98–5.106, pp. 136–138.
"""
@component function SoilWaterEquilibrium(; name = :SoilWaterEquilibrium)

    @constants begin
        one_m = 1.0, [description = "Unit length for nondimensionalization", unit = u"m"]
        ψ_min = -1.0e5, [description = "Minimum matric potential (converted from -1e8 mm)", unit = u"m"]
    end

    @parameters begin
        θ_sat, [description = "Porosity", unit = u"m^3/m^3"]
        ψ_sat, [description = "Saturated matric potential", unit = u"m"]
        B, [description = "Clapp-Hornberger exponent (dimensionless)"]
        z_h_upper, [description = "Upper interface depth", unit = u"m"]
        z_h_lower, [description = "Lower interface depth", unit = u"m"]
        z_v, [description = "Water table depth", unit = u"m"]
    end

    @variables begin
        θ_E(t), [description = "Layer-averaged equilibrium volumetric water content (Eq. 5.101)", unit = u"m^3/m^3"]
        ψ_E(t), [description = "Equilibrium matric potential at node (Eq. 5.106)", unit = u"m"]
    end

    eqs = [
        # Eq. 5.101: Layer-averaged equilibrium water content
        # θ̄_{E,i} = (θ_sat * ψ_sat) / ((z_{h,i} - z_{h,i-1})(1 - 1/B))
        #   × [((ψ_sat - z_v + z_{h,i})/ψ_sat)^{1-1/B} - ((ψ_sat - z_v + z_{h,i-1})/ψ_sat)^{1-1/B}]
        # Simplified: clamp to [0, θ_sat] per Eq. 5.105
        θ_E ~ max(0.0 * θ_sat, min(θ_sat,
            (θ_sat * ψ_sat) / ((z_h_lower - z_h_upper) * (1.0 - 1.0 / B)) *
            (((ψ_sat - z_v + z_h_lower) / ψ_sat)^(1.0 - 1.0 / B) -
             ((ψ_sat - z_v + z_h_upper) / ψ_sat)^(1.0 - 1.0 / B))
        )),

        # Eq. 5.106: Equilibrium matric potential
        ψ_E ~ max(ψ_sat * (max(θ_E, 0.01 * θ_sat) / θ_sat)^(-B) * one_m / one_m, ψ_min),
    ]

    return System(eqs, t; name)
end


"""
    RichardsEquation(; name=:RichardsEquation, N=10)

Implements Richards' equation (Eq. 5.59) for vertical soil water movement,
discretized into `N` layers using the method of lines.

The equation `∂θ/∂t = -∂q/∂z - Q` is semi-discretized into an N-layer ODE system
where each layer has a volumetric water content `θ_liq[i]` as a state variable.
Soil hydraulic properties follow Clapp and Hornberger (1978) as parameterized in
Eqs. 5.69–5.74. Soil water flux between layers follows the modified Darcy's law
of Zeng and Decker (2009) (Eq. 5.88). Equilibrium moisture profiles follow
Eqs. 5.101 and 5.106.

Flux sign convention: positive upward (matching `SoilWaterFlux`).
- Top boundary: `q_0 = -q_infl` (infiltration is downward).
- Bottom boundary: `q_N = -q_bottom` (drainage is downward).

Frozen soil effects are not included (f_frz = 0).

**Reference**: Oleson et al. (2010), Chapter 5, Eqs. 5.59, 5.69–5.74, 5.76,
5.87–5.88, 5.98–5.106, pp. 127–138.
"""
@component function RichardsEquation(; name = :RichardsEquation, N = 10)

    @constants begin
        # Paper gives k_coeff = 0.0070556 mm/s; multiply by 0.001 to get m/s
        k_coeff = 0.0070556e-3, [description = "Coefficient for saturated hydraulic conductivity (Eq. 5.70), converted from mm/s to m/s", unit = u"m/s"]
        # Paper gives ψ_coeff = -10.0 mm; multiply by 0.001 to get m
        ψ_coeff = -10.0e-3, [description = "Coefficient for saturated matric potential (Eq. 5.74), converted from mm to m", unit = u"m"]
        ψ_min = -1.0e5, [description = "Minimum soil matric potential (converted from -1e8 mm to m)", unit = u"m"]
        one_m = 1.0, [description = "Unit length for nondimensionalization", unit = u"m"]
        one_mps = 1.0, [description = "Unit hydraulic conductivity for nondimensionalization", unit = u"m/s"]
    end

    @parameters begin
        pct_sand, [description = "Percent sand content (dimensionless)"]
        pct_clay, [description = "Percent clay content (dimensionless)"]
        z_v, [description = "Water table depth", unit = u"m"]
        q_infl, [description = "Infiltration rate at surface (positive downward)", unit = u"m/s"]
        q_bottom, [description = "Bottom boundary flux (positive downward, 0 for zero-flux)", unit = u"m/s"]
        e[1:N], [description = "Evapotranspiration sink per layer", unit = u"m/s"]
        z_node[1:N], [description = "Node depth of each layer", unit = u"m"]
        Δz_layer[1:N], [description = "Layer thickness", unit = u"m"]
        z_interface[1:(N + 1)], [description = "Interface depth (N+1 interfaces for N layers)", unit = u"m"]
    end

    @variables begin
        # State variables (N ODEs)
        θ_liq[1:N](t), [description = "Volumetric liquid water content per layer", unit = u"m^3/m^3"]
        # Soil properties from texture (4 algebraic)
        k_sat(t), [description = "Saturated hydraulic conductivity (Eq. 5.70)", unit = u"m/s"]
        θ_sat(t), [description = "Porosity (Eq. 5.71)", unit = u"m^3/m^3"]
        B_CH(t), [description = "Clapp-Hornberger exponent (Eq. 5.72) (dimensionless)"]
        ψ_sat(t), [description = "Saturated matric potential (Eq. 5.74)", unit = u"m"]
        # Per-layer algebraic variables
        ψ[1:N](t), [description = "Matric potential per layer (Eq. 5.73)", unit = u"m"]
        θ_E[1:N](t), [description = "Equilibrium volumetric water content per layer (Eq. 5.101)", unit = u"m^3/m^3"]
        ψ_E[1:N](t), [description = "Equilibrium matric potential per layer (Eq. 5.106)", unit = u"m"]
        # Interface variables
        k_interface[1:(N - 1)](t), [description = "Hydraulic conductivity at internal interfaces (Eq. 5.69)", unit = u"m/s"]
        q[1:(N - 1)](t), [description = "Water flux at internal interfaces (Eq. 5.88, positive upward)", unit = u"m/s"]
    end

    eqs = Equation[]

    # ===== Soil properties from texture (4 equations) =====
    # Eq. 5.70: Saturated hydraulic conductivity
    push!(eqs, k_sat ~ k_coeff * 10.0^(-0.884 + 0.0153 * pct_sand))
    # Eq. 5.71: Porosity
    push!(eqs, θ_sat ~ 0.489 - 0.00126 * pct_sand)
    # Eq. 5.72: Clapp-Hornberger exponent
    push!(eqs, B_CH ~ 2.91 + 0.159 * pct_clay)
    # Eq. 5.74: Saturated matric potential
    push!(eqs, ψ_sat ~ ψ_coeff * 10.0^(1.88 - 0.0131 * pct_sand))

    # ===== Per-layer matric potential (N equations) =====
    for i in 1:N
        # Eq. 5.73: ψ_i = ψ_sat * (θ_i/θ_sat)^{-B}, clamped to [ψ_min, 0]
        push!(eqs, ψ[i] ~ max(ψ_sat * (max(θ_liq[i], 0.01 * θ_sat) / θ_sat)^(-B_CH) * one_m / one_m, ψ_min))
    end

    # ===== Per-layer equilibrium θ_E (N equations) =====
    for i in 1:N
        # Eq. 5.101: Layer-averaged equilibrium water content
        push!(eqs, θ_E[i] ~ max(0.0 * θ_sat, min(θ_sat,
            (θ_sat * ψ_sat) / ((z_interface[i + 1] - z_interface[i]) * (1.0 - 1.0 / B_CH)) *
            (((ψ_sat - z_v + z_interface[i + 1]) / ψ_sat)^(1.0 - 1.0 / B_CH) -
             ((ψ_sat - z_v + z_interface[i]) / ψ_sat)^(1.0 - 1.0 / B_CH))
        )))
    end

    # ===== Per-layer equilibrium ψ_E (N equations) =====
    for i in 1:N
        # Eq. 5.106: Equilibrium matric potential
        push!(eqs, ψ_E[i] ~ max(ψ_sat * (max(θ_E[i], 0.01 * θ_sat) / θ_sat)^(-B_CH) * one_m / one_m, ψ_min))
    end

    # ===== Interface hydraulic conductivity (N-1 equations) =====
    for j in 1:(N - 1)
        # Eq. 5.69: k = k_sat * (θ/θ_sat)^{2B+3}, averaged at interface
        S_j = max(θ_liq[j], 0.01 * θ_sat) / θ_sat
        S_jp1 = max(θ_liq[j + 1], 0.01 * θ_sat) / θ_sat
        S_avg = 0.5 * (S_j + S_jp1)
        push!(eqs, k_interface[j] ~ k_sat * S_avg^(2.0 * B_CH + 3.0) * one_mps / one_mps)
    end

    # ===== Interface flux (N-1 equations) =====
    for j in 1:(N - 1)
        # Eq. 5.88: Soil water flux (positive upward)
        push!(eqs, q[j] ~ -k_interface[j] * ((ψ[j] - ψ[j + 1]) + (ψ_E[j + 1] - ψ_E[j])) / (z_node[j + 1] - z_node[j]))
    end

    # ===== Conservation ODEs (N equations) =====
    # Eq. 5.76: Δz[i] * dθ/dt = q_below - q_above - e[i]
    # With BCs: q_above for layer 1 = -q_infl, q_below for layer N = -q_bottom
    for i in 1:N
        if i == 1
            # Top layer: upper boundary is surface infiltration
            push!(eqs, D(θ_liq[i]) ~ (q[1] + q_infl - e[i]) / Δz_layer[i])
        elseif i == N
            # Bottom layer: lower boundary is bottom flux
            push!(eqs, D(θ_liq[i]) ~ (-q_bottom - q[N - 1] - e[i]) / Δz_layer[i])
        else
            # Interior layer
            push!(eqs, D(θ_liq[i]) ~ (q[i] - q[i - 1] - e[i]) / Δz_layer[i])
        end
    end

    return System(eqs, t; name)
end


"""
    GroundwaterDrainage(; name=:GroundwaterDrainage)

Computes sub-surface drainage for the pervious road based on the SIMTOP scheme,
following Section 5.4 of Oleson et al. (2010).

Drainage is computed as an exponential function of water table depth (Eq. 5.140),
modified by a frozen soil impermeable fraction (Eq. 5.141). Drainage is restricted
in frozen soils where ice content exceeds liquid water content (Eq. 5.139).

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.4, Eqs. 5.138–5.141, pp. 141–142.
"""
@component function GroundwaterDrainage(; name = :GroundwaterDrainage)

    @constants begin
        f_drai = 2.5, [description = "Drainage decay factor (p. 142)", unit = u"m^-1"]
        q_drai_max = 5.5e-3, [description = "Maximum drainage rate at surface (p. 142)", unit = u"kg/(m^2*s)"]
        α = 3.0, [description = "Scale-dependent parameter (Eq. 5.141) (dimensionless)"]
    end

    @parameters begin
        z_v, [description = "Water table depth", unit = u"m"]
        f_ice_weighted, [description = "Ice-weighted impermeable fraction (dimensionless)"]
    end

    @variables begin
        f_imp(t), [description = "Impermeable fraction for drainage (Eq. 5.141) (dimensionless)"]
        q_drai(t), [description = "Sub-surface drainage (Eq. 5.140)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.141: Impermeable fraction
        f_imp ~ max((exp(-α * (1.0 - f_ice_weighted)) - exp(-α)) / (1.0 - exp(-α)), 0.0),

        # Eq. 5.140: Modified drainage with frozen soil fraction
        q_drai ~ (1.0 - f_imp) * q_drai_max * exp(-f_drai * z_v),
    ]

    return System(eqs, t; name)
end


"""
    WaterTableDepth(; name=:WaterTableDepth)

Computes the water table depth from aquifer water storage, following
Section 5.4 of Oleson et al. (2010).

When the water table is below the soil column, the depth is computed from the
aquifer water storage (Eq. 5.146). Water table depth is constrained to
0.05 ≤ z_v ≤ 80 m.

Aquifer water W_a is in kg/m² (equivalent to mm of water column height).
Conversion: W_a [kg/m²] / ρ_liq [kg/m³] gives meters of water.

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.4, Eq. 5.146, p. 143.
"""
@component function WaterTableDepth(; name = :WaterTableDepth)

    @constants begin
        S_y = 0.2, [description = "Specific yield (Eq. 5.146) (dimensionless)"]
        z_v_min = 0.05, [description = "Minimum water table depth", unit = u"m"]
        z_v_max = 80.0, [description = "Maximum water table depth", unit = u"m"]
        z_offset = 25.0, [description = "Offset for water table calculation", unit = u"m"]
        ρ_liq = 1000.0, [description = "Density of liquid water (Table 1.4)", unit = u"kg/m^3"]
    end

    @parameters begin
        z_h_bottom, [description = "Depth of bottom soil interface", unit = u"m"]
        W_a, [description = "Water stored in unconfined aquifer", unit = u"kg/m^2"]
    end

    @variables begin
        z_v(t), [description = "Water table depth (Eq. 5.146)", unit = u"m"]
    end

    eqs = [
        # Eq. 5.146: z_v = z_{h,N_levsoi} + 25 - W_a / (10^3 * S_y)
        # Original: W_a in mm, so W_a/1000 gives m. With W_a in kg/m²,
        # W_a/ρ_liq gives m of water, then divide by S_y for water table depth.
        z_v ~ max(z_v_min, min(z_v_max,
            z_h_bottom + z_offset - W_a / (ρ_liq * S_y))),
    ]

    return System(eqs, t; name)
end


"""
    AquiferWaterBalance(; name=:AquiferWaterBalance)

Computes the update of aquifer water storage, following Section 5.4 of
Oleson et al. (2010).

For the case when the water table is below the soil column, the aquifer water
W_a (kg/m²) is updated from the balance of recharge and drainage (Eq. 5.143).
Since 1 mm water = 1 kg/m², we use kg/m² for storage and kg/(m²·s) for rates,
which ensures unit consistency.

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.4, Eqs. 5.142–5.147, pp. 142–144.
"""
@component function AquiferWaterBalance(; name = :AquiferWaterBalance)

    @parameters begin
        q_recharge, [description = "Recharge rate to aquifer", unit = u"kg/(m^2*s)"]
        q_drai, [description = "Sub-surface drainage rate", unit = u"kg/(m^2*s)"]
    end

    @variables begin
        dW_a(t), [description = "Rate of change of aquifer water (Eq. 5.143)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.143: dW_a/dt = q_recharge - q_drai
        # Both sides in kg/(m²·s)
        dW_a ~ q_recharge - q_drai,
    ]

    return System(eqs, t; name)
end


"""
    SnowCappingRunoff(; name=:SnowCappingRunoff)

Computes runoff from snow-capping when the snow water equivalent exceeds
1000 kg/m², following Section 5.5 of Oleson et al. (2010).

For snow-capped surfaces, precipitation and dew are routed to solid and liquid
runoff terms (Eqs. 5.153–5.154).

**Reference**: Oleson et al. (2010), Chapter 5, Section 5.5, Eqs. 5.153–5.154, p. 146.
"""
@component function SnowCappingRunoff(; name = :SnowCappingRunoff)

    @parameters begin
        q_grnd_ice, [description = "Solid precipitation reaching surface", unit = u"kg/(m^2*s)"]
        q_grnd_liq, [description = "Liquid precipitation reaching surface", unit = u"kg/(m^2*s)"]
        q_frost, [description = "Frost rate", unit = u"kg/(m^2*s)"]
        q_dew, [description = "Dew rate", unit = u"kg/(m^2*s)"]
    end

    @variables begin
        q_snwcp_ice(t), [description = "Solid runoff from snow capping (Eq. 5.153)", unit = u"kg/(m^2*s)"]
        q_snwcp_liq(t), [description = "Liquid runoff from snow capping (Eq. 5.154)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.153: Solid runoff from snow capping
        q_snwcp_ice ~ q_grnd_ice + q_frost,

        # Eq. 5.154: Liquid runoff from snow capping
        q_snwcp_liq ~ q_grnd_liq + q_dew,
    ]

    return System(eqs, t; name)
end


"""
    SurfaceLayerUpdate(; name=:SurfaceLayerUpdate)

Computes surface layer ice and liquid water content updates from dew, frost,
and sublimation, following Section 5.4 of Oleson et al. (2010).

These updates apply to the top layer of all surfaces (roof, pervious road,
impervious road).

**Reference**: Oleson et al. (2010), Chapter 5, Eqs. 5.150–5.152, p. 146.
"""
@component function SurfaceLayerUpdate(; name = :SurfaceLayerUpdate)

    @parameters begin
        q_sdew, [description = "Surface dew rate", unit = u"kg/(m^2*s)"]
        q_frost, [description = "Frost rate", unit = u"kg/(m^2*s)"]
        q_subl, [description = "Sublimation rate", unit = u"kg/(m^2*s)"]
    end

    @variables begin
        dw_liq_1(t), [description = "Rate of change of liquid water from dew (Eq. 5.150)", unit = u"kg/(m^2*s)"]
        dw_ice_1_frost(t), [description = "Rate of change of ice from frost (Eq. 5.151)", unit = u"kg/(m^2*s)"]
        dw_ice_1_subl(t), [description = "Rate of change of ice from sublimation (Eq. 5.152)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.150: Surface dew update (rate form)
        dw_liq_1 ~ q_sdew,

        # Eq. 5.151: Surface frost update (rate form)
        dw_ice_1_frost ~ q_frost,

        # Eq. 5.152: Surface sublimation update (rate form)
        dw_ice_1_subl ~ -q_subl,
    ]

    return System(eqs, t; name)
end


"""
    PerviousRoadWaterBalance(; name=:PerviousRoadWaterBalance)

Describes the overall water balance equation for the pervious road column,
following Eq. 5.1 of Oleson et al. (2010).

The pervious road has a full soil column with infiltration, surface runoff,
sub-surface drainage, redistribution within the soil, and groundwater interactions.

**Reference**: Oleson et al. (2010), Chapter 5, Eq. 5.1, p. 110.
"""
@component function PerviousRoadWaterBalance(; name = :PerviousRoadWaterBalance)

    @parameters begin
        q_rain, [description = "Liquid precipitation rate", unit = u"kg/(m^2*s)"]
        q_sno, [description = "Solid precipitation rate", unit = u"kg/(m^2*s)"]
        E_prvrd, [description = "Total evaporation (pervious road)", unit = u"kg/(m^2*s)"]
        q_over, [description = "Surface runoff", unit = u"kg/(m^2*s)"]
        q_drai, [description = "Sub-surface drainage", unit = u"kg/(m^2*s)"]
        q_rgwl, [description = "Liquid runoff from groundwater", unit = u"kg/(m^2*s)"]
        q_snwcp_ice, [description = "Solid runoff from snow capping", unit = u"kg/(m^2*s)"]
    end

    @variables begin
        water_input(t), [description = "Net water input rate (Eq. 5.1)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.1: Water balance of pervious road
        # ΔW_sno + Σ(Δw_liq,i + Δw_ice,i) + ΔW_a = (q_rain + q_sno - E_prvrd - q_over - q_drai - q_rgwl - q_snwcp,ice)Δt
        water_input ~ q_rain + q_sno - E_prvrd - q_over - q_drai - q_rgwl - q_snwcp_ice,
    ]

    return System(eqs, t; name)
end


"""
    ImperviousWaterBalance(; name=:ImperviousWaterBalance)

Describes the water balance for the roof and impervious road surfaces,
following Eqs. 5.2–5.3 of Oleson et al. (2010).

These surfaces have limited water storage capacity (1 kg/m² ponding limit),
no sub-surface drainage, and intercept both liquid and solid precipitation.

**Reference**: Oleson et al. (2010), Chapter 5, Eqs. 5.2–5.5, pp. 112.
"""
@component function ImperviousWaterBalance(; name = :ImperviousWaterBalance)

    @constants begin
        w_pond_max = 1.0, [description = "Maximum ponding limit (p. 112)", unit = u"kg/m^2"]
        zero_kgpm2s = 0.0, [description = "Zero flux", unit = u"kg/(m^2*s)"]
    end

    @parameters begin
        q_rain, [description = "Liquid precipitation rate (Eq. 5.4)", unit = u"kg/(m^2*s)"]
        q_sno, [description = "Solid precipitation rate (Eq. 5.5)", unit = u"kg/(m^2*s)"]
        E_surface, [description = "Total evaporation from surface", unit = u"kg/(m^2*s)"]
        q_rgwl, [description = "Liquid runoff from groundwater", unit = u"kg/(m^2*s)"]
        q_snwcp_ice, [description = "Solid runoff from snow capping", unit = u"kg/(m^2*s)"]
    end

    @variables begin
        water_input(t), [description = "Net water input rate (Eqs. 5.2-5.3)", unit = u"kg/(m^2*s)"]
        q_grnd_liq(t), [description = "Liquid precipitation reaching surface (Eq. 5.4)", unit = u"kg/(m^2*s)"]
        q_grnd_ice(t), [description = "Solid precipitation reaching surface (Eq. 5.5)", unit = u"kg/(m^2*s)"]
    end

    eqs = [
        # Eq. 5.4: Liquid precipitation reaching surface
        q_grnd_liq ~ q_rain,

        # Eq. 5.5: Solid precipitation reaching surface
        q_grnd_ice ~ q_sno,

        # Eqs. 5.2-5.3: Water balance for roof/impervious road
        water_input ~ q_rain + q_sno - E_surface - q_rgwl - q_snwcp_ice,
    ]

    return System(eqs, t; name)
end
