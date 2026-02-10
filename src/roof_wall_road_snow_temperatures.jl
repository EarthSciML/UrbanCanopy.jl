export SoilThermalProperties, SnowThermalProperties, UrbanSurfaceThermalProperties,
    InterfaceThermalConductivity, HeatFlux,
    SurfaceEnergyFlux, BuildingTemperature, WasteHeatAirConditioning,
    PhaseChangeEnergy, SnowLayerGeometry,
    WasteHeatAllocation, AdjustedLayerThickness, HeatingCoolingFlux,
    PhaseChangeAdjustment, SnowMeltNoLayers,
    UniformGrid, ExponentialGrid,
    FreezingPointDepression, SnowSoilBlendedHeatCapacity,
    TotalPhaseChangeEnergy, LayerPhaseChangeEnergy

"""
    SnowLayerGeometry(; name=:SnowLayerGeometry)

Computes the snow layer structure for a given snow depth, following Section 4.1
of Oleson et al. (2010). The snow pack has up to 5 layers, with thickness
assignments depending on the total snow depth z_sno (pp. 92–93).

This component takes the total snow depth and number of snow layers as inputs
and computes the layer thicknesses, node depths, and interface depths
(Eqs. 4.9–4.10).

Note: The actual number of snow layers (snl) is determined by the snow depth
according to the rules on pp. 92–93. This component computes geometry for a
single snow layer (the simplest case), since multi-layer snow packing logic
is procedural and typically handled outside the equation system.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, pp. 92–93.
"""
@component function SnowLayerGeometry(; name = :SnowLayerGeometry)

    @parameters begin
        z_sno, [description = "Total snow depth", unit = u"m"]
    end

    @variables begin
        Δz_snow_0(t), [description = "Thickness of bottom snow layer (layer 0)", unit = u"m"]
        z_node_snow_0(t), [description = "Node depth of bottom snow layer", unit = u"m"]
        z_interface_snow_neg1(t), [description = "Interface depth above bottom snow layer", unit = u"m"]
    end

    @constants begin
        zero_m = 0.0, [description = "Zero depth reference", unit = u"m"]
    end

    eqs = [
        # For single snow layer (snl = -1): Δz_0 = z_sno (pp. 92)
        Δz_snow_0 ~ z_sno,

        # Eq. 4.9: z_i = z_{h,i} - 0.5 * Δz_i; z_{h,0} = 0
        z_node_snow_0 ~ zero_m - 0.5 * Δz_snow_0,

        # Eq. 4.10: z_{h,i} = z_{h,i+1} - Δz_{i+1}; z_{h,0} = 0
        z_interface_snow_neg1 ~ zero_m - Δz_snow_0,
    ]

    return System(eqs, t; name)
end


"""
    SoilThermalProperties(; name=:SoilThermalProperties)

Computes the thermal conductivity and volumetric heat capacity for soil layers
used in pervious road columns, following Section 4.3 of Oleson et al. (2010).

Thermal conductivity follows Farouki (1981) using the Kersten number method
(Eqs. 4.77–4.82). Volumetric heat capacity follows de Vries (1963) (Eqs. 4.85–4.86).

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.3, pp. 107–108.
"""
@component function SoilThermalProperties(; name = :SoilThermalProperties)

    # Physical constants (Table 1.4)
    @constants begin
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
        ρ_liq = 1000.0, [description = "Density of liquid water (Table 1.4)", unit = u"kg/m^3"]
        ρ_ice = 917.0, [description = "Density of ice (Table 1.4)", unit = u"kg/m^3"]
        C_liq = 4188.0, [description = "Specific heat of liquid water (Table 1.4)", unit = u"J/(kg*K)"]
        C_ice = 2117.27, [description = "Specific heat of ice (Table 1.4)", unit = u"J/(kg*K)"]
        λ_liq_const = 0.6, [description = "Thermal conductivity of liquid water (Table 1.4)", unit = u"W/(m*K)"]
        λ_ice_const = 2.29, [description = "Thermal conductivity of ice (Table 1.4)", unit = u"W/(m*K)"]
    end

    # Reference constants for non-integer exponents and unit handling
    @constants begin
        one_Wm⁻¹K⁻¹ = 1.0, [description = "Unit thermal conductivity for non-dimensionalization", unit = u"W/(m*K)"]
        one_kgm⁻³ = 1.0, [description = "Unit density for non-dimensionalization", unit = u"kg/m^3"]
        one_Jm⁻³K⁻¹ = 1.0, [description = "Unit volumetric heat capacity for non-dimensionalization", unit = u"J/(m^3*K)"]
    end

    # Reference constants for empirical formula unit tracking (Eq. 4.79)
    @constants begin
        λ_s_sand_coeff = 8.8, [description = "Empirical coefficient for sand in Eq. 4.79 (W m⁻¹ K⁻¹ per %)", unit = u"W/(m*K)"]
        λ_s_clay_coeff = 2.92, [description = "Empirical coefficient for clay in Eq. 4.79 (W m⁻¹ K⁻¹ per %)", unit = u"W/(m*K)"]
        ρ_mineral = 2700.0, [description = "Mineral soil particle density", unit = u"kg/m^3"]
        λ_dry_numerator_coeff = 0.135, [description = "Empirical numerator coefficient in Eq. 4.80 (W m² kg⁻¹ K⁻¹)", unit = u"W*m^2/(kg*K)"]
        λ_dry_numerator_offset = 64.7, [description = "Empirical numerator offset in Eq. 4.80", unit = u"W/(m*K)"]
        λ_dry_denom_offset = 2700.0, [description = "Empirical denominator offset in Eq. 4.80 (dimensionless)"]
        λ_dry_denom_coeff = 0.947, [description = "Empirical denominator coefficient in Eq. 4.80 (dimensionless)"]
        c_s_sand_coeff = 2.128e6, [description = "Empirical coefficient for sand in Eq. 4.86 (J m⁻³ K⁻¹ per %)", unit = u"J/(m^3*K)"]
        c_s_clay_coeff = 2.385e6, [description = "Empirical coefficient for clay in Eq. 4.86 (J m⁻³ K⁻¹ per %)", unit = u"J/(m^3*K)"]
    end

    @parameters begin
        pct_sand, [description = "Percent sand (dimensionless)"]
        pct_clay, [description = "Percent clay (dimensionless)"]
        θ_sat, [description = "Porosity / volumetric water content at saturation (dimensionless)"]
        θ_liq, [description = "Volumetric liquid water content (dimensionless)"]
        T_i, [description = "Temperature of layer", unit = u"K"]
        w_ice, [description = "Ice mass in layer", unit = u"kg/m^2"]
        w_liq, [description = "Liquid water mass in layer", unit = u"kg/m^2"]
        Δz, [description = "Layer thickness", unit = u"m"]
    end

    @variables begin
        λ_s(t), [description = "Thermal conductivity of soil solids (Eq. 4.79)", unit = u"W/(m*K)"]
        λ_sat(t), [description = "Saturated thermal conductivity (Eq. 4.78)", unit = u"W/(m*K)"]
        λ_dry(t), [description = "Dry thermal conductivity (Eq. 4.80)", unit = u"W/(m*K)"]
        ρ_d(t), [description = "Bulk density of dry soil", unit = u"kg/m^3"]
        S_r(t), [description = "Degree of saturation (Eq. 4.82) (dimensionless)"]
        K_e(t), [description = "Kersten number (Eq. 4.81) (dimensionless)"]
        λ_soil(t), [description = "Soil thermal conductivity (Eq. 4.77)", unit = u"W/(m*K)"]
        c_s(t), [description = "Heat capacity of soil solids (Eq. 4.86)", unit = u"J/(m^3*K)"]
        c_soil(t), [description = "Volumetric heat capacity of soil (Eq. 4.85)", unit = u"J/(m^3*K)"]
    end

    eqs = [
        # Eq. 4.79: λ_{s,i} = (8.80 * %sand + 2.92 * %clay) / (%sand + %clay)
        λ_s ~ (λ_s_sand_coeff * pct_sand + λ_s_clay_coeff * pct_clay) / (pct_sand + pct_clay),

        # Eq. 4.78: Saturated thermal conductivity
        # For T_i ≥ T_f: λ_sat = λ_s^(1-θ_sat) * λ_liq^θ_sat
        # For T_i < T_f: λ_sat = λ_s^(1-θ_sat) * λ_liq^θ_liq * λ_ice^(θ_sat - θ_liq)
        λ_sat ~ ifelse(
            T_i ≥ T_f,
            (λ_s / one_Wm⁻¹K⁻¹)^(1 - θ_sat) * (λ_liq_const / one_Wm⁻¹K⁻¹)^θ_sat * one_Wm⁻¹K⁻¹,
            (λ_s / one_Wm⁻¹K⁻¹)^(1 - θ_sat) * (λ_liq_const / one_Wm⁻¹K⁻¹)^θ_liq * (λ_ice_const / one_Wm⁻¹K⁻¹)^(θ_sat - θ_liq) * one_Wm⁻¹K⁻¹
        ),

        # Bulk density: ρ_d = 2700 * (1 - θ_sat) (Section 4.3)
        ρ_d ~ ρ_mineral * (1 - θ_sat),

        # Eq. 4.80: λ_dry = (0.135 * ρ_d + 64.7) / (2700 - 0.947 * ρ_d)
        # Denominator is dimensionless (ρ_d divided by reference density)
        λ_dry ~ (λ_dry_numerator_coeff * ρ_d + λ_dry_numerator_offset) / (λ_dry_denom_offset - λ_dry_denom_coeff * ρ_d / one_kgm⁻³),

        # Eq. 4.82: S_r = (w_liq / (ρ_liq * Δz) + w_ice / (ρ_ice * Δz)) / θ_sat
        S_r ~ (w_liq / (ρ_liq * Δz) + w_ice / (ρ_ice * Δz)) / θ_sat,

        # Eq. 4.81: Kersten number (unfrozen: K_e = log(S_r) + 1 ≥ 0; frozen: K_e = S_r)
        K_e ~ ifelse(
            T_i ≥ T_f,
            max(log(max(S_r, 1.0e-10)) + 1.0, 0.0),
            S_r
        ),

        # Eq. 4.77: λ = K_e * λ_sat + (1 - K_e) * λ_dry for S_r > 1.0e-7
        λ_soil ~ ifelse(
            S_r > 1.0e-7,
            K_e * λ_sat + (1 - K_e) * λ_dry,
            λ_dry
        ),

        # Eq. 4.86: c_s = (2.128 * %sand + 2.385 * %clay) / (%sand + %clay) × 10^6
        c_s ~ (c_s_sand_coeff * pct_sand + c_s_clay_coeff * pct_clay) / (pct_sand + pct_clay),

        # Eq. 4.85: c_i = c_s * (1 - θ_sat) + (w_ice / Δz) * C_ice + (w_liq / Δz) * C_liq
        c_soil ~ c_s * (1 - θ_sat) + (w_ice / Δz) * C_ice + (w_liq / Δz) * C_liq,
    ]

    return System(eqs, t; name)
end


"""
    SnowThermalProperties(; name=:SnowThermalProperties)

Computes the thermal conductivity and volumetric heat capacity for snow layers,
following Section 4.3 of Oleson et al. (2010).

Thermal conductivity follows Jordan (1991) (Eq. 4.83). Volumetric heat capacity
depends on the ice and liquid water content (Eq. 4.87).

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.3, pp. 108–109.
"""
@component function SnowThermalProperties(; name = :SnowThermalProperties)

    # Physical constants (Table 1.4)
    @constants begin
        C_liq = 4188.0, [description = "Specific heat of liquid water (Table 1.4)", unit = u"J/(kg*K)"]
        C_ice = 2117.27, [description = "Specific heat of ice (Table 1.4)", unit = u"J/(kg*K)"]
        λ_air = 0.023, [description = "Thermal conductivity of air (Table 1.4)", unit = u"W/(m*K)"]
        λ_ice_const = 2.29, [description = "Thermal conductivity of ice (Table 1.4)", unit = u"W/(m*K)"]
    end

    @constants begin
        one_kgm⁻³ = 1.0, [description = "Unit density for non-dimensionalization", unit = u"kg/m^3"]
        jordan_linear_coeff = 7.75e-5, [description = "Jordan (1991) linear density coefficient in Eq. 4.83", unit = u"m^3/kg"]
        jordan_quad_coeff = 1.105e-6, [description = "Jordan (1991) quadratic density coefficient in Eq. 4.83", unit = u"m^6/kg^2"]
    end

    @parameters begin
        w_ice, [description = "Ice mass in snow layer", unit = u"kg/m^2"]
        w_liq, [description = "Liquid water mass in snow layer", unit = u"kg/m^2"]
        Δz, [description = "Snow layer thickness", unit = u"m"]
    end

    @variables begin
        ρ_sno(t), [description = "Snow bulk density (Eq. 4.84)", unit = u"kg/m^3"]
        λ_snow(t), [description = "Snow thermal conductivity (Eq. 4.83)", unit = u"W/(m*K)"]
        c_snow(t), [description = "Snow volumetric heat capacity (Eq. 4.87)", unit = u"J/(m^3*K)"]
    end

    eqs = [
        # Eq. 4.84: ρ_sno = (w_ice + w_liq) / Δz
        ρ_sno ~ (w_ice + w_liq) / Δz,

        # Eq. 4.83: λ = λ_air + (7.75e-5 * ρ_sno + 1.105e-6 * ρ_sno²) * (λ_ice - λ_air)
        λ_snow ~ λ_air + (jordan_linear_coeff * ρ_sno + jordan_quad_coeff * ρ_sno^2) * (λ_ice_const - λ_air),

        # Eq. 4.87: c = (w_ice / Δz) * C_ice + (w_liq / Δz) * C_liq
        c_snow ~ (w_ice / Δz) * C_ice + (w_liq / Δz) * C_liq,
    ]

    return System(eqs, t; name)
end


"""
    UrbanSurfaceThermalProperties(; name=:UrbanSurfaceThermalProperties)

Represents the prescribed thermal conductivity and heat capacity for roof, wall,
and impervious road layers, following Section 4.3 of Oleson et al. (2010).

For roofs, walls, and impervious road layers (i = 1, ..., N_imprvrd), the thermal
conductivity and heat capacity are specified by the surface dataset (Table 1.3).
This component takes those prescribed values as parameters.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.3, pp. 106.
"""
@component function UrbanSurfaceThermalProperties(; name = :UrbanSurfaceThermalProperties)

    @parameters begin
        λ_prescribed, [description = "Prescribed thermal conductivity (Table 1.3)", unit = u"W/(m*K)"]
        c_prescribed, [description = "Prescribed volumetric heat capacity (Table 1.3)", unit = u"J/(m^3*K)"]
    end

    @variables begin
        λ_surf(t), [description = "Surface thermal conductivity", unit = u"W/(m*K)"]
        c_surf(t), [description = "Surface volumetric heat capacity", unit = u"J/(m^3*K)"]
    end

    eqs = [
        λ_surf ~ λ_prescribed,
        c_surf ~ c_prescribed,
    ]

    return System(eqs, t; name)
end


"""
    InterfaceThermalConductivity(; name=:InterfaceThermalConductivity)

Computes the thermal conductivity at the interface between two adjacent layers
using the harmonic mean formula (Eq. 4.12) from Section 4.1 of Oleson et al. (2010).

The interface conductivity is derived from continuity of heat flux at the interface
(Eq. 4.13), resulting in a weighted harmonic mean of the two adjacent layer
conductivities.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, Eq. 4.12, pp. 94.
"""
@component function InterfaceThermalConductivity(; name = :InterfaceThermalConductivity)

    @parameters begin
        λ_i, [description = "Thermal conductivity of layer i", unit = u"W/(m*K)"]
        λ_ip1, [description = "Thermal conductivity of layer i+1", unit = u"W/(m*K)"]
        z_i, [description = "Node depth of layer i", unit = u"m"]
        z_ip1, [description = "Node depth of layer i+1", unit = u"m"]
        z_h, [description = "Interface depth between layers i and i+1", unit = u"m"]
    end

    @variables begin
        λ_interface(t), [description = "Thermal conductivity at interface (Eq. 4.12)", unit = u"W/(m*K)"]
    end

    eqs = [
        # Eq. 4.12: λ[z_{h,i}] = λ_i * λ_{i+1} * (z_{i+1} - z_i) /
        #           [λ_i * (z_{i+1} - z_{h,i}) + λ_{i+1} * (z_{h,i} - z_i)]
        λ_interface ~ λ_i * λ_ip1 * (z_ip1 - z_i) / (λ_i * (z_ip1 - z_h) + λ_ip1 * (z_h - z_i)),
    ]

    return System(eqs, t; name)
end


"""
    HeatFlux(; name=:HeatFlux)

Computes the heat flux from layer i to layer i+1 using Fourier's law discretized
across the interface (Eq. 4.11) from Section 4.1 of Oleson et al. (2010).

The flux is defined as positive upwards: F_i = -λ[z_{h,i}] * (T_i - T_{i+1}) / (z_{i+1} - z_i).

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, Eq. 4.11, pp. 94.
"""
@component function HeatFlux(; name = :HeatFlux)

    @parameters begin
        λ_interface, [description = "Thermal conductivity at the interface", unit = u"W/(m*K)"]
        T_i, [description = "Temperature of layer i", unit = u"K"]
        T_ip1, [description = "Temperature of layer i+1", unit = u"K"]
        z_i, [description = "Node depth of layer i", unit = u"m"]
        z_ip1, [description = "Node depth of layer i+1", unit = u"m"]
    end

    @variables begin
        F(t), [description = "Heat flux from layer i to i+1 (Eq. 4.11)", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.11: F_i = -λ[z_{h,i}] * (T_i - T_{i+1}) / (z_{i+1} - z_i)
        F ~ -λ_interface * (T_i - T_ip1) / (z_ip1 - z_i),
    ]

    return System(eqs, t; name)
end


"""
    SurfaceEnergyFlux(; name=:SurfaceEnergyFlux)

Computes the net heat flux into each urban surface and its derivative with respect
to surface temperature, following Eqs. 4.26, 4.28, and 4.29 of Oleson et al. (2010).

The surface heat flux h is the sum of absorbed solar radiation, net longwave radiation,
sensible and latent heat fluxes, and waste heat / air conditioning terms. The
derivative ∂h/∂T is used in the linearized Crank-Nicholson scheme (Eq. 4.19).

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.26–4.29, pp. 97–98.
"""
@component function SurfaceEnergyFlux(; name = :SurfaceEnergyFlux)

    # Physical constants
    @constants begin
        σ_SB = 5.67e-8, [description = "Stefan-Boltzmann constant (Table 1.4)", unit = u"W/(m^2*K^4)"]
        one_K = 1.0, [description = "Unit temperature for non-dimensionalization", unit = u"K"]
    end

    @parameters begin
        S_g, [description = "Absorbed solar radiation", unit = u"W/m^2"]
        L_g, [description = "Net longwave radiation (downward positive)", unit = u"W/m^2"]
        H_g, [description = "Sensible heat flux", unit = u"W/m^2"]
        λE_g, [description = "Latent heat flux", unit = u"W/m^2"]
        H_wasteheat_g, [description = "Waste heat flux applied to this surface", unit = u"W/m^2"]
        H_aircond_g, [description = "Air conditioning heat applied to this surface", unit = u"W/m^2"]
        ε_g, [description = "Surface emissivity (dimensionless)"]
        T_g_n, [description = "Surface temperature at time n", unit = u"K"]
        dH_g_dT, [description = "Derivative of sensible heat flux w.r.t. T_g", unit = u"W/(m^2*K)"]
        dλE_g_dT, [description = "Derivative of latent heat flux w.r.t. T_g", unit = u"W/(m^2*K)"]
    end

    @variables begin
        h(t), [description = "Net heat flux into the surface (Eq. 4.26)", unit = u"W/m^2"]
        dh_dT(t), [description = "Partial derivative of h w.r.t. T_g (Eq. 4.28)", unit = u"W/(m^2*K)"]
        dL_g_dT(t), [description = "Partial derivative of net longwave radiation (Eq. 4.29)", unit = u"W/(m^2*K)"]
    end

    eqs = [
        # Eq. 4.26: h = S̄_g - L̄_g - H_g - λE_g + H_wasteheat + H_aircond
        h ~ S_g - L_g - H_g - λE_g + H_wasteheat_g + H_aircond_g,

        # Eq. 4.29: ∂L̄_g/∂T_g = 4 * ε_g * σ * (T_g^n)^3
        dL_g_dT ~ 4 * ε_g * σ_SB * (T_g_n / one_K)^3 * one_K^3,

        # Eq. 4.28: ∂h/∂T_g = -∂L̄_g/∂T_g - ∂H_g/∂T_g - ∂(λE_g)/∂T_g
        dh_dT ~ -dL_g_dT - dH_g_dT - dλE_g_dT,
    ]

    return System(eqs, t; name)
end


"""
    BuildingTemperature(; name=:BuildingTemperature)

Computes the internal building temperature from a weighted combination of
inner layer wall and roof temperatures, following Eqs. 4.37–4.38 of
Oleson et al. (2010).

The building temperature T_{iB} is constrained between prescribed minimum
and maximum values (T_{iB,min} and T_{iB,max}).

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.37–4.38, pp. 99.
"""
@component function BuildingTemperature(; name = :BuildingTemperature)

    @parameters begin
        T_inner_shdwall, [description = "Inner layer temperature of shaded wall", unit = u"K"]
        T_inner_sunwall, [description = "Inner layer temperature of sunlit wall", unit = u"K"]
        T_inner_roof, [description = "Inner layer temperature of roof", unit = u"K"]
        H_canyon, [description = "Building/canyon height", unit = u"m"]
        W_roof, [description = "Roof fraction (dimensionless)"]
        H_W, [description = "Height-to-width ratio (dimensionless)"]
    end

    @variables begin
        L_roof(t), [description = "Roof length in infinite canyon (Eq. 4.38)", unit = u"m"]
        T_iB_unclamped(t), [description = "Unclamped internal building temperature (Eq. 4.37)", unit = u"K"]
    end

    eqs = [
        # Eq. 4.38: L = (H / (H/W)) * (W_roof / (1 - W_roof))
        L_roof ~ (H_canyon / H_W) * (W_roof / (1 - W_roof)),

        # Eq. 4.37: T_iB = [H * (T_shdwall + T_sunwall) + L_roof * T_roof] / (2H + L_roof)
        T_iB_unclamped ~ (H_canyon * (T_inner_shdwall + T_inner_sunwall) + L_roof * T_inner_roof) / (2 * H_canyon + L_roof),
    ]

    return System(eqs, t; name)
end


"""
    WasteHeatAirConditioning(; name=:WasteHeatAirConditioning)

Computes waste heat from space heating/air conditioning and heat removed by
air conditioning, following Eqs. 4.55–4.56 of Oleson et al. (2010).

Waste heat is distributed to the pervious and impervious road surfaces. Air
conditioning heat removal equals the cooling flux (Eq. 4.56).

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.51–4.56, pp. 101–102.
"""
@component function WasteHeatAirConditioning(; name = :WasteHeatAirConditioning)

    @constants begin
        f_heat = 1.0 / 0.75, [description = "Heating efficiency factor (Eq. 4.55) (dimensionless)"]
        f_cool = 1.0 / 0.25, [description = "Cooling efficiency factor (Eq. 4.55) (dimensionless)"]
        H_wasteheat_max = 100.0, [description = "Maximum waste heat limit (Eq. 4.55)", unit = u"W/m^2"]
    end

    @parameters begin
        F_heat_roof, [description = "Heating flux for roof", unit = u"W/m^2"]
        F_cool_roof, [description = "Cooling flux for roof", unit = u"W/m^2"]
        F_heat_sunwall, [description = "Heating flux for sunlit wall", unit = u"W/m^2"]
        F_cool_sunwall, [description = "Cooling flux for sunlit wall", unit = u"W/m^2"]
        F_heat_shdwall, [description = "Heating flux for shaded wall", unit = u"W/m^2"]
        F_cool_shdwall, [description = "Cooling flux for shaded wall", unit = u"W/m^2"]
        W_roof, [description = "Roof fraction (dimensionless)"]
        H_W, [description = "Height-to-width ratio (dimensionless)"]
    end

    @variables begin
        H_wasteheat_unclamped(t), [description = "Unclamped total waste heat (Eq. 4.55)", unit = u"W/m^2"]
        H_wasteheat(t), [description = "Total waste heat (clamped to max) (Eq. 4.55)", unit = u"W/m^2"]
        H_aircond(t), [description = "Heat removed by air conditioning (Eq. 4.56)", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.55: Total waste heat
        H_wasteheat_unclamped ~ W_roof * (f_heat * F_heat_roof + f_cool * F_cool_roof) +
            (1 - W_roof) * H_W * (
            f_heat * F_heat_sunwall + f_cool * F_cool_sunwall +
                f_heat * F_heat_shdwall + f_cool * F_cool_shdwall
        ),

        # Clamped to H_wasteheat_max
        H_wasteheat ~ min(H_wasteheat_unclamped, H_wasteheat_max),

        # Eq. 4.56: H_aircond = F_cool (sum of all cooling fluxes)
        H_aircond ~ F_cool_roof + F_cool_sunwall + F_cool_shdwall,
    ]

    return System(eqs, t; name)
end


"""
    PhaseChangeEnergy(; name=:PhaseChangeEnergy)

Computes the energy excess or deficit for phase change assessment in a layer,
following Eq. 4.59 of Oleson et al. (2010).

For the top layer (i = snl+1), the energy includes the surface heat flux and
its derivative. For interior layers (i > snl+1), only the heat fluxes and
thermal storage are included.

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.59–4.65, pp. 103–104.
"""
@component function PhaseChangeEnergy(; name = :PhaseChangeEnergy, layer_type = :interior)

    @constants begin
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
        L_f = 3.337e5, [description = "Latent heat of fusion (Table 1.4)", unit = u"J/kg"]
        α_CN = 0.5, [description = "Crank-Nicholson weight (dimensionless)"]
    end

    @parameters begin
        c_i, [description = "Volumetric heat capacity", unit = u"J/(m^3*K)"]
        Δz_i, [description = "Layer thickness", unit = u"m"]
        Δt, [description = "Time step", unit = u"s"]
        T_i_n, [description = "Temperature at time n", unit = u"K"]
        F_i_n, [description = "Heat flux below at time n", unit = u"W/m^2"]
        F_i_np1, [description = "Heat flux below at time n+1", unit = u"W/m^2"]
    end

    @variables begin
        H_excess(t), [description = "Energy excess/deficit for phase change (Eq. 4.59)", unit = u"W/m^2"]
    end

    if layer_type == :top
        @parameters begin
            h_n, [description = "Surface heat flux at time n", unit = u"W/m^2"]
            dh_dT, [description = "Derivative of surface heat flux", unit = u"W/(m^2*K)"]
        end

        eqs = [
            # Eq. 4.59 (i = snl+1):
            # H_i = h + (∂h/∂T)(T_f - T_i^n) + α*F_i^n + (1-α)*F_i^{n+1}
            #       - (c_i*Δz_i/Δt)(T_f - T_i^n)
            H_excess ~ h_n + dh_dT * (T_f - T_i_n) + α_CN * F_i_n + (1 - α_CN) * F_i_np1 -
                (c_i * Δz_i / Δt) * (T_f - T_i_n),
        ]
    else
        @parameters begin
            F_im1_n, [description = "Heat flux above at time n", unit = u"W/m^2"]
            F_im1_np1, [description = "Heat flux above at time n+1", unit = u"W/m^2"]
        end

        eqs = [
            # Eq. 4.59 (i = snl+2,...,N_levgrnd):
            # H_i = α*(F_i^n - F_{i-1}^n) + (1-α)*(F_i^{n+1} - F_{i-1}^{n+1})
            #       - (c_i*Δz_i/Δt)(T_f - T_i^n)
            H_excess ~ α_CN * (F_i_n - F_im1_n) + (1 - α_CN) * (F_i_np1 - F_im1_np1) -
                (c_i * Δz_i / Δt) * (T_f - T_i_n),
        ]
    end

    return System(eqs, t; name)
end


"""
    WasteHeatAllocation(; name=:WasteHeatAllocation)

Allocates total waste heat and air conditioning heat to individual urban surfaces,
following Eq. 4.27 of Oleson et al. (2010).

Waste heat and air conditioning are applied ONLY to the pervious and impervious road
surfaces (divided by `1 - W_roof`). Walls and roof receive zero waste heat and zero
air conditioning.

**Reference**: Oleson et al. (2010), Chapter 4, Eq. 4.27, pp. 97.
"""
@component function WasteHeatAllocation(; name = :WasteHeatAllocation)

    @parameters begin
        H_wasteheat, [description = "Total waste heat (Eq. 4.55)", unit = u"W/m^2"]
        H_aircond, [description = "Total heat removed by air conditioning (Eq. 4.56)", unit = u"W/m^2"]
        W_roof, [description = "Roof fraction (dimensionless)"]
    end

    @variables begin
        H_wasteheat_prvrd(t), [description = "Waste heat to pervious road (Eq. 4.27)", unit = u"W/m^2"]
        H_wasteheat_imprvrd(t), [description = "Waste heat to impervious road (Eq. 4.27)", unit = u"W/m^2"]
        H_aircond_prvrd(t), [description = "Air conditioning heat to pervious road (Eq. 4.27)", unit = u"W/m^2"]
        H_aircond_imprvrd(t), [description = "Air conditioning heat to impervious road (Eq. 4.27)", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.27: H_{wasteheat,prvrd} = H_{wasteheat,imprvrd} = H_wasteheat / (1 - W_roof)
        H_wasteheat_prvrd ~ H_wasteheat / (1 - W_roof),
        H_wasteheat_imprvrd ~ H_wasteheat / (1 - W_roof),

        # Eq. 4.27: H_{aircond,prvrd} = H_{aircond,imprvrd} = H_aircond / (1 - W_roof)
        H_aircond_prvrd ~ H_aircond / (1 - W_roof),
        H_aircond_imprvrd ~ H_aircond / (1 - W_roof),
    ]

    return System(eqs, t; name)
end


"""
    AdjustedLayerThickness(; name=:AdjustedLayerThickness)

Computes the adjusted top layer thickness for pervious and impervious road surfaces,
following Eq. 4.30 of Oleson et al. (2010).

The adjustment uses a tunable parameter `c_a = 0.34` to compensate for the difference
between layer-averaged and surface temperature in the numerical scheme.

**Reference**: Oleson et al. (2010), Chapter 4, Eq. 4.30, pp. 98.
"""
@component function AdjustedLayerThickness(; name = :AdjustedLayerThickness)

    @constants begin
        c_a = 0.34, [description = "Tunable parameter from Z.-L. Yang (1998) (dimensionless)"]
    end

    @parameters begin
        z_i, [description = "Node depth of top layer i", unit = u"m"]
        z_h_im1, [description = "Interface depth above layer i (z_{h,i-1})", unit = u"m"]
        z_ip1, [description = "Node depth of layer i+1", unit = u"m"]
    end

    @variables begin
        Δz_star(t), [description = "Adjusted top layer thickness (Eq. 4.30)", unit = u"m"]
    end

    eqs = [
        # Eq. 4.30: Δz_{i*} = 0.5 * [z_i - z_{h,i-1} + c_a*(z_{i+1} - z_{h,i-1})]
        Δz_star ~ 0.5 * (z_i - z_h_im1 + c_a * (z_ip1 - z_h_im1)),
    ]

    return System(eqs, t; name)
end


"""
    HeatingCoolingFlux(; name=:HeatingCoolingFlux)

Computes the heating and cooling fluxes applied to roofs, sunlit walls, and shaded walls,
following Eqs. 4.51–4.54 of Oleson et al. (2010).

Heating is applied when the internal building temperature `T_iB` falls below the
prescribed minimum `T_{iB,min}`. Cooling is applied when `T_iB` exceeds the prescribed
maximum `T_{iB,max}`.

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.51–4.54, pp. 101.
"""
@component function HeatingCoolingFlux(; name = :HeatingCoolingFlux)

    @constants begin
        α_CN = 0.5, [description = "Crank-Nicholson weight (dimensionless)"]
        zero_flux = 0.0, [description = "Zero flux reference", unit = u"W/m^2"]
    end

    @parameters begin
        T_iB, [description = "Internal building temperature (Eq. 4.37)", unit = u"K"]
        T_iB_min, [description = "Minimum prescribed building temperature", unit = u"K"]
        T_iB_max, [description = "Maximum prescribed building temperature", unit = u"K"]
        F_bottom_n, [description = "Bottom heat flux at time n (Eq. 4.53)", unit = u"W/m^2"]
        F_bottom_np1, [description = "Bottom heat flux at time n+1 (Eq. 4.54)", unit = u"W/m^2"]
    end

    @variables begin
        F_combined(t), [description = "Combined bottom flux (α*F^n + (1-α)*F^{n+1})", unit = u"W/m^2"]
        F_heat(t), [description = "Heating flux (Eq. 4.51)", unit = u"W/m^2"]
        F_cool(t), [description = "Cooling flux (Eq. 4.52)", unit = u"W/m^2"]
    end

    eqs = [
        # Combined flux: α*F^n + (1-α)*F^{n+1}
        F_combined ~ α_CN * F_bottom_n + (1 - α_CN) * F_bottom_np1,

        # Eq. 4.51: F_heat = |combined flux| if T_iB < T_min, else 0
        F_heat ~ ifelse(T_iB < T_iB_min, abs(F_combined), zero_flux),

        # Eq. 4.52: F_cool = |combined flux| if T_iB > T_max, else 0
        F_cool ~ ifelse(T_iB > T_iB_max, abs(F_combined), zero_flux),
    ]

    return System(eqs, t; name)
end


"""
    PhaseChangeAdjustment(; name=:PhaseChangeAdjustment)

Performs the ice/liquid water mass adjustment and temperature correction after
phase change, following Eqs. 4.60–4.65 of Oleson et al. (2010).

Given the excess/deficit energy `H_i` from Eq. 4.59, this component computes
the melt/freeze amount `H_m`, adjusts ice mass, conserves liquid water, computes
residual energy `H_{i*}`, and corrects the temperature.

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.60–4.65, pp. 103–104.
"""
@component function PhaseChangeAdjustment(; name = :PhaseChangeAdjustment, layer_type = :interior)

    @constants begin
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
        L_f = 3.337e5, [description = "Latent heat of fusion (Table 1.4)", unit = u"J/kg"]
        zero_kgm2 = 0.0, [description = "Zero mass reference", unit = u"kg/m^2"]
    end

    @parameters begin
        H_i, [description = "Excess/deficit energy for phase change (Eq. 4.59)", unit = u"W/m^2"]
        w_ice_n, [description = "Ice mass at time n", unit = u"kg/m^2"]
        w_liq_n, [description = "Liquid water mass at time n", unit = u"kg/m^2"]
        Δt, [description = "Time step", unit = u"s"]
        c_i, [description = "Volumetric heat capacity", unit = u"J/(m^3*K)"]
        Δz_i, [description = "Layer thickness", unit = u"m"]
    end

    @variables begin
        H_m(t), [description = "Mass change potential H_i*Δt/L_f (Eq. 4.60)", unit = u"kg/m^2"]
        w_ice_np1(t), [description = "Ice mass at time n+1 (Eqs. 4.60-4.62)", unit = u"kg/m^2"]
        w_liq_np1(t), [description = "Liquid water mass at time n+1 (Eq. 4.63)", unit = u"kg/m^2"]
        H_residual(t), [description = "Residual energy after phase change (Eq. 4.64)", unit = u"W/m^2"]
        T_np1(t), [description = "Corrected temperature (Eq. 4.65)", unit = u"K"]
    end

    if layer_type == :top
        @parameters begin
            dh_dT, [description = "Derivative of surface heat flux", unit = u"W/(m^2*K)"]
        end

        eqs = [
            # H_m = H_i * Δt / L_f
            H_m ~ H_i * Δt / L_f,

            # Eq. 4.60 (melting) / Eq. 4.61 (freezing): adjust ice mass
            # For melting (H_m > 0): w_ice^{n+1} = max(w_ice^n - H_m, 0)
            # For freezing (H_m < 0): w_ice^{n+1} = min(w_liq^n + w_ice^n, w_ice^n - H_m)
            w_ice_np1 ~ ifelse(
                H_m > zero_kgm2,
                max(w_ice_n - H_m, zero_kgm2),
                min(w_liq_n + w_ice_n, w_ice_n - H_m)
            ),

            # Eq. 4.63: w_liq^{n+1} = max(w_liq^n + w_ice^n - w_ice^{n+1}, 0)
            w_liq_np1 ~ max(w_liq_n + w_ice_n - w_ice_np1, zero_kgm2),

            # Eq. 4.64: H_{i*} = H_i - L_f * (w_ice^n - w_ice^{n+1}) / Δt
            H_residual ~ H_i - L_f * (w_ice_n - w_ice_np1) / Δt,

            # Eq. 4.65 (top layer): T^{n+1} = T_f + (Δt/(c*Δz)) * H_{i*} / (1 - Δt/(c*Δz) * ∂h/∂T)
            T_np1 ~ T_f + (Δt / (c_i * Δz_i)) * H_residual / (1 - (Δt / (c_i * Δz_i)) * dh_dT),
        ]
    else
        eqs = [
            H_m ~ H_i * Δt / L_f,

            w_ice_np1 ~ ifelse(
                H_m > zero_kgm2,
                max(w_ice_n - H_m, zero_kgm2),
                min(w_liq_n + w_ice_n, w_ice_n - H_m)
            ),

            w_liq_np1 ~ max(w_liq_n + w_ice_n - w_ice_np1, zero_kgm2),

            H_residual ~ H_i - L_f * (w_ice_n - w_ice_np1) / Δt,

            # Eq. 4.65 (interior layers): T^{n+1} = T_f + (Δt/(c*Δz)) * H_{i*}
            T_np1 ~ T_f + (Δt / (c_i * Δz_i)) * H_residual,
        ]
    end

    return System(eqs, t; name)
end


"""
    SnowMeltNoLayers(; name=:SnowMeltNoLayers)

Handles the special case of snow melt when snow is present (W_sno > 0) but there
are no explicit snow layers (snl = 0), following Eqs. 4.66–4.71 of Oleson et al. (2010).

When the snow mass is too small for explicit snow layers, snow melt is computed
from the excess energy in the top soil layer. Snow mass and depth are reduced
proportionally.

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.66–4.71, pp. 105.
"""
@component function SnowMeltNoLayers(; name = :SnowMeltNoLayers)

    @constants begin
        L_f = 3.337e5, [description = "Latent heat of fusion (Table 1.4)", unit = u"J/kg"]
        zero_kgm2 = 0.0, [description = "Zero mass reference", unit = u"kg/m^2"]
        zero_m = 0.0, [description = "Zero depth reference", unit = u"m"]
    end

    @parameters begin
        H_1, [description = "Excess energy in top soil layer", unit = u"W/m^2"]
        W_sno_n, [description = "Snow mass at time n", unit = u"kg/m^2"]
        z_sno_n, [description = "Snow depth at time n", unit = u"m"]
        Δt, [description = "Time step", unit = u"s"]
    end

    @variables begin
        W_sno_np1(t), [description = "Snow mass at time n+1 (Eq. 4.66)", unit = u"kg/m^2"]
        z_sno_np1(t), [description = "Snow depth at time n+1 (Eq. 4.67)", unit = u"m"]
        H_residual(t), [description = "Residual energy after snow melt (Eq. 4.68)", unit = u"W/m^2"]
        M_1S(t), [description = "Snow melt rate (Eq. 4.70)", unit = u"kg/(m^2*s)"]
        E_p1S(t), [description = "Phase change energy (Eq. 4.71)", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.66: W_sno^{n+1} = max(W_sno^n - H_1*Δt/L_f, 0)
        W_sno_np1 ~ max(W_sno_n - H_1 * Δt / L_f, zero_kgm2),

        # Eq. 4.67: z_sno^{n+1} = (W_sno^{n+1} / W_sno^n) * z_sno^n
        z_sno_np1 ~ ifelse(
            W_sno_n > zero_kgm2,
            (W_sno_np1 / W_sno_n) * z_sno_n,
            zero_m
        ),

        # Eq. 4.68: H_{1*} = H_1 - L_f * (W_sno^n - W_sno^{n+1}) / Δt
        H_residual ~ H_1 - L_f * (W_sno_n - W_sno_np1) / Δt,

        # Eq. 4.70: M_{1S} = max((W_sno^n - W_sno^{n+1}) / Δt, 0)
        M_1S ~ max((W_sno_n - W_sno_np1) / Δt, zero_kgm2 / Δt),

        # Eq. 4.71: E_{p,1S} = L_f * M_{1S}
        E_p1S ~ L_f * M_1S,
    ]

    return System(eqs, t; name)
end


"""
    UniformGrid(; name=:UniformGrid, N=15)

Computes the uniform grid discretization for roofs and walls, following Eqs. 4.5–4.7
of Oleson et al. (2010). Roofs and walls are discretized into `N` layers of equal
thickness, with node depths at layer midpoints and interface depths between layers.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, Eqs. 4.5–4.7, pp. 91.
"""
@component function UniformGrid(; name = :UniformGrid, N = 15)

    @parameters begin
        Δz_total, [description = "Total thickness of the roof or wall (Table 1.3)", unit = u"m"]
    end

    # Note: z_interface has N+1 entries (indices 1 to N+1) corresponding to
    # paper indices i=0,...,N. z_interface[1] = z_{h,0}, z_interface[N+1] = z_{h,N}.
    @variables begin
        z_node[1:N](t), [description = "Node depth of layer i (Eq. 4.5)", unit = u"m"]
        Δz_layer[1:N](t), [description = "Layer thickness (Eq. 4.6)", unit = u"m"]
        z_interface[1:(N + 1)](t), [description = "Interface depth (Eq. 4.7); index j maps to paper i=j-1", unit = u"m"]
    end

    eqs = Equation[]

    for i in 1:N
        # Eq. 4.5: z_i = (i - 0.5) * (Δz / N_levgrnd)
        push!(eqs, z_node[i] ~ (i - 0.5) * (Δz_total / N))
    end

    for i in 1:N
        if i == 1
            # Eq. 4.6: Δz_1 = 0.5 * (z_1 + z_2)
            push!(eqs, Δz_layer[i] ~ 0.5 * (z_node[1] + z_node[2]))
        elseif i == N
            # Eq. 4.6: Δz_N = z_N - z_{N-1}
            push!(eqs, Δz_layer[i] ~ z_node[N] - z_node[N - 1])
        else
            # Eq. 4.6: Δz_i = 0.5 * (z_{i+1} - z_{i-1})
            push!(eqs, Δz_layer[i] ~ 0.5 * (z_node[i + 1] - z_node[i - 1]))
        end
    end

    # Eq. 4.7: Interface depths (shifted by +1 for 1-based indexing)
    # z_{h,0} = 0  →  z_interface[1]
    push!(eqs, z_interface[1] ~ 0 * Δz_total)
    for i in 1:(N - 1)
        # z_{h,i} = 0.5 * (z_i + z_{i+1})  →  z_interface[i+1]
        push!(eqs, z_interface[i + 1] ~ 0.5 * (z_node[i] + z_node[i + 1]))
    end
    # z_{h,N} = z_N + 0.5 * Δz_N  →  z_interface[N+1]
    push!(eqs, z_interface[N + 1] ~ z_node[N] + 0.5 * Δz_layer[N])

    return System(eqs, t; name)
end


"""
    ExponentialGrid(; name=:ExponentialGrid, N=15)

Computes the exponential grid discretization for pervious and impervious road columns,
following Eq. 4.8 and Eqs. 4.6–4.7 of Oleson et al. (2010). The scaling factor
`f_s = 0.025` produces finer resolution near the surface, which is needed because
roads overlie real soil.

Layer thicknesses and interface depths are computed from Eqs. 4.6 and 4.7.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, Eqs. 4.6–4.8, pp. 91–92.
"""
@component function ExponentialGrid(; name = :ExponentialGrid, N = 15)

    @constants begin
        f_s = 0.025, [description = "Scaling factor for road grid (Eq. 4.8)", unit = u"m"]
    end

    # Note: z_interface has N+1 entries (indices 1 to N+1) corresponding to
    # paper indices i=0,...,N. z_interface[1] = z_{h,0}, z_interface[N+1] = z_{h,N}.
    @variables begin
        z_node[1:N](t), [description = "Node depth of layer i (Eq. 4.8)", unit = u"m"]
        Δz_layer[1:N](t), [description = "Layer thickness (Eq. 4.6)", unit = u"m"]
        z_interface[1:(N + 1)](t), [description = "Interface depth (Eq. 4.7); index j maps to paper i=j-1", unit = u"m"]
    end

    eqs = Equation[]

    for i in 1:N
        # Eq. 4.8: z_i = f_s * {exp[0.5*(i - 0.5)] - 1}
        push!(eqs, z_node[i] ~ f_s * (exp(0.5 * (i - 0.5)) - 1))
    end

    for i in 1:N
        if i == 1
            # Eq. 4.6: Δz_1 = 0.5 * (z_1 + z_2)
            push!(eqs, Δz_layer[i] ~ 0.5 * (z_node[1] + z_node[2]))
        elseif i == N
            # Eq. 4.6: Δz_N = z_N - z_{N-1}
            push!(eqs, Δz_layer[i] ~ z_node[N] - z_node[N - 1])
        else
            # Eq. 4.6: Δz_i = 0.5 * (z_{i+1} - z_{i-1})
            push!(eqs, Δz_layer[i] ~ 0.5 * (z_node[i + 1] - z_node[i - 1]))
        end
    end

    # Eq. 4.7: Interface depths (shifted by +1 for 1-based indexing)
    # z_{h,0} = 0  →  z_interface[1]
    push!(eqs, z_interface[1] ~ 0 * f_s)
    for i in 1:(N - 1)
        # z_{h,i} = 0.5 * (z_i + z_{i+1})  →  z_interface[i+1]
        push!(eqs, z_interface[i + 1] ~ 0.5 * (z_node[i] + z_node[i + 1]))
    end
    # z_{h,N} = z_N + 0.5 * Δz_N  →  z_interface[N+1]
    push!(eqs, z_interface[N + 1] ~ z_node[N] + 0.5 * Δz_layer[N])

    return System(eqs, t; name)
end


"""
    FreezingPointDepression(; name=:FreezingPointDepression)

Computes the maximum liquid water content in a soil layer when the temperature is
below the freezing point, following Eq. 4.58 of Oleson et al. (2010).

The concept of supercooled soil water from Niu and Yang (2006) is used, where liquid
water coexists with ice below freezing through a freezing point depression equation.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.2, Eq. 4.58, pp. 103.
"""
@component function FreezingPointDepression(; name = :FreezingPointDepression)

    @constants begin
        T_f = 273.15, [description = "Freezing temperature of water (Table 1.4)", unit = u"K"]
        L_f = 3.337e5, [description = "Latent heat of fusion (Table 1.4)", unit = u"J/kg"]
        g = 9.80616, [description = "Gravitational acceleration (Table 1.4)", unit = u"m/s^2"]
        one_ref = 1.0, [description = "Dimensionless reference for non-dimensionalization (dimensionless)"]
    end

    @parameters begin
        T_i, [description = "Layer temperature", unit = u"K"]
        θ_sat, [description = "Porosity (dimensionless)"]
        Δz, [description = "Layer thickness", unit = u"m"]
        ψ_sat, [description = "Saturated matric potential (converted to SI)", unit = u"m"]
        B_i, [description = "Clapp and Hornberger exponent (dimensionless)"]
        ρ_liq, [description = "Density of liquid water", unit = u"kg/m^3"]
    end

    @variables begin
        w_liq_max(t), [description = "Maximum liquid water when T < T_f (Eq. 4.58)", unit = u"kg/m^2"]
    end

    eqs = [
        # Eq. 4.58: w_{liq,max,i} = Δz_i * θ_{sat,i} * ρ_liq *
        #   [10^3 * L_f * (T_f - T_i) / (g * T_i * ψ_{sat,i})]^{-1/B_i}
        # The paper uses ψ_sat in mm. Since our ψ_sat is already in SI (m),
        # the 10^3 conversion factor is absorbed. ψ_sat is negative by
        # convention (soil matric potential), so we negate it to keep the
        # base positive for the fractional exponent.
        # Dimensional check (dimensionless inside brackets):
        #   L_f [J/kg] * ΔT [K] / (g [m/s^2] * T [K] * |ψ| [m])
        #   = [m^2/s^2] / [m^2/s^2] = dimensionless
        w_liq_max ~ ifelse(
            T_i < T_f,
            Δz * θ_sat * ρ_liq * (L_f * (T_f - T_i) / (g * T_i * (-ψ_sat)) * one_ref)^(-1 / B_i),
            Δz * θ_sat * ρ_liq
        ),
    ]

    return System(eqs, t; name)
end


"""
    SnowSoilBlendedHeatCapacity(; name=:SnowSoilBlendedHeatCapacity)

Computes the blended heat capacity of the top soil layer when snow is present but
there are no explicit snow layers (snl = 0), following Eq. 4.88 of Oleson et al. (2010).

The heat capacity is the soil heat capacity (from Eq. 4.85) plus an additional term
from the snow mass distributed over the top layer.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.3, Eq. 4.88, pp. 109.
"""
@component function SnowSoilBlendedHeatCapacity(; name = :SnowSoilBlendedHeatCapacity)

    @constants begin
        C_ice = 2117.27, [description = "Specific heat of ice (Table 1.4)", unit = u"J/(kg*K)"]
    end

    @parameters begin
        c_soil, [description = "Soil volumetric heat capacity from Eq. 4.85", unit = u"J/(m^3*K)"]
        W_sno, [description = "Snow mass", unit = u"kg/m^2"]
        Δz, [description = "Top layer thickness", unit = u"m"]
    end

    @variables begin
        c_blended(t), [description = "Blended heat capacity (Eq. 4.88)", unit = u"J/(m^3*K)"]
    end

    eqs = [
        # Eq. 4.88: c_i = c_i* + C_ice * W_sno / Δz_i
        c_blended ~ c_soil + C_ice * W_sno / Δz,
    ]

    return System(eqs, t; name)
end


"""
    LayerPhaseChangeEnergy(; name=:LayerPhaseChangeEnergy)

Computes the phase change energy for a single layer, following Eq. 4.73 of
Oleson et al. (2010). This is used to sum contributions across all layers to
obtain the total phase change energy (Eq. 4.72).

**Reference**: Oleson et al. (2010), Chapter 4, Eq. 4.73, pp. 106.
"""
@component function LayerPhaseChangeEnergy(; name = :LayerPhaseChangeEnergy)

    @constants begin
        L_f = 3.337e5, [description = "Latent heat of fusion (Table 1.4)", unit = u"J/kg"]
    end

    @parameters begin
        w_ice_n, [description = "Ice mass at time n", unit = u"kg/m^2"]
        w_ice_np1, [description = "Ice mass at time n+1", unit = u"kg/m^2"]
        Δt, [description = "Time step", unit = u"s"]
    end

    @variables begin
        E_p(t), [description = "Phase change energy for layer (Eq. 4.73)", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.73: E_{p,i} = L_f * (w_{ice,i}^n - w_{ice,i}^{n+1}) / Δt
        E_p ~ L_f * (w_ice_n - w_ice_np1) / Δt,
    ]

    return System(eqs, t; name)
end


"""
    TotalPhaseChangeEnergy(; name=:TotalPhaseChangeEnergy)

Computes the total phase change energy for the column, following Eq. 4.72 of
Oleson et al. (2010). The total is the sum of the surface snow melt energy (Eq. 4.71)
and the per-layer phase change energies (Eq. 4.73).

**Reference**: Oleson et al. (2010), Chapter 4, Eqs. 4.72–4.73, pp. 105–106.
"""
@component function TotalPhaseChangeEnergy(; name = :TotalPhaseChangeEnergy)

    @parameters begin
        E_p1S, [description = "Surface snow melt energy (Eq. 4.71)", unit = u"W/m^2"]
        E_p_layers, [description = "Sum of per-layer phase change energies (Eq. 4.73)", unit = u"W/m^2"]
    end

    @variables begin
        E_p_total(t), [description = "Total phase change energy (Eq. 4.72)", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.72: E_p = E_{p,1S} + Σ E_{p,i}
        E_p_total ~ E_p1S + E_p_layers,
    ]

    return System(eqs, t; name)
end
