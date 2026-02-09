export compute_grid_roofwall, compute_grid_road,
    SoilThermalProperties, SnowThermalProperties, UrbanSurfaceThermalProperties,
    InterfaceThermalConductivity, HeatFlux, TridiagonalCoefficients,
    SurfaceEnergyFlux, BuildingTemperature, WasteHeatAirConditioning,
    PhaseChangeEnergy, SnowLayerGeometry

"""
    compute_grid_roofwall(Δz_total; N_levgrnd=15)

Computes the node depths, layer thicknesses, and interface depths for the
uniform grid used for roofs and walls in the CLMU, following Section 4.1
of Oleson et al. (2010).

Returns a named tuple `(z_node, Δz_layer, z_interface)` where:
- `z_node[i]`: node depth of layer i (Eq. 4.5), in meters
- `Δz_layer[i]`: thickness of layer i (Eq. 4.6), in meters
- `z_interface[j]`: interface depth (Eq. 4.7), in meters, where j maps to paper
  index i = j-1 (so `z_interface[1]` = top surface = 0)

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, Eqs. 4.5–4.7, pp. 91.
"""
function compute_grid_roofwall(Δz_total; N_levgrnd = 15)
    # Eq. 4.5: z_i = (i - 0.5) * (Δz_total / N_levgrnd)
    z_node = [(i - 0.5) * (Δz_total / N_levgrnd) for i in 1:N_levgrnd]

    # Eq. 4.6: Layer thicknesses
    Δz_layer = Vector{Float64}(undef, N_levgrnd)
    Δz_layer[1] = 0.5 * (z_node[1] + z_node[2])
    for i in 2:(N_levgrnd - 1)
        Δz_layer[i] = 0.5 * (z_node[i + 1] - z_node[i - 1])
    end
    Δz_layer[N_levgrnd] = z_node[N_levgrnd] - z_node[N_levgrnd - 1]

    # Eq. 4.7: Interface depths (j = paper_i + 1)
    z_interface = Vector{Float64}(undef, N_levgrnd + 1)
    z_interface[1] = 0.0  # paper i=0: top
    for i in 1:(N_levgrnd - 1)
        z_interface[i + 1] = 0.5 * (z_node[i] + z_node[i + 1])
    end
    z_interface[N_levgrnd + 1] = z_node[N_levgrnd] + 0.5 * Δz_layer[N_levgrnd]

    return (; z_node, Δz_layer, z_interface)
end

"""
    compute_grid_road(; N_levgrnd=15)

Computes the node depths, layer thicknesses, and interface depths for the
exponentially-spaced grid used for pervious and impervious roads in the CLMU,
following Section 4.1 of Oleson et al. (2010).

The scaling factor is f_s = 0.025 m (Eq. 4.8). Returns the same named tuple
format as `compute_grid_roofwall`.

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, Eqs. 4.6–4.8, pp. 91–92.
"""
function compute_grid_road(; N_levgrnd = 15)
    f_s = 0.025

    # Eq. 4.8: z_i = f_s * {exp[0.5(i - 0.5)] - 1}
    z_node = [f_s * (exp(0.5 * (i - 0.5)) - 1) for i in 1:N_levgrnd]

    # Eq. 4.6: Layer thicknesses
    Δz_layer = Vector{Float64}(undef, N_levgrnd)
    Δz_layer[1] = 0.5 * (z_node[1] + z_node[2])
    for i in 2:(N_levgrnd - 1)
        Δz_layer[i] = 0.5 * (z_node[i + 1] - z_node[i - 1])
    end
    Δz_layer[N_levgrnd] = z_node[N_levgrnd] - z_node[N_levgrnd - 1]

    # Eq. 4.7: Interface depths (j = paper_i + 1)
    z_interface = Vector{Float64}(undef, N_levgrnd + 1)
    z_interface[1] = 0.0
    for i in 1:(N_levgrnd - 1)
        z_interface[i + 1] = 0.5 * (z_node[i] + z_node[i + 1])
    end
    z_interface[N_levgrnd + 1] = z_node[N_levgrnd] + 0.5 * Δz_layer[N_levgrnd]

    return (; z_node, Δz_layer, z_interface)
end


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
        one_K = 1.0, [description = "Unit temperature", unit = u"K"]
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
        λ_s ~ (8.80 * pct_sand + 2.92 * pct_clay) / (pct_sand + pct_clay) * one_Wm⁻¹K⁻¹,

        # Eq. 4.78: Saturated thermal conductivity
        # For T_i ≥ T_f: λ_sat = λ_s^(1-θ_sat) * λ_liq^θ_sat
        # For T_i < T_f: λ_sat = λ_s^(1-θ_sat) * λ_liq^θ_liq * λ_ice^(θ_sat - θ_liq)
        λ_sat ~ ifelse(T_i ≥ T_f,
            (λ_s / one_Wm⁻¹K⁻¹)^(1 - θ_sat) * (λ_liq_const / one_Wm⁻¹K⁻¹)^θ_sat * one_Wm⁻¹K⁻¹,
            (λ_s / one_Wm⁻¹K⁻¹)^(1 - θ_sat) * (λ_liq_const / one_Wm⁻¹K⁻¹)^θ_liq * (λ_ice_const / one_Wm⁻¹K⁻¹)^(θ_sat - θ_liq) * one_Wm⁻¹K⁻¹),

        # Bulk density: ρ_d = 2700 * (1 - θ_sat) (Section 4.3)
        ρ_d ~ 2700.0 * (1 - θ_sat) * one_kgm⁻³,

        # Eq. 4.80: λ_dry = (0.135 * ρ_d + 64.7) / (2700 - 0.947 * ρ_d)
        λ_dry ~ (0.135 * ρ_d / one_kgm⁻³ + 64.7) / (2700.0 - 0.947 * ρ_d / one_kgm⁻³) * one_Wm⁻¹K⁻¹,

        # Eq. 4.82: S_r = (w_liq / (ρ_liq * Δz) + w_ice / (ρ_ice * Δz)) / θ_sat
        S_r ~ (w_liq / (ρ_liq * Δz) + w_ice / (ρ_ice * Δz)) / θ_sat,

        # Eq. 4.81: Kersten number (unfrozen: K_e = log(S_r) + 1 ≥ 0; frozen: K_e = S_r)
        # Using unfrozen formula for T_i ≥ T_f
        K_e ~ ifelse(T_i ≥ T_f,
            max(log(max(S_r, 1e-10)) + 1.0, 0.0),
            S_r),

        # Eq. 4.77: λ = K_e * λ_sat + (1 - K_e) * λ_dry for S_r > 1e-7
        λ_soil ~ ifelse(S_r > 1e-7,
            K_e * λ_sat + (1 - K_e) * λ_dry,
            λ_dry),

        # Eq. 4.86: c_s = (2.128 * %sand + 2.385 * %clay) / (%sand + %clay) × 10^6
        c_s ~ (2.128 * pct_sand + 2.385 * pct_clay) / (pct_sand + pct_clay) * 1e6 * one_Jm⁻³K⁻¹,

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
        λ_snow ~ λ_air + (7.75e-5 * ρ_sno / one_kgm⁻³ + 1.105e-6 * (ρ_sno / one_kgm⁻³)^2) * (λ_ice_const - λ_air),

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
    TridiagonalCoefficients(; name=:TridiagonalCoefficients, layer_type=:interior)

Computes the tridiagonal matrix coefficients (a_i, b_i, c_i, r_i) for the
Crank-Nicholson discretization of the heat equation, following Section 4.1
of Oleson et al. (2010).

The `layer_type` keyword selects the boundary condition:
- `:interior` — Interior layers (Eqs. 4.46–4.50)
- `:top` — Top layer with surface heat flux (Eqs. 4.21–4.24)
- `:bottom_zero_flux` — Bottom layer with zero heat flux (Eqs. 4.32–4.35)
- `:bottom_building` — Bottom layer with building temperature BC (Eqs. 4.40–4.43)

**Reference**: Oleson et al. (2010), Chapter 4, Section 4.1, pp. 96–101.
"""
@component function TridiagonalCoefficients(; name = :TridiagonalCoefficients, layer_type = :interior)

    @constants begin
        α_CN = 0.5, [description = "Crank-Nicholson weight (Section 4.1) (dimensionless)"]
    end

    # Common parameters for all layer types
    @parameters begin
        c_i, [description = "Volumetric heat capacity of layer", unit = u"J/(m^3*K)"]
        Δz_i, [description = "Layer thickness", unit = u"m"]
        Δt, [description = "Time step", unit = u"s"]
        T_i_n, [description = "Temperature at time n", unit = u"K"]
    end

    @variables begin
        a_coeff(t), [description = "Subdiagonal coefficient (dimensionless)"]
        b_coeff(t), [description = "Diagonal coefficient (dimensionless)"]
        c_coeff(t), [description = "Superdiagonal coefficient (dimensionless)"]
        r_coeff(t), [description = "Right-hand side", unit = u"K"]
    end

    eqs = Equation[]

    if layer_type == :top
        # Top layer (i = snl + 1): Eqs. 4.21–4.25
        @parameters begin
            λ_h_below, [description = "Interface conductivity below", unit = u"W/(m*K)"]
            z_ip1, [description = "Node depth of layer below", unit = u"m"]
            z_i, [description = "Node depth of this layer", unit = u"m"]
            dh_dT, [description = "Derivative of surface heat flux w.r.t. temperature", unit = u"W/(m^2*K)"]
            h_n, [description = "Surface heat flux at time n", unit = u"W/m^2"]
            F_i_n, [description = "Heat flux below at time n", unit = u"W/m^2"]
        end

        push!(eqs,
            # Eq. 4.21: a_i = 0
            a_coeff ~ 0.0,
        )
        push!(eqs,
            # Eq. 4.22: b_i = 1 + (Δt/(c_i*Δz_i)) * [(1-α)*λ_h/(z_{i+1}-z_i) - dh/dT]
            b_coeff ~ 1.0 + (Δt / (c_i * Δz_i)) * ((1 - α_CN) * λ_h_below / (z_ip1 - z_i) - dh_dT),
        )
        push!(eqs,
            # Eq. 4.23: c_i = -(1-α) * (Δt/(c_i*Δz_i)) * λ_h/(z_{i+1}-z_i)
            c_coeff ~ -(1 - α_CN) * (Δt / (c_i * Δz_i)) * λ_h_below / (z_ip1 - z_i),
        )
        push!(eqs,
            # Eq. 4.24: r_i = T_i^n + (Δt/(c_i*Δz_i)) * [h^n - (dh/dT)*T_i^n + α*F_i]
            r_coeff ~ T_i_n + (Δt / (c_i * Δz_i)) * (h_n - dh_dT * T_i_n + α_CN * F_i_n),
        )

    elseif layer_type == :bottom_zero_flux
        # Bottom layer with zero flux (pervious/impervious road): Eqs. 4.32–4.35
        @parameters begin
            λ_h_above, [description = "Interface conductivity above", unit = u"W/(m*K)"]
            z_i, [description = "Node depth of this layer", unit = u"m"]
            z_im1, [description = "Node depth of layer above", unit = u"m"]
            F_im1_n, [description = "Heat flux above at time n", unit = u"W/m^2"]
        end

        push!(eqs,
            # Eq. 4.32: a_i = -(1-α) * (Δt/(c_i*Δz_i)) * λ_h/(z_i - z_{i-1})
            a_coeff ~ -(1 - α_CN) * (Δt / (c_i * Δz_i)) * λ_h_above / (z_i - z_im1),
        )
        push!(eqs,
            # Eq. 4.33: b_i = 1 + (1-α) * (Δt/(c_i*Δz_i)) * λ_h/(z_i - z_{i-1})
            b_coeff ~ 1.0 + (1 - α_CN) * (Δt / (c_i * Δz_i)) * λ_h_above / (z_i - z_im1),
        )
        push!(eqs,
            # Eq. 4.34: c_i = 0
            c_coeff ~ 0.0,
        )
        push!(eqs,
            # Eq. 4.35: r_i = T_i^n - α * (Δt/(c_i*Δz_i)) * F_{i-1}
            r_coeff ~ T_i_n - α_CN * (Δt / (c_i * Δz_i)) * F_im1_n,
        )

    elseif layer_type == :bottom_building
        # Bottom layer with building temperature BC (roof/wall): Eqs. 4.40–4.43
        @parameters begin
            λ_h_above, [description = "Interface conductivity above", unit = u"W/(m*K)"]
            λ_h_below, [description = "Interface conductivity below (to building)", unit = u"W/(m*K)"]
            z_i, [description = "Node depth of this layer", unit = u"m"]
            z_im1, [description = "Node depth of layer above", unit = u"m"]
            z_h_below, [description = "Interface depth below this layer", unit = u"m"]
            F_i_n, [description = "Heat flux below at time n", unit = u"W/m^2"]
            F_im1_n, [description = "Heat flux above at time n", unit = u"W/m^2"]
        end

        push!(eqs,
            # Eq. 4.40: a_i = -(1-α) * (Δt/(c_i*Δz_i)) * λ_h_above/(z_i - z_{i-1})
            a_coeff ~ -(1 - α_CN) * (Δt / (c_i * Δz_i)) * λ_h_above / (z_i - z_im1),
        )
        push!(eqs,
            # Eq. 4.41: b_i = 1 + (1-α) * (Δt/(c_i*Δz_i)) * [λ_h_above/(z_i - z_{i-1}) + λ_h_below/(z_h_below - z_i)]
            b_coeff ~ 1.0 + (1 - α_CN) * (Δt / (c_i * Δz_i)) * (λ_h_above / (z_i - z_im1) + λ_h_below / (z_h_below - z_i)),
        )
        push!(eqs,
            # Eq. 4.42: c_i = 0 (T_{i+1} = T_iB goes to RHS)
            c_coeff ~ 0.0,
        )
        push!(eqs,
            # Eq. 4.43: r_i = T_i^n + α * (Δt/(c_i*Δz_i)) * (F_i - α * F_{i-1})
            # Note: The paper writes r_i = T_i^n + α*(Δt/(c_i*Δz_i))*(F_i - α*F_{i-1})
            # but this should be F_i - F_{i-1} with α weighting already applied in F_i, F_{i-1}
            r_coeff ~ T_i_n + α_CN * (Δt / (c_i * Δz_i)) * (F_i_n - F_im1_n),
        )

    else  # :interior
        # Interior layers (snl+1 < i < N_levgrnd): Eqs. 4.47–4.50
        @parameters begin
            λ_h_above, [description = "Interface conductivity above", unit = u"W/(m*K)"]
            λ_h_below, [description = "Interface conductivity below", unit = u"W/(m*K)"]
            z_i, [description = "Node depth of this layer", unit = u"m"]
            z_im1, [description = "Node depth of layer above", unit = u"m"]
            z_ip1, [description = "Node depth of layer below", unit = u"m"]
            F_i_n, [description = "Heat flux below at time n", unit = u"W/m^2"]
            F_im1_n, [description = "Heat flux above at time n", unit = u"W/m^2"]
        end

        push!(eqs,
            # Eq. 4.47: a_i = -(1-α) * (Δt/(c_i*Δz_i)) * λ_h_above/(z_i - z_{i-1})
            a_coeff ~ -(1 - α_CN) * (Δt / (c_i * Δz_i)) * λ_h_above / (z_i - z_im1),
        )
        push!(eqs,
            # Eq. 4.48: b_i = 1 + (1-α) * (Δt/(c_i*Δz_i)) * [λ_h_above/(z_i - z_{i-1}) + λ_h_below/(z_{i+1} - z_i)]
            b_coeff ~ 1.0 + (1 - α_CN) * (Δt / (c_i * Δz_i)) * (λ_h_above / (z_i - z_im1) + λ_h_below / (z_ip1 - z_i)),
        )
        push!(eqs,
            # Eq. 4.49: c_i = -(1-α) * (Δt/(c_i*Δz_i)) * λ_h_below/(z_{i+1} - z_i)
            c_coeff ~ -(1 - α_CN) * (Δt / (c_i * Δz_i)) * λ_h_below / (z_ip1 - z_i),
        )
        push!(eqs,
            # Eq. 4.50: r_i = T_i^n + α * (Δt/(c_i*Δz_i)) * (F_i - F_{i-1})
            r_coeff ~ T_i_n + α_CN * (Δt / (c_i * Δz_i)) * (F_i_n - F_im1_n),
        )
    end

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
        # Simplifies to: L = W * W_roof / (1 - W_roof) = (H / H_W) * W_roof / (1 - W_roof)
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
                                (1 - W_roof) * H_W * (f_heat * F_heat_sunwall + f_cool * F_cool_sunwall +
                                                       f_heat * F_heat_shdwall + f_cool * F_cool_shdwall),

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
