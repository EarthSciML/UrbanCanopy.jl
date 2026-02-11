# Hydrology

## Overview

Chapter 5 of the CLMU technical note describes the hydrological processes
simulated for the pervious road, roof, and impervious road surfaces. The
hydrology for the pervious road follows CLM4 for bare soil surfaces and includes:
- Snow accumulation and melt, water transfer between snow layers
- Infiltration, evaporation, surface runoff
- Sub-surface drainage, redistribution within the soil column
- Groundwater discharge and recharge

The roof and impervious road are hydrologically inactive except for their
limited capacity to intercept, store (1 kg m``^{-2}``), and evaporate liquid
precipitation and snow. Sunlit and shaded walls are hydrologically inactive.

The implementation is split into modular components covering: snow density,
ice/water content in snow, snow compaction, snow layer combination, soil
hydraulic properties, surface runoff and infiltration, soil water flux,
equilibrium soil moisture, groundwater drainage, water table depth, aquifer
water balance, snow capping runoff, surface layer updates, and water balance
equations for pervious and impervious surfaces.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp. Chapter 5: Hydrology (pp. 110-146).

```@docs
SnowDensity
SnowIceContent
SnowWaterContent
SnowCompaction
SnowLayerCombination
SoilHydraulicProperties
SurfaceRunoffInfiltration
SoilWaterFlux
SoilWaterEquilibrium
GroundwaterDrainage
WaterTableDepth
AquiferWaterBalance
SnowCappingRunoff
SurfaceLayerUpdate
PerviousRoadWaterBalance
ImperviousWaterBalance
```

## Implementation

The hydrology module implements the equations from Chapter 5 as 16 modular
ModelingToolkit components. Each component can be used independently or
composed into larger systems.

### Snow Density

```@example hydrology
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using UrbanCanopy

sys = SnowDensity()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys)
```

### Soil Hydraulic Properties

```@example hydrology
sys_soil = SoilHydraulicProperties()

vars = unknowns(sys_soil)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_soil)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_soil)
```

### Groundwater Drainage

```@example hydrology
sys_gw = GroundwaterDrainage()

vars = unknowns(sys_gw)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_gw)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_gw)
```

### Water Table Depth

```@example hydrology
sys_wt = WaterTableDepth()

vars = unknowns(sys_wt)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_wt)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_wt)
```

## Analysis

The following figures illustrate the key parametric relationships implemented
in the hydrology module, using the component functions directly.

```@example hydrology
using OrdinaryDiffEqDefault
using Plots

nothing # hide
```

### Snow Density vs Temperature (Eq. 5.10)

Bulk density of newly fallen snow as a function of atmospheric temperature,
following Anderson (1976). The density is piecewise: 50 kg/m³ at very cold
temperatures (below ``T_f - 15``), increasing as ``50 + 1.7(T_{atm} - T_f + 15)^{1.5}``
for intermediate temperatures, and capped at approximately 169 kg/m³ above
``T_f + 2`` K (cf. Eq. 5.10, p. 115).

```@example hydrology
sys_sd = SnowDensity()
compiled_sd = mtkcompile(sys_sd)

T_f = 273.15
temps = range(250.0, 280.0, length=200)
densities = Float64[]

prob = ODEProblem(compiled_sd, [compiled_sd.T_atm => 250.0], (0.0, 1.0))
for T in temps
    p = remake(prob; p = [compiled_sd.T_atm => T])
    sol = solve(p)
    push!(densities, sol[compiled_sd.ρ_sno][end])
end

p = plot(temps, densities, linewidth=2, label="ρ_sno",
    xlabel="Atmospheric Temperature (K)",
    ylabel="Snow Density (kg/m³)",
    title="New Snow Density vs Temperature (Eq. 5.10)")
vline!(p, [T_f - 15], linestyle=:dash, color=:gray, label="T_f - 15 K")
vline!(p, [T_f + 2], linestyle=:dash, color=:red, label="T_f + 2 K")
p
```

Below ``T_f - 15`` K (258.15 K), the snow density is constant at 50 kg/m³.
Between ``T_f - 15`` and ``T_f + 2`` K, density increases nonlinearly following
the 1.5-power law. Above ``T_f + 2`` K, density is capped at 50 + 1.7(17)^1.5 ≈ 169 kg/m³.

### Soil Hydraulic Properties vs Sand Content (Eqs. 5.70-5.74)

Soil hydraulic properties as a function of sand and clay content, using the
Clapp and Hornberger (1978) and Cosby et al. (1984) pedotransfer functions.
These are shown for a saturated soil (``\theta = \theta_{sat}``).

```@example hydrology
sys_hp = SoilHydraulicProperties()
compiled_hp = mtkcompile(sys_hp)

sand_range = range(5.0, 95.0, length=100)
ksat_vals = Float64[]
θsat_vals = Float64[]
ψsat_vals = Float64[]
B_vals = Float64[]

pct_clay = 20.0
prob_hp = ODEProblem(compiled_hp,
    [compiled_hp.pct_sand => 50.0, compiled_hp.pct_clay => pct_clay,
     compiled_hp.θ_i => 0.4],
    (0.0, 1.0))

for sand in sand_range
    θ_sat_val = 0.489 - 0.00126 * sand
    p = remake(prob_hp; p = [compiled_hp.pct_sand => sand,
                              compiled_hp.pct_clay => pct_clay,
                              compiled_hp.θ_i => θ_sat_val])
    sol = solve(p)
    push!(ksat_vals, sol[compiled_hp.k_sat][end])
    push!(θsat_vals, sol[compiled_hp.θ_sat][end])
    push!(ψsat_vals, abs(sol[compiled_hp.ψ_sat][end]))
    push!(B_vals, sol[compiled_hp.B][end])
end

p1 = plot(sand_range, ksat_vals .* 1e3, linewidth=2, label=false, yscale=:log10,
    xlabel="Sand Content (%)", ylabel="k_sat (mm/s)",
    title="Saturated Hydraulic Conductivity (Eq. 5.70)")

p2 = plot(sand_range, θsat_vals, linewidth=2, label=false,
    xlabel="Sand Content (%)", ylabel="θ_sat (m³/m³)",
    title="Porosity (Eq. 5.71)")

p3 = plot(sand_range, ψsat_vals .* 1e3, linewidth=2, label=false, yscale=:log10,
    xlabel="Sand Content (%)", ylabel="|ψ_sat| (mm)",
    title="Saturated Matric Potential (Eq. 5.74)")

clay_range = range(5.0, 60.0, length=100)
B_clay = 2.91 .+ 0.159 .* clay_range
p4 = plot(clay_range, B_clay, linewidth=2, label=false,
    xlabel="Clay Content (%)", ylabel="B",
    title="Clapp-Hornberger Exponent (Eq. 5.72)")

p = plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 600))
p
```

Sandy soils have higher saturated hydraulic conductivity (more permeable), lower
porosity, and weaker matric potential (less capillary suction). The
Clapp-Hornberger exponent B increases linearly with clay content, reflecting
the broader pore-size distribution of clay-rich soils. These relationships
correspond to Eqs. 5.70-5.74 of the CLMU technical note.

### Groundwater Drainage vs Water Table Depth (Eq. 5.140)

Sub-surface drainage as an exponential function of water table depth, modified
by the frozen soil impermeable fraction (cf. Eq. 5.140, p. 142).

```@example hydrology
sys_dr = GroundwaterDrainage()
compiled_dr = mtkcompile(sys_dr)

zv_range = range(0.0, 5.0, length=200)
prob_dr = ODEProblem(compiled_dr,
    [compiled_dr.z_v => 0.0, compiled_dr.f_ice_weighted => 0.0],
    (0.0, 1.0))

f_ice_values = [0.0, 0.25, 0.5, 0.75, 1.0]
p = plot(xlabel="Water Table Depth (m)",
    ylabel="Drainage Rate (kg m⁻² s⁻¹)",
    title="Sub-surface Drainage vs Water Table Depth (Eq. 5.140)",
    legend=:topright)

for f_ice in f_ice_values
    drain_vals = Float64[]
    for zv in zv_range
        pr = remake(prob_dr; p = [compiled_dr.z_v => zv,
                                   compiled_dr.f_ice_weighted => f_ice])
        sol = solve(pr)
        push!(drain_vals, sol[compiled_dr.q_drai][end])
    end
    plot!(p, zv_range, drain_vals, linewidth=2,
        label="f_ice = $(f_ice)")
end
p
```

Drainage decreases exponentially with water table depth (decay factor
``f_{drai} = 2.5`` m``^{-1}``), reaching the maximum rate of 5.5×10``^{-3}``
kg m``^{-2}`` s``^{-1}`` at the surface. Frozen soil (higher ``f_{ice}``)
reduces drainage through the impermeable fraction ``f_{imp}`` (Eq. 5.141);
fully frozen soil (``f_{ice} = 1``) eliminates drainage entirely.

### Water Table Depth vs Aquifer Storage (Eq. 5.146)

Water table depth as a function of aquifer water storage, for a soil column
with bottom interface at 3.5 m depth (cf. Eq. 5.146, p. 143).

```@example hydrology
sys_wt = WaterTableDepth()
compiled_wt = mtkcompile(sys_wt)

wa_range = range(0.0, 6000.0, length=200)
zv_vals = Float64[]

prob_wt = ODEProblem(compiled_wt,
    [compiled_wt.z_h_bottom => 3.5, compiled_wt.W_a => 0.0],
    (0.0, 1.0))

for wa in wa_range
    p = remake(prob_wt; p = [compiled_wt.W_a => wa])
    sol = solve(p)
    push!(zv_vals, sol[compiled_wt.z_v][end])
end

p = plot(wa_range, zv_vals, linewidth=2, label=false,
    xlabel="Aquifer Water Storage W_a (kg/m²)",
    ylabel="Water Table Depth z_v (m)",
    title="Water Table Depth vs Aquifer Storage (Eq. 5.146)",
    yflip=true)
hline!(p, [0.05], linestyle=:dash, color=:red, label="z_v min (0.05 m)")
hline!(p, [80.0], linestyle=:dash, color=:blue, label="z_v max (80 m)")
p
```

The water table depth decreases (rises) linearly with increasing aquifer water
storage according to ``z_v = z_{h,bottom} + 25 - W_a / (\rho_{liq} \cdot S_y)``,
with specific yield ``S_y = 0.2``. The water table depth is clamped between
0.05 m and 80 m. For the configuration shown (``z_{h,bottom} = 3.5`` m),
at ``W_a = 0`` the water table is at 28.5 m depth, and it reaches the minimum
depth of 0.05 m when ``W_a \approx 5690`` kg/m².

### Snow Layer Combination - Enthalpy Conservation (Eq. 5.41)

Temperature of a combined snow layer from two layers with different
temperatures, demonstrating enthalpy conservation (Eqs. 5.38-5.42, p. 122).

```@example hydrology
sys_comb = SnowLayerCombination()
compiled_comb = mtkcompile(sys_comb)

T_f = 273.15
T2_range = range(245.0, 273.0, length=100)
Tc_vals = Float64[]

prob_comb = ODEProblem(compiled_comb,
    [compiled_comb.Δz_1 => 0.05, compiled_comb.Δz_2 => 0.1,
     compiled_comb.w_liq_1 => 0.5, compiled_comb.w_liq_2 => 1.0,
     compiled_comb.w_ice_1 => 3.0, compiled_comb.w_ice_2 => 5.0,
     compiled_comb.T_1 => 270.0, compiled_comb.T_2 => 260.0],
    (0.0, 1.0))

for T2 in T2_range
    p = remake(prob_comb; p = [compiled_comb.T_2 => T2])
    sol = solve(p)
    push!(Tc_vals, sol[compiled_comb.T_c][end])
end

p = plot(T2_range, Tc_vals, linewidth=2, label="T_combined",
    xlabel="Layer 2 Temperature (K)",
    ylabel="Combined Temperature (K)",
    title="Snow Layer Combination Temperature (Eq. 5.41)")
hline!(p, [270.0], linestyle=:dash, color=:gray, label="T_1 = 270 K")
plot!(p, T2_range, T2_range, linestyle=:dot, color=:gray, label="T_2")
p
```

The combined temperature lies between the two layer temperatures, weighted by
the thermal mass (ice and liquid heat capacities) of each layer. Because layer 2
has more mass (5 kg/m² ice + 1 kg/m² liquid vs 3 + 0.5 for layer 1), the
combined temperature is closer to ``T_2`` than ``T_1``.
