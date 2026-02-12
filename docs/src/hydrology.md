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
RichardsEquation
GroundwaterDrainage
WaterTableDepth
AquiferWaterBalance
SnowCappingRunoff
SurfaceLayerUpdate
PerviousRoadWaterBalance
ImperviousWaterBalance
ImperviousRunoff
InterfaceHydraulicConductivity
SoilWaterContentCalc
```

## Implementation

The hydrology module implements the equations from Chapter 5 as 17 modular
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
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys)
```

### Snow Ice Content

```@example hydrology
sys_ice = SnowIceContent()

vars = unknowns(sys_ice)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_ice)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_ice)
```

### Snow Water Content

```@example hydrology
sys_swc = SnowWaterContent()

vars = unknowns(sys_swc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_swc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_swc)
```

### Snow Compaction

```@example hydrology
sys_sc = SnowCompaction()

vars = unknowns(sys_sc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_sc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_sc)
```

### Snow Layer Combination

```@example hydrology
sys_slc = SnowLayerCombination()

vars = unknowns(sys_slc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_slc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_slc)
```

### Soil Hydraulic Properties

```@example hydrology
sys_soil = SoilHydraulicProperties()

vars = unknowns(sys_soil)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_soil)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
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
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_gw)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
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
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_wt)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_wt)
```

### Richards' Equation (PDE Soil Water Movement)

Richards' equation is implemented as a PDE using MethodOfLines.jl for automatic
spatial discretization. The function returns a named tuple containing the
discretized ODE problem and symbolic variables.

```@example hydrology
using MethodOfLines, DomainSets

result_re = RichardsEquation(;
    N_layers = 10, Δz_total = 2.0,
    pct_sand = 50.0, pct_clay = 20.0,
    θ_top_val = 0.3, θ_bottom_val = 0.3,
    θ_init_val = 0.3,
)
println("Soil properties for 50% sand, 20% clay:")
println("  k_sat = ", result_re.k_sat_val, " m/s")
println("  θ_sat = ", result_re.θ_sat_val, " m³/m³")
println("  B     = ", result_re.B_val)
println("  ψ_sat = ", result_re.ψ_sat_val, " m")
```

### Surface Runoff and Infiltration

```@example hydrology
sys_sri = SurfaceRunoffInfiltration()

vars = unknowns(sys_sri)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_sri)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_sri)
```

### Soil Water Flux

```@example hydrology
sys_swf = SoilWaterFlux()

vars = unknowns(sys_swf)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_swf)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_swf)
```

### Soil Water Equilibrium

```@example hydrology
sys_swe = SoilWaterEquilibrium()

vars = unknowns(sys_swe)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_swe)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_swe)
```

### Aquifer Water Balance

```@example hydrology
sys_awb = AquiferWaterBalance()

vars = unknowns(sys_awb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_awb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_awb)
```

### Snow Capping Runoff

```@example hydrology
sys_scr = SnowCappingRunoff()

vars = unknowns(sys_scr)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_scr)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_scr)
```

### Surface Layer Update

```@example hydrology
sys_slu = SurfaceLayerUpdate()

vars = unknowns(sys_slu)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_slu)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_slu)
```

### Pervious Road Water Balance

```@example hydrology
sys_prwb = PerviousRoadWaterBalance()

vars = unknowns(sys_prwb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_prwb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_prwb)
```

### Impervious Water Balance

```@example hydrology
sys_iwb = ImperviousWaterBalance()

vars = unknowns(sys_iwb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_iwb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_iwb)
```

### Impervious Runoff

```@example hydrology
sys_ir = ImperviousRunoff()

vars = unknowns(sys_ir)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_ir)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_ir)
```

### Interface Hydraulic Conductivity

```@example hydrology
sys_ihc = InterfaceHydraulicConductivity()

vars = unknowns(sys_ihc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_ihc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_ihc)
```

### Soil Water Content

```@example hydrology
sys_swcc = SoilWaterContentCalc()

vars = unknowns(sys_swcc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example hydrology
params = parameters(sys_swcc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [DynamicQuantities.dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example hydrology
equations(sys_swcc)
```

### Table 5.1: Snow Layer Thickness Bounds

The minimum and maximum thicknesses for the five snow layers used in the model
(Table 5.1, p. 121):

| Layer | ``\Delta z_{min}`` (m) | ``(\Delta z_{max})_l`` (m) | ``(\Delta z_{max})_u`` (m) |
|:-----:|:----------------------:|:--------------------------:|:--------------------------:|
|   1 (top)   | 0.010 | 0.03 | 0.02 |
|   2   | 0.015 | 0.07 | 0.05 |
|   3   | 0.025 | 0.18 | 0.11 |
|   4   | 0.055 | 0.41 | 0.23 |
|   5 (bottom)  | 0.115 | --- | --- |

Where ``(\Delta z_{max})_l`` and ``(\Delta z_{max})_u`` are the maximum thickness
at the lower bound (when number of layers equals the layer index) and upper bound
(when number of layers exceeds the layer index), respectively.

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

### Richards' Equation - Wetting Front Propagation (Eq. 5.59)

Wetting front propagating downward through a soil column, driven by a wet
Dirichlet boundary at the top and dry initial conditions. The soil column is
discretized using MethodOfLines.jl.

```@example hydrology
using MethodOfLines, DomainSets

θ_sat_val = 0.489 - 0.00126 * 50.0
Δz_total = 2.0
N_layers = 15
θ_dry = 0.1 * θ_sat_val
θ_wet = 0.8 * θ_sat_val

result_infl = RichardsEquation(;
    N_layers = N_layers, Δz_total = Δz_total,
    pct_sand = 50.0, pct_clay = 20.0,
    θ_top_val = θ_wet,
    θ_bottom_val = θ_dry,
    θ_init_val = θ_dry,
)

tspan = (0.0, 50000.0)
prob_infl = remake(result_infl.prob; tspan = tspan)
sol_infl = solve(prob_infl, saveat = 10000.0)
θ_mat_infl = sol_infl[result_infl.θ(result_infl.t_pde, result_infl.z_var)]

# Spatial grid points from the discretization
dz = Δz_total / N_layers
z_grid = range(dz / 2, Δz_total - dz / 2, length = size(θ_mat_infl, 2))

p = plot(xlabel = "Volumetric Water Content (m³/m³)", ylabel = "Depth (m)",
    title = "Wetting Front Propagation (Eq. 5.59)",
    yflip = true, legend = :bottomleft)
for (j, t_val) in enumerate(sol_infl.t)
    plot!(p, θ_mat_infl[j, :], collect(z_grid), linewidth = 2,
        label = "t = $(Int(t_val)) s", marker = :circle, markersize = 3)
end
vline!(p, [θ_sat_val], linestyle = :dash, color = :gray, label = "θ_sat")
p
```

The wetting front propagates downward over time as the wet top boundary drives
water into the dry soil column through both capillary diffusion and gravity.
The front moves faster through sandy soils (higher ``k_{sat}``) and slows as
the soil approaches saturation.

### Richards' Equation - Gravity Drainage (Eq. 5.68)

Starting from a moderately wet uniform initial condition with equal Dirichlet
boundary conditions at top and bottom, the soil column evolves under the combined
effects of nonlinear diffusion and gravity drainage.

```@example hydrology
θ_mid = 0.4 * θ_sat_val

result_grav = RichardsEquation(;
    N_layers = N_layers, Δz_total = Δz_total,
    pct_sand = 50.0, pct_clay = 20.0,
    θ_top_val = θ_mid,
    θ_bottom_val = θ_mid,
    θ_init_val = θ_mid,
)

tspan_grav = (0.0, 200000.0)
prob_grav = remake(result_grav.prob; tspan = tspan_grav)
sol_grav = solve(prob_grav, saveat = 40000.0)
θ_mat_grav = sol_grav[result_grav.θ(result_grav.t_pde, result_grav.z_var)]

p = plot(xlabel = "Volumetric Water Content (m³/m³)", ylabel = "Depth (m)",
    title = "Gravity Drainage (Eq. 5.68)",
    yflip = true, legend = :bottomleft)
for (j, t_val) in enumerate(sol_grav.t)
    plot!(p, θ_mat_grav[j, :], collect(z_grid), linewidth = 2,
        label = "t = $(Int(t_val)) s", marker = :circle, markersize = 3)
end
vline!(p, [θ_sat_val], linestyle = :dash, color = :gray, label = "θ_sat")
p
```

With equal boundary conditions and gravity, the system evolves toward a steady
state where the diffusive flux balances the gravitational flux. The rate of
redistribution depends on the hydraulic conductivity, which is a strong function
of water content through the Clapp-Hornberger relation (Eq. 5.69).

### Interface Hydraulic Conductivity vs Saturation (Eq. 5.69)

Hydraulic conductivity at a layer interface as a function of saturation degree
for different Clapp-Hornberger exponents B, showing the strong nonlinearity
(cf. Eq. 5.69, p. 131). The frozen fraction is set to zero.

```@example hydrology
sys_ihc = InterfaceHydraulicConductivity()
compiled_ihc = mtkcompile(sys_ihc)

k_sat = 1.0e-5  # m/s
θ_sat = 0.4     # m³/m³
sat_range = range(0.05, 1.0, length = 100)

prob_ihc = ODEProblem(compiled_ihc, [
    compiled_ihc.k_sat_h => k_sat, compiled_ihc.θ_upper => 0.1,
    compiled_ihc.θ_lower => 0.1,
    compiled_ihc.θ_sat_upper => θ_sat, compiled_ihc.θ_sat_lower => θ_sat,
    compiled_ihc.B_i => 6.0, compiled_ihc.f_frz_upper => 0.0,
    compiled_ihc.f_frz_lower => 0.0,
], (0.0, 1.0))

p = plot(xlabel = "Saturation Degree (θ/θ_sat)",
    ylabel = "k/k_sat",
    title = "Interface Hydraulic Conductivity (Eq. 5.69)",
    yscale = :log10, legend = :bottomright)

for B in [3.0, 6.0, 9.0, 12.0]
    k_vals = Float64[]
    for s in sat_range
        θ_val = s * θ_sat
        pr = remake(prob_ihc; p = [
            compiled_ihc.θ_upper => θ_val, compiled_ihc.θ_lower => θ_val,
            compiled_ihc.B_i => B,
        ])
        sol = solve(pr)
        push!(k_vals, sol[compiled_ihc.k_h][end] / k_sat)
    end
    plot!(p, collect(sat_range), k_vals, linewidth = 2, label = "B = $(Int(B))")
end
p
```

Higher Clapp-Hornberger exponents (clay-rich soils) produce steeper conductivity
curves, reflecting the narrower pore-size distribution and stronger moisture
dependence of hydraulic conductivity. At full saturation (``\theta/\theta_{sat} = 1``),
the conductivity equals ``k_{sat}`` regardless of B.

### Impervious Surface Runoff vs Precipitation (Eq. 5.47)

Runoff from impervious surfaces (roof and impervious road) as a function of
precipitation rate, for different initial surface water contents. The ponding
limit is 1 kg/m² (cf. Eqs. 5.47-5.48, pp. 124-125).

```@example hydrology
sys_ir = ImperviousRunoff()
compiled_ir = mtkcompile(sys_ir)

q_range = range(0.0, 5.0, length = 200)

prob_ir = ODEProblem(compiled_ir, [
    compiled_ir.w_liq_1 => 0.5, compiled_ir.q_liq_0 => 0.0,
    compiled_ir.q_seva => 0.0, compiled_ir.has_snow => 0.0,
], (0.0, 1.0))

p = plot(xlabel = "Liquid Precipitation Rate (kg m⁻² s⁻¹)",
    ylabel = "Surface Runoff (kg m⁻² s⁻¹)",
    title = "Impervious Surface Runoff (Eq. 5.47)",
    legend = :topleft)

for w0 in [0.0, 0.5, 1.0, 1.5]
    runoff_vals = Float64[]
    for q in q_range
        pr = remake(prob_ir; p = [compiled_ir.w_liq_1 => w0, compiled_ir.q_liq_0 => q])
        sol = solve(pr)
        push!(runoff_vals, sol[compiled_ir.q_over][end])
    end
    plot!(p, collect(q_range), runoff_vals, linewidth = 2,
        label = "w_liq = $(w0) kg/m²")
end
plot!(p, q_range, q_range, linestyle = :dash, color = :gray, label = "1:1 line")
p
```

When surface water exceeds the ponding limit (1 kg/m²), any additional input
becomes runoff immediately. Higher initial surface water content (``w_{liq,1}``)
shifts the runoff curve to the left, reflecting less available storage.
The dashed 1:1 line represents the theoretical maximum when all precipitation
becomes runoff.
