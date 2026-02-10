# Roof, Wall, Road, Snow Temperatures

## Overview

Chapter 4 of the CLMU technical note describes the computation of temperatures for
roof, sunlit wall, shaded wall, pervious road, and impervious road columns with
optional snow overlays. Heat conduction through each surface column is governed by
the 1D heat equation (Eq. 4.4):

``c \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}\!\left[\lambda \frac{\partial T}{\partial z}\right]``

which is discretized automatically using [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl).

The implementation provides modular components for:
1. **PDE heat conduction**: Automatic spatial discretization of the heat equation
   for roofs/walls (uniform grid, Eq. 4.5) and roads (exponential grid, Eq. 4.8)
   using MethodOfLines.jl
2. **Snow layer geometry**: Node depths, thicknesses, and interfaces for up to 5
   snow layers (Eqs. 4.9--4.10)
3. **Thermal properties**: Soil conductivity via Farouki (1981) / Kersten number
   (Eqs. 4.77--4.82), snow conductivity via Jordan (1991) (Eq. 4.83), and
   volumetric heat capacities from de Vries (1963) (Eqs. 4.85--4.87)
4. **Interface conductivity**: Harmonic mean across adjacent layers (Eq. 4.12)
5. **Heat flux**: Fourier's law discretized across interfaces (Eq. 4.11)
6. **Surface energy flux**: Net surface heat flux and its temperature derivative
   for the implicit scheme (Eqs. 4.26--4.29)
7. **Building temperature**: Weighted wall/roof average for interior building
   temperature (Eqs. 4.37--4.38)
8. **Waste heat and air conditioning**: HVAC waste heat distribution and air
   conditioning heat removal (Eqs. 4.55--4.56)
9. **Waste heat allocation**: Distribution of waste heat and AC to pervious and
   impervious road surfaces (Eq. 4.27)
10. **Adjusted layer thickness**: Road top layer thickness adjustment for
    numerical accuracy (Eq. 4.30)
11. **Heating/cooling flux**: HVAC heating and cooling fluxes based on building
    temperature vs prescribed limits (Eqs. 4.51--4.54)
12. **Phase change energy**: Energy excess/deficit for freezing/thawing assessment
    (Eq. 4.59)
13. **Phase change adjustment**: Ice/liquid mass adjustment and temperature
    correction after phase change (Eqs. 4.60--4.65)
14. **Snow melt without layers**: Snow melt when snow is present but has no
    explicit layers (Eqs. 4.66--4.71)
15. **Grid discretization**: Uniform grid for roofs/walls and exponential grid
    for roads (Eqs. 4.5--4.8)
16. **Freezing point depression**: Maximum liquid water content below freezing
    via supercooled soil water (Eq. 4.58)
17. **Snow-soil blended heat capacity**: Top layer heat capacity blending when
    snow is present but has no explicit layers (Eq. 4.88)
18. **Phase change energy**: Per-layer and total phase change energy for
    freezing/thawing assessment (Eqs. 4.72--4.73)

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp. Chapter 4: Roof, Wall, Road, Snow Temperatures (pp. 91--109).

```@docs
SnowLayerGeometry
SoilThermalProperties
SnowThermalProperties
UrbanSurfaceThermalProperties
InterfaceThermalConductivity
HeatFlux
SurfaceEnergyFlux
BuildingTemperature
WasteHeatAirConditioning
WasteHeatAllocation
AdjustedLayerThickness
HeatingCoolingFlux
PhaseChangeEnergy
PhaseChangeAdjustment
SnowMeltNoLayers
UniformGrid
ExponentialGrid
FreezingPointDepression
SnowSoilBlendedHeatCapacity
LayerPhaseChangeEnergy
TotalPhaseChangeEnergy
RoofWallHeatConduction
RoadHeatConduction
```

## Implementation

### Thermal Properties Components

#### Soil Thermal Properties

```@example ch4_temps
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Statistics
using UrbanCanopy

sys_soil = SoilThermalProperties()

vars = unknowns(sys_soil)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
params = parameters(sys_soil)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example ch4_temps
eqs = equations(sys_soil)
```

#### Snow Thermal Properties

```@example ch4_temps
sys_snow = SnowThermalProperties()

vars = unknowns(sys_snow)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_snow)
```

#### Waste Heat Allocation

```@example ch4_temps
sys_wasteheat = WasteHeatAllocation()

vars = unknowns(sys_wasteheat)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_wasteheat)
```

#### Adjusted Layer Thickness

```@example ch4_temps
sys_alt = AdjustedLayerThickness()

vars = unknowns(sys_alt)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_alt)
```

#### Heating/Cooling Flux

```@example ch4_temps
sys_hcf = HeatingCoolingFlux()

vars = unknowns(sys_hcf)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_hcf)
```

#### Phase Change Adjustment

```@example ch4_temps
sys_pca = PhaseChangeAdjustment(; layer_type=:interior)

vars = unknowns(sys_pca)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_pca)
```

#### Snow Melt Without Layers

```@example ch4_temps
sys_sml = SnowMeltNoLayers()

vars = unknowns(sys_sml)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_sml)
```

#### Grid Discretization

```@example ch4_temps
sys_ug = UniformGrid(N=5)

vars = unknowns(sys_ug)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_ug)
```

```@example ch4_temps
sys_eg = ExponentialGrid(N=5)

vars = unknowns(sys_eg)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_eg)
```

#### Freezing Point Depression

```@example ch4_temps
sys_fpd = FreezingPointDepression()

vars = unknowns(sys_fpd)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
params = parameters(sys_fpd)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example ch4_temps
eqs = equations(sys_fpd)
```

#### Snow-Soil Blended Heat Capacity

```@example ch4_temps
sys_ssbhc = SnowSoilBlendedHeatCapacity()

vars = unknowns(sys_ssbhc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_ssbhc)
```

#### Phase Change Energy

```@example ch4_temps
sys_lpe = LayerPhaseChangeEnergy()

vars = unknowns(sys_lpe)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_lpe)
```

```@example ch4_temps
sys_tpce = TotalPhaseChangeEnergy()

vars = unknowns(sys_tpce)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example ch4_temps
eqs = equations(sys_tpce)
```

## Analysis

### Soil Thermal Conductivity vs Saturation (Eqs. 4.77--4.82)

The Farouki (1981) thermal conductivity parameterization uses the Kersten number
``K_e`` to interpolate between dry (``\lambda_{dry}``) and saturated
(``\lambda_{sat}``) conductivity. For unfrozen soil, ``K_e = \max(\log S_r + 1, 0)``
where ``S_r`` is the degree of saturation; for frozen soil, ``K_e = S_r``.

```@example ch4_temps
using OrdinaryDiffEqDefault

sys = SoilThermalProperties()
compiled = mtkcompile(sys)

function solve_soil(compiled, params)
    prob = ODEProblem(compiled, params, (0.0, 1.0))
    solve(prob)
end

function soil_params(compiled; pct_sand=60.0, pct_clay=20.0, θ_sat=0.4, θ_liq=0.0, T_i=293.15, w_liq=0.0, w_ice=0.0, Δz=0.5)
    Dict(
        compiled.pct_sand => pct_sand,
        compiled.pct_clay => pct_clay,
        compiled.θ_sat => θ_sat,
        compiled.θ_liq => θ_liq,
        compiled.T_i => T_i,
        compiled.w_ice => w_ice,
        compiled.w_liq => w_liq,
        compiled.Δz => Δz,
    )
end

# Sweep saturation for unfrozen and frozen conditions
S_r_range = 0.01:0.01:1.0
λ_unfrozen = Float64[]
λ_frozen = Float64[]
K_e_unfrozen = Float64[]
K_e_frozen = Float64[]

Δz = 0.5
θ_sat = 0.4
ρ_liq = 1000.0
ρ_ice = 917.0

for Sr in S_r_range
    # Unfrozen: w_liq = S_r * θ_sat * ρ_liq * Δz, θ_liq = S_r * θ_sat
    w_liq = Sr * θ_sat * ρ_liq * Δz
    p = soil_params(compiled; w_liq=w_liq, w_ice=0.0, θ_liq=Sr*θ_sat, T_i=293.15)
    sol = solve_soil(compiled, p)
    push!(λ_unfrozen, sol[compiled.λ_soil][end])
    push!(K_e_unfrozen, sol[compiled.K_e][end])

    # Frozen: w_ice = S_r * θ_sat * ρ_ice * Δz, θ_liq = 0 (all frozen)
    w_ice = Sr * θ_sat * ρ_ice * Δz
    p = soil_params(compiled; w_liq=0.0, w_ice=w_ice, θ_liq=0.0, T_i=263.15)
    sol = solve_soil(compiled, p)
    push!(λ_frozen, sol[compiled.λ_soil][end])
    push!(K_e_frozen, sol[compiled.K_e][end])
end

using Plots

p1 = plot(S_r_range, K_e_unfrozen, label="Unfrozen (log S_r + 1)", linewidth=2,
    xlabel="Degree of Saturation S_r", ylabel="Kersten Number K_e",
    title="Kersten Number (Eq. 4.81)", legend=:topleft)
plot!(p1, S_r_range, K_e_frozen, label="Frozen (K_e = S_r)", linewidth=2)

p2 = plot(S_r_range, λ_unfrozen, label="Unfrozen (T = 293 K)", linewidth=2,
    xlabel="Degree of Saturation S_r", ylabel="λ (W m⁻¹ K⁻¹)",
    title="Soil Thermal Conductivity (Eq. 4.77)", legend=:topleft)
plot!(p2, S_r_range, λ_frozen, label="Frozen (T = 263 K)", linewidth=2)

p = plot(p1, p2, layout=(2, 1), size=(700, 700))
p
```

The unfrozen Kersten number shows a logarithmic dependence on saturation,
reaching zero around ``S_r \approx 0.37`` (i.e., ``e^{-1}``), below which the
conductivity equals ``\lambda_{dry}``. The frozen Kersten number is simply linear
in saturation. At full saturation, the frozen conductivity exceeds the unfrozen
value because ice has a higher thermal conductivity than liquid water.

### Soil Thermal Conductivity by Soil Type

Different soil compositions (sand/clay fractions) yield different solid-phase
conductivities ``\lambda_s`` (Eq. 4.79) and saturated conductivities ``\lambda_{sat}``
(Eq. 4.78). Sandy soils conduct heat more efficiently than clay soils.

```@example ch4_temps
soil_types = [
    ("Sand (90/5)", 90.0, 5.0),
    ("Sandy loam (60/15)", 60.0, 15.0),
    ("Loam (40/25)", 40.0, 25.0),
    ("Clay loam (30/35)", 30.0, 35.0),
    ("Clay (10/50)", 10.0, 50.0),
]

S_r_range_soils = 0.01:0.02:1.0
p = plot(xlabel="Degree of Saturation S_r", ylabel="λ (W m⁻¹ K⁻¹)",
    title="Soil Conductivity by Type (Eqs. 4.77-4.80)", legend=:topleft)

for (name, sand, clay) in soil_types
    vals = Float64[]
    for Sr in S_r_range_soils
        w_liq = Sr * 0.4 * 1000.0 * 0.5
        pr = soil_params(compiled; pct_sand=sand, pct_clay=clay, w_liq=w_liq, θ_liq=Sr*0.4)
        sol = solve_soil(compiled, pr)
        push!(vals, sol[compiled.λ_soil][end])
    end
    plot!(p, S_r_range_soils, vals, label=name, linewidth=2)
end
p
```

### Snow Thermal Conductivity vs Density (Eq. 4.83)

Snow thermal conductivity follows the Jordan (1991) parameterization, which is
a quadratic function of snow bulk density ``\rho_{sno}`` (Eq. 4.83). As snow
compacts, the increased contact area between ice grains raises the effective
conductivity from near-air values (fresh snow) towards ice-like values (dense
glacial ice).

```@example ch4_temps
sys_sn = SnowThermalProperties()
compiled_sn = mtkcompile(sys_sn)

ρ_range = 50.0:10.0:700.0
λ_vals = Float64[]
c_vals = Float64[]
Δz_sn = 0.2

for ρ in ρ_range
    w_ice = ρ * Δz_sn
    p = Dict(compiled_sn.w_ice => w_ice, compiled_sn.w_liq => 0.0, compiled_sn.Δz => Δz_sn)
    prob = ODEProblem(compiled_sn, p, (0.0, 1.0))
    sol = solve(prob)
    push!(λ_vals, sol[compiled_sn.λ_snow][end])
    push!(c_vals, sol[compiled_sn.c_snow][end])
end

p1 = plot(ρ_range, λ_vals, linewidth=2, label="λ_snow (Jordan 1991)",
    xlabel="Snow Density (kg m⁻³)", ylabel="λ (W m⁻¹ K⁻¹)",
    title="Snow Thermal Conductivity (Eq. 4.83)", legend=:topleft)
hline!(p1, [0.023], label="λ_air", linestyle=:dash, color=:gray)
hline!(p1, [2.29], label="λ_ice", linestyle=:dash, color=:blue)

p2 = plot(ρ_range, c_vals ./ 1e6, linewidth=2, label="c_snow",
    xlabel="Snow Density (kg m⁻³)", ylabel="c (MJ m⁻³ K⁻¹)",
    title="Snow Heat Capacity (Eq. 4.87)", legend=:topleft, color=:red)

p = plot(p1, p2, layout=(2, 1), size=(700, 700))
p
```

Fresh snow (``\rho \approx 50`` kg m⁻³) has a conductivity only slightly above air,
making it an effective insulator. As snow compacts beyond 400 kg m⁻³, the
quadratic term dominates and conductivity increases rapidly. The heat capacity
scales linearly with density since it depends on the mass of ice per unit volume.

### Interface Thermal Conductivity (Eq. 4.12)

The interface thermal conductivity is a harmonic mean weighted by the distances
from each node to the interface. When the interface is at the midpoint, this
reduces to the standard harmonic mean of the two layer conductivities.

```@example ch4_temps
sys_itc = InterfaceThermalConductivity()
compiled_itc = mtkcompile(sys_itc)

λ_2_range = 0.1:0.1:10.0
λ_int_midpoint = Float64[]
λ_int_offset = Float64[]

for λ2 in λ_2_range
    # Midpoint interface
    p = Dict(compiled_itc.λ_i => 1.0, compiled_itc.λ_ip1 => λ2,
        compiled_itc.z_i => 0.0, compiled_itc.z_ip1 => 1.0, compiled_itc.z_h => 0.5)
    prob = ODEProblem(compiled_itc, p, (0.0, 1.0))
    sol = solve(prob)
    push!(λ_int_midpoint, sol[compiled_itc.λ_interface][end])

    # Interface at 1/4 point
    p2 = Dict(compiled_itc.λ_i => 1.0, compiled_itc.λ_ip1 => λ2,
        compiled_itc.z_i => 0.0, compiled_itc.z_ip1 => 1.0, compiled_itc.z_h => 0.25)
    prob2 = ODEProblem(compiled_itc, p2, (0.0, 1.0))
    sol2 = solve(prob2)
    push!(λ_int_offset, sol2[compiled_itc.λ_interface][end])
end

p = plot(λ_2_range, λ_int_midpoint, label="Midpoint (z_h = 0.5)", linewidth=2,
    xlabel="λ₂ (W m⁻¹ K⁻¹)", ylabel="λ_interface (W m⁻¹ K⁻¹)",
    title="Interface Conductivity (Eq. 4.12), λ₁ = 1", legend=:topleft)
plot!(p, λ_2_range, λ_int_offset, label="Offset (z_h = 0.25)", linewidth=2)
plot!(p, λ_2_range, λ_2_range, label="Arithmetic mean", linestyle=:dash, color=:gray)
plot!(p, λ_2_range, 2 .* λ_2_range ./ (1 .+ λ_2_range), label="Harmonic mean", linestyle=:dot, color=:red)
p
```

The interface conductivity always lies below the arithmetic mean, reflecting
the physical principle that heat flow is limited by the less conductive layer.
When the interface is closer to layer 1 (offset case), the result is weighted
more heavily toward ``\lambda_1``.

### Building Temperature (Eqs. 4.37--4.38)

The internal building temperature is a weighted average of the inner wall and
roof temperatures, with weights proportional to the surface areas exposed to
the building interior. The roof contribution depends on the roof fraction
``W_{roof}`` and the canyon height-to-width ratio.

```@example ch4_temps
sys_bt = BuildingTemperature()
compiled_bt = mtkcompile(sys_bt)

W_roof_range = 0.1:0.02:0.9
T_iB_vals = Float64[]

for wr in W_roof_range
    p = Dict(
        compiled_bt.H_canyon => 10.0, compiled_bt.H_W => 1.0, compiled_bt.W_roof => wr,
        compiled_bt.T_inner_shdwall => 293.0, compiled_bt.T_inner_sunwall => 295.0,
        compiled_bt.T_inner_roof => 290.0,
    )
    prob = ODEProblem(compiled_bt, p, (0.0, 1.0))
    sol = solve(prob)
    push!(T_iB_vals, sol[compiled_bt.T_iB_unclamped][end])
end

p = plot(W_roof_range, T_iB_vals, linewidth=2, label="T_iB",
    xlabel="Roof Fraction W_roof", ylabel="Building Temperature (K)",
    title="Internal Building Temperature (Eqs. 4.37-4.38)", legend=:topright)
hline!(p, [293.0], label="T_shdwall", linestyle=:dash, color=:blue)
hline!(p, [295.0], label="T_sunwall", linestyle=:dash, color=:red)
hline!(p, [290.0], label="T_roof", linestyle=:dash, color=:green)
p
```

As the roof fraction increases, the building temperature shifts from the
wall-dominated average towards the roof temperature. At small roof fractions,
the two wall surfaces dominate and ``T_{iB}`` lies between ``T_{sunwall}`` and
``T_{shdwall}``.

### Grid Discretization (Eqs. 4.5--4.8)

Roofs and walls use a uniform grid (Eq. 4.5) while roads use an exponential
grid (Eq. 4.8) that provides finer resolution near the surface. The exponential
scaling factor ``f_s = 0.025`` m concentrates layers where the diurnal thermal
wave has the largest amplitude.

```@example ch4_temps
using OrdinaryDiffEqDefault

N = 10
sys_ug = UniformGrid(N=N)
compiled_ug = mtkcompile(sys_ug)
prob_ug = ODEProblem(compiled_ug, [compiled_ug.Δz_total => 0.5], (0.0, 1.0))
sol_ug = solve(prob_ug)

sys_eg = ExponentialGrid(N=N)
compiled_eg = mtkcompile(sys_eg)
prob_eg = ODEProblem(compiled_eg, [], (0.0, 1.0))
sol_eg = solve(prob_eg)

obs_ug = Dict(string(o.lhs) => o.lhs for o in observed(compiled_ug))
obs_eg = Dict(string(o.lhs) => o.lhs for o in observed(compiled_eg))

z_uniform = [sol_ug[obs_ug["z_node[$i](t)"]][end] for i in 1:N]
z_exp = [sol_eg[obs_eg["z_node[$i](t)"]][end] for i in 1:N]

p = plot(1:N, z_uniform, marker=:circle, label="Uniform (roof/wall)", linewidth=2,
    xlabel="Layer Index", ylabel="Node Depth (m)",
    title="Grid Node Depths (Eqs. 4.5, 4.8)", legend=:topleft)
plot!(p, 1:N, z_exp, marker=:square, label="Exponential (road)", linewidth=2)
p
```

The uniform grid spaces nodes equally, while the exponential grid clusters
nodes near the surface (small indices) and spaces them progressively farther
apart at depth, providing better resolution where surface forcing varies most.

### Freezing Point Depression (Eq. 4.58)

Below the freezing point, liquid water coexists with ice in soil through a
freezing point depression mechanism. The maximum liquid water content decreases
with temperature according to the matric potential relationship.

```@example ch4_temps
sys_fpd = FreezingPointDepression()
compiled_fpd = mtkcompile(sys_fpd)

T_range = 253.15:0.5:273.15
w_liq_vals = Float64[]

for T in T_range
    p = Dict(
        compiled_fpd.T_i => T,
        compiled_fpd.θ_sat => 0.4,
        compiled_fpd.Δz => 0.5,
        compiled_fpd.ψ_sat => -0.1,
        compiled_fpd.B_i => 5.0,
        compiled_fpd.ρ_liq => 1000.0,
    )
    prob = ODEProblem(compiled_fpd, p, (0.0, 1.0))
    sol = solve(prob)
    push!(w_liq_vals, sol[compiled_fpd.w_liq_max][end])
end

p = plot(T_range .- 273.15, w_liq_vals, linewidth=2, label="w_liq_max",
    xlabel="Temperature (°C)", ylabel="Max Liquid Water (kg m⁻²)",
    title="Freezing Point Depression (Eq. 4.58)", legend=:topleft)
vline!(p, [0.0], label="T_f", linestyle=:dash, color=:gray)
p
```

At the freezing point (0°C), the maximum liquid water equals the fully
saturated value ``\Delta z \cdot \theta_{sat} \cdot \rho_{liq}``. As
temperature decreases below freezing, the available liquid water decreases
rapidly following the Clapp and Hornberger (1978) relationship.

### Roof/Wall Heat Conduction (Eq. 4.4)

The 1D heat equation for roof and wall surfaces is discretized using
MethodOfLines.jl. The following example demonstrates the steady-state
temperature profile with a constant surface heat flux of 50 W/m² at the
top and a fixed building temperature of 290 K at the bottom.

At steady state (``\partial T / \partial t = 0``), the heat equation reduces to
``d^2T/dz^2 = 0``, giving a linear temperature profile ``T(z) = T_{iB} + h(L-z)/\lambda``.

```@example ch4_temps
using MethodOfLines, DomainSets

L = 0.3
λ_val = 1.5
h_val = 50.0
T_iB = 290.0

result = RoofWallHeatConduction(;
    Δz_total = L, N_layers = 30,
    λ_val = λ_val, c_val = 2.0e6,
    h_top = h_val, T_bottom = T_iB,
)

# Solve for long time to reach steady state
prob2 = remake(result.prob; tspan = (0.0, 100000.0))
sol = solve(prob2)

T_mat = sol[result.T(result.t_pde, result.z_var)]
z_disc = sol[result.z_var]
T_final = T_mat[end, :]

# Analytical solution
z_analytical = range(0, L, length=100)
T_analytical = T_iB .+ h_val .* (L .- z_analytical) ./ λ_val

p = plot(z_analytical, T_analytical, label="Analytical", linewidth=2, linestyle=:dash,
    xlabel="Depth z (m)", ylabel="Temperature (K)",
    title="Roof/Wall Steady-State Temperature Profile", legend=:topright)
scatter!(p, z_disc, T_final, label="MethodOfLines", markersize=4)
p
```

### Road Heat Conduction with Zero-Flux Bottom (Eq. 4.4)

Road surfaces use the same heat equation but with a zero heat flux (insulated)
boundary condition at the bottom. With zero flux at both boundaries, the
temperature remains uniform at its initial value, demonstrating energy
conservation.

```@example ch4_temps
T_init = 288.15
result_road = RoadHeatConduction(;
    N_layers = 15,
    λ_val = 1.5, c_val = 2.0e6,
    h_top = 0.0, T_init_val = T_init,
)

sol_road = solve(result_road.prob, saveat = 0.1)
T_mat_road = sol_road[result_road.T(result_road.t_pde, result_road.z_var)]

t_disc = sol_road[result_road.t_pde]
avg_T = [mean(T_mat_road[i, :]) for i in 1:length(t_disc)]

using Statistics
p = plot(t_disc, avg_T, linewidth=2, label="Mean temperature",
    xlabel="Time (s)", ylabel="Temperature (K)",
    title="Road Energy Conservation (zero flux BCs)", legend=:topright)
hline!(p, [T_init], label="Initial T", linestyle=:dash, color=:red)
p
```

With insulated (zero-flux) boundary conditions on both sides, the average
temperature remains constant, confirming energy conservation in the
discretized system.
