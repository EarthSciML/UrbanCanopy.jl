# Urban Canopy Model

## Overview

The `UrbanCanopyModel` is the top-level composed diagnostic model for the Community Land
Model Urban (CLMU) parameterization. It wires together four subsystem components
([`OfflineModeForcing`](@ref), [`CLMUAtmosphere`](@ref), [`UrbanRadiation`](@ref), and
[`HeatMomentumFluxes`](@ref)) into a single system that can be driven by prescribing
meteorological forcing and surface conditions.

In this diagnostic configuration, surface temperatures (roof, walls, roads) are
prescribed as parameters rather than computed prognostically. The model computes
urban canopy layer (UCL) air temperature, sensible and latent heat fluxes, and
radiation budgets.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp.

```@docs
UrbanCanopyModel
```

## Implementation

### Creating the System

```@example ucm
using ModelingToolkit, UrbanCanopy, DynamicQuantities, DataFrames, Symbolics

sys = UrbanCanopyModel()
nothing # hide
```

### Top-Level Parameters

```@example ucm
params = parameters(sys)
parent_params = filter(p -> !occursin("₊", string(Symbolics.tosymbol(p, escape=false))), params)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in parent_params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in parent_params],
    :Description => [ModelingToolkit.getdescription(p) for p in parent_params]
)
```

### Subsystem Variables

```@example ucm
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Equations

```@example ucm
eqs = equations(sys)
```

## Analysis

### UCL Air Temperature vs Atmospheric Temperature

This example shows how the urban canopy layer air temperature (T\_ac) varies
with atmospheric temperature for different canyon height-to-width ratios.
Surfaces warmer than the atmosphere create the urban heat island effect.

```@example ucm
using OrdinaryDiffEqDefault, Plots

compiled = mtkcompile(sys)

function solve_ucm(compiled, params)
    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)
    return sol
end

# Baseline parameters
function make_params(compiled; T_atm_val=300.0, H_W_val=1.0, S_atm_val=500.0)
    Dict(
        compiled.S_atm => S_atm_val,
        compiled.T_atm => T_atm_val,
        compiled.P_atm => 101325.0,
        compiled.q_atm => 0.01,
        compiled.W_atm => 4.0,
        compiled.P_precip => 0.0,
        compiled.μ_zen => 0.7854,
        compiled.H_W => H_W_val,
        compiled.W_roof => 0.25,
        compiled.H_canyon => 10.0,
        compiled.f_prvrd => 0.5,
        compiled.H_w_est => 5.0,
        compiled.z_0 => 0.5,
        compiled.z_d => 5.0,
        compiled.α_roof_dir_vis => 0.2,
        compiled.α_roof_dir_nir => 0.2,
        compiled.α_roof_dif_vis => 0.2,
        compiled.α_roof_dif_nir => 0.2,
        compiled.α_road_dir_vis => 0.08,
        compiled.α_road_dir_nir => 0.08,
        compiled.α_road_dif_vis => 0.08,
        compiled.α_road_dif_nir => 0.08,
        compiled.α_wall_dir_vis => 0.25,
        compiled.α_wall_dir_nir => 0.25,
        compiled.α_wall_dif_vis => 0.25,
        compiled.α_wall_dif_nir => 0.25,
        compiled.ε_roof => 0.9,
        compiled.ε_road => 0.95,
        compiled.ε_wall => 0.9,
        compiled.T_g_roof => 310.0,
        compiled.T_g_prvrd => 305.0,
        compiled.T_g_imprvrd => 308.0,
        compiled.T_g_sunwall => 312.0,
        compiled.T_g_shdwall => 303.0,
        compiled.q_g_roof => 0.015,
        compiled.q_g_prvrd => 0.012,
        compiled.q_g_imprvrd => 0.014,
        compiled.f_wet_roof => 0.0,
        compiled.f_wet_imprvrd => 0.0,
        compiled.ζ_in => 0.1,
    )
end

T_atm_range = 285.0:1.0:310.0
H_W_vals = [0.5, 1.0, 2.0]

p = plot(xlabel="Atmospheric Temperature (K)", ylabel="UCL Air Temperature T_ac (K)",
         title="UCL Temperature vs Atmospheric Temperature",
         legend=:topleft)
plot!(p, T_atm_range, T_atm_range, label="T_ac = T_atm (1:1 line)", linestyle=:dash, color=:black)

for hw in H_W_vals
    T_ac_vals = Float64[]
    for T_atm_val in T_atm_range
        params = make_params(compiled; T_atm_val=T_atm_val, H_W_val=hw)
        sol = solve_ucm(compiled, params)
        push!(T_ac_vals, sol[compiled.flux.T_ac][end])
    end
    plot!(p, T_atm_range, T_ac_vals, label="H/W = $hw", linewidth=2)
end

p
```

### Diurnal Cycle Heat Map of UCL Air Temperature

This example creates a spatial heat map showing how T\_ac varies over a grid
as solar radiation changes through the day. The grid represents spatial variation
in atmospheric temperature (latitude axis) and urban density (longitude axis,
via H/W ratio).

```@example ucm
T_atm_grid = 290.0:2.0:310.0  # Temperature gradient (rows)
H_W_grid = 0.3:0.1:2.0        # Urban density gradient (columns)

# Simulate a diurnal cycle by varying S_atm (solar radiation proxy)
S_atm_vals = [0.0, 100.0, 300.0, 500.0, 700.0, 500.0, 300.0, 100.0, 0.0]
hour_labels = ["0:00", "3:00", "6:00", "9:00", "12:00", "15:00", "18:00", "21:00", "24:00"]

anim = @animate for (i, S_val) in enumerate(S_atm_vals)
    T_ac_grid = zeros(length(T_atm_grid), length(H_W_grid))
    for (ri, T_val) in enumerate(T_atm_grid)
        for (ci, hw_val) in enumerate(H_W_grid)
            params = make_params(compiled; T_atm_val=T_val, H_W_val=hw_val, S_atm_val=max(S_val, 1.0))
            sol = solve_ucm(compiled, params)
            T_ac_grid[ri, ci] = sol[compiled.flux.T_ac][end]
        end
    end
    heatmap(collect(H_W_grid), collect(T_atm_grid), T_ac_grid,
            xlabel="Canyon H/W Ratio", ylabel="Atmospheric Temperature (K)",
            title="UCL Air Temperature at $(hour_labels[i])",
            colorbar_title="T_ac (K)", clims=(290.0, 315.0),
            color=:thermal)
end

gif(anim, fps=2)
```

### GEOS-FP Coupled System

The `UrbanCanopyModel` can be coupled to GEOS-FP reanalysis data via the
`EarthSciDataExt` package extension. This coupling automatically maps GEOS-FP
surface meteorological fields to the model's forcing parameters and computes
derived quantities (Monin-Obukhov stability, solar zenith angle) from GEOS-FP
friction velocity, sensible heat flux, latitude, longitude, and time.

The following example creates a coupled system over a domain spanning part of
the continental United States and inspects the coupling equations.

```@example ucm
using EarthSciData, EarthSciMLBase, Dates

domain = DomainInfo(
    DateTime(2016, 5, 1),
    DateTime(2016, 5, 4);
    lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
    latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
    levrange = 1:15,
)

csys = couple(
    UrbanCanopyModel(),
    GEOSFP("4x5", domain),
    domain,
)

coupled_sys = convert(System, csys)
nothing # hide
```

#### Coupling Equations

The coupling creates observed equations that map GEOS-FP meteorological
variables to the `UrbanCanopyModel` forcing parameters. Below we extract
the top-level coupling equations showing how each GEOS-FP field drives the
urban canopy model.

```@example ucm
obs = observed(coupled_sys)

# Extract the top-level coupling equations (UrbanCanopyModel parameters driven by GEOS-FP)
coupling_terms = ["S_atm", "T_atm", "P_atm", "q_atm", "W_atm", "P_precip", "ζ_in", "μ_zen"]
coupling_eqs = filter(obs) do eq
    lhs_str = string(eq.lhs)
    occursin("UrbanCanopyModel", lhs_str) &&
        count("₊", lhs_str) == 1 &&
        any(term -> occursin(term, lhs_str), coupling_terms)
end

for eq in coupling_eqs
    println(eq)
end
```

#### Coupled System Variables

```@example ucm
coupled_vars = unknowns(coupled_sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in coupled_vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in coupled_vars],
    :Description => [ModelingToolkit.getdescription(v) for v in coupled_vars]
)
```
