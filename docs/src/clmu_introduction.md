# CLMU Introduction: Atmospheric Interface

## Overview

The Community Land Model Urban (CLMU) parameterization provides a physically-based
representation of urban surfaces for use within the Community Land Model (CLM4).
Chapter 1 of the technical note introduces the model structure and defines the
atmospheric interface: the forcing variables provided by the atmospheric model,
the diagnostic equations for air density and vapor pressure, and the reference
height computation.

The urban surface is represented using the canyon concept (Oke, 1987), where the
considerable complexity of the urban surface is reduced to a single urban canyon
consisting of a canyon floor of width W bordered by two facing buildings of height H.
The urban canyon is divided into five columns: roof, sunlit wall, shaded wall,
pervious road, and impervious road.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp.

```@docs
CLMUAtmosphere
```

## Implementation

The `CLMUAtmosphere` component implements three diagnostic equations from Chapter 1:

1. **Atmospheric vapor pressure** from specific humidity and pressure
2. **Air density** from pressure, vapor pressure, and temperature
3. **Reference height** for flux computations

### State Variables

```@example clmu_intro
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using UrbanCanopy

sys = CLMUAtmosphere()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example clmu_intro
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example clmu_intro
eqs = equations(sys)
```

## Analysis

### Air Density vs. Temperature

The air density equation from Chapter 1 (p. 17) computes atmospheric density as a
function of pressure, temperature, and vapor pressure. For dry air, this reduces to
the ideal gas law: ``\rho = P / (R_{da} T)``. The following figure shows how air
density varies with temperature at standard pressure for both dry and moist air.

```@example clmu_intro
using OrdinaryDiffEqDefault
using Plots

compiled = mtkcompile(sys)

T_range = 250.0:5.0:320.0
ρ_dry = Float64[]
ρ_moist = Float64[]

for T_val in T_range
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.q_atm => 0.0,
            compiled.P_atm => 101325.0,
            compiled.T_atm => T_val,
            compiled.z_prime_atm => 30.0,
            compiled.z_0 => 1.0,
            compiled.z_d => 5.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    push!(ρ_dry, sol[compiled.ρ_atm][end])

    prob_m = ODEProblem(
        compiled,
        Dict(
            compiled.q_atm => 0.015,
            compiled.P_atm => 101325.0,
            compiled.T_atm => T_val,
            compiled.z_prime_atm => 30.0,
            compiled.z_0 => 1.0,
            compiled.z_d => 5.0,
        ),
        (0.0, 1.0),
    )
    sol_m = solve(prob_m)
    push!(ρ_moist, sol_m[compiled.ρ_atm][end])
end

p = plot(T_range, ρ_dry, label="Dry air (q=0)", xlabel="Temperature (K)",
    ylabel="Air density (kg/m³)", title="Air Density vs Temperature at P=101325 Pa",
    linewidth=2, legend=:topright)
plot!(p, T_range, ρ_moist, label="Moist air (q=0.015 kg/kg)", linewidth=2, linestyle=:dash)
p
```

At all temperatures, moist air is less dense than dry air at the same pressure, consistent
with the lower molecular weight of water vapor (18.016 kg/kmol) compared to dry air
(28.966 kg/kmol).

### Vapor Pressure vs. Specific Humidity

The vapor pressure equation relates the atmospheric vapor pressure to specific humidity
and total pressure. This relationship is approximately linear for small specific humidities
but becomes nonlinear at higher moisture content.

```@example clmu_intro
q_range = 0.0:0.001:0.03
e_vals = Float64[]

for q_val in q_range
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.q_atm => q_val,
            compiled.P_atm => 101325.0,
            compiled.T_atm => 293.15,
            compiled.z_prime_atm => 30.0,
            compiled.z_0 => 1.0,
            compiled.z_d => 5.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    push!(e_vals, sol[compiled.e_atm][end])
end

p = plot(q_range .* 1000, e_vals, xlabel="Specific humidity (g/kg)",
    ylabel="Vapor pressure (Pa)",
    title="Atmospheric Vapor Pressure vs Specific Humidity\nat P=101325 Pa",
    linewidth=2, legend=false)
p
```
