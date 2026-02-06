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
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example clmu_intro
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example clmu_intro
eqs = equations(sys)
```

### Table 1.1: Atmospheric Input to Urban Model

The atmospheric model provides the following forcing variables to the urban model at each time step (Table 1.1, p. 16). The reference height ``z'_{atm}`` is the height above the surface defined as the roughness length ``z_0`` plus displacement height ``z_d``, so the reference height used for flux computations is ``z_{atm} = z'_{atm} + z_0 + z_d``.

| Variable | Symbol | Units |
|----------|--------|-------|
| Reference height | ``z'_{atm}`` | m |
| Zonal wind at ``z_{atm}`` | ``u_{atm}`` | m/s |
| Meridional wind at ``z_{atm}`` | ``v_{atm}`` | m/s |
| Potential temperature | ``\theta_{atm}`` | K |
| Specific humidity at ``z_{atm}`` | ``q_{atm}`` | kg/kg |
| Pressure at ``z_{atm}`` | ``P_{atm}`` | Pa |
| Temperature at ``z_{atm}`` | ``T_{atm}`` | K |
| Incident longwave radiation | ``L_{atm}\downarrow`` | W/m² |
| Liquid precipitation | ``q_{rain}`` | kg/(m²·s) |
| Solid precipitation | ``q_{sno}`` | kg/(m²·s) |
| Direct beam visible solar radiation | ``S_{atm}\downarrow^{\mu}_{vis}`` | W/m² |
| Direct beam near-infrared solar radiation | ``S_{atm}\downarrow^{\mu}_{nir}`` | W/m² |
| Diffuse visible solar radiation | ``S_{atm}\downarrow_{vis}`` | W/m² |
| Diffuse near-infrared solar radiation | ``S_{atm}\downarrow_{nir}`` | W/m² |

### Table 1.4: Physical Constants

Physical constants used by the CLMU, shared by all components of CCSM (Table 1.4, p. 25).

| Constant | Symbol | Value | Units |
|----------|--------|-------|-------|
| Boltzmann constant | ``\kappa`` | 1.38065 × 10⁻²³ | J/K |
| Avogadro's number | ``N_A`` | 6.02214 × 10²⁶ | molecule/kmol |
| Universal gas constant | ``R_{gas}`` | ``N_A \kappa`` | J/(K·kmol) |
| Molecular weight of dry air | ``MW_{da}`` | 28.966 | kg/kmol |
| Dry air gas constant | ``R_{da}`` | ``R_{gas}/MW_{da}`` | J/(K·kg) |
| Molecular weight of water vapor | ``MW_{wv}`` | 18.016 | kg/kmol |
| Water vapor gas constant | ``R_{wv}`` | ``R_{gas}/MW_{wv}`` | J/(K·kg) |

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
