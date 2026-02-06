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

### Table 1.2: Urban Model Output to Atmospheric Model

The urban model provides the following area-averaged fluxes to the atmospheric model (Table 1.2, p. 18).

| Variable | Symbol | Units |
|----------|--------|-------|
| Latent heat flux | ``\lambda E`` | W/m² |
| Sensible heat flux | ``H`` | W/m² |
| Water vapor flux | ``E`` | mm/s |
| Zonal momentum flux | ``\tau_x`` | kg/(m·s²) |
| Meridional momentum flux | ``\tau_y`` | kg/(m·s²) |
| Emitted longwave radiation | ``L\uparrow`` | W/m² |
| Direct beam visible albedo | ``I\uparrow^{\mu}_{vis}`` | - |
| Direct beam near-infrared albedo | ``I\uparrow^{\mu}_{nir}`` | - |
| Diffuse visible albedo | ``I\uparrow_{vis}`` | - |
| Diffuse near-infrared albedo | ``I\uparrow_{nir}`` | - |
| Absorbed solar radiation | ``\vec{S}`` | W/m² |
| Radiative temperature | ``T_{rad}`` | K |
| Temperature at 2 meter height | ``T_{2m}`` | K |
| Specific humidity at 2 meter height | ``q_{2m}`` | kg/kg |
| Snow water equivalent | ``W_{sno}`` | m |
| Aerodynamic resistance | ``r_{am}`` | s/m |
| Friction velocity | ``u_*`` | m/s |

Note: ``\lambda`` is either the latent heat of vaporization ``\lambda_{vap}`` or latent heat of sublimation ``\lambda_{sub}`` (Table 1.4) depending on the thermal state of surface water on the roof, pervious and impervious road. Dust flux and net ecosystem exchange are set to zero for urban areas.

### Table 1.3: Input Data Required for the Urban Model

Required input data for urban landunits are listed in Table 1.3 (p. 23). This data is provided by the surface dataset at the required spatial resolution.

| Data | Symbol | Units |
|------|--------|-------|
| Percent urban | - | % |
| Canyon height to width ratio | ``H/W`` | - |
| Roof fraction | ``W_{roof}`` | - |
| Pervious road fraction | ``f_{prvrd}`` | - |
| Emissivity of roof | ``\varepsilon_{roof}`` | - |
| Emissivity of impervious road | ``\varepsilon_{imprvrd}`` | - |
| Emissivity of pervious road | ``\varepsilon_{prvrd}`` | - |
| Emissivity of walls | ``\varepsilon_{wall}`` | - |
| Building height | ``H`` | m |
| Roof albedo (visible direct/diffuse, near-infrared direct/diffuse) | ``\alpha_{roof}`` | - |
| Wall albedo (visible direct/diffuse, near-infrared direct/diffuse) | ``\alpha_{walls}`` | - |
| Impervious road albedo (visible direct/diffuse, near-infrared direct/diffuse) | ``\alpha_{imprvrd}`` | - |
| Pervious road albedo (visible direct/diffuse, near-infrared direct/diffuse) | ``\alpha_{prvrd}`` | - |
| Roof thermal conductivity | ``\lambda_{roof,i}`` | W/(m·K) |
| Wall thermal conductivity | ``\lambda_{wall,i}`` | W/(m·K) |
| Impervious road thermal conductivity | ``\lambda_{imprvrd,i}`` | W/(m·K) |
| Pervious road thermal conductivity | ``\lambda_{prvrd,i}`` | W/(m·K) |
| Roof volumetric heat capacity | ``c_{roof,i}`` | J/(m³·K) |
| Wall volumetric heat capacity | ``c_{wall,i}`` | J/(m³·K) |
| Impervious road volumetric heat capacity | ``c_{imprvrd,i}`` | J/(m³·K) |
| Pervious road volumetric heat capacity | ``c_{prvrd,i}`` | J/(m³·K) |
| Maximum interior building temperature | ``T_{iB,max}`` | K |
| Minimum interior building temperature | ``T_{iB,min}`` | K |
| Height of wind source in canyon | ``H_w`` | m |
| Number of impervious road layers | ``N_{imprvrd}`` | - |
| Wall thickness | ``\Delta z_{wall}`` | m |
| Roof thickness | ``\Delta z_{roof}`` | m |
| Percent sand, percent clay of pervious road (soil) | %sand, %clay | % |
| Grid cell latitude and longitude | ``\phi, \theta`` | degrees |

### Table 1.4: Physical Constants

Physical constants used by the CLMU, shared by all components of CCSM (Table 1.4, p. 25).

| Constant | Symbol | Value | Units |
|----------|--------|-------|-------|
| Pi | ``\pi`` | 3.14159265358979323846 | - |
| Acceleration of gravity | ``g`` | 9.80616 | m/s² |
| Standard pressure | ``P_{std}`` | 101325 | Pa |
| Stefan-Boltzmann constant | ``\sigma`` | 5.67 × 10⁻⁸ | W/(m²·K⁴) |
| Boltzmann constant | ``\kappa`` | 1.38065 × 10⁻²³ | J/K |
| Avogadro's number | ``N_A`` | 6.02214 × 10²⁶ | molecule/kmol |
| Universal gas constant | ``R_{gas}`` | ``N_A \kappa`` | J/(K·kmol) |
| Molecular weight of dry air | ``MW_{da}`` | 28.966 | kg/kmol |
| Dry air gas constant | ``R_{da}`` | ``R_{gas}/MW_{da}`` | J/(K·kg) |
| Molecular weight of water vapor | ``MW_{wv}`` | 18.016 | kg/kmol |
| Water vapor gas constant | ``R_{wv}`` | ``R_{gas}/MW_{wv}`` | J/(K·kg) |
| Von Karman constant | ``k`` | 0.4 | - |
| Freezing temperature of fresh water | ``T_f`` | 273.15 | K |
| Density of liquid water | ``\rho_{liq}`` | 1000 | kg/m³ |
| Density of ice | ``\rho_{ice}`` | 917 | kg/m³ |
| Specific heat capacity of dry air | ``C_p`` | 1.00464 × 10³ | J/(kg·K) |
| Specific heat capacity of water | ``C_{liq}`` | 4.188 × 10³ | J/(kg·K) |
| Specific heat capacity of ice | ``C_{ice}`` | 2.11727 × 10³ | J/(kg·K) |
| Latent heat of vaporization | ``\lambda_{vap}`` | 2.501 × 10⁶ | J/kg |
| Latent heat of fusion | ``L_f`` | 3.337 × 10⁵ | J/kg |
| Latent heat of sublimation | ``\lambda_{sub}`` | ``\lambda_{vap} + L_f`` | J/kg |
| Thermal conductivity of water | ``\lambda_{liq}`` | 0.6 | W/(m·K) |
| Thermal conductivity of ice | ``\lambda_{ice}`` | 2.29 | W/(m·K) |
| Thermal conductivity of air | ``\lambda_{air}`` | 0.023 | W/(m·K) |
| Radius of the earth | ``R_e`` | 6.37122 × 10⁶ | m |

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
