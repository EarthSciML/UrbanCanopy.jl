# Offline Mode: Atmospheric Forcing

## Overview

In offline mode (uncoupled from an atmospheric model), the Community Land Model Urban
(CLMU) requires atmospheric forcing data from observed datasets. Chapter 6 of the
technical note describes how raw forcing data is processed into the specific variables
needed by the urban model.

The `OfflineModeForcing` component implements the algebraic relationships for:
- Partitioning total solar radiation into direct/diffuse and visible/near-infrared
  components (Eqs. 6.2-6.8) using polynomial fits derived from CAM output
- Computing atmospheric downwelling longwave radiation from vapor pressure and
  temperature using the Idso (1981) formula (Eqs. 6.9-6.10)
- Partitioning precipitation into rain and snow based on temperature (Eqs. 6.11-6.13)
- Deriving wind components, potential temperature, and reference height (p. 149)

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp.

```@docs
OfflineModeForcing
```

## Implementation

### State Variables

```@example offline_mode
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using UrbanCanopy

sys = OfflineModeForcing()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example offline_mode
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example offline_mode
eqs = equations(sys)
```

## Analysis

### Solar Radiation Partitioning (Eqs. 6.2-6.8)

The total incident solar radiation is split into four components: direct beam visible,
direct beam near-infrared, diffuse visible, and diffuse near-infrared. The visible
fraction ``\alpha = 0.5`` (Eq. 6.6), and the direct-to-total ratios ``R_{vis}`` and
``R_{nir}`` are polynomial functions of the solar flux (Eqs. 6.7-6.8), clamped to
[0.01, 0.99].

```@example offline_mode
using OrdinaryDiffEqDefault
using Plots

compiled = mtkcompile(sys)

S_range = 0.0:10.0:1200.0
dir_vis = Float64[]
dir_nir = Float64[]
dif_vis = Float64[]
dif_nir = Float64[]

for S_val in S_range
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => S_val,
            compiled.T_atm => 293.15,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.01,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    push!(dir_vis, sol[compiled.S_atm_dir_vis][end])
    push!(dir_nir, sol[compiled.S_atm_dir_nir][end])
    push!(dif_vis, sol[compiled.S_atm_dif_vis][end])
    push!(dif_nir, sol[compiled.S_atm_dif_nir][end])
end

p = plot(S_range, dir_vis, label="Direct visible", xlabel="Total solar radiation (W/m²)",
    ylabel="Component radiation (W/m²)",
    title="Solar Radiation Partitioning (Eqs. 6.2-6.5)",
    linewidth=2, legend=:topleft)
plot!(p, S_range, dir_nir, label="Direct NIR", linewidth=2)
plot!(p, S_range, dif_vis, label="Diffuse visible", linewidth=2, linestyle=:dash)
plot!(p, S_range, dif_nir, label="Diffuse NIR", linewidth=2, linestyle=:dash)
p
```

At low solar radiation (overcast), the diffuse fraction dominates. As total solar
radiation increases (clearer skies), the direct beam fraction increases, consistent
with less atmospheric scattering under clear conditions.

### Direct Fraction vs. Solar Flux (Eqs. 6.7-6.8)

The polynomial fits for ``R_{vis}`` and ``R_{nir}`` show how the direct-to-total
radiation ratio varies with solar intensity.

```@example offline_mode
R_vis_vals = Float64[]
R_nir_vals = Float64[]

for S_val in S_range
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => S_val,
            compiled.T_atm => 293.15,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.01,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    push!(R_vis_vals, sol[compiled.R_vis][end])
    push!(R_nir_vals, sol[compiled.R_nir][end])
end

p = plot(S_range, R_vis_vals, label="R_vis (visible)", xlabel="Total solar radiation (W/m²)",
    ylabel="Direct fraction (dimensionless)",
    title="Direct-to-Total Radiation Ratio (Eqs. 6.7-6.8)",
    linewidth=2, legend=:right, ylims=(0, 1))
plot!(p, S_range, R_nir_vals, label="R_nir (NIR)", linewidth=2)
hline!(p, [0.01, 0.99], label="Clamp bounds", linestyle=:dot, color=:gray)
p
```

The near-infrared direct fraction is consistently higher than the visible direct
fraction, reflecting the greater atmospheric scattering of shorter (visible)
wavelengths relative to longer (NIR) wavelengths.

### Longwave Radiation vs. Temperature (Eq. 6.9)

The Idso (1981) formula computes atmospheric downwelling longwave radiation from
vapor pressure and temperature. The dominant dependence is on ``T_{atm}^4`` through
the Stefan-Boltzmann law, modulated by an emissivity factor that increases with
atmospheric moisture.

```@example offline_mode
T_range = 250.0:2.0:320.0
L_dry = Float64[]
L_moist = Float64[]

for T_val in T_range
    prob_d = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 0.0,
            compiled.T_atm => T_val,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.001,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol_d = solve(prob_d)
    push!(L_dry, sol_d[compiled.L_atm_down][end])

    prob_m = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 0.0,
            compiled.T_atm => T_val,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.02,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol_m = solve(prob_m)
    push!(L_moist, sol_m[compiled.L_atm_down][end])
end

p = plot(T_range, L_dry, label="Dry (q=0.001 kg/kg)", xlabel="Temperature (K)",
    ylabel="Longwave radiation (W/m²)",
    title="Atmospheric Downwelling Longwave (Eq. 6.9, Idso 1981)",
    linewidth=2, legend=:topleft)
plot!(p, T_range, L_moist, label="Moist (q=0.02 kg/kg)", linewidth=2, linestyle=:dash)
p
```

Higher atmospheric moisture content increases the effective emissivity of the
atmosphere, resulting in greater downwelling longwave radiation at all temperatures.

### Precipitation Phase Partitioning (Eqs. 6.11-6.13)

The rain fraction ``f_P`` varies linearly from 0 (all snow) at or below the freezing
point ``T_f = 273.15`` K to 1 (all rain) at ``T_f + 2`` K. This simple linear ramp
avoids an abrupt transition between rain and snow.

```@example offline_mode
T_precip_range = 268.0:0.1:280.0
f_P_vals = Float64[]
rain_vals = Float64[]
snow_vals = Float64[]
P_val = 1e-3  # 1 mm/s

for T_val in T_precip_range
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 0.0,
            compiled.T_atm => T_val,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.005,
            compiled.W_atm => 3.0,
            compiled.P_precip => P_val,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    push!(f_P_vals, sol[compiled.f_P][end])
    push!(rain_vals, sol[compiled.q_rain][end])
    push!(snow_vals, sol[compiled.q_snow][end])
end

p = plot(T_precip_range .- 273.15, f_P_vals,
    label="Rain fraction f_P", xlabel="Temperature (°C)",
    ylabel="Fraction / Rate (mm/s)",
    title="Precipitation Phase Partitioning (Eqs. 6.11-6.13)",
    linewidth=2, legend=:right)
plot!(p, T_precip_range .- 273.15, rain_vals .* 1000,
    label="Rain rate (mm/s)", linewidth=2, linestyle=:dash)
plot!(p, T_precip_range .- 273.15, snow_vals .* 1000,
    label="Snow rate (mm/s)", linewidth=2, linestyle=:dot)
vline!(p, [0.0, 2.0], label="Transition range", linestyle=:dot, color=:gray)
p
```

The transition zone spans 0°C to 2°C. Below 0°C, all precipitation falls as snow;
above 2°C, all precipitation falls as rain. Between these limits, a mixture of rain
and snow occurs with the rain fraction ``f_P = 0.5(T_{atm} - T_f)`` (Eq. 6.13).

### Energy Conservation Check

A key validation of the solar radiation partitioning is that the sum of all four
components equals the total input solar radiation for all values of ``S_{atm}``.

```@example offline_mode
total_check = Float64[]

for S_val in S_range
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => S_val,
            compiled.T_atm => 293.15,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.01,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    total = sol[compiled.S_atm_dir_vis][end] + sol[compiled.S_atm_dir_nir][end] +
            sol[compiled.S_atm_dif_vis][end] + sol[compiled.S_atm_dif_nir][end]
    push!(total_check, total - S_val)
end

p = plot(S_range, total_check, xlabel="Total solar radiation (W/m²)",
    ylabel="Residual (W/m²)",
    title="Energy Conservation: Sum of Components - Total",
    linewidth=2, legend=false)
p
```

The residual is effectively zero (machine precision) across the full range of solar
radiation values, confirming exact energy conservation in the partitioning scheme.
