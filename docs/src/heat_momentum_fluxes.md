# Heat and Momentum Fluxes

## Overview

Chapter 3 of the CLMU technical note describes the computation of heat and momentum
fluxes for urban surfaces using Monin-Obukhov similarity theory. These fluxes
drive the surface energy balance and determine the exchange of energy and momentum
between the urban surface and the atmospheric boundary layer.

The `HeatMomentumFluxes` component computes:
1. **Roughness length and displacement height** from canyon geometry using the
   MacDonald et al. (1998) formulation (Eqs. 3.55-3.58)
2. **Stability correction functions** ``\psi_m`` and ``\psi_h`` for four regimes:
   very unstable (``\zeta < -1.574``), unstable (``-1.574 \leq \zeta < 0``),
   stable (``0 \leq \zeta \leq 1``), and very stable (``\zeta > 1``)
   (Eqs. 3.31-3.46)
3. **Convective velocity** ``U_c`` for unstable conditions using the convective
   velocity scale ``w_*`` (Eqs. 3.25, 3.29-3.30)
4. **Friction velocity** ``u_*`` and aerodynamic resistances ``r_{am}``,
   ``r_{ah}``, ``r_{aw}`` (Eqs. 3.26, 3.65-3.67)
5. **Canyon wind speed** for three flow regimes: isolated (``H/W < 0.5``),
   wake interference (``0.5 \leq H/W < 1``), and skimming (``H/W \geq 1``)
   (Eqs. 3.59-3.62)
6. **Surface resistance** for heat and moisture exchange (Eq. 3.68)
7. **UCL air temperature and humidity** as conductance-weighted averages
   (Eqs. 3.75, 3.93)
8. **Sensible heat fluxes** for all five surfaces (roof, pervious road,
   impervious road, sunlit wall, shaded wall) and the area-weighted total
   (Eqs. 3.69-3.74)
9. **Water vapor fluxes** for roof, pervious road, and impervious road
   (walls have zero evaporation) (Eqs. 3.76-3.81)
10. **Momentum fluxes** ``\tau_x`` and ``\tau_y`` (Eqs. 3.6-3.7)
11. **Partial derivatives** ``\partial H / \partial T_g`` for all five surfaces,
    needed for the implicit soil temperature calculation (Eqs. 3.95-3.99)

The Monin-Obukhov stability parameter ``\zeta`` is taken as an input parameter,
since in the original CLM implementation it is computed through an iterative
procedure (Section 3.2.3). This avoids an algebraic loop and allows the system
to be solved as a directed algebraic system.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp. Chapter 3: Heat and Momentum Fluxes (pp. 61-89).

```@docs
HeatMomentumFluxes
```

## Implementation

The `HeatMomentumFluxes` component implements the heat and momentum flux
parameterization with 48 equations and 48 unknowns.

### State Variables

```@example heat_mom
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using UrbanCanopy

sys = HeatMomentumFluxes()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example heat_mom
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example heat_mom
eqs = equations(sys)
```

## Analysis

The following figures explore the behavior of the heat and momentum flux system
across different canyon geometries and atmospheric stability conditions.

```@example heat_mom
using OrdinaryDiffEqDefault
using Plots

compiled = mtkcompile(sys)

function solve_hmf(compiled, params)
    prob = ODEProblem(compiled, params, (0.0, 1.0))
    solve(prob)
end

function base_hmf_params(compiled; H_W=1.0, ζ_in=0.1)
    Dict(
        compiled.H_W => H_W,
        compiled.H_canyon => 10.0,
        compiled.W_roof => 0.25,
        compiled.f_prvrd => 0.5,
        compiled.H_w_est => 5.0,
        compiled.u_atm => 3.0,
        compiled.v_atm => 2.0,
        compiled.θ_atm => 300.0,
        compiled.q_atm => 0.01,
        compiled.ρ_atm => 1.2,
        compiled.z_atm_m => 30.0,
        compiled.z_atm_h => 30.0,
        compiled.z_atm_w => 30.0,
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
        compiled.ζ_in => ζ_in,
    )
end
nothing # hide
```

### Roughness Parameters vs Canyon H/W Ratio

Displacement height ``d`` and roughness length ``z_{0m}`` as functions of the
canyon height-to-width ratio, computed using the MacDonald et al. (1998) formulation
(Eqs. 3.55-3.57). The displacement height approaches the building height for
deep canyons, while the roughness length peaks at intermediate H/W and decreases
for very deep canyons.

```@example heat_mom
hw_range = 0.1:0.05:5.0
d_vals = Float64[]
z0m_vals = Float64[]

for hw in hw_range
    p = base_hmf_params(compiled; H_W=hw)
    sol = solve_hmf(compiled, p)
    push!(d_vals, sol[compiled.d_canopy][end])
    push!(z0m_vals, sol[compiled.z_0m_canopy][end])
end

p1 = plot(hw_range, d_vals, label="d (displacement height)",
    xlabel="Canyon H/W Ratio", ylabel="Height (m)",
    title="Roughness Parameters vs H/W (Eqs. 3.55-3.57)",
    linewidth=2, legend=:topleft)
plot!(p1, hw_range, z0m_vals, label="z₀ₘ (roughness length)", linewidth=2)
hline!(p1, [10.0], label="H_canyon = 10 m", linestyle=:dash, color=:gray)
p1
```

### Stability Corrections and Friction Velocity

The stability correction functions ``\psi_m`` and ``\psi_h`` modify the
log-law wind and temperature profiles. Positive ``\psi`` (unstable conditions)
enhances turbulent transfer; negative ``\psi`` (stable conditions) suppresses it.
The friction velocity responds accordingly, being higher in unstable conditions
when turbulent mixing is enhanced.

```@example heat_mom
ζ_range = -3.0:0.05:3.0
ψ_m_vals = Float64[]
ψ_h_vals = Float64[]
ustar_vals = Float64[]

for ζ in ζ_range
    if abs(ζ) < 0.001
        ζ = sign(ζ) == -1 ? -0.001 : 0.001
    end
    p = base_hmf_params(compiled; ζ_in=ζ)
    sol = solve_hmf(compiled, p)
    push!(ψ_m_vals, sol[compiled.ψ_m_atm][end])
    push!(ψ_h_vals, sol[compiled.ψ_h_atm][end])
    push!(ustar_vals, sol[compiled.u_star][end])
end

p1 = plot(ζ_range, ψ_m_vals, label="ψ_m", linewidth=2,
    xlabel="ζ = (z-d)/L", ylabel="ψ",
    title="Stability Corrections (Eqs. 3.31-3.46)",
    legend=:topright)
plot!(p1, ζ_range, ψ_h_vals, label="ψ_h", linewidth=2)
vline!(p1, [0.0], linestyle=:dash, color=:gray, label="Neutral")
vline!(p1, [-1.574], linestyle=:dot, color=:red, label="ζ_m transition")
vline!(p1, [-0.465], linestyle=:dot, color=:blue, label="ζ_h transition")

p2 = plot(ζ_range, ustar_vals, label="u*", linewidth=2, color=:green,
    xlabel="ζ = (z-d)/L", ylabel="u* (m/s)",
    title="Friction Velocity vs Stability",
    legend=:topright)
vline!(p2, [0.0], linestyle=:dash, color=:gray, label="Neutral")

p = plot(p1, p2, layout=(2, 1), size=(700, 700))
p
```

### Canyon Wind Speed and Flow Regimes

The canyon wind speed exhibits three flow regimes depending on the canyon
height-to-width ratio: isolated roughness flow for ``H/W < 0.5``, wake
interference for ``0.5 \leq H/W < 1``, and skimming flow for ``H/W \geq 1``
(Eqs. 3.60-3.62). Deeper canyons trap wind less effectively, reducing the
canyon wind speed.

```@example heat_mom
hw_range_wind = 0.1:0.05:5.0
Ucan_vals = Float64[]
Uac_vals = Float64[]
Vr_vals = Float64[]

for hw in hw_range_wind
    p = base_hmf_params(compiled; H_W=hw)
    sol = solve_hmf(compiled, p)
    push!(Ucan_vals, sol[compiled.U_can][end])
    push!(Uac_vals, sol[compiled.U_ac][end])
    push!(Vr_vals, sol[compiled.V_r][end])
end

p = plot(hw_range_wind, Ucan_vals, label="U_can (horizontal)", linewidth=2,
    xlabel="Canyon H/W Ratio", ylabel="Wind Speed (m/s)",
    title="Canyon Wind Speed vs H/W (Eqs. 3.59-3.62)",
    legend=:topright)
plot!(p, hw_range_wind, Uac_vals, label="U_ac (total canyon)", linewidth=2)
hline!(p, [sqrt(3.0^2 + 2.0^2)], label="V_r (reference wind)", linestyle=:dash, color=:gray)
vline!(p, [0.5], linestyle=:dot, color=:red, label="Isolated/Wake boundary")
vline!(p, [1.0], linestyle=:dot, color=:blue, label="Wake/Skimming boundary")
p
```

### Sensible Heat Fluxes vs Canyon H/W

Sensible heat fluxes from each urban surface as a function of canyon
height-to-width ratio (Eqs. 3.69-3.74). Positive values indicate upward
heat flux (surface warming the atmosphere). The UCL temperature is a weighted
average of all surface temperatures and the atmospheric temperature.

```@example heat_mom
hw_range_H = 0.1:0.1:5.0
H_roof_vals = Float64[]
H_prvrd_vals = Float64[]
H_imprvrd_vals = Float64[]
H_sunwall_vals = Float64[]
H_shdwall_vals = Float64[]
H_total_vals = Float64[]
T_ac_vals = Float64[]

for hw in hw_range_H
    p = base_hmf_params(compiled; H_W=hw)
    sol = solve_hmf(compiled, p)
    push!(H_roof_vals, sol[compiled.H_roof][end])
    push!(H_prvrd_vals, sol[compiled.H_prvrd][end])
    push!(H_imprvrd_vals, sol[compiled.H_imprvrd][end])
    push!(H_sunwall_vals, sol[compiled.H_sunwall][end])
    push!(H_shdwall_vals, sol[compiled.H_shdwall][end])
    push!(H_total_vals, sol[compiled.H_total][end])
    push!(T_ac_vals, sol[compiled.T_ac][end])
end

p1 = plot(hw_range_H, H_roof_vals, label="H_roof", linewidth=2,
    xlabel="Canyon H/W Ratio", ylabel="Sensible Heat Flux (W m⁻²)",
    title="Surface Sensible Heat Fluxes (Eqs. 3.69-3.74)",
    legend=:right)
plot!(p1, hw_range_H, H_prvrd_vals, label="H_prvrd", linewidth=2)
plot!(p1, hw_range_H, H_imprvrd_vals, label="H_imprvrd", linewidth=2)
plot!(p1, hw_range_H, H_sunwall_vals, label="H_sunwall", linewidth=2)
plot!(p1, hw_range_H, H_shdwall_vals, label="H_shdwall", linewidth=2)
plot!(p1, hw_range_H, H_total_vals, label="H_total", linewidth=2, linestyle=:dash, color=:black)

p2 = plot(hw_range_H, T_ac_vals, label="T_ac", linewidth=2, color=:red,
    xlabel="Canyon H/W Ratio", ylabel="Temperature (K)",
    title="UCL Air Temperature (Eq. 3.75)",
    legend=:topright)
hline!(p2, [300.0], label="θ_atm", linestyle=:dash, color=:gray)

p = plot(p1, p2, layout=(2, 1), size=(700, 700))
p
```

The UCL temperature shifts towards the surface temperatures as the canyon
deepens and surface-to-air conductances dominate over the atmospheric conductance.
Individual surface heat fluxes adjust as the UCL temperature changes with H/W.

### Aerodynamic and Surface Resistances vs Stability

Aerodynamic resistance decreases under unstable conditions (enhanced turbulent
mixing) and increases under stable conditions (suppressed mixing). The surface
resistance depends on the canyon wind speed through forced convection (Eq. 3.68).

```@example heat_mom
ζ_range_r = -2.0:0.1:2.0
r_am_vals = Float64[]
r_ah_vals = Float64[]
r_s_vals = Float64[]

for ζ in ζ_range_r
    if abs(ζ) < 0.001
        ζ = sign(ζ) == -1 ? -0.001 : 0.001
    end
    p = base_hmf_params(compiled; ζ_in=ζ)
    sol = solve_hmf(compiled, p)
    push!(r_am_vals, sol[compiled.r_am][end])
    push!(r_ah_vals, sol[compiled.r_ah][end])
    push!(r_s_vals, sol[compiled.r_s_u][end])
end

p = plot(ζ_range_r, r_am_vals, label="r_am (momentum)", linewidth=2,
    xlabel="ζ = (z-d)/L", ylabel="Resistance (s/m)",
    title="Resistances vs Stability (Eqs. 3.65-3.68)",
    legend=:topleft, yscale=:log10)
plot!(p, ζ_range_r, r_ah_vals, label="r_ah (heat)", linewidth=2)
plot!(p, ζ_range_r, r_s_vals, label="r_s (surface)", linewidth=2)
vline!(p, [0.0], linestyle=:dash, color=:gray, label="Neutral")
p
```

The aerodynamic resistances ``r_{am}`` and ``r_{ah}`` vary strongly with
stability, while the surface resistance ``r_s`` shows a weaker dependence
because it is controlled primarily by the canyon wind speed rather than
directly by the stability corrections.

### Convective Velocity vs Stability

The convective velocity ``U_c`` (Eqs. 3.29-3.30) augments the effective wind
speed ``V_a`` under unstable conditions. ``U_c = \beta w_*`` where ``w_*`` is
the convective velocity scale computed from the friction velocity and the
Monin-Obukhov length. Under stable conditions ``U_c = 0`` and ``V_a`` reduces
to the horizontal wind speed.

```@example heat_mom
ζ_range_uc = -3.0:0.05:1.0
Uc_vals = Float64[]
Va_vals = Float64[]
Vr_vals_stab = Float64[]

for ζ in ζ_range_uc
    if abs(ζ) < 0.001
        ζ = sign(ζ) == -1 ? -0.001 : 0.001
    end
    p = base_hmf_params(compiled; ζ_in=ζ)
    sol = solve_hmf(compiled, p)
    push!(Uc_vals, sol[compiled.U_c][end])
    push!(Va_vals, sol[compiled.V_a][end])
    push!(Vr_vals_stab, sol[compiled.V_r][end])
end

p1 = plot(ζ_range_uc, Uc_vals, label="U_c (convective velocity)", linewidth=2,
    xlabel="ζ = (z-d)/L", ylabel="Speed (m/s)",
    title="Convective Velocity and Effective Wind Speed (Eqs. 3.25, 3.29-3.30)",
    legend=:topright)
plot!(p1, ζ_range_uc, Va_vals, label="V_a (effective wind)", linewidth=2)
plot!(p1, ζ_range_uc, Vr_vals_stab, label="V_r (horizontal wind)", linewidth=2, linestyle=:dash)
vline!(p1, [0.0], linestyle=:dash, color=:gray, label="Neutral")
p1
```

Under strongly unstable conditions, the convective velocity can exceed the
mean horizontal wind, significantly increasing the effective wind speed and
turbulent exchange.

### Partial Derivatives of Sensible Heat Flux

The partial derivatives ``\partial H / \partial T_g`` (Eqs. 3.95-3.99) quantify
the sensitivity of each surface's sensible heat flux to changes in its ground
temperature. These are needed for the implicit coupling with the soil temperature
equation (Chapter 4). Larger values indicate stronger thermal coupling between
the surface and the atmosphere.

```@example heat_mom
hw_range_dH = 0.1:0.1:5.0
dH_roof_vals = Float64[]
dH_prvrd_vals = Float64[]
dH_imprvrd_vals = Float64[]
dH_sunwall_vals = Float64[]
dH_shdwall_vals = Float64[]

for hw in hw_range_dH
    p = base_hmf_params(compiled; H_W=hw)
    sol = solve_hmf(compiled, p)
    push!(dH_roof_vals, sol[compiled.dH_roof_dT][end])
    push!(dH_prvrd_vals, sol[compiled.dH_prvrd_dT][end])
    push!(dH_imprvrd_vals, sol[compiled.dH_imprvrd_dT][end])
    push!(dH_sunwall_vals, sol[compiled.dH_sunwall_dT][end])
    push!(dH_shdwall_vals, sol[compiled.dH_shdwall_dT][end])
end

p = plot(hw_range_dH, dH_roof_vals, label="∂H_roof/∂T", linewidth=2,
    xlabel="Canyon H/W Ratio", ylabel="∂H/∂T (W m⁻² K⁻¹)",
    title="Sensible Heat Flux Derivatives (Eqs. 3.95-3.99)",
    legend=:right)
plot!(p, hw_range_dH, dH_prvrd_vals, label="∂H_prvrd/∂T", linewidth=2)
plot!(p, hw_range_dH, dH_imprvrd_vals, label="∂H_imprvrd/∂T", linewidth=2)
plot!(p, hw_range_dH, dH_sunwall_vals, label="∂H_sunwall/∂T", linewidth=2)
plot!(p, hw_range_dH, dH_shdwall_vals, label="∂H_shdwall/∂T", linewidth=2)
p
```

The partial derivatives reflect the thermal coupling strength of each surface.
Surfaces with larger area fractions (e.g., roof) have stronger coupling, while
the wall derivatives increase with H/W as wall area grows relative to the
canyon floor.
