# Albedos and Radiative Fluxes

## Overview

Chapter 2 of the CLMU technical note describes the computation of albedos and
radiative fluxes for urban surfaces. The effects of canyon geometry on the
radiation balance are a key driver of urban-rural energy balance differences
(Oke et al. 1991). Shadowing of urban surfaces affects incident radiation and
thus temperature, and multiple reflections of radiation between canyon surfaces
must be accounted for (Harman et al. 2004).

The `UrbanRadiation` component computes:
1. **View factors** between canyon surfaces (road, walls, sky) as functions of the
   height-to-width ratio ``H/W`` (Eqs. 2.22-2.28)
2. **Incident direct beam solar radiation** on walls and road, integrated over
   canyon orientations (Eqs. 2.11-2.17)
3. **Incident diffuse solar radiation** on canyon surfaces (Eqs. 2.29-2.32)
4. **Absorbed and reflected solar radiation** with closed-form multi-reflection
   solution (Eqs. 2.34-2.95)
5. **Longwave radiation** including emission, reflection, and multi-reflection
   within the canyon (Eqs. 2.96-2.175)

The iterative multi-reflection scheme from the paper (Eqs. 2.63-2.85 for solar,
2.138-2.170 for longwave) is solved in closed form using the matrix geometric
series: ``S_{absorbed} = (I - A)(I - TA)^{-1} S_0``, where ``A`` is the albedo
matrix, ``T`` is the view-factor transport matrix, and ``S_0`` is the initial
incident radiation vector. The 3x3 matrix inverse is computed analytically.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp. Chapter 2: Albedos and Radiative Fluxes (pp. 26-60).

```@docs
UrbanRadiation
```

## Implementation

The `UrbanRadiation` component implements the radiation balance for the urban
canyon system with 62 equations, 62 variables, and 31 parameters.

### State Variables

```@example albedo_rad
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using UrbanCanopy

sys = UrbanRadiation()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example albedo_rad
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example albedo_rad
eqs = equations(sys)
```

## Analysis

The following figures reproduce the key results from Chapter 2 of the CLMU
technical note using the implementation above.

```@example albedo_rad
using OrdinaryDiffEqDefault
using Plots

compiled = mtkcompile(sys)

function solve_urban(compiled, params)
    prob = ODEProblem(compiled, params, (0.0, 1.0))
    solve(prob)
end

function base_params(compiled; H_W=1.0, μ_zen=0.5236)
    Dict(
        compiled.H_W => H_W,
        compiled.W_roof => 0.25,
        compiled.S_atm_dir_vis => 200.0,
        compiled.S_atm_dir_nir => 200.0,
        compiled.S_atm_dif_vis => 100.0,
        compiled.S_atm_dif_nir => 100.0,
        compiled.L_atm_down => 340.0,
        compiled.μ_zen => μ_zen,
        compiled.α_roof_dir_vis => 0.15,
        compiled.α_roof_dir_nir => 0.15,
        compiled.α_roof_dif_vis => 0.15,
        compiled.α_roof_dif_nir => 0.15,
        compiled.α_road_dir_vis => 0.08,
        compiled.α_road_dir_nir => 0.08,
        compiled.α_road_dif_vis => 0.08,
        compiled.α_road_dif_nir => 0.08,
        compiled.α_wall_dir_vis => 0.25,
        compiled.α_wall_dir_nir => 0.25,
        compiled.α_wall_dif_vis => 0.25,
        compiled.α_wall_dif_nir => 0.25,
        compiled.ε_roof => 0.90,
        compiled.ε_road => 0.95,
        compiled.ε_wall => 0.85,
        compiled.T_roof => 292.16,
        compiled.T_road => 292.16,
        compiled.T_sunwall => 292.16,
        compiled.T_shdwall => 292.16,
    )
end
nothing # hide
```

### Figure 2.4: View Factors vs Canyon H/W Ratio

View factors as a function of canyon height to width ratio (cf. Figure 2.4, p. 37).
``\Psi_{road-sky}`` is the fraction of radiation reaching the sky from the road,
``\Psi_{road-wall}`` is the fraction reaching the wall from the road,
``\Psi_{wall-sky}`` and ``\Psi_{wall-road}`` are the fractions reaching the sky
and road from the wall, and ``\Psi_{wall-wall}`` is the fraction reaching the
opposite wall.

```@example albedo_rad
hw_range = 10 .^ range(-2, 1, length=200)

ψ_road_sky = Float64[]
ψ_road_wall = Float64[]
ψ_wall_sky = Float64[]
ψ_wall_road = Float64[]
ψ_wall_wall = Float64[]

for hw in hw_range
    p = base_params(compiled; H_W=hw)
    sol = solve_urban(compiled, p)
    push!(ψ_road_sky, sol[compiled.Ψ_road_sky][end])
    push!(ψ_road_wall, sol[compiled.Ψ_road_wall][end])
    push!(ψ_wall_sky, sol[compiled.Ψ_wall_sky][end])
    push!(ψ_wall_road, sol[compiled.Ψ_wall_road][end])
    push!(ψ_wall_wall, sol[compiled.Ψ_wall_wall][end])
end

p = plot(hw_range, ψ_road_sky, label="Ψ_road-sky", xscale=:log10,
    xlabel="Canyon H/W Ratio", ylabel="View Factor",
    title="View Factors vs Canyon H/W Ratio (Figure 2.4)",
    linewidth=2, legend=:right, ylims=(0, 1))
plot!(p, hw_range, ψ_road_wall, label="Ψ_road-wall", linewidth=2)
plot!(p, hw_range, ψ_wall_sky, label="Ψ_wall-sky, Ψ_wall-road", linewidth=2)
plot!(p, hw_range, ψ_wall_wall, label="Ψ_wall-wall", linewidth=2)
p
```

At low H/W ratios the road-sky view factor is close to one (flat surface behavior),
while at high H/W ratios most radiative interactions occur between the walls and
between the walls and road.

### Figure 2.5: Absorbed Solar Radiation vs H/W Ratio

Solar radiation absorbed by urban surfaces for solar zenith angles of 30 degrees
(top) and 60 degrees (bottom) (cf. Figure 2.5, p. 47). The atmospheric solar
radiation is ``S_{atm}\downarrow^{\mu} = 400`` and ``S_{atm}\downarrow = 200``
W m``^{-2}``. Albedos: roof = 0.15, road = 0.08, wall = 0.25. Note that wall
fluxes are per unit wall area; the canyon total converts them to per unit ground
area using the H/W ratio.

```@example albedo_rad
hw_range_lin = 0.1:0.1:10.0

function compute_absorbed_solar(compiled, hw_range, μ_zen)
    road = Float64[]
    sunwall = Float64[]
    shdwall = Float64[]
    roof = Float64[]
    canyon = Float64[]

    for hw in hw_range
        p = base_params(compiled; H_W=hw, μ_zen=μ_zen)
        p[compiled.S_atm_dir_vis] = 200.0
        p[compiled.S_atm_dir_nir] = 200.0
        p[compiled.S_atm_dif_vis] = 100.0
        p[compiled.S_atm_dif_nir] = 100.0
        sol = solve_urban(compiled, p)
        push!(road, sol[compiled.S_net_road][end])
        push!(sunwall, sol[compiled.S_net_sunwall][end])
        push!(shdwall, sol[compiled.S_net_shdwall][end])
        push!(roof, sol[compiled.S_net_roof][end])
        push!(canyon, sol[compiled.S_net_uc][end])
    end
    return road, sunwall, shdwall, roof, canyon
end

road30, sw30, sh30, roof30, can30 = compute_absorbed_solar(compiled, hw_range_lin, deg2rad(30))
road60, sw60, sh60, roof60, can60 = compute_absorbed_solar(compiled, hw_range_lin, deg2rad(60))

p1 = plot(hw_range_lin, can30, label="Canyon", linewidth=2,
    xlabel="Canyon H/W Ratio",
    ylabel="Absorbed Solar Radiation (W m⁻²)",
    title="Solar Zenith Angle = 30°",
    legend=:right, ylims=(0, 700))
plot!(p1, hw_range_lin, roof30, label="Roof", linewidth=2)
plot!(p1, hw_range_lin, road30, label="Road", linewidth=2)
plot!(p1, hw_range_lin, sw30, label="Sunlit Wall", linewidth=2)
plot!(p1, hw_range_lin, sh30, label="Shaded Wall", linewidth=2)

p2 = plot(hw_range_lin, can60, label="Canyon", linewidth=2,
    xlabel="Canyon H/W Ratio",
    ylabel="Absorbed Solar Radiation (W m⁻²)",
    title="Solar Zenith Angle = 60°",
    legend=:right, ylims=(0, 700))
plot!(p2, hw_range_lin, roof60, label="Roof", linewidth=2)
plot!(p2, hw_range_lin, road60, label="Road", linewidth=2)
plot!(p2, hw_range_lin, sw60, label="Sunlit Wall", linewidth=2)
plot!(p2, hw_range_lin, sh60, label="Shaded Wall", linewidth=2)

p = plot(p1, p2, layout=(2, 1), size=(700, 800))
p
```

The absorbed solar radiation for the roof is independent of H/W and solar zenith
angle. At both zenith angles, absorbed road radiation decreases rapidly with
increasing H/W as buildings shade the road. The canyon total increases slowly
with H/W because radiation trapping in deeper canyons increases total absorption.

### Figure 2.6: Canyon Albedo vs H/W Ratio

Direct beam and diffuse canyon albedo (excluding roof) as a function of H/W ratio
for solar zenith angles from 0 to 85 degrees (cf. Figure 2.6, p. 48). Albedos:
road = 0.08, wall = 0.25.

```@example albedo_rad
hw_range_alb = 0.05:0.05:3.0
zenith_angles = 0:5:85

p = plot(xlabel="Canyon H/W Ratio", ylabel="Canyon Albedo",
    title="Canyon Albedo vs H/W Ratio (Figure 2.6)",
    legend=false, ylims=(0, 0.14))

# Direct beam albedo for each zenith angle (solid lines)
for (i, zen_deg) in enumerate(zenith_angles)
    zen_rad = deg2rad(zen_deg)
    # Avoid exactly 0 or 90 degrees
    if zen_rad < 0.001
        zen_rad = 0.001
    end
    alb_dir = Float64[]
    for hw in hw_range_alb
        params = base_params(compiled; H_W=hw, μ_zen=zen_rad)
        sol = solve_urban(compiled, params)
        # Average VIS and NIR direct beam albedo
        a_vis = sol[compiled.α_uc_dir_vis][end]
        a_nir = sol[compiled.α_uc_dir_nir][end]
        push!(alb_dir, (a_vis + a_nir) / 2)
    end
    plot!(p, hw_range_alb, alb_dir, linewidth=1.5, color=:black, linestyle=:solid)
end

# Diffuse albedo (dashed line) - independent of zenith angle
alb_dif = Float64[]
for hw in hw_range_alb
    params = base_params(compiled; H_W=hw)
    sol = solve_urban(compiled, params)
    a_vis = sol[compiled.α_uc_dif_vis][end]
    a_nir = sol[compiled.α_uc_dif_nir][end]
    push!(alb_dif, (a_vis + a_nir) / 2)
end
plot!(p, hw_range_alb, alb_dif, linewidth=2, color=:black, linestyle=:dash)

# Add annotations
annotate!(p, 0.6, 0.12, text("Horizon", 8))
annotate!(p, 0.4, 0.035, text("Zenith", 8))
annotate!(p, 2.5, 0.09, text("Diffuse", 8))
annotate!(p, 2.5, 0.04, text("Direct Beam", 8))
p
```

The canyon albedo generally decreases with H/W as more radiation is trapped and
absorbed within deeper canyons. The trapping of solar radiation is less effective
at larger solar zenith angles (near horizon), where the albedo can increase at
small H/W because the higher-albedo walls dominate the radiative exchange.

### Figure 2.7: Net Longwave Radiation vs H/W Ratio

Net longwave radiation (positive toward the atmosphere) for urban surfaces with
two different emissivity configurations (cf. Figure 2.7, p. 59). The atmospheric
longwave radiation is ``L_{atm}\downarrow = 340`` W m``^{-2}`` and the
temperature of each surface is 292.16 K. Wall fluxes are per unit wall area.

```@example albedo_rad
hw_range_lw = 0.1:0.1:10.0

function compute_net_lw(compiled, hw_range; ε_roof, ε_road, ε_wall)
    road_lw = Float64[]
    wall_lw = Float64[]
    roof_lw = Float64[]
    canyon_lw = Float64[]

    for hw in hw_range
        p = base_params(compiled; H_W=hw)
        p[compiled.ε_roof] = ε_roof
        p[compiled.ε_road] = ε_road
        p[compiled.ε_wall] = ε_wall
        sol = solve_urban(compiled, p)
        push!(road_lw, sol[compiled.L_net_road][end])
        push!(wall_lw, sol[compiled.L_net_sunwall][end])  # walls equal when T equal
        push!(roof_lw, sol[compiled.L_net_roof][end])
        push!(canyon_lw, sol[compiled.L_net_uc][end])
    end
    return road_lw, wall_lw, roof_lw, canyon_lw
end

road1, wall1, roof1, can1 = compute_net_lw(compiled, hw_range_lw;
    ε_roof=0.95, ε_road=0.95, ε_wall=0.95)
road2, wall2, roof2, can2 = compute_net_lw(compiled, hw_range_lw;
    ε_roof=0.90, ε_road=0.94, ε_wall=0.85)

p1 = plot(hw_range_lw, can1, label="Canyon", linewidth=2,
    xlabel="Canyon H/W Ratio",
    ylabel="Net Longwave Radiation (W m⁻²)",
    title="Emissivity: Roof=0.95, Road=0.95, Wall=0.95",
    legend=:right, ylims=(0, 80))
plot!(p1, hw_range_lw, roof1, label="Roof", linewidth=2)
plot!(p1, hw_range_lw, road1, label="Road", linewidth=2)
plot!(p1, hw_range_lw, wall1, label="Walls", linewidth=2)

p2 = plot(hw_range_lw, can2, label="Canyon", linewidth=2,
    xlabel="Canyon H/W Ratio",
    ylabel="Net Longwave Radiation (W m⁻²)",
    title="Emissivity: Roof=0.90, Road=0.94, Wall=0.85",
    legend=:right, ylims=(0, 80))
plot!(p2, hw_range_lw, roof2, label="Roof", linewidth=2)
plot!(p2, hw_range_lw, road2, label="Road", linewidth=2)
plot!(p2, hw_range_lw, wall2, label="Walls", linewidth=2)

p = plot(p1, p2, layout=(2, 1), size=(700, 800))
p
```

The net longwave radiation for the roof is independent of H/W and increases with
higher emissivity. The net longwave for road and walls decreases rapidly with
increasing H/W as more radiation is trapped within the canyon. The walls have
lower net longwave than the road because their sky view factors are smaller. The
canyon net longwave (sum of road and wall fluxes converted to per unit ground
area) increases slowly with H/W because of the larger surface area of the walls.
