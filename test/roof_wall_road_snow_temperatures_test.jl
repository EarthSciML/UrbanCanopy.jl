@testsnippet TempSetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy
end

# ========================================================================
# Snow Layer Geometry Tests
# ========================================================================

@testitem "SnowLayerGeometry - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowLayerGeometry()

    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowLayerGeometry - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowLayerGeometry()
    compiled = mtkcompile(sys)

    z_sno_val = 0.02  # single layer case (0.01 ≤ z_sno ≤ 0.03)

    prob = ODEProblem(compiled, [compiled.z_sno => z_sno_val], (0.0, 1.0))
    sol = solve(prob)

    # Layer 0 thickness = z_sno
    @test sol[compiled.Δz_snow_0][end] ≈ z_sno_val rtol = 1.0e-10

    # Node at midpoint of layer: z = 0 - 0.5 * Δz (negative = above surface)
    @test sol[compiled.z_node_snow_0][end] ≈ -0.5 * z_sno_val rtol = 1.0e-10

    # Top interface: z_{h,-1} = 0 - Δz_0
    @test sol[compiled.z_interface_snow_neg1][end] ≈ -z_sno_val rtol = 1.0e-10
end

# ========================================================================
# Soil Thermal Properties Tests
# ========================================================================

@testitem "SoilThermalProperties - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SoilThermalProperties()

    @test length(equations(sys)) == 9
    @test length(unknowns(sys)) == 9

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SoilThermalProperties - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SoilThermalProperties()
    compiled = mtkcompile(sys)

    # Typical sandy soil: 60% sand, 20% clay
    params = Dict(
        compiled.pct_sand => 60.0,
        compiled.pct_clay => 20.0,
        compiled.θ_sat => 0.4,
        compiled.θ_liq => 0.1,  # w_liq/(ρ_liq*Δz) = 50/(1000*0.5)
        compiled.T_i => 293.15,
        compiled.w_ice => 0.0,
        compiled.w_liq => 50.0,
        compiled.Δz => 0.5,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.79: λ_s = (8.80 * 60 + 2.92 * 20) / (60 + 20) = 7.33
    @test sol[compiled.λ_s][end] ≈ 7.33 rtol = 1.0e-6

    # Eq. 4.86: c_s
    expected_cs = (2.128 * 60 + 2.385 * 20) / (60 + 20) * 1.0e6
    @test sol[compiled.c_s][end] ≈ expected_cs rtol = 1.0e-6

    # Bulk density: ρ_d = 2700 * (1 - 0.4) = 1620
    @test sol[compiled.ρ_d][end] ≈ 1620.0 rtol = 1.0e-6

    # Eq. 4.80: λ_dry
    expected_λ_dry = (0.135 * 1620 + 64.7) / (2700 - 0.947 * 1620)
    @test sol[compiled.λ_dry][end] ≈ expected_λ_dry rtol = 1.0e-6

    # Eq. 4.82: S_r = (50/(1000*0.5) + 0/(917*0.5)) / 0.4 = 0.25
    @test sol[compiled.S_r][end] ≈ 0.25 rtol = 1.0e-6

    # Thermal conductivity: physical bounds
    @test sol[compiled.λ_soil][end] > 0
    @test sol[compiled.λ_soil][end] ≥ sol[compiled.λ_dry][end] - 1.0e-10

    # Heat capacity: positive
    @test sol[compiled.c_soil][end] > 0
end

@testitem "SoilThermalProperties - Dry Soil Limit" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SoilThermalProperties()
    compiled = mtkcompile(sys)

    params = Dict(
        compiled.pct_sand => 50.0,
        compiled.pct_clay => 30.0,
        compiled.θ_sat => 0.4,
        compiled.θ_liq => 0.0,
        compiled.T_i => 293.15,
        compiled.w_ice => 0.0,
        compiled.w_liq => 0.0,
        compiled.Δz => 0.5,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # For S_r = 0, λ = λ_dry (Eq. 4.77)
    @test sol[compiled.S_r][end] ≈ 0.0 atol = 1.0e-15
    @test sol[compiled.λ_soil][end] ≈ sol[compiled.λ_dry][end] rtol = 1.0e-6

    # Dry soil heat capacity: c = c_s * (1 - θ_sat) only
    expected_c = sol[compiled.c_s][end] * (1 - 0.4)
    @test sol[compiled.c_soil][end] ≈ expected_c rtol = 1.0e-6
end

@testitem "SoilThermalProperties - Frozen Kersten" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SoilThermalProperties()
    compiled = mtkcompile(sys)

    params = Dict(
        compiled.pct_sand => 60.0,
        compiled.pct_clay => 20.0,
        compiled.θ_sat => 0.4,
        compiled.θ_liq => 0.0,  # Fully frozen, no liquid water
        compiled.T_i => 263.15,
        compiled.w_ice => 100.0,
        compiled.w_liq => 0.0,
        compiled.Δz => 0.5,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.81: For T < T_f, K_e = S_r
    S_r = (100.0 / (917.0 * 0.5)) / 0.4
    @test sol[compiled.K_e][end] ≈ S_r rtol = 1.0e-6
end

# ========================================================================
# Snow Thermal Properties Tests
# ========================================================================

@testitem "SnowThermalProperties - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowThermalProperties()

    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowThermalProperties - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowThermalProperties()
    compiled = mtkcompile(sys)

    Δz_val = 0.1
    w_ice_val = 10.0
    w_liq_val = 0.0

    params = Dict(
        compiled.w_ice => w_ice_val,
        compiled.w_liq => w_liq_val,
        compiled.Δz => Δz_val,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.84: ρ_sno = (w_ice + w_liq) / Δz = 100
    @test sol[compiled.ρ_sno][end] ≈ 100.0 rtol = 1.0e-10

    # Eq. 4.83: λ = λ_air + (7.75e-5 * ρ + 1.105e-6 * ρ²) * (λ_ice - λ_air)
    λ_air = 0.023
    λ_ice = 2.29
    ρ = 100.0
    expected_λ = λ_air + (7.75e-5 * ρ + 1.105e-6 * ρ^2) * (λ_ice - λ_air)
    @test sol[compiled.λ_snow][end] ≈ expected_λ rtol = 1.0e-6

    # Eq. 4.87: c = (w_ice/Δz) * C_ice
    C_ice = 2117.27
    expected_c = (w_ice_val / Δz_val) * C_ice
    @test sol[compiled.c_snow][end] ≈ expected_c rtol = 1.0e-6

    # Snow conductivity between air and ice values
    @test sol[compiled.λ_snow][end] > λ_air
    @test sol[compiled.λ_snow][end] < λ_ice
end

@testitem "SnowThermalProperties - Density Monotonicity" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowThermalProperties()
    compiled = mtkcompile(sys)

    params_low = Dict(compiled.w_ice => 5.0, compiled.w_liq => 0.0, compiled.Δz => 0.1)
    params_high = Dict(compiled.w_ice => 30.0, compiled.w_liq => 0.0, compiled.Δz => 0.1)

    prob_low = ODEProblem(compiled, params_low, (0.0, 1.0))
    prob_high = ODEProblem(compiled, params_high, (0.0, 1.0))
    sol_low = solve(prob_low)
    sol_high = solve(prob_high)

    @test sol_high[compiled.λ_snow][end] > sol_low[compiled.λ_snow][end]
    @test sol_high[compiled.c_snow][end] > sol_low[compiled.c_snow][end]
end

# ========================================================================
# Interface Thermal Conductivity Tests
# ========================================================================

@testitem "InterfaceThermalConductivity - Harmonic Mean" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = InterfaceThermalConductivity()
    compiled = mtkcompile(sys)

    # Equal conductivities: interface = same
    params = Dict(
        compiled.λ_i => 1.0, compiled.λ_ip1 => 1.0,
        compiled.z_i => 0.0, compiled.z_ip1 => 1.0, compiled.z_h => 0.5,
    )
    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)
    @test sol[compiled.λ_interface][end] ≈ 1.0 rtol = 1.0e-10

    # Different conductivities: harmonic mean at midpoint
    params2 = Dict(
        compiled.λ_i => 1.0, compiled.λ_ip1 => 3.0,
        compiled.z_i => 0.0, compiled.z_ip1 => 1.0, compiled.z_h => 0.5,
    )
    prob2 = ODEProblem(compiled, params2, (0.0, 1.0))
    sol2 = solve(prob2)

    # Eq. 4.12: λ = 1*3*1 / (1*0.5 + 3*0.5) = 3/2 = 1.5
    @test sol2[compiled.λ_interface][end] ≈ 1.5 rtol = 1.0e-10
end

# ========================================================================
# Heat Flux Tests
# ========================================================================

@testitem "HeatFlux - Direction and Magnitude" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = HeatFlux()
    compiled = mtkcompile(sys)

    # F = -λ * (T_i - T_{i+1}) / (z_{i+1} - z_i)
    params = Dict(
        compiled.λ_interface => 1.0,
        compiled.T_i => 300.0, compiled.T_ip1 => 290.0,
        compiled.z_i => 0.05, compiled.z_ip1 => 0.15,
    )
    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # F = -1.0 * (300 - 290) / (0.15 - 0.05) = -100 W/m²
    @test sol[compiled.F][end] ≈ -100.0 rtol = 1.0e-10
    @test sol[compiled.F][end] < 0  # downward: heat flows from hot to cold

    # Equal temperatures: zero flux
    params_eq = Dict(
        compiled.λ_interface => 1.0,
        compiled.T_i => 300.0, compiled.T_ip1 => 300.0,
        compiled.z_i => 0.05, compiled.z_ip1 => 0.15,
    )
    prob_eq = ODEProblem(compiled, params_eq, (0.0, 1.0))
    sol_eq = solve(prob_eq)
    @test sol_eq[compiled.F][end] ≈ 0.0 atol = 1.0e-12
end

# ========================================================================
# Surface Energy Flux Tests
# ========================================================================

@testitem "SurfaceEnergyFlux - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SurfaceEnergyFlux()
    compiled = mtkcompile(sys)

    params = Dict(
        compiled.S_g => 200.0, compiled.L_g => 50.0,
        compiled.H_g => 30.0, compiled.λE_g => 20.0,
        compiled.H_wasteheat_g => 10.0, compiled.H_aircond_g => 5.0,
        compiled.ε_g => 0.95, compiled.T_g_n => 300.0,
        compiled.dH_g_dT => 5.0, compiled.dλE_g_dT => 2.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.26
    expected_h = 200.0 - 50.0 - 30.0 - 20.0 + 10.0 + 5.0
    @test sol[compiled.h][end] ≈ expected_h rtol = 1.0e-10

    # Eq. 4.29
    σ = 5.67e-8
    expected_dLdT = 4 * 0.95 * σ * 300.0^3
    @test sol[compiled.dL_g_dT][end] ≈ expected_dLdT rtol = 1.0e-6

    # Eq. 4.28
    expected_dhdT = -expected_dLdT - 5.0 - 2.0
    @test sol[compiled.dh_dT][end] ≈ expected_dhdT rtol = 1.0e-6
end

# ========================================================================
# Building Temperature Tests
# ========================================================================

@testitem "BuildingTemperature - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = BuildingTemperature()
    compiled = mtkcompile(sys)

    H = 10.0
    H_W = 1.0
    W_roof = 0.3
    T_shd = 293.0
    T_sun = 295.0
    T_roof_val = 290.0

    params = Dict(
        compiled.H_canyon => H, compiled.H_W => H_W,
        compiled.W_roof => W_roof, compiled.T_inner_shdwall => T_shd,
        compiled.T_inner_sunwall => T_sun, compiled.T_inner_roof => T_roof_val,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.38: L = (H / H_W) * (W_roof / (1 - W_roof))
    expected_L = (H / H_W) * (W_roof / (1 - W_roof))
    @test sol[compiled.L_roof][end] ≈ expected_L rtol = 1.0e-10

    # Eq. 4.37
    expected_T = (H * (T_shd + T_sun) + expected_L * T_roof_val) / (2 * H + expected_L)
    @test sol[compiled.T_iB_unclamped][end] ≈ expected_T rtol = 1.0e-10

    # T_iB should be between min and max of input temperatures
    @test sol[compiled.T_iB_unclamped][end] ≥ min(T_shd, T_sun, T_roof_val) - 1.0e-10
    @test sol[compiled.T_iB_unclamped][end] ≤ max(T_shd, T_sun, T_roof_val) + 1.0e-10
end

@testitem "BuildingTemperature - Equal Temperatures" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = BuildingTemperature()
    compiled = mtkcompile(sys)

    T_uniform = 295.0
    params = Dict(
        compiled.H_canyon => 10.0, compiled.H_W => 1.0, compiled.W_roof => 0.3,
        compiled.T_inner_shdwall => T_uniform,
        compiled.T_inner_sunwall => T_uniform,
        compiled.T_inner_roof => T_uniform,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)
    @test sol[compiled.T_iB_unclamped][end] ≈ T_uniform rtol = 1.0e-10
end

# ========================================================================
# Waste Heat / Air Conditioning Tests
# ========================================================================

@testitem "WasteHeatAirConditioning - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = WasteHeatAirConditioning()
    compiled = mtkcompile(sys)

    params = Dict(
        compiled.F_heat_roof => 10.0, compiled.F_cool_roof => 5.0,
        compiled.F_heat_sunwall => 8.0, compiled.F_cool_sunwall => 3.0,
        compiled.F_heat_shdwall => 7.0, compiled.F_cool_shdwall => 2.0,
        compiled.W_roof => 0.3, compiled.H_W => 1.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    f_heat = 1 / 0.75
    f_cool = 1 / 0.25

    expected_unclamped = 0.3 * (f_heat * 10.0 + f_cool * 5.0) +
        0.7 * 1.0 * (f_heat * 8.0 + f_cool * 3.0 + f_heat * 7.0 + f_cool * 2.0)
    @test sol[compiled.H_wasteheat_unclamped][end] ≈ expected_unclamped rtol = 1.0e-6
    @test sol[compiled.H_wasteheat][end] ≈ min(expected_unclamped, 100.0) rtol = 1.0e-6
    @test sol[compiled.H_aircond][end] ≈ 5.0 + 3.0 + 2.0 rtol = 1.0e-10
end

@testitem "WasteHeatAirConditioning - Max Clamping" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = WasteHeatAirConditioning()
    compiled = mtkcompile(sys)

    params = Dict(
        compiled.F_heat_roof => 100.0, compiled.F_cool_roof => 50.0,
        compiled.F_heat_sunwall => 80.0, compiled.F_cool_sunwall => 30.0,
        compiled.F_heat_shdwall => 70.0, compiled.F_cool_shdwall => 20.0,
        compiled.W_roof => 0.3, compiled.H_W => 1.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)
    @test sol[compiled.H_wasteheat][end] ≈ 100.0 rtol = 1.0e-6
end

# ========================================================================
# Phase Change Energy Tests
# ========================================================================

@testitem "PhaseChangeEnergy Interior" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeEnergy(; layer_type = :interior)
    compiled = mtkcompile(sys)

    T_f = 273.15
    α = 0.5
    c_i = 2.0e6
    Δz = 0.1
    Δt = 3600.0
    T_n = 274.0
    F_i_n = -10.0
    F_im1_n = -5.0
    F_i_np1 = -12.0
    F_im1_np1 = -6.0

    params = Dict(
        compiled.c_i => c_i, compiled.Δz_i => Δz, compiled.Δt => Δt,
        compiled.T_i_n => T_n,
        compiled.F_i_n => F_i_n, compiled.F_i_np1 => F_i_np1,
        compiled.F_im1_n => F_im1_n, compiled.F_im1_np1 => F_im1_np1,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    expected = α * (F_i_n - F_im1_n) + (1 - α) * (F_i_np1 - F_im1_np1) -
        (c_i * Δz / Δt) * (T_f - T_n)
    @test sol[compiled.H_excess][end] ≈ expected rtol = 1.0e-10
end

@testitem "PhaseChangeEnergy Top" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeEnergy(; layer_type = :top)
    compiled = mtkcompile(sys)

    T_f = 273.15
    α = 0.5
    c_i = 2.0e6
    Δz = 0.01
    Δt = 3600.0
    T_n = 272.0
    h_n = 100.0
    dh_dT = -15.0
    F_i_n = -5.0
    F_i_np1 = -6.0

    params = Dict(
        compiled.c_i => c_i, compiled.Δz_i => Δz, compiled.Δt => Δt,
        compiled.T_i_n => T_n,
        compiled.h_n => h_n, compiled.dh_dT => dh_dT,
        compiled.F_i_n => F_i_n, compiled.F_i_np1 => F_i_np1,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    expected = h_n + dh_dT * (T_f - T_n) + α * F_i_n + (1 - α) * F_i_np1 -
        (c_i * Δz / Δt) * (T_f - T_n)
    @test sol[compiled.H_excess][end] ≈ expected rtol = 1.0e-10
end

# ========================================================================
# Waste Heat Allocation Tests
# ========================================================================

@testitem "WasteHeatAllocation - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = WasteHeatAllocation()

    @test length(equations(sys)) == 4
    @test length(unknowns(sys)) == 4

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "WasteHeatAllocation - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = WasteHeatAllocation()
    compiled = mtkcompile(sys)

    H_wh = 80.0
    H_ac = 20.0
    W_roof = 0.3

    params = Dict(
        compiled.H_wasteheat => H_wh,
        compiled.H_aircond => H_ac,
        compiled.W_roof => W_roof,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.27: waste heat divided by (1 - W_roof)
    expected_wh = H_wh / (1 - W_roof)
    @test sol[compiled.H_wasteheat_prvrd][end] ≈ expected_wh rtol = 1.0e-10
    @test sol[compiled.H_wasteheat_imprvrd][end] ≈ expected_wh rtol = 1.0e-10

    expected_ac = H_ac / (1 - W_roof)
    @test sol[compiled.H_aircond_prvrd][end] ≈ expected_ac rtol = 1.0e-10
    @test sol[compiled.H_aircond_imprvrd][end] ≈ expected_ac rtol = 1.0e-10
end

@testitem "WasteHeatAllocation - Equal Distribution" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = WasteHeatAllocation()
    compiled = mtkcompile(sys)

    # Both road types get the same allocation
    params = Dict(
        compiled.H_wasteheat => 50.0,
        compiled.H_aircond => 10.0,
        compiled.W_roof => 0.25,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.H_wasteheat_prvrd][end] ≈ sol[compiled.H_wasteheat_imprvrd][end] rtol = 1.0e-10
    @test sol[compiled.H_aircond_prvrd][end] ≈ sol[compiled.H_aircond_imprvrd][end] rtol = 1.0e-10
end

# ========================================================================
# Adjusted Layer Thickness Tests
# ========================================================================

@testitem "AdjustedLayerThickness - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = AdjustedLayerThickness()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "AdjustedLayerThickness - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = AdjustedLayerThickness()
    compiled = mtkcompile(sys)

    z_i = 0.05
    z_h_im1 = 0.0
    z_ip1 = 0.15
    c_a = 0.34

    params = Dict(
        compiled.z_i => z_i,
        compiled.z_h_im1 => z_h_im1,
        compiled.z_ip1 => z_ip1,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.30: Δz* = 0.5 * [z_i - z_{h,i-1} + c_a*(z_{i+1} - z_{h,i-1})]
    expected = 0.5 * (z_i - z_h_im1 + c_a * (z_ip1 - z_h_im1))
    @test sol[compiled.Δz_star][end] ≈ expected rtol = 1.0e-10
end

@testitem "AdjustedLayerThickness - Positive Output" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = AdjustedLayerThickness()
    compiled = mtkcompile(sys)

    # Typical road layer geometry
    params = Dict(
        compiled.z_i => 0.1,
        compiled.z_h_im1 => 0.0,
        compiled.z_ip1 => 0.3,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.Δz_star][end] > 0
end

# ========================================================================
# Heating/Cooling Flux Tests
# ========================================================================

@testitem "HeatingCoolingFlux - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = HeatingCoolingFlux()

    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "HeatingCoolingFlux - Heating Active" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = HeatingCoolingFlux()
    compiled = mtkcompile(sys)

    # T_iB < T_min => heating active
    params = Dict(
        compiled.T_iB => 285.0,
        compiled.T_iB_min => 290.0,
        compiled.T_iB_max => 300.0,
        compiled.F_bottom_n => -20.0,
        compiled.F_bottom_np1 => -25.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Combined flux: 0.5*(-20) + 0.5*(-25) = -22.5
    expected_combined = 0.5 * (-20.0) + 0.5 * (-25.0)
    @test sol[compiled.F_combined][end] ≈ expected_combined rtol = 1.0e-10

    # Heating is active: F_heat = |F_combined|
    @test sol[compiled.F_heat][end] ≈ abs(expected_combined) rtol = 1.0e-6

    # Cooling is inactive
    @test sol[compiled.F_cool][end] ≈ 0.0 atol = 1.0e-10
end

@testitem "HeatingCoolingFlux - Cooling Active" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = HeatingCoolingFlux()
    compiled = mtkcompile(sys)

    # T_iB > T_max => cooling active
    params = Dict(
        compiled.T_iB => 305.0,
        compiled.T_iB_min => 290.0,
        compiled.T_iB_max => 300.0,
        compiled.F_bottom_n => 30.0,
        compiled.F_bottom_np1 => 35.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    expected_combined = 0.5 * 30.0 + 0.5 * 35.0
    @test sol[compiled.F_cool][end] ≈ abs(expected_combined) rtol = 1.0e-6
    @test sol[compiled.F_heat][end] ≈ 0.0 atol = 1.0e-10
end

@testitem "HeatingCoolingFlux - No HVAC Needed" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = HeatingCoolingFlux()
    compiled = mtkcompile(sys)

    # T_min < T_iB < T_max => no heating or cooling
    params = Dict(
        compiled.T_iB => 295.0,
        compiled.T_iB_min => 290.0,
        compiled.T_iB_max => 300.0,
        compiled.F_bottom_n => 10.0,
        compiled.F_bottom_np1 => 12.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.F_heat][end] ≈ 0.0 atol = 1.0e-10
    @test sol[compiled.F_cool][end] ≈ 0.0 atol = 1.0e-10
end

# ========================================================================
# Phase Change Adjustment Tests
# ========================================================================

@testitem "PhaseChangeAdjustment Interior - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeAdjustment(; layer_type = :interior)

    @test length(equations(sys)) == 5
    @test length(unknowns(sys)) == 5

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "PhaseChangeAdjustment Interior - Melting" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeAdjustment(; layer_type = :interior)
    compiled = mtkcompile(sys)

    L_f = 3.337e5
    T_f = 273.15
    Δt = 3600.0
    c_i = 2.0e6
    Δz_i = 0.1

    # Positive H_i means energy available for melting
    H_i = 50.0  # W/m²
    w_ice_n = 10.0  # kg/m²
    w_liq_n = 5.0   # kg/m²

    params = Dict(
        compiled.H_i => H_i,
        compiled.w_ice_n => w_ice_n,
        compiled.w_liq_n => w_liq_n,
        compiled.Δt => Δt,
        compiled.c_i => c_i,
        compiled.Δz_i => Δz_i,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # H_m = H_i * Δt / L_f
    expected_Hm = H_i * Δt / L_f
    @test sol[compiled.H_m][end] ≈ expected_Hm rtol = 1.0e-10

    # Ice decreases by melting
    @test sol[compiled.w_ice_np1][end] < w_ice_n

    # Mass conservation: w_ice + w_liq = const
    total_before = w_ice_n + w_liq_n
    total_after = sol[compiled.w_ice_np1][end] + sol[compiled.w_liq_np1][end]
    @test total_after ≈ total_before rtol = 1.0e-10
end

@testitem "PhaseChangeAdjustment Interior - Freezing" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeAdjustment(; layer_type = :interior)
    compiled = mtkcompile(sys)

    L_f = 3.337e5

    # Negative H_i means energy deficit => freezing
    H_i = -50.0  # W/m²
    w_ice_n = 5.0   # kg/m²
    w_liq_n = 10.0  # kg/m²

    params = Dict(
        compiled.H_i => H_i,
        compiled.w_ice_n => w_ice_n,
        compiled.w_liq_n => w_liq_n,
        compiled.Δt => 3600.0,
        compiled.c_i => 2.0e6,
        compiled.Δz_i => 0.1,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Ice increases by freezing
    @test sol[compiled.w_ice_np1][end] > w_ice_n

    # Mass conservation
    total_before = w_ice_n + w_liq_n
    total_after = sol[compiled.w_ice_np1][end] + sol[compiled.w_liq_np1][end]
    @test total_after ≈ total_before rtol = 1.0e-10
end

@testitem "PhaseChangeAdjustment Top - Temperature Correction" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeAdjustment(; layer_type = :top)
    compiled = mtkcompile(sys)

    T_f = 273.15
    L_f = 3.337e5
    Δt = 3600.0
    c_i = 2.0e6
    Δz_i = 0.1
    dh_dT = -10.0

    # Small melting event
    H_i = 10.0
    w_ice_n = 20.0
    w_liq_n = 5.0

    params = Dict(
        compiled.H_i => H_i,
        compiled.w_ice_n => w_ice_n,
        compiled.w_liq_n => w_liq_n,
        compiled.Δt => Δt,
        compiled.c_i => c_i,
        compiled.Δz_i => Δz_i,
        compiled.dh_dT => dh_dT,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.65 (top): T = T_f + (Δt/(c*Δz)) * H_res / (1 - (Δt/(c*Δz)) * dh/dT)
    H_res = sol[compiled.H_residual][end]
    expected_T = T_f + (Δt / (c_i * Δz_i)) * H_res / (1 - (Δt / (c_i * Δz_i)) * dh_dT)
    @test sol[compiled.T_np1][end] ≈ expected_T rtol = 1.0e-10
end

# ========================================================================
# Snow Melt No Layers Tests
# ========================================================================

@testitem "SnowMeltNoLayers - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowMeltNoLayers()

    @test length(equations(sys)) == 5
    @test length(unknowns(sys)) == 5

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowMeltNoLayers - Partial Melt" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowMeltNoLayers()
    compiled = mtkcompile(sys)

    L_f = 3.337e5
    H_1 = 50.0    # Excess energy
    W_sno = 10.0   # Snow mass
    z_sno = 0.05   # Snow depth
    Δt = 3600.0

    params = Dict(
        compiled.H_1 => H_1,
        compiled.W_sno_n => W_sno,
        compiled.z_sno_n => z_sno,
        compiled.Δt => Δt,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.66: W_sno^{n+1} = max(W_sno - H_1*Δt/L_f, 0)
    expected_W = max(W_sno - H_1 * Δt / L_f, 0)
    @test sol[compiled.W_sno_np1][end] ≈ expected_W rtol = 1.0e-10

    # Snow decreased
    @test sol[compiled.W_sno_np1][end] < W_sno

    # Eq. 4.67: z_sno^{n+1} = (W^{n+1}/W^n) * z^n
    expected_z = (expected_W / W_sno) * z_sno
    @test sol[compiled.z_sno_np1][end] ≈ expected_z rtol = 1.0e-10

    # Melt rate positive
    @test sol[compiled.M_1S][end] > 0

    # Eq. 4.71: E_p = L_f * M_1S
    @test sol[compiled.E_p1S][end] ≈ L_f * sol[compiled.M_1S][end] rtol = 1.0e-10
end

@testitem "SnowMeltNoLayers - Complete Melt" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowMeltNoLayers()
    compiled = mtkcompile(sys)

    L_f = 3.337e5
    # Large excess energy that would melt all snow
    H_1 = 500.0
    W_sno = 0.5    # Small snow mass
    z_sno = 0.002
    Δt = 3600.0

    params = Dict(
        compiled.H_1 => H_1,
        compiled.W_sno_n => W_sno,
        compiled.z_sno_n => z_sno,
        compiled.Δt => Δt,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # All snow melts
    @test sol[compiled.W_sno_np1][end] ≈ 0.0 atol = 1.0e-10
    @test sol[compiled.z_sno_np1][end] ≈ 0.0 atol = 1.0e-10

    # Residual energy is positive (excess after melting all snow)
    @test sol[compiled.H_residual][end] > 0
end

@testitem "SnowMeltNoLayers - Zero Excess Energy" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowMeltNoLayers()
    compiled = mtkcompile(sys)

    # No excess energy => no melting
    params = Dict(
        compiled.H_1 => 0.0,
        compiled.W_sno_n => 5.0,
        compiled.z_sno_n => 0.03,
        compiled.Δt => 3600.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Snow unchanged
    @test sol[compiled.W_sno_np1][end] ≈ 5.0 rtol = 1.0e-10
    @test sol[compiled.z_sno_np1][end] ≈ 0.03 rtol = 1.0e-10
    @test sol[compiled.M_1S][end] ≈ 0.0 atol = 1.0e-10
end

# ========================================================================
# Urban Surface Thermal Properties Tests
# ========================================================================

@testitem "UrbanSurfaceThermalProperties - Passthrough" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = UrbanSurfaceThermalProperties()
    compiled = mtkcompile(sys)

    params = Dict(compiled.λ_prescribed => 1.5, compiled.c_prescribed => 2.5e6)

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.λ_surf][end] ≈ 1.5 rtol = 1.0e-10
    @test sol[compiled.c_surf][end] ≈ 2.5e6 rtol = 1.0e-10
end

# ========================================================================
# Uniform Grid Tests
# ========================================================================

@testitem "UniformGrid - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = UniformGrid(N = 3)

    @test length(equations(sys)) == 10  # 3 nodes + 3 thicknesses + 4 interfaces
    @test length(unknowns(sys)) == 10

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "UniformGrid - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    N = 3
    sys = UniformGrid(N = N)
    compiled = mtkcompile(sys)

    Δz_total = 0.3  # 30 cm total thickness

    prob = ODEProblem(compiled, [compiled.Δz_total => Δz_total], (0.0, 1.0))
    sol = solve(prob)

    # Access observed variables through the observed list
    obs = Dict(string(o.lhs) => o.lhs for o in observed(compiled))

    # Eq. 4.5: z_i = (i - 0.5) * (Δz / N)
    for i in 1:N
        expected_z = (i - 0.5) * (Δz_total / N)
        @test sol[obs["z_node[$i](t)"]][end] ≈ expected_z rtol = 1.0e-10
    end

    # All layer thicknesses should be equal for uniform grid = Δz_total / N
    expected_Δz = Δz_total / N
    for i in 1:N
        @test sol[obs["Δz_layer[$i](t)"]][end] ≈ expected_Δz rtol = 1.0e-6
    end

    # z_interface[1] = z_{h,0} = 0 (top surface)
    @test sol[obs["z_interface[1](t)"]][end] ≈ 0.0 atol = 1.0e-12

    # z_interface[N+1] = z_{h,N} = Δz_total (bottom)
    @test sol[obs["z_interface[$(N + 1)](t)"]][end] ≈ Δz_total rtol = 1.0e-6
end

@testitem "UniformGrid - Node at Midpoint" setup = [TempSetup] tags = [:ch4_temps] begin
    N = 5
    sys = UniformGrid(N = N)
    compiled = mtkcompile(sys)

    prob = ODEProblem(compiled, [compiled.Δz_total => 1.0], (0.0, 1.0))
    sol = solve(prob)

    obs = Dict(string(o.lhs) => o.lhs for o in observed(compiled))

    # Each node should be at the midpoint of its layer
    for i in 1:N
        node = sol[obs["z_node[$i](t)"]][end]
        z_top = sol[obs["z_interface[$i](t)"]][end]
        z_bot = sol[obs["z_interface[$(i + 1)](t)"]][end]
        midpoint = 0.5 * (z_top + z_bot)
        @test node ≈ midpoint rtol = 1.0e-6
    end
end

# ========================================================================
# Exponential Grid Tests
# ========================================================================

@testitem "ExponentialGrid - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = ExponentialGrid(N = 3)

    @test length(equations(sys)) == 10
    @test length(unknowns(sys)) == 10

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "ExponentialGrid - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    N = 5
    sys = ExponentialGrid(N = N)
    compiled = mtkcompile(sys)

    prob = ODEProblem(compiled, [], (0.0, 1.0))
    sol = solve(prob)

    # Access observed variables through the observed list
    obs = Dict(string(o.lhs) => o.lhs for o in observed(compiled))

    f_s = 0.025

    # Eq. 4.8: z_i = f_s * (exp(0.5*(i-0.5)) - 1)
    for i in 1:N
        expected_z = f_s * (exp(0.5 * (i - 0.5)) - 1)
        @test sol[obs["z_node[$i](t)"]][end] ≈ expected_z rtol = 1.0e-10
    end

    # z_interface[1] = z_{h,0} = 0
    @test sol[obs["z_interface[1](t)"]][end] ≈ 0.0 atol = 1.0e-12

    # Layer thicknesses should increase with depth (exponential spacing)
    for i in 2:N
        @test sol[obs["Δz_layer[$i](t)"]][end] > sol[obs["Δz_layer[$(i - 1)](t)"]][end]
    end
end

@testitem "ExponentialGrid - Increasing Depth" setup = [TempSetup] tags = [:ch4_temps] begin
    N = 15
    sys = ExponentialGrid(N = N)
    compiled = mtkcompile(sys)

    prob = ODEProblem(compiled, [], (0.0, 1.0))
    sol = solve(prob)

    # Access observed variables through the observed list
    obs = Dict(string(o.lhs) => o.lhs for o in observed(compiled))

    # Node depths should be monotonically increasing
    for i in 2:N
        @test sol[obs["z_node[$i](t)"]][end] > sol[obs["z_node[$(i - 1)](t)"]][end]
    end

    # Interface depths should be monotonically increasing
    for i in 2:(N + 1)
        @test sol[obs["z_interface[$i](t)"]][end] > sol[obs["z_interface[$(i - 1)](t)"]][end]
    end
end

# ========================================================================
# Freezing Point Depression Tests
# ========================================================================

@testitem "FreezingPointDepression - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = FreezingPointDepression()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "FreezingPointDepression - Above Freezing" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = FreezingPointDepression()
    compiled = mtkcompile(sys)

    # Above freezing: w_liq_max = Δz * θ_sat * ρ_liq (fully saturated)
    params = Dict(
        compiled.T_i => 280.0,  # Above T_f = 273.15
        compiled.θ_sat => 0.4,
        compiled.Δz => 0.5,
        compiled.ψ_sat => -0.1,  # Negative matric potential in m
        compiled.B_i => 5.0,
        compiled.ρ_liq => 1000.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    expected = 0.5 * 0.4 * 1000.0  # Δz * θ_sat * ρ_liq = 200 kg/m²
    @test sol[compiled.w_liq_max][end] ≈ expected rtol = 1.0e-10
end

@testitem "FreezingPointDepression - Below Freezing" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = FreezingPointDepression()
    compiled = mtkcompile(sys)

    T_f = 273.15
    T_i = 268.15  # 5 K below freezing
    θ_sat = 0.4
    Δz = 0.5
    ψ_sat = -0.1  # m (negative in SI)
    B_i = 5.0
    ρ_liq = 1000.0
    L_f = 3.337e5
    g = 9.80616

    params = Dict(
        compiled.T_i => T_i,
        compiled.θ_sat => θ_sat,
        compiled.Δz => Δz,
        compiled.ψ_sat => ψ_sat,
        compiled.B_i => B_i,
        compiled.ρ_liq => ρ_liq,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Below freezing: w_liq_max < Δz * θ_sat * ρ_liq
    max_possible = Δz * θ_sat * ρ_liq
    @test sol[compiled.w_liq_max][end] < max_possible
    @test sol[compiled.w_liq_max][end] > 0
end

@testitem "FreezingPointDepression - Colder Means Less Liquid" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = FreezingPointDepression()
    compiled = mtkcompile(sys)

    base_params = Dict(
        compiled.θ_sat => 0.4,
        compiled.Δz => 0.5,
        compiled.ψ_sat => -0.1,
        compiled.B_i => 5.0,
        compiled.ρ_liq => 1000.0,
    )

    # Slightly below freezing
    params1 = merge(base_params, Dict(compiled.T_i => 272.15))
    prob1 = ODEProblem(compiled, params1, (0.0, 1.0))
    sol1 = solve(prob1)

    # Much below freezing
    params2 = merge(base_params, Dict(compiled.T_i => 263.15))
    prob2 = ODEProblem(compiled, params2, (0.0, 1.0))
    sol2 = solve(prob2)

    # Colder temperature => less liquid water
    @test sol2[compiled.w_liq_max][end] < sol1[compiled.w_liq_max][end]
end

# ========================================================================
# Snow-Soil Blended Heat Capacity Tests
# ========================================================================

@testitem "SnowSoilBlendedHeatCapacity - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowSoilBlendedHeatCapacity()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowSoilBlendedHeatCapacity - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowSoilBlendedHeatCapacity()
    compiled = mtkcompile(sys)

    c_soil = 2.0e6
    W_sno = 5.0
    Δz = 0.1
    C_ice = 2117.27

    params = Dict(
        compiled.c_soil => c_soil,
        compiled.W_sno => W_sno,
        compiled.Δz => Δz,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.88: c_i = c_soil + C_ice * W_sno / Δz
    expected = c_soil + C_ice * W_sno / Δz
    @test sol[compiled.c_blended][end] ≈ expected rtol = 1.0e-10
end

@testitem "SnowSoilBlendedHeatCapacity - No Snow" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = SnowSoilBlendedHeatCapacity()
    compiled = mtkcompile(sys)

    c_soil = 2.0e6
    params = Dict(
        compiled.c_soil => c_soil,
        compiled.W_sno => 0.0,
        compiled.Δz => 0.1,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # With no snow, blended = soil only
    @test sol[compiled.c_blended][end] ≈ c_soil rtol = 1.0e-10
end

# ========================================================================
# Layer Phase Change Energy Tests
# ========================================================================

@testitem "LayerPhaseChangeEnergy - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = LayerPhaseChangeEnergy()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "LayerPhaseChangeEnergy - Melting" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = LayerPhaseChangeEnergy()
    compiled = mtkcompile(sys)

    L_f = 3.337e5
    w_ice_n = 10.0
    w_ice_np1 = 8.0  # 2 kg/m² melted
    Δt = 3600.0

    params = Dict(
        compiled.w_ice_n => w_ice_n,
        compiled.w_ice_np1 => w_ice_np1,
        compiled.Δt => Δt,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.73: E_p = L_f * (w_ice_n - w_ice_np1) / Δt
    expected = L_f * (w_ice_n - w_ice_np1) / Δt
    @test sol[compiled.E_p][end] ≈ expected rtol = 1.0e-10
    @test sol[compiled.E_p][end] > 0  # Positive for melting
end

@testitem "LayerPhaseChangeEnergy - No Change" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = LayerPhaseChangeEnergy()
    compiled = mtkcompile(sys)

    params = Dict(
        compiled.w_ice_n => 10.0,
        compiled.w_ice_np1 => 10.0,  # No change
        compiled.Δt => 3600.0,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.E_p][end] ≈ 0.0 atol = 1.0e-12
end

# ========================================================================
# Total Phase Change Energy Tests
# ========================================================================

@testitem "TotalPhaseChangeEnergy - Structural" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = TotalPhaseChangeEnergy()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "TotalPhaseChangeEnergy - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = TotalPhaseChangeEnergy()
    compiled = mtkcompile(sys)

    E_p1S = 50.0
    E_p_layers = 30.0

    params = Dict(
        compiled.E_p1S => E_p1S,
        compiled.E_p_layers => E_p_layers,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.72: E_p = E_{p,1S} + Σ E_{p,i}
    @test sol[compiled.E_p_total][end] ≈ E_p1S + E_p_layers rtol = 1.0e-10
end

# ========================================================================
# PDE Heat Conduction Tests (MethodOfLines.jl)
# ========================================================================

@testitem "RoofWallHeatConduction - Steady State" setup = [TempSetup] tags = [:ch4_temps] begin
    using MethodOfLines, DomainSets

    # Steady-state test: zero surface flux with constant T_bottom
    # should yield uniform temperature
    T_bottom = 293.15
    result = RoofWallHeatConduction(;
        Δz_total = 0.3, N_layers = 15,
        λ_val = 1.0, c_val = 2.0e6,
        h_top = 0.0, T_bottom = T_bottom,
    )

    sol = solve(result.prob, saveat = 0.1)
    T_mat = sol[result.T(result.t_pde, result.z_var)]

    # At steady state with zero flux at top and T_iB at bottom,
    # the temperature should be T_bottom everywhere
    T_final = T_mat[end, :]
    for Tv in T_final
        @test Tv ≈ T_bottom rtol = 1.0e-3
    end
end

@testitem "RoofWallHeatConduction - Heat Flux Response" setup = [TempSetup] tags = [:ch4_temps] begin
    using MethodOfLines, DomainSets

    # Apply constant positive heat flux at top
    T_bottom = 288.15
    h_top = 100.0  # 100 W/m²
    result = RoofWallHeatConduction(;
        Δz_total = 0.3, N_layers = 15,
        λ_val = 1.0, c_val = 2.0e6,
        h_top = h_top, T_bottom = T_bottom,
    )

    sol = solve(result.prob, saveat = 0.1)
    T_mat = sol[result.T(result.t_pde, result.z_var)]

    # Top surface should be warmer than bottom after heating
    T_final = T_mat[end, :]
    @test T_final[1] > T_bottom  # top layer heated above bottom temperature

    # Temperature should be monotonically decreasing from top to bottom
    for i in 1:(length(T_final) - 1)
        @test T_final[i] ≥ T_final[i + 1] - 1.0e-6
    end
end

@testitem "RoofWallHeatConduction - Analytical Steady State with Flux" setup = [TempSetup] tags = [:ch4_temps] begin
    using MethodOfLines, DomainSets

    # Analytical solution for steady state with constant flux at top:
    # At steady state: ∂T/∂t = 0, so d²T/dz² = 0 -> T(z) = a*z + b
    # BC: -λ*dT/dz|_{z=0} = h -> -λ*a = h -> a = -h/λ
    # BC: T(L) = T_iB -> T_iB = a*L + b -> b = T_iB - a*L = T_iB + h*L/λ
    # T(z) = T_iB + h*(L-z)/λ
    # T(0) = T_iB + h*L/λ

    L = 0.3
    λ = 1.5
    h = 50.0
    T_iB = 290.0

    result = RoofWallHeatConduction(;
        Δz_total = L, N_layers = 30,  # Use finer grid for better accuracy
        λ_val = λ, c_val = 2.0e6,
        h_top = h, T_bottom = T_iB,
    )

    # Solve for longer time to reach steady state
    tspan = (0.0, 100000.0)
    prob2 = remake(result.prob; tspan = tspan)
    sol = solve(prob2)
    T_mat = sol[result.T(result.t_pde, result.z_var)]

    T_final = T_mat[end, :]

    # Surface temperature should match analytical solution
    T_surface_analytical = T_iB + h * L / λ
    @test T_final[1] ≈ T_surface_analytical rtol = 0.05

    # Bottom temperature should be T_iB
    @test T_final[end] ≈ T_iB rtol = 1.0e-3
end

@testitem "RoadHeatConduction - Steady State Zero Flux" setup = [TempSetup] tags = [:ch4_temps] begin
    using MethodOfLines, DomainSets

    # With zero flux at both boundaries, temperature should remain at initial value
    T_init = 288.15
    result = RoadHeatConduction(;
        N_layers = 15,
        λ_val = 1.5, c_val = 2.0e6,
        h_top = 0.0, T_init_val = T_init,
    )

    sol = solve(result.prob, saveat = 0.1)
    T_mat = sol[result.T(result.t_pde, result.z_var)]

    # All temperatures should remain at initial value (zero flux everywhere)
    T_final = T_mat[end, :]
    for Tv in T_final
        @test Tv ≈ T_init rtol = 1.0e-3
    end
end

@testitem "RoadHeatConduction - Heat Flux Response" setup = [TempSetup] tags = [:ch4_temps] begin
    using MethodOfLines, DomainSets

    # Apply heat flux at top, zero flux at bottom
    T_init = 288.15
    h_top = 100.0
    result = RoadHeatConduction(;
        N_layers = 15,
        λ_val = 1.5, c_val = 2.0e6,
        h_top = h_top, T_init_val = T_init,
    )

    sol = solve(result.prob, saveat = 0.1)
    T_mat = sol[result.T(result.t_pde, result.z_var)]

    # With positive heat flux at top, all temperatures should increase
    T_final = T_mat[end, :]
    @test T_final[1] > T_init  # top gets heated

    # Surface should be warmest
    @test T_final[1] ≥ T_final[end] - 1.0e-6
end

@testitem "RoadHeatConduction - Energy Conservation" setup = [TempSetup] tags = [:ch4_temps] begin
    using MethodOfLines, DomainSets

    # With zero flux at both boundaries, total energy should be conserved
    T_init = 288.15
    result = RoadHeatConduction(;
        N_layers = 15,
        λ_val = 1.5, c_val = 2.0e6,
        h_top = 0.0, T_init_val = T_init,
    )

    sol = solve(result.prob, saveat = [0.0, 1.0])
    T_mat = sol[result.T(result.t_pde, result.z_var)]

    # With zero flux, average temperature should be constant
    T_initial = T_mat[1, :]
    T_final = T_mat[end, :]
    @test sum(T_initial) ≈ sum(T_final) rtol = 1.0e-3
end
