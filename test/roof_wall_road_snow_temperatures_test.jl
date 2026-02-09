@testsnippet TempSetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy
end

# ========================================================================
# Grid Discretization Tests
# ========================================================================

@testitem "compute_grid_roofwall - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    Δz_total = 0.3
    N = 15
    grid = compute_grid_roofwall(Δz_total; N_levgrnd = N)

    # Eq. 4.5: z_i = (i - 0.5) * (Δz / N)
    for i in 1:N
        expected = (i - 0.5) * (Δz_total / N)
        @test grid.z_node[i] ≈ expected rtol = 1e-10
    end

    # Eq. 4.6: For uniform grid, all interior layer thicknesses should be Δz/N
    spacing = Δz_total / N
    for i in 2:(N - 1)
        @test grid.Δz_layer[i] ≈ spacing rtol = 1e-10
    end

    # Eq. 4.7: Top interface should be at 0
    @test grid.z_interface[1] ≈ 0.0 atol = 1e-15

    # Sum of layer thicknesses should equal total thickness
    @test sum(grid.Δz_layer) ≈ Δz_total rtol = 1e-10

    # Bottom interface should be at total thickness
    @test grid.z_interface[N + 1] ≈ Δz_total rtol = 1e-10

    # Node depths should be monotonically increasing
    for i in 1:(N - 1)
        @test grid.z_node[i + 1] > grid.z_node[i]
    end
end

@testitem "compute_grid_road - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    grid = compute_grid_road()
    f_s = 0.025
    N = 15

    # Eq. 4.8: z_i = f_s * {exp[0.5(i - 0.5)] - 1}
    for i in 1:N
        expected = f_s * (exp(0.5 * (i - 0.5)) - 1)
        @test grid.z_node[i] ≈ expected rtol = 1e-10
    end

    # Exponential grid: node depths monotonically increasing
    for i in 1:(N - 1)
        @test grid.z_node[i + 1] > grid.z_node[i]
    end

    # Layer thicknesses should increase (exponential spacing)
    for i in 1:(N - 1)
        @test grid.Δz_layer[i + 1] > grid.Δz_layer[i]
    end

    # Top interface at 0
    @test grid.z_interface[1] ≈ 0.0 atol = 1e-15

    # Sum of layer thicknesses should approximately equal bottom interface depth
    @test sum(grid.Δz_layer) ≈ grid.z_interface[N + 1] rtol = 1e-10
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
    @test sol[compiled.Δz_snow_0][end] ≈ z_sno_val rtol = 1e-10

    # Node at midpoint of layer: z = 0 - 0.5 * Δz (negative = above surface)
    @test sol[compiled.z_node_snow_0][end] ≈ -0.5 * z_sno_val rtol = 1e-10

    # Top interface: z_{h,-1} = 0 - Δz_0
    @test sol[compiled.z_interface_snow_neg1][end] ≈ -z_sno_val rtol = 1e-10
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
    @test sol[compiled.λ_s][end] ≈ 7.33 rtol = 1e-6

    # Eq. 4.86: c_s
    expected_cs = (2.128 * 60 + 2.385 * 20) / (60 + 20) * 1e6
    @test sol[compiled.c_s][end] ≈ expected_cs rtol = 1e-6

    # Bulk density: ρ_d = 2700 * (1 - 0.4) = 1620
    @test sol[compiled.ρ_d][end] ≈ 1620.0 rtol = 1e-6

    # Eq. 4.80: λ_dry
    expected_λ_dry = (0.135 * 1620 + 64.7) / (2700 - 0.947 * 1620)
    @test sol[compiled.λ_dry][end] ≈ expected_λ_dry rtol = 1e-6

    # Eq. 4.82: S_r = (50/(1000*0.5) + 0/(917*0.5)) / 0.4 = 0.25
    @test sol[compiled.S_r][end] ≈ 0.25 rtol = 1e-6

    # Thermal conductivity: physical bounds
    @test sol[compiled.λ_soil][end] > 0
    @test sol[compiled.λ_soil][end] ≥ sol[compiled.λ_dry][end] - 1e-10

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
    @test sol[compiled.S_r][end] ≈ 0.0 atol = 1e-15
    @test sol[compiled.λ_soil][end] ≈ sol[compiled.λ_dry][end] rtol = 1e-6

    # Dry soil heat capacity: c = c_s * (1 - θ_sat) only
    expected_c = sol[compiled.c_s][end] * (1 - 0.4)
    @test sol[compiled.c_soil][end] ≈ expected_c rtol = 1e-6
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
    @test sol[compiled.K_e][end] ≈ S_r rtol = 1e-6
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
    @test sol[compiled.ρ_sno][end] ≈ 100.0 rtol = 1e-10

    # Eq. 4.83: λ = λ_air + (7.75e-5 * ρ + 1.105e-6 * ρ²) * (λ_ice - λ_air)
    λ_air = 0.023
    λ_ice = 2.29
    ρ = 100.0
    expected_λ = λ_air + (7.75e-5 * ρ + 1.105e-6 * ρ^2) * (λ_ice - λ_air)
    @test sol[compiled.λ_snow][end] ≈ expected_λ rtol = 1e-6

    # Eq. 4.87: c = (w_ice/Δz) * C_ice
    C_ice = 2117.27
    expected_c = (w_ice_val / Δz_val) * C_ice
    @test sol[compiled.c_snow][end] ≈ expected_c rtol = 1e-6

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
    @test sol[compiled.λ_interface][end] ≈ 1.0 rtol = 1e-10

    # Different conductivities: harmonic mean at midpoint
    params2 = Dict(
        compiled.λ_i => 1.0, compiled.λ_ip1 => 3.0,
        compiled.z_i => 0.0, compiled.z_ip1 => 1.0, compiled.z_h => 0.5,
    )
    prob2 = ODEProblem(compiled, params2, (0.0, 1.0))
    sol2 = solve(prob2)

    # Eq. 4.12: λ = 1*3*1 / (1*0.5 + 3*0.5) = 3/2 = 1.5
    @test sol2[compiled.λ_interface][end] ≈ 1.5 rtol = 1e-10
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
    @test sol[compiled.F][end] ≈ -100.0 rtol = 1e-10
    @test sol[compiled.F][end] < 0  # downward: heat flows from hot to cold

    # Equal temperatures: zero flux
    params_eq = Dict(
        compiled.λ_interface => 1.0,
        compiled.T_i => 300.0, compiled.T_ip1 => 300.0,
        compiled.z_i => 0.05, compiled.z_ip1 => 0.15,
    )
    prob_eq = ODEProblem(compiled, params_eq, (0.0, 1.0))
    sol_eq = solve(prob_eq)
    @test sol_eq[compiled.F][end] ≈ 0.0 atol = 1e-12
end

# ========================================================================
# Tridiagonal Coefficients Tests
# ========================================================================

@testitem "TridiagonalCoefficients Interior - Equation Verification" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = TridiagonalCoefficients(; layer_type = :interior)
    compiled = mtkcompile(sys)

    α = 0.5
    c_i = 2e6
    Δz = 0.02
    Δt = 3600.0
    T_n = 290.0
    λ_above = 1.5
    λ_below = 1.5
    z_im1 = 0.01
    z_i = 0.03
    z_ip1 = 0.05
    F_i_n = -λ_below * (T_n - 285.0) / (z_ip1 - z_i)
    F_im1_n = -λ_above * (295.0 - T_n) / (z_i - z_im1)

    params = Dict(
        compiled.c_i => c_i, compiled.Δz_i => Δz, compiled.Δt => Δt,
        compiled.T_i_n => T_n, compiled.λ_h_above => λ_above,
        compiled.λ_h_below => λ_below, compiled.z_im1 => z_im1,
        compiled.z_i => z_i, compiled.z_ip1 => z_ip1,
        compiled.F_i_n => F_i_n, compiled.F_im1_n => F_im1_n,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.47
    expected_a = -(1 - α) * (Δt / (c_i * Δz)) * λ_above / (z_i - z_im1)
    @test sol[compiled.a_coeff][end] ≈ expected_a rtol = 1e-10

    # Eq. 4.48
    expected_b = 1.0 + (1 - α) * (Δt / (c_i * Δz)) * (λ_above / (z_i - z_im1) + λ_below / (z_ip1 - z_i))
    @test sol[compiled.b_coeff][end] ≈ expected_b rtol = 1e-10

    # Eq. 4.49
    expected_c = -(1 - α) * (Δt / (c_i * Δz)) * λ_below / (z_ip1 - z_i)
    @test sol[compiled.c_coeff][end] ≈ expected_c rtol = 1e-10

    # Eq. 4.50
    expected_r = T_n + α * (Δt / (c_i * Δz)) * (F_i_n - F_im1_n)
    @test sol[compiled.r_coeff][end] ≈ expected_r rtol = 1e-10

    # Diagonal dominance: |b| ≥ |a| + |c|
    @test abs(sol[compiled.b_coeff][end]) ≥ abs(sol[compiled.a_coeff][end]) + abs(sol[compiled.c_coeff][end]) - 1e-10
end

@testitem "TridiagonalCoefficients Top" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = TridiagonalCoefficients(; layer_type = :top)
    compiled = mtkcompile(sys)

    α = 0.5
    c_i = 2e6
    Δz = 0.01
    Δt = 3600.0
    T_n = 295.0
    λ_below = 1.0
    z_i = 0.005
    z_ip1 = 0.015
    dh_dT_val = -15.0
    h_n_val = 50.0
    F_i_n = -20.0

    params = Dict(
        compiled.c_i => c_i, compiled.Δz_i => Δz, compiled.Δt => Δt,
        compiled.T_i_n => T_n, compiled.λ_h_below => λ_below,
        compiled.z_i => z_i, compiled.z_ip1 => z_ip1,
        compiled.dh_dT => dh_dT_val, compiled.h_n => h_n_val,
        compiled.F_i_n => F_i_n,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    # Eq. 4.21
    @test sol[compiled.a_coeff][end] ≈ 0.0 atol = 1e-15

    # Eq. 4.22
    expected_b = 1.0 + (Δt / (c_i * Δz)) * ((1 - α) * λ_below / (z_ip1 - z_i) - dh_dT_val)
    @test sol[compiled.b_coeff][end] ≈ expected_b rtol = 1e-10

    # Eq. 4.23
    expected_c = -(1 - α) * (Δt / (c_i * Δz)) * λ_below / (z_ip1 - z_i)
    @test sol[compiled.c_coeff][end] ≈ expected_c rtol = 1e-10

    # Eq. 4.24
    expected_r = T_n + (Δt / (c_i * Δz)) * (h_n_val - dh_dT_val * T_n + α * F_i_n)
    @test sol[compiled.r_coeff][end] ≈ expected_r rtol = 1e-10
end

@testitem "TridiagonalCoefficients BottomZeroFlux" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = TridiagonalCoefficients(; layer_type = :bottom_zero_flux)
    compiled = mtkcompile(sys)

    α = 0.5
    c_i = 2e6
    Δz = 0.5
    Δt = 3600.0
    T_n = 285.0
    λ_above = 2.0
    z_im1 = 2.0
    z_i = 3.0
    F_im1_n = -5.0

    params = Dict(
        compiled.c_i => c_i, compiled.Δz_i => Δz, compiled.Δt => Δt,
        compiled.T_i_n => T_n, compiled.λ_h_above => λ_above,
        compiled.z_im1 => z_im1, compiled.z_i => z_i,
        compiled.F_im1_n => F_im1_n,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.c_coeff][end] ≈ 0.0 atol = 1e-15
    expected_a = -(1 - α) * (Δt / (c_i * Δz)) * λ_above / (z_i - z_im1)
    @test sol[compiled.a_coeff][end] ≈ expected_a rtol = 1e-10
    expected_b = 1.0 + (1 - α) * (Δt / (c_i * Δz)) * λ_above / (z_i - z_im1)
    @test sol[compiled.b_coeff][end] ≈ expected_b rtol = 1e-10
end

@testitem "TridiagonalCoefficients BottomBuilding" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = TridiagonalCoefficients(; layer_type = :bottom_building)
    compiled = mtkcompile(sys)

    α = 0.5
    c_i = 2e6
    Δz = 0.02
    Δt = 3600.0
    T_n = 292.0
    λ_above = 1.0
    λ_below_val = 1.0
    z_im1 = 0.27
    z_i = 0.29
    z_h_below_val = 0.30
    F_i_n = -10.0
    F_im1_n = -5.0

    params = Dict(
        compiled.c_i => c_i, compiled.Δz_i => Δz, compiled.Δt => Δt,
        compiled.T_i_n => T_n, compiled.λ_h_above => λ_above,
        compiled.λ_h_below => λ_below_val, compiled.z_im1 => z_im1,
        compiled.z_i => z_i, compiled.z_h_below => z_h_below_val,
        compiled.F_i_n => F_i_n, compiled.F_im1_n => F_im1_n,
    )

    prob = ODEProblem(compiled, params, (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.c_coeff][end] ≈ 0.0 atol = 1e-15
    expected_a = -(1 - α) * (Δt / (c_i * Δz)) * λ_above / (z_i - z_im1)
    @test sol[compiled.a_coeff][end] ≈ expected_a rtol = 1e-10
    expected_b = 1.0 + (1 - α) * (Δt / (c_i * Δz)) * (λ_above / (z_i - z_im1) + λ_below_val / (z_h_below_val - z_i))
    @test sol[compiled.b_coeff][end] ≈ expected_b rtol = 1e-10
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
    @test sol[compiled.h][end] ≈ expected_h rtol = 1e-10

    # Eq. 4.29
    σ = 5.67e-8
    expected_dLdT = 4 * 0.95 * σ * 300.0^3
    @test sol[compiled.dL_g_dT][end] ≈ expected_dLdT rtol = 1e-6

    # Eq. 4.28
    expected_dhdT = -expected_dLdT - 5.0 - 2.0
    @test sol[compiled.dh_dT][end] ≈ expected_dhdT rtol = 1e-6
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
    @test sol[compiled.L_roof][end] ≈ expected_L rtol = 1e-10

    # Eq. 4.37
    expected_T = (H * (T_shd + T_sun) + expected_L * T_roof_val) / (2 * H + expected_L)
    @test sol[compiled.T_iB_unclamped][end] ≈ expected_T rtol = 1e-10

    # T_iB should be between min and max of input temperatures
    @test sol[compiled.T_iB_unclamped][end] ≥ min(T_shd, T_sun, T_roof_val) - 1e-10
    @test sol[compiled.T_iB_unclamped][end] ≤ max(T_shd, T_sun, T_roof_val) + 1e-10
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
    @test sol[compiled.T_iB_unclamped][end] ≈ T_uniform rtol = 1e-10
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
    @test sol[compiled.H_wasteheat_unclamped][end] ≈ expected_unclamped rtol = 1e-6
    @test sol[compiled.H_wasteheat][end] ≈ min(expected_unclamped, 100.0) rtol = 1e-6
    @test sol[compiled.H_aircond][end] ≈ 5.0 + 3.0 + 2.0 rtol = 1e-10
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
    @test sol[compiled.H_wasteheat][end] ≈ 100.0 rtol = 1e-6
end

# ========================================================================
# Phase Change Energy Tests
# ========================================================================

@testitem "PhaseChangeEnergy Interior" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeEnergy(; layer_type = :interior)
    compiled = mtkcompile(sys)

    T_f = 273.15
    α = 0.5
    c_i = 2e6
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
    @test sol[compiled.H_excess][end] ≈ expected rtol = 1e-10
end

@testitem "PhaseChangeEnergy Top" setup = [TempSetup] tags = [:ch4_temps] begin
    sys = PhaseChangeEnergy(; layer_type = :top)
    compiled = mtkcompile(sys)

    T_f = 273.15
    α = 0.5
    c_i = 2e6
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
    @test sol[compiled.H_excess][end] ≈ expected rtol = 1e-10
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

    @test sol[compiled.λ_surf][end] ≈ 1.5 rtol = 1e-10
    @test sol[compiled.c_surf][end] ≈ 2.5e6 rtol = 1e-10
end
