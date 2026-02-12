@testsnippet OfflineModeSetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy
end

@testitem "Structural Verification" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()

    # Verify equation count: 18 equations total
    @test length(equations(sys)) == 18

    # Verify unknown count: 18 unknowns
    @test length(unknowns(sys)) == 18

    # Verify key unknown names
    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    for name in [
            "R_vis_raw", "R_nir_raw", "R_vis", "R_nir",
            "S_atm_dir_vis", "S_atm_dir_nir", "S_atm_dif_vis", "S_atm_dif_nir",
            "e_atm", "L_atm_down",
            "f_P_raw", "f_P", "q_rain", "q_snow",
            "u_atm", "v_atm", "θ_atm", "z_prime_atm",
        ]
        @test name in unk_names
    end

    # Verify key parameter names
    param_names =
        Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    for name in ["S_atm", "T_atm", "P_atm", "q_atm", "W_atm", "P_precip"]
        @test name in param_names
    end

    # Verify system compiles
    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "Unit Verification" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()

    # Check units of unknowns
    expected_units = Dict(
        "S_atm_dir_vis" => u"W/m^2",
        "S_atm_dir_nir" => u"W/m^2",
        "S_atm_dif_vis" => u"W/m^2",
        "S_atm_dif_nir" => u"W/m^2",
        "e_atm" => u"Pa",
        "L_atm_down" => u"W/m^2",
        "q_rain" => u"kg/(m^2*s)",
        "q_snow" => u"kg/(m^2*s)",
        "u_atm" => u"m/s",
        "v_atm" => u"m/s",
        "θ_atm" => u"K",
        "z_prime_atm" => u"m",
    )

    for v in unknowns(sys)
        name = string(Symbolics.tosymbol(v, escape = false))
        if haskey(expected_units, name)
            @test ModelingToolkit.get_unit(v) == expected_units[name]
        end
    end

    # Check units of parameters
    param_units = Dict(
        "S_atm" => u"W/m^2",
        "T_atm" => u"K",
        "P_atm" => u"Pa",
        "q_atm" => u"kg/kg",
        "W_atm" => u"m/s",
        "P_precip" => u"kg/(m^2*s)",
    )

    for p in parameters(sys)
        name = string(Symbolics.tosymbol(p, escape = false))
        if haskey(param_units, name)
            @test ModelingToolkit.get_unit(p) == param_units[name]
        end
    end
end

@testitem "Eq. 6.7-6.8 - Direct Fraction Polynomials" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    # Standard conditions: S_atm = 500 W/m², T_atm = 293.15 K, P_atm = 101325 Pa
    # α = 0.5, so α*S_atm = 250 W/m² for visible, (1-α)*S_atm = 250 W/m² for NIR
    # R_vis = 0.17639 + 0.00380*250 + (-9.0039e-6)*250^2 + 8.1351e-9*250^3
    #       = 0.17639 + 0.95 + (-0.562744) + 0.127111 = 0.691257
    α = 0.5
    S_val = 500.0
    vis_flux = α * S_val  # 250
    nir_flux = (1 - α) * S_val  # 250

    R_vis_expected = 0.17639 + 0.0038 * vis_flux +
        (-9.0039e-6) * vis_flux^2 + 8.1351e-9 * vis_flux^3
    R_nir_expected = 0.29548 + 0.00504 * nir_flux +
        (-1.4957e-5) * nir_flux^2 + 1.4881e-8 * nir_flux^3

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

    @test sol[compiled.R_vis][end] ≈ R_vis_expected rtol = 1.0e-10
    @test sol[compiled.R_nir][end] ≈ R_nir_expected rtol = 1.0e-10

    # Verify R values are within physical bounds
    @test 0.01 ≤ sol[compiled.R_vis][end] ≤ 0.99
    @test 0.01 ≤ sol[compiled.R_nir][end] ≤ 0.99
end

@testitem "Eq. 6.7-6.8 - Direct Fraction Clamping" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    # At zero solar radiation, the polynomial gives R_vis = a_0 = 0.17639
    # which is within bounds; test clamping by using very high S_atm
    # where the polynomial may go outside [0.01, 0.99]
    prob_zero = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 0.0,
            compiled.T_atm => 293.15,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.01,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol_zero = solve(prob_zero)

    # At S_atm = 0: R_vis_raw = a_0 = 0.17639, R_nir_raw = b_0 = 0.29548
    @test sol_zero[compiled.R_vis][end] ≈ 0.17639 rtol = 1.0e-10
    @test sol_zero[compiled.R_nir][end] ≈ 0.29548 rtol = 1.0e-10

    # Very high solar radiation - test that clamping works
    prob_high = remake(prob_zero; p = [compiled.S_atm => 2000.0])
    sol_high = solve(prob_high)
    @test 0.01 ≤ sol_high[compiled.R_vis][end] ≤ 0.99
    @test 0.01 ≤ sol_high[compiled.R_nir][end] ≤ 0.99
end

@testitem "Eqs. 6.2-6.5 - Solar Radiation Partitioning" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    S_val = 800.0  # W/m²

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

    dir_vis = sol[compiled.S_atm_dir_vis][end]
    dir_nir = sol[compiled.S_atm_dir_nir][end]
    dif_vis = sol[compiled.S_atm_dif_vis][end]
    dif_nir = sol[compiled.S_atm_dif_nir][end]

    # Energy conservation: sum of all 4 components should equal S_atm
    @test dir_vis + dir_nir + dif_vis + dif_nir ≈ S_val rtol = 1.0e-10

    # Visible components should sum to α * S_atm = 0.5 * S_atm
    @test dir_vis + dif_vis ≈ 0.5 * S_val rtol = 1.0e-10

    # NIR components should sum to (1-α) * S_atm = 0.5 * S_atm
    @test dir_nir + dif_nir ≈ 0.5 * S_val rtol = 1.0e-10

    # All components should be non-negative
    @test dir_vis ≥ 0
    @test dir_nir ≥ 0
    @test dif_vis ≥ 0
    @test dif_nir ≥ 0
end

@testitem "Eqs. 6.2-6.5 - Zero Solar Radiation" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    # When S_atm = 0, all radiation components should be 0
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 0.0,
            compiled.T_atm => 293.15,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.01,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    @test sol[compiled.S_atm_dir_vis][end] ≈ 0.0 atol = 1.0e-15
    @test sol[compiled.S_atm_dir_nir][end] ≈ 0.0 atol = 1.0e-15
    @test sol[compiled.S_atm_dif_vis][end] ≈ 0.0 atol = 1.0e-15
    @test sol[compiled.S_atm_dif_nir][end] ≈ 0.0 atol = 1.0e-15
end

@testitem "Eq. 6.10 - Vapor Pressure" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    # e_atm = P_atm * q_atm / (0.622 + 0.378 * q_atm)
    # At q_atm = 0.01, P_atm = 101325:
    # e_atm = 101325 * 0.01 / (0.622 + 0.378 * 0.01)
    #       = 1013.25 / 0.62578 = 1618.86... (using exact MW ratios)
    q = 0.01
    P = 101325.0
    ε = 0.018016 / 0.028966  # exact ratio
    e_expected = P * q / (ε + (1 - ε) * q)

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 500.0,
            compiled.T_atm => 293.15,
            compiled.P_atm => P,
            compiled.q_atm => q,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    @test sol[compiled.e_atm][end] ≈ e_expected rtol = 1.0e-10

    # Dry air: q_atm = 0 → e_atm = 0
    prob_dry = remake(prob; p = [compiled.q_atm => 0.0])
    sol_dry = solve(prob_dry)
    @test sol_dry[compiled.e_atm][end] ≈ 0.0 atol = 1.0e-15
end

@testitem "Eq. 6.9 - Longwave Radiation" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    # L_atm↓ = [0.70 + 5.95e-5 * 0.01 * e_atm * exp(1500/T_atm)] * σ * T_atm^4
    σ = 5.67e-8
    T = 293.15
    q = 0.01
    P = 101325.0
    ε = 0.018016 / 0.028966
    e = P * q / (ε + (1 - ε) * q)

    L_expected = (0.7 + 5.95e-5 * 0.01 * e * exp(1500.0 / T)) * σ * T^4

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 500.0,
            compiled.T_atm => T,
            compiled.P_atm => P,
            compiled.q_atm => q,
            compiled.W_atm => 5.0,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    @test sol[compiled.L_atm_down][end] ≈ L_expected rtol = 1.0e-8

    # Longwave radiation should always be positive
    @test sol[compiled.L_atm_down][end] > 0

    # Longwave should increase with temperature (Stefan-Boltzmann)
    prob_hot = remake(prob; p = [compiled.T_atm => 313.15])
    sol_hot = solve(prob_hot)
    @test sol_hot[compiled.L_atm_down][end] > sol[compiled.L_atm_down][end]

    # Longwave should increase with humidity (more water vapor → more emission)
    prob_moist = remake(prob; p = [compiled.q_atm => 0.02])
    sol_moist = solve(prob_moist)
    @test sol_moist[compiled.L_atm_down][end] > sol[compiled.L_atm_down][end]
end

@testitem "Eqs. 6.11-6.13 - Precipitation Partitioning" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    P_val = 1.0e-3  # 1 mm/s precipitation rate

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 0.0,
            compiled.T_atm => 274.15,  # 1 K above freezing → f_P = 0.5
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.005,
            compiled.W_atm => 3.0,
            compiled.P_precip => P_val,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    # Conservation: q_rain + q_snow = P_precip
    @test sol[compiled.q_rain][end] + sol[compiled.q_snow][end] ≈ P_val rtol = 1.0e-10

    # At T = T_f + 1 K: f_P = 0.5 * 1 = 0.5
    @test sol[compiled.f_P][end] ≈ 0.5 rtol = 1.0e-10
    @test sol[compiled.q_rain][end] ≈ 0.5 * P_val rtol = 1.0e-10
    @test sol[compiled.q_snow][end] ≈ 0.5 * P_val rtol = 1.0e-10
end

@testitem "Eq. 6.13 - Rain Fraction Limiting Cases" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    P_val = 1.0e-3
    base_params = Dict(
        compiled.S_atm => 0.0,
        compiled.P_atm => 101325.0,
        compiled.q_atm => 0.005,
        compiled.W_atm => 3.0,
        compiled.P_precip => P_val,
    )

    # Well below freezing: T_atm = 263.15 K (T_f - 10) → f_P = clamp(0.5*(-10), 0, 1) = 0 → all snow
    prob_cold = ODEProblem(compiled, merge(base_params, Dict(compiled.T_atm => 263.15)), (0.0, 1.0))
    sol_cold = solve(prob_cold)
    @test sol_cold[compiled.f_P][end] ≈ 0.0 atol = 1.0e-15
    @test sol_cold[compiled.q_snow][end] ≈ P_val rtol = 1.0e-10
    @test sol_cold[compiled.q_rain][end] ≈ 0.0 atol = 1.0e-15

    # At freezing: T_atm = T_f = 273.15 K → f_P = clamp(0, 0, 1) = 0 → all snow
    prob_freeze = ODEProblem(compiled, merge(base_params, Dict(compiled.T_atm => 273.15)), (0.0, 1.0))
    sol_freeze = solve(prob_freeze)
    @test sol_freeze[compiled.f_P][end] ≈ 0.0 atol = 1.0e-10
    @test sol_freeze[compiled.q_snow][end] ≈ P_val rtol = 1.0e-10

    # Well above freezing: T_atm = 283.15 K (T_f + 10) → f_P = clamp(0.5*10, 0, 1) = 1 → all rain
    prob_warm = ODEProblem(compiled, merge(base_params, Dict(compiled.T_atm => 283.15)), (0.0, 1.0))
    sol_warm = solve(prob_warm)
    @test sol_warm[compiled.f_P][end] ≈ 1.0 atol = 1.0e-10
    @test sol_warm[compiled.q_rain][end] ≈ P_val rtol = 1.0e-10
    @test sol_warm[compiled.q_snow][end] ≈ 0.0 atol = 1.0e-15

    # At T_f + 2 K → f_P = clamp(0.5*2, 0, 1) = 1 → just at the transition
    prob_trans = ODEProblem(compiled, merge(base_params, Dict(compiled.T_atm => 275.15)), (0.0, 1.0))
    sol_trans = solve(prob_trans)
    @test sol_trans[compiled.f_P][end] ≈ 1.0 atol = 1.0e-10
end

@testitem "Wind Decomposition and Derived Variables" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    W_val = 10.0
    T_val = 300.0

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.S_atm => 500.0,
            compiled.T_atm => T_val,
            compiled.P_atm => 101325.0,
            compiled.q_atm => 0.01,
            compiled.W_atm => W_val,
            compiled.P_precip => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    # u_atm = v_atm = W_atm / √2
    expected_uv = W_val / sqrt(2.0)
    @test sol[compiled.u_atm][end] ≈ expected_uv rtol = 1.0e-10
    @test sol[compiled.v_atm][end] ≈ expected_uv rtol = 1.0e-10

    # u² + v² should equal W²
    @test sol[compiled.u_atm][end]^2 + sol[compiled.v_atm][end]^2 ≈ W_val^2 rtol = 1.0e-10

    # θ_atm = T_atm
    @test sol[compiled.θ_atm][end] ≈ T_val rtol = 1.0e-10

    # z_prime_atm = 30 m
    @test sol[compiled.z_prime_atm][end] ≈ 30.0 rtol = 1.0e-10
end

@testitem "Qualitative Properties" setup = [OfflineModeSetup] tags = [:offline_mode] begin
    sys = OfflineModeForcing()
    compiled = mtkcompile(sys)

    base_params = Dict(
        compiled.S_atm => 500.0,
        compiled.T_atm => 293.15,
        compiled.P_atm => 101325.0,
        compiled.q_atm => 0.01,
        compiled.W_atm => 5.0,
        compiled.P_precip => 1.0e-4,
    )

    prob = ODEProblem(compiled, base_params, (0.0, 1.0))
    sol = solve(prob)

    # All radiation components should be non-negative
    @test sol[compiled.S_atm_dir_vis][end] ≥ 0
    @test sol[compiled.S_atm_dir_nir][end] ≥ 0
    @test sol[compiled.S_atm_dif_vis][end] ≥ 0
    @test sol[compiled.S_atm_dif_nir][end] ≥ 0
    @test sol[compiled.L_atm_down][end] > 0

    # Precipitation rates should be non-negative
    @test sol[compiled.q_rain][end] ≥ 0
    @test sol[compiled.q_snow][end] ≥ 0

    # Wind components should be non-negative (since W_atm > 0)
    @test sol[compiled.u_atm][end] > 0
    @test sol[compiled.v_atm][end] > 0

    # Higher solar radiation should give more direct beam radiation
    prob_high_S = remake(prob; p = [compiled.S_atm => 900.0])
    sol_high_S = solve(prob_high_S)
    @test sol_high_S[compiled.S_atm_dir_vis][end] > sol[compiled.S_atm_dir_vis][end]
    @test sol_high_S[compiled.S_atm_dir_nir][end] > sol[compiled.S_atm_dir_nir][end]

    # Direct fraction should generally increase with solar intensity (more clear-sky)
    # This is a qualitative property of the polynomial fits
    @test sol_high_S[compiled.R_vis][end] > sol[compiled.R_vis][end] ||
        sol_high_S[compiled.R_vis][end] ≈ sol[compiled.R_vis][end]
end
