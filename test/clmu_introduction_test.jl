@testsnippet CLMUIntroSetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy
end

@testitem "Structural Verification" setup = [CLMUIntroSetup] tags = [:clmu_intro] begin
    sys = CLMUAtmosphere()

    # Verify equation count
    @test length(equations(sys)) == 3

    # Verify unknown count (e_atm, ρ_atm, z_atm)
    @test length(unknowns(sys)) == 3

    # Verify unknown names
    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "e_atm" in unk_names
    @test "ρ_atm" in unk_names
    @test "z_atm" in unk_names

    # Verify key parameter names exist
    param_names =
        Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "P_atm" in param_names
    @test "T_atm" in param_names
    @test "q_atm" in param_names
    @test "z_prime_atm" in param_names
    @test "z_0" in param_names
    @test "z_d" in param_names

    # Verify system compiles
    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "Unit Verification" setup = [CLMUIntroSetup] tags = [:clmu_intro] begin
    sys = CLMUAtmosphere()

    # Verify units of unknowns
    for v in unknowns(sys)
        name = string(Symbolics.tosymbol(v, escape = false))
        u = ModelingToolkit.get_unit(v)
        if name == "e_atm"
            @test u == u"Pa"
        elseif name == "ρ_atm"
            @test u == u"kg/m^3"
        elseif name == "z_atm"
            @test u == u"m"
        end
    end

    # Verify units of key parameters
    for p in parameters(sys)
        name = string(Symbolics.tosymbol(p, escape = false))
        u = ModelingToolkit.get_unit(p)
        if name == "P_atm"
            @test u == u"Pa"
        elseif name == "T_atm"
            @test u == u"K"
        elseif name == "q_atm"
            @test u == u"kg/kg"
        elseif name == "z_prime_atm"
            @test u == u"m"
        end
    end
end

@testitem "Equation Verification - Vapor Pressure" setup = [CLMUIntroSetup] tags = [:clmu_intro] begin
    sys = CLMUAtmosphere()
    compiled = mtkcompile(sys)

    # Test vapor pressure equation: e_atm = q_atm * P_atm / (0.622 + 0.378 * q_atm)
    # At q_atm = 0.01 kg/kg, P_atm = 101325 Pa:
    q_test = 0.01
    P_test = 101325.0
    ε = 0.018016 / 0.028966  # MW_wv / MW_da ≈ 0.62197
    e_expected = q_test * P_test / (ε + (1 - ε) * q_test)

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.q_atm => q_test,
            compiled.P_atm => P_test,
            compiled.T_atm => 293.15,
            compiled.z_prime_atm => 30.0,
            compiled.z_0 => 1.0,
            compiled.z_d => 5.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    @test sol[compiled.e_atm][end] ≈ e_expected rtol = 1.0e-10

    # Additional check: dry air (q_atm = 0) should give e_atm = 0
    prob_dry = remake(prob; p = [compiled.q_atm => 0.0])
    sol_dry = solve(prob_dry)
    @test sol_dry[compiled.e_atm][end] ≈ 0.0 atol = 1.0e-15
end

@testitem "Equation Verification - Air Density" setup = [CLMUIntroSetup] tags = [:clmu_intro] begin
    sys = CLMUAtmosphere()
    compiled = mtkcompile(sys)

    # Test air density: ρ_atm = (P_atm - 0.378 * e_atm) / (R_da * T_atm)
    # For dry air (q=0, e=0) at standard conditions:
    # ρ = P / (R_da * T) = 101325 / (287.042 * 293.15) ≈ 1.204 kg/m³
    R_gas = 6.02214e23 * 1.38065e-23
    R_da = R_gas / 0.028966
    T_test = 293.15
    P_test = 101325.0
    ρ_expected_dry = P_test / (R_da * T_test)

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.q_atm => 0.0,
            compiled.P_atm => P_test,
            compiled.T_atm => T_test,
            compiled.z_prime_atm => 30.0,
            compiled.z_0 => 1.0,
            compiled.z_d => 5.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    @test sol[compiled.ρ_atm][end] ≈ ρ_expected_dry rtol = 1.0e-6

    # Test with moist air (q_atm = 0.01)
    q_test = 0.01
    ε = 0.018016 / 0.028966
    e_test = q_test * P_test / (ε + (1 - ε) * q_test)
    ρ_expected_moist = (P_test - (1 - ε) * e_test) / (R_da * T_test)

    prob_moist = remake(prob; p = [compiled.q_atm => q_test])
    sol_moist = solve(prob_moist)
    @test sol_moist[compiled.ρ_atm][end] ≈ ρ_expected_moist rtol = 1.0e-6

    # Moist air should be less dense than dry air at same T and P
    @test sol_moist[compiled.ρ_atm][end] < sol[compiled.ρ_atm][end]
end

@testitem "Equation Verification - Reference Height" setup = [CLMUIntroSetup] tags = [:clmu_intro] begin
    sys = CLMUAtmosphere()
    compiled = mtkcompile(sys)

    # z_atm = z'_atm + z_0 + z_d (Footnote 1, p.16)
    z_prime = 30.0
    z_0_val = 1.0
    z_d_val = 5.0

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.q_atm => 0.01,
            compiled.P_atm => 101325.0,
            compiled.T_atm => 293.15,
            compiled.z_prime_atm => z_prime,
            compiled.z_0 => z_0_val,
            compiled.z_d => z_d_val,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)
    @test sol[compiled.z_atm][end] ≈ z_prime + z_0_val + z_d_val rtol = 1.0e-10
end

@testitem "Qualitative Properties" setup = [CLMUIntroSetup] tags = [:clmu_intro] begin
    sys = CLMUAtmosphere()
    compiled = mtkcompile(sys)

    base_params = Dict(
        compiled.q_atm => 0.01,
        compiled.P_atm => 101325.0,
        compiled.T_atm => 293.15,
        compiled.z_prime_atm => 30.0,
        compiled.z_0 => 1.0,
        compiled.z_d => 5.0,
    )

    prob = ODEProblem(compiled, base_params, (0.0, 1.0))
    sol = solve(prob)

    # All outputs should be positive
    @test sol[compiled.e_atm][end] > 0
    @test sol[compiled.ρ_atm][end] > 0
    @test sol[compiled.z_atm][end] > 0

    # Higher pressure at same temperature should give higher density
    prob_high_P = remake(prob; p = [compiled.P_atm => 110000.0])
    sol_high_P = solve(prob_high_P)
    @test sol_high_P[compiled.ρ_atm][end] > sol[compiled.ρ_atm][end]

    # Higher temperature at same pressure should give lower density
    prob_high_T = remake(prob; p = [compiled.T_atm => 313.15])
    sol_high_T = solve(prob_high_T)
    @test sol_high_T[compiled.ρ_atm][end] < sol[compiled.ρ_atm][end]

    # Vapor pressure should increase with specific humidity
    prob_high_q = remake(prob; p = [compiled.q_atm => 0.02])
    sol_high_q = solve(prob_high_q)
    @test sol_high_q[compiled.e_atm][end] > sol[compiled.e_atm][end]
end
