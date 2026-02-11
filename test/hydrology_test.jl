@testsnippet HydrologySetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy
end

# ============================================================================
# Structural Verification Tests
# ============================================================================

@testitem "SnowDensity - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowDensity()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "ρ_sno" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "T_atm" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowIceContent - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowIceContent()

    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "dz_sno_dt" in unk_names
    @test "q_ice_top" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "q_grnd_ice" in param_names
    @test "q_frost" in param_names
    @test "q_subl" in param_names
    @test "ρ_sno" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowWaterContent - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowWaterContent()

    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "θ_ice" in unk_names
    @test "θ_liq" in unk_names
    @test "q_liq_out" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "w_liq" in param_names
    @test "w_ice" in param_names
    @test "Δz" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowCompaction - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowCompaction()

    @test length(equations(sys)) == 7
    @test length(unknowns(sys)) == 7

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "C_R1" in unk_names
    @test "C_R2" in unk_names
    @test "C_R3" in unk_names
    @test "C_R" in unk_names
    @test "η" in unk_names
    @test "c_1" in unk_names
    @test "c_2" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "T_i" in param_names
    @test "w_ice" in param_names
    @test "w_liq" in param_names
    @test "Δz" in param_names
    @test "P_s" in param_names
    @test "f_ice_n" in param_names
    @test "f_ice_n1" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowLayerCombination - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowLayerCombination()

    @test length(equations(sys)) == 6
    @test length(unknowns(sys)) == 6

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "Δz_c" in unk_names
    @test "w_liq_c" in unk_names
    @test "w_ice_c" in unk_names
    @test "T_c" in unk_names
    @test "h_1" in unk_names
    @test "h_2" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "Δz_1" in param_names
    @test "Δz_2" in param_names
    @test "w_liq_1" in param_names
    @test "w_liq_2" in param_names
    @test "w_ice_1" in param_names
    @test "w_ice_2" in param_names
    @test "T_1" in param_names
    @test "T_2" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SoilHydraulicProperties - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SoilHydraulicProperties()

    @test length(equations(sys)) == 5
    @test length(unknowns(sys)) == 5

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "k_sat" in unk_names
    @test "θ_sat" in unk_names
    @test "B" in unk_names
    @test "ψ_sat" in unk_names
    @test "ψ" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "pct_sand" in param_names
    @test "pct_clay" in param_names
    @test "θ_i" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SurfaceRunoffInfiltration - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SurfaceRunoffInfiltration()

    @test length(equations(sys)) == 7
    @test length(unknowns(sys)) == 7

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "f_frz_1" in unk_names
    @test "f_sat" in unk_names
    @test "s" in unk_names
    @test "v" in unk_names
    @test "q_infl_max" in unk_names
    @test "q_over" in unk_names
    @test "q_infl" in unk_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SoilWaterFlux - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SoilWaterFlux()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "q" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "k_interface" in param_names
    @test "ψ_upper" in param_names
    @test "ψ_lower" in param_names
    @test "ψ_E_upper" in param_names
    @test "ψ_E_lower" in param_names
    @test "z_upper" in param_names
    @test "z_lower" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SoilWaterEquilibrium - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SoilWaterEquilibrium()

    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "θ_E" in unk_names
    @test "ψ_E" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "θ_sat" in param_names
    @test "ψ_sat" in param_names
    @test "B" in param_names
    @test "z_h_upper" in param_names
    @test "z_h_lower" in param_names
    @test "z_v" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "GroundwaterDrainage - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = GroundwaterDrainage()

    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "f_imp" in unk_names
    @test "q_drai" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "z_v" in param_names
    @test "f_ice_weighted" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "WaterTableDepth - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = WaterTableDepth()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "z_v" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "z_h_bottom" in param_names
    @test "W_a" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "AquiferWaterBalance - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = AquiferWaterBalance()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "dW_a" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "q_recharge" in param_names
    @test "q_drai" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SnowCappingRunoff - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowCappingRunoff()

    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "q_snwcp_ice" in unk_names
    @test "q_snwcp_liq" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "q_grnd_ice" in param_names
    @test "q_grnd_liq" in param_names
    @test "q_frost" in param_names
    @test "q_dew" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "SurfaceLayerUpdate - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SurfaceLayerUpdate()

    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "dw_liq_1" in unk_names
    @test "dw_ice_1_frost" in unk_names
    @test "dw_ice_1_subl" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "q_sdew" in param_names
    @test "q_frost" in param_names
    @test "q_subl" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "PerviousRoadWaterBalance - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = PerviousRoadWaterBalance()

    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "water_input" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "q_rain" in param_names
    @test "q_sno" in param_names
    @test "E_prvrd" in param_names
    @test "q_over" in param_names
    @test "q_drai" in param_names
    @test "q_rgwl" in param_names
    @test "q_snwcp_ice" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "ImperviousWaterBalance - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = ImperviousWaterBalance()

    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    @test "water_input" in unk_names
    @test "q_grnd_liq" in unk_names
    @test "q_grnd_ice" in unk_names

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "q_rain" in param_names
    @test "q_sno" in param_names
    @test "E_surface" in param_names
    @test "q_over" in param_names
    @test "q_rgwl" in param_names
    @test "q_snwcp_ice" in param_names

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

# ============================================================================
# Equation Verification Tests
# ============================================================================

@testitem "SnowDensity - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowDensity()
    compiled = mtkcompile(sys)

    # Test warm case: T_atm = T_f + 2 = 275.15 K => rho_sno = 50 + 1.7*(17)^1.5
    T_f = 273.15
    ρ_expected_warm = 50.0 + 1.7 * (17.0)^1.5

    prob = ODEProblem(compiled, [compiled.T_atm => T_f + 2.0], (0.0, 1.0))
    sol = solve(prob)
    @test sol[compiled.ρ_sno][end] ≈ ρ_expected_warm rtol = 1.0e-6

    # Test cold case: T_atm = 250 K (< T_f - 15 = 258.15) => rho_sno = 50
    prob_cold = remake(prob; p = [compiled.T_atm => 250.0])
    sol_cold = solve(prob_cold)
    @test sol_cold[compiled.ρ_sno][end] ≈ 50.0 rtol = 1.0e-6

    # Test mid-range: T_atm = 268.15 K (T_f - 5) => rho_sno = 50 + 1.7*(10)^1.5
    ρ_expected_mid = 50.0 + 1.7 * (10.0)^1.5
    prob_mid = remake(prob; p = [compiled.T_atm => 268.15])
    sol_mid = solve(prob_mid)
    @test sol_mid[compiled.ρ_sno][end] ≈ ρ_expected_mid rtol = 1.0e-6

    # Test at exact lower threshold: T_atm = T_f - 15 = 258.15 K => rho_sno = 50
    # At the boundary T_atm = T_f + T_thresh_low, the condition is T_atm > T_f - 15,
    # so T_atm = T_f - 15 exactly should give rho_sno = 50 (the else branch)
    prob_boundary = remake(prob; p = [compiled.T_atm => 258.15])
    sol_boundary = solve(prob_boundary)
    @test sol_boundary[compiled.ρ_sno][end] ≈ 50.0 rtol = 1.0e-6

    # Test just above T_f + 2: should give ρ_warm
    prob_above = remake(prob; p = [compiled.T_atm => 276.0])
    sol_above = solve(prob_above)
    @test sol_above[compiled.ρ_sno][end] ≈ ρ_expected_warm rtol = 1.0e-6
end

@testitem "SnowIceContent - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowIceContent()
    compiled = mtkcompile(sys)

    # Test dz_sno_dt = q_grnd_ice / ρ_sno
    # For q_grnd_ice = 0.001 kg/(m^2*s), ρ_sno = 100 kg/m^3:
    # dz_sno_dt = 0.001 / 100 = 1e-5 m/s
    prob = ODEProblem(
        compiled,
        [
            compiled.q_grnd_ice => 0.001, compiled.ρ_sno => 100.0,
            compiled.q_frost => 0.0002, compiled.q_subl => 0.0001,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)
    @test sol[compiled.dz_sno_dt][end] ≈ 1.0e-5 rtol = 1.0e-6

    # Test q_ice_top = q_grnd_ice + (q_frost - q_subl)
    # = 0.001 + (0.0002 - 0.0001) = 0.0011
    @test sol[compiled.q_ice_top][end] ≈ 0.0011 rtol = 1.0e-6
end

@testitem "SnowWaterContent - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowWaterContent()
    compiled = mtkcompile(sys)

    # Test with known values:
    # w_ice = 5.0 kg/m^2, w_liq = 2.0 kg/m^2, Δz = 0.1 m
    # θ_ice = w_ice / (Δz * ρ_ice) = 5 / (0.1 * 917) = 0.05453
    # θ_liq = w_liq / (Δz * ρ_liq) = 2 / (0.1 * 1000) = 0.02
    # q_liq_out = max(ρ_liq * (θ_liq - S_r * (1 - θ_ice)) * Δz / 1s, 0)
    #           = max(1000 * (0.02 - 0.033 * (1 - 0.05453)) * 0.1, 0)
    #           = max(1000 * (0.02 - 0.031201) * 0.1, 0)
    #           = max(-1.1201, 0) = 0
    prob = ODEProblem(
        compiled,
        [compiled.w_ice => 5.0, compiled.w_liq => 2.0, compiled.Δz => 0.1],
        (0.0, 1.0)
    )
    sol = solve(prob)

    θ_ice_expected = 5.0 / (0.1 * 917.0)
    θ_liq_expected = 2.0 / (0.1 * 1000.0)
    @test sol[compiled.θ_ice][end] ≈ θ_ice_expected rtol = 1.0e-6
    @test sol[compiled.θ_liq][end] ≈ θ_liq_expected rtol = 1.0e-6
    @test sol[compiled.q_liq_out][end] ≈ 0.0 atol = 1.0e-10

    # Test with higher water content where flow occurs:
    # w_liq = 10.0 kg/m^2, w_ice = 1.0 kg/m^2, Δz = 0.1 m
    # θ_ice = 1 / (0.1 * 917) = 0.01091
    # θ_liq = 10 / (0.1 * 1000) = 0.1
    # q_liq_out = max(1000 * (0.1 - 0.033 * (1 - 0.01091)) * 0.1, 0)
    #           = max(1000 * (0.1 - 0.03264) * 0.1, 0)
    #           = max(6.736, 0) = 6.736
    prob_wet = remake(prob; p = [compiled.w_liq => 10.0, compiled.w_ice => 1.0, compiled.Δz => 0.1])
    sol_wet = solve(prob_wet)

    θ_ice_wet = 1.0 / (0.1 * 917.0)
    θ_liq_wet = 10.0 / (0.1 * 1000.0)
    q_expected = 1000.0 * (θ_liq_wet - 0.033 * (1.0 - θ_ice_wet)) * 0.1 / 1.0
    @test sol_wet[compiled.θ_ice][end] ≈ θ_ice_wet rtol = 1.0e-6
    @test sol_wet[compiled.θ_liq][end] ≈ θ_liq_wet rtol = 1.0e-6
    @test sol_wet[compiled.q_liq_out][end] ≈ q_expected rtol = 1.0e-4
end

@testitem "SoilHydraulicProperties - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SoilHydraulicProperties()
    compiled = mtkcompile(sys)

    # Reference values for %sand = 50, %clay = 20
    # k_sat = 0.0070556e-3 * 10^(-0.884 + 0.0153*50) m/s
    pct_sand = 50.0
    pct_clay = 20.0
    k_sat_expected = 0.0070556e-3 * 10.0^(-0.884 + 0.0153 * pct_sand)
    θ_sat_expected = 0.489 - 0.00126 * pct_sand  # = 0.426
    B_expected = 2.91 + 0.159 * pct_clay  # = 6.09
    ψ_sat_expected = -10.0e-3 * 10.0^(1.88 - 0.0131 * pct_sand)  # ≈ -0.1679 m

    # Use θ_i = θ_sat for initial test (saturated soil)
    prob = ODEProblem(
        compiled,
        [
            compiled.pct_sand => pct_sand, compiled.pct_clay => pct_clay,
            compiled.θ_i => θ_sat_expected,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)

    @test sol[compiled.k_sat][end] ≈ k_sat_expected rtol = 1.0e-6
    @test sol[compiled.θ_sat][end] ≈ θ_sat_expected rtol = 1.0e-6
    @test sol[compiled.B][end] ≈ B_expected rtol = 1.0e-6
    @test sol[compiled.ψ_sat][end] ≈ ψ_sat_expected rtol = 1.0e-6

    # At saturation (θ_i = θ_sat), ψ should equal ψ_sat
    @test sol[compiled.ψ][end] ≈ ψ_sat_expected rtol = 1.0e-4
end

@testitem "WaterTableDepth - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = WaterTableDepth()
    compiled = mtkcompile(sys)

    # For W_a = 4800 kg/m^2, z_h_bottom = 3.5 m:
    # z_v = 3.5 + 25 - 4800 / (1000 * 0.2) = 3.5 + 25 - 24 = 4.5 m
    prob = ODEProblem(
        compiled,
        [compiled.z_h_bottom => 3.5, compiled.W_a => 4800.0],
        (0.0, 1.0)
    )
    sol = solve(prob)
    @test sol[compiled.z_v][end] ≈ 4.5 rtol = 1.0e-6

    # Test clamping at minimum (z_v_min = 0.05 m)
    # Large W_a should push z_v below minimum
    prob_min = remake(prob; p = [compiled.z_h_bottom => 3.5, compiled.W_a => 6000.0])
    sol_min = solve(prob_min)
    # z_v = 3.5 + 25 - 6000/(1000*0.2) = 28.5 - 30 = -1.5 -> clamped to 0.05
    @test sol_min[compiled.z_v][end] ≈ 0.05 rtol = 1.0e-6

    # Test clamping at maximum (z_v_max = 80 m)
    # Very small W_a should push z_v to maximum
    prob_max = remake(prob; p = [compiled.z_h_bottom => 3.5, compiled.W_a => 0.0])
    sol_max = solve(prob_max)
    # z_v = 3.5 + 25 - 0 = 28.5 m (within bounds)
    @test sol_max[compiled.z_v][end] ≈ 28.5 rtol = 1.0e-6
end

@testitem "GroundwaterDrainage - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = GroundwaterDrainage()
    compiled = mtkcompile(sys)

    # For z_v = 1.0 m, f_ice_weighted = 0 (no frozen soil):
    # f_imp = max((exp(-3*(1-0)) - exp(-3)) / (1 - exp(-3)), 0)
    #       = max((exp(-3) - exp(-3)) / (1 - exp(-3)), 0) = 0
    # q_drai = (1-0) * 5.5e-3 * exp(-2.5 * 1.0)
    q_drai_expected = 5.5e-3 * exp(-2.5 * 1.0)

    prob = ODEProblem(
        compiled,
        [compiled.z_v => 1.0, compiled.f_ice_weighted => 0.0],
        (0.0, 1.0)
    )
    sol = solve(prob)

    @test sol[compiled.f_imp][end] ≈ 0.0 atol = 1.0e-10
    @test sol[compiled.q_drai][end] ≈ q_drai_expected rtol = 1.0e-6

    # Test with fully frozen soil: f_ice_weighted = 1.0
    # f_imp = max((exp(-3*(1-1)) - exp(-3)) / (1 - exp(-3)), 0)
    #       = max((1 - exp(-3)) / (1 - exp(-3)), 0) = 1.0
    # q_drai = (1-1) * ... = 0
    prob_frozen = remake(prob; p = [compiled.z_v => 1.0, compiled.f_ice_weighted => 1.0])
    sol_frozen = solve(prob_frozen)
    @test sol_frozen[compiled.f_imp][end] ≈ 1.0 rtol = 1.0e-6
    @test sol_frozen[compiled.q_drai][end] ≈ 0.0 atol = 1.0e-10

    # Test at z_v = 0 m:
    # q_drai = 5.5e-3 * exp(0) = 5.5e-3
    prob_surface = remake(prob; p = [compiled.z_v => 0.0, compiled.f_ice_weighted => 0.0])
    sol_surface = solve(prob_surface)
    @test sol_surface[compiled.q_drai][end] ≈ 5.5e-3 rtol = 1.0e-6
end

@testitem "SnowCappingRunoff - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowCappingRunoff()
    compiled = mtkcompile(sys)

    # q_snwcp_ice = q_grnd_ice + q_frost = 0.001 + 0.0005 = 0.0015
    # q_snwcp_liq = q_grnd_liq + q_dew = 0.002 + 0.0003 = 0.0023
    prob = ODEProblem(
        compiled,
        [
            compiled.q_grnd_ice => 0.001, compiled.q_frost => 0.0005,
            compiled.q_grnd_liq => 0.002, compiled.q_dew => 0.0003,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)

    @test sol[compiled.q_snwcp_ice][end] ≈ 0.0015 rtol = 1.0e-6
    @test sol[compiled.q_snwcp_liq][end] ≈ 0.0023 rtol = 1.0e-6
end

@testitem "SnowLayerCombination - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowLayerCombination()
    compiled = mtkcompile(sys)

    T_f = 273.15
    C_ice = 2117.27
    C_liq = 4188.0
    L_f = 3.337e5

    # Two layers with known properties
    Δz1 = 0.05; Δz2 = 0.1
    w_liq1 = 0.5; w_liq2 = 1.0
    w_ice1 = 3.0; w_ice2 = 5.0
    T1 = 270.0; T2 = 265.0

    prob = ODEProblem(
        compiled,
        [
            compiled.Δz_1 => Δz1, compiled.Δz_2 => Δz2,
            compiled.w_liq_1 => w_liq1, compiled.w_liq_2 => w_liq2,
            compiled.w_ice_1 => w_ice1, compiled.w_ice_2 => w_ice2,
            compiled.T_1 => T1, compiled.T_2 => T2,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)

    # Combined thickness
    @test sol[compiled.Δz_c][end] ≈ Δz1 + Δz2 rtol = 1.0e-6

    # Combined water/ice
    @test sol[compiled.w_liq_c][end] ≈ w_liq1 + w_liq2 rtol = 1.0e-6
    @test sol[compiled.w_ice_c][end] ≈ w_ice1 + w_ice2 rtol = 1.0e-6

    # Enthalpy calculation
    h1_expected = (C_ice * w_ice1 + C_liq * w_liq1) * (T1 - T_f) + L_f * w_liq1
    h2_expected = (C_ice * w_ice2 + C_liq * w_liq2) * (T2 - T_f) + L_f * w_liq2
    @test sol[compiled.h_1][end] ≈ h1_expected rtol = 1.0e-6
    @test sol[compiled.h_2][end] ≈ h2_expected rtol = 1.0e-6

    # Combined temperature
    w_liq_c = w_liq1 + w_liq2
    w_ice_c = w_ice1 + w_ice2
    T_c_expected = T_f + (h1_expected + h2_expected - L_f * w_liq_c) / (C_ice * w_ice_c + C_liq * w_liq_c)
    @test sol[compiled.T_c][end] ≈ T_c_expected rtol = 1.0e-6

    # Verify enthalpy conservation: h_combined should equal h_1 + h_2
    h_c = (C_ice * w_ice_c + C_liq * w_liq_c) * (sol[compiled.T_c][end] - T_f) + L_f * w_liq_c
    @test h_c ≈ h1_expected + h2_expected rtol = 1.0e-6
end

@testitem "SnowCompaction - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowCompaction()
    compiled = mtkcompile(sys)

    T_f = 273.15
    c_3 = 2.777e-6  # s^-1
    c_4 = 0.04  # K^-1
    η_0 = 9.0e5  # kg*s/m^2
    c_5 = 0.08  # K^-1
    c_6 = 0.023  # m^3/kg

    # Test at T = T_f (no temperature offset), low density, no liquid
    # c_1 = 1 (density <= 100 kg/m^3), c_2 = 1 (no liquid)
    # C_R1 = -c_3 * 1 * 1 * exp(-c_4 * 0) = -c_3
    T_i = T_f
    w_ice = 5.0  # kg/m^2
    w_liq = 0.0  # kg/m^2
    Δz = 0.1     # m (w_ice/Δz = 50 kg/m^3 < 100)
    P_s = 10.0   # kg/m^2
    f_ice_n = 1.0
    f_ice_n1 = 1.0

    prob = ODEProblem(
        compiled,
        [
            compiled.T_i => T_i, compiled.w_ice => w_ice, compiled.w_liq => w_liq,
            compiled.Δz => Δz, compiled.P_s => P_s,
            compiled.f_ice_n => f_ice_n, compiled.f_ice_n1 => f_ice_n1,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)

    # C_R1 at T_f with c_1=1, c_2=1: C_R1 = -c_3
    @test sol[compiled.C_R1][end] ≈ -c_3 rtol = 1.0e-6
    @test sol[compiled.c_1][end] ≈ 1.0 rtol = 1.0e-6
    @test sol[compiled.c_2][end] ≈ 1.0 rtol = 1.0e-6

    # η = η_0 * exp(c_5*(T_f-T_f) + c_6*50) = η_0 * exp(c_6 * 50)
    η_expected = η_0 * exp(c_5 * 0.0 + c_6 * (w_ice / Δz))
    @test sol[compiled.η][end] ≈ η_expected rtol = 1.0e-6

    # C_R2 = -P_s / η
    C_R2_expected = -P_s / η_expected
    @test sol[compiled.C_R2][end] ≈ C_R2_expected rtol = 1.0e-6

    # C_R3 = -1 * max(0, (1.0-1.0)/1.0) = 0 (no melting)
    @test sol[compiled.C_R3][end] ≈ 0.0 atol = 1.0e-10

    # Total: C_R = C_R1 + C_R2 + C_R3
    C_R_expected = -c_3 + C_R2_expected + 0.0
    @test sol[compiled.C_R][end] ≈ C_R_expected rtol = 1.0e-6

    # Test high-density case: w_ice/Δz > 100 kg/m^3
    # c_1 = exp(-0.046*(150-100)) = exp(-2.3)
    w_ice_dense = 15.0  # w_ice/Δz = 150 kg/m^3
    prob_dense = remake(prob; p = [compiled.w_ice => w_ice_dense])
    sol_dense = solve(prob_dense)
    c_1_expected = exp(-0.046 * (w_ice_dense / Δz - 100.0))
    @test sol_dense[compiled.c_1][end] ≈ c_1_expected rtol = 1.0e-6

    # Test with liquid water: c_2 = 2 when w_liq/Δz > 0.01
    w_liq_wet = 0.01  # w_liq/Δz = 0.1 > 0.01
    prob_wet = remake(prob; p = [compiled.w_liq => w_liq_wet])
    sol_wet = solve(prob_wet)
    @test sol_wet[compiled.c_2][end] ≈ 2.0 rtol = 1.0e-6
end

@testitem "SoilWaterFlux - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SoilWaterFlux()
    compiled = mtkcompile(sys)

    # q = -k_interface * ((ψ_upper - ψ_lower) + (ψ_E_lower - ψ_E_upper)) / (z_lower - z_upper)
    # With ψ_upper = -0.5 m, ψ_lower = -1.0 m, ψ_E = same, z_upper = 0.1 m, z_lower = 0.3 m
    # q = -1e-5 * ((-0.5 - (-1.0)) + (-0.8 - (-0.3))) / (0.3 - 0.1)
    #   = -1e-5 * (0.5 + (-0.5)) / 0.2
    #   = -1e-5 * 0 / 0.2 = 0
    prob = ODEProblem(
        compiled,
        [
            compiled.k_interface => 1.0e-5, compiled.ψ_upper => -0.5, compiled.ψ_lower => -1.0,
            compiled.ψ_E_upper => -0.3, compiled.ψ_E_lower => -0.8,
            compiled.z_upper => 0.1, compiled.z_lower => 0.3,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)
    @test sol[compiled.q][end] ≈ 0.0 atol = 1.0e-15

    # Non-zero flux test:
    # ψ_upper = -0.2, ψ_lower = -0.5, ψ_E_upper = -0.1, ψ_E_lower = -0.1
    # q = -1e-5 * ((-0.2 - (-0.5)) + (-0.1 - (-0.1))) / (0.3 - 0.1)
    #   = -1e-5 * (0.3 + 0) / 0.2 = -1.5e-5
    prob2 = remake(
        prob; p = [
            compiled.k_interface => 1.0e-5,
            compiled.ψ_upper => -0.2, compiled.ψ_lower => -0.5,
            compiled.ψ_E_upper => -0.1, compiled.ψ_E_lower => -0.1,
            compiled.z_upper => 0.1, compiled.z_lower => 0.3,
        ]
    )
    sol2 = solve(prob2)
    @test sol2[compiled.q][end] ≈ -1.5e-5 rtol = 1.0e-6
end

@testitem "AquiferWaterBalance - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = AquiferWaterBalance()
    compiled = mtkcompile(sys)

    # dW_a = q_recharge - q_drai
    prob = ODEProblem(
        compiled,
        [compiled.q_recharge => 0.003, compiled.q_drai => 0.001],
        (0.0, 1.0)
    )
    sol = solve(prob)
    @test sol[compiled.dW_a][end] ≈ 0.002 rtol = 1.0e-6

    # Zero balance
    prob_zero = remake(prob; p = [compiled.q_recharge => 0.001, compiled.q_drai => 0.001])
    sol_zero = solve(prob_zero)
    @test sol_zero[compiled.dW_a][end] ≈ 0.0 atol = 1.0e-15
end

@testitem "SurfaceLayerUpdate - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SurfaceLayerUpdate()
    compiled = mtkcompile(sys)

    prob = ODEProblem(
        compiled,
        [compiled.q_sdew => 0.0005, compiled.q_frost => 0.0003, compiled.q_subl => 0.0002],
        (0.0, 1.0)
    )
    sol = solve(prob)

    # dw_liq_1 = q_sdew = 0.0005
    @test sol[compiled.dw_liq_1][end] ≈ 0.0005 rtol = 1.0e-6
    # dw_ice_1_frost = q_frost = 0.0003
    @test sol[compiled.dw_ice_1_frost][end] ≈ 0.0003 rtol = 1.0e-6
    # dw_ice_1_subl = -q_subl = -0.0002
    @test sol[compiled.dw_ice_1_subl][end] ≈ -0.0002 rtol = 1.0e-6
end

@testitem "PerviousRoadWaterBalance - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = PerviousRoadWaterBalance()
    compiled = mtkcompile(sys)

    # water_input = q_rain + q_sno - E_prvrd - q_over - q_drai - q_rgwl - q_snwcp_ice
    prob = ODEProblem(
        compiled,
        [
            compiled.q_rain => 0.01, compiled.q_sno => 0.005,
            compiled.E_prvrd => 0.002, compiled.q_over => 0.003,
            compiled.q_drai => 0.001, compiled.q_rgwl => 0.0005,
            compiled.q_snwcp_ice => 0.0001,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)
    expected = 0.01 + 0.005 - 0.002 - 0.003 - 0.001 - 0.0005 - 0.0001
    @test sol[compiled.water_input][end] ≈ expected rtol = 1.0e-6
end

@testitem "ImperviousWaterBalance - Equation Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = ImperviousWaterBalance()
    compiled = mtkcompile(sys)

    prob = ODEProblem(
        compiled,
        [
            compiled.q_rain => 0.01, compiled.q_sno => 0.005,
            compiled.E_surface => 0.002, compiled.q_over => 0.003,
            compiled.q_rgwl => 0.0005, compiled.q_snwcp_ice => 0.0001,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)

    # q_grnd_liq = q_rain
    @test sol[compiled.q_grnd_liq][end] ≈ 0.01 rtol = 1.0e-6
    # q_grnd_ice = q_sno
    @test sol[compiled.q_grnd_ice][end] ≈ 0.005 rtol = 1.0e-6
    # water_input = q_rain + q_sno - E_surface - q_over - q_rgwl - q_snwcp_ice (Eqs. 5.2-5.3)
    expected = 0.01 + 0.005 - 0.002 - 0.003 - 0.0005 - 0.0001
    @test sol[compiled.water_input][end] ≈ expected rtol = 1.0e-6
end

# ============================================================================
# Qualitative Behavior Tests
# ============================================================================

@testitem "SnowDensity - Qualitative Behavior" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowDensity()
    compiled = mtkcompile(sys)

    T_f = 273.15

    # Snow density should increase with temperature between T_f-15 and T_f+2
    prob = ODEProblem(compiled, [compiled.T_atm => 260.0], (0.0, 1.0))
    temperatures = [260.0, 264.0, 268.0, 272.0]
    densities = Float64[]
    for T in temperatures
        p = remake(prob; p = [compiled.T_atm => T])
        s = solve(p)
        push!(densities, s[compiled.ρ_sno][end])
    end
    # Density should be monotonically increasing
    for i in 2:length(densities)
        @test densities[i] > densities[i - 1]
    end

    # All densities should be at least 50 kg/m^3
    for ρ in densities
        @test ρ >= 50.0
    end

    # Very cold temperature should give exactly 50
    prob_cold = remake(prob; p = [compiled.T_atm => 240.0])
    sol_cold = solve(prob_cold)
    @test sol_cold[compiled.ρ_sno][end] ≈ 50.0 rtol = 1.0e-6
end

@testitem "GroundwaterDrainage - Qualitative Behavior" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = GroundwaterDrainage()
    compiled = mtkcompile(sys)

    # Drainage should decrease exponentially with water table depth
    prob = ODEProblem(compiled, [compiled.z_v => 0.5, compiled.f_ice_weighted => 0.0], (0.0, 1.0))
    depths = [0.5, 1.0, 2.0, 3.0, 5.0]
    drainages = Float64[]
    for z in depths
        p = remake(prob; p = [compiled.z_v => z])
        s = solve(p)
        push!(drainages, s[compiled.q_drai][end])
    end
    # Drainage should be monotonically decreasing with depth
    for i in 2:length(drainages)
        @test drainages[i] < drainages[i - 1]
    end
    # All drainage should be positive
    for q in drainages
        @test q > 0.0
    end

    # Higher frozen fraction should reduce drainage
    prob_partial = remake(prob; p = [compiled.z_v => 1.0, compiled.f_ice_weighted => 0.5])
    sol_partial = solve(prob_partial)
    prob_unfrozen = remake(prob; p = [compiled.z_v => 1.0, compiled.f_ice_weighted => 0.0])
    sol_unfrozen = solve(prob_unfrozen)
    @test sol_partial[compiled.q_drai][end] < sol_unfrozen[compiled.q_drai][end]
end

@testitem "SnowCompaction - Qualitative Behavior" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowCompaction()
    compiled = mtkcompile(sys)

    T_f = 273.15

    # All compaction rates should be negative (compression)
    prob = ODEProblem(
        compiled,
        [
            compiled.T_i => 268.0, compiled.w_ice => 5.0, compiled.w_liq => 0.0,
            compiled.Δz => 0.1, compiled.P_s => 10.0,
            compiled.f_ice_n => 1.0, compiled.f_ice_n1 => 1.0,
        ],
        (0.0, 1.0)
    )
    sol = solve(prob)
    @test sol[compiled.C_R1][end] < 0.0
    @test sol[compiled.C_R2][end] < 0.0
    @test sol[compiled.C_R][end] < 0.0

    # Compaction rate magnitude should increase with temperature (closer to T_f)
    # because exp(-c_4*(T_f-T_i)) increases as T_i -> T_f
    temps = [255.0, 260.0, 265.0, 270.0]
    cr1_values = Float64[]
    for T in temps
        p = remake(prob; p = [compiled.T_i => T])
        s = solve(p)
        push!(cr1_values, abs(s[compiled.C_R1][end]))
    end
    # |C_R1| should increase with temperature
    for i in 2:length(cr1_values)
        @test cr1_values[i] > cr1_values[i - 1]
    end

    # Higher overburden pressure -> larger |C_R2|
    prob_low_P = remake(prob; p = [compiled.P_s => 5.0])
    sol_low_P = solve(prob_low_P)
    prob_high_P = remake(prob; p = [compiled.P_s => 20.0])
    sol_high_P = solve(prob_high_P)
    @test abs(sol_high_P[compiled.C_R2][end]) > abs(sol_low_P[compiled.C_R2][end])
end

@testitem "WaterTableDepth - Qualitative Behavior" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = WaterTableDepth()
    compiled = mtkcompile(sys)

    # More aquifer water should give shallower water table
    prob = ODEProblem(compiled, [compiled.z_h_bottom => 3.5, compiled.W_a => 3000.0], (0.0, 1.0))
    wa_values = [2000.0, 3000.0, 4000.0, 5000.0]
    zv_values = Float64[]
    for wa in wa_values
        p = remake(prob; p = [compiled.W_a => wa])
        s = solve(p)
        push!(zv_values, s[compiled.z_v][end])
    end
    # z_v should decrease (shallower) as W_a increases
    for i in 2:length(zv_values)
        @test zv_values[i] <= zv_values[i - 1]
    end
end

@testitem "SnowWaterContent - Qualitative Behavior" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SnowWaterContent()
    compiled = mtkcompile(sys)

    # More liquid water should give more outflow
    prob = ODEProblem(compiled, [compiled.w_liq => 5.0, compiled.w_ice => 1.0, compiled.Δz => 0.1], (0.0, 1.0))

    w_liq_values = [5.0, 10.0, 20.0, 50.0]
    outflows = Float64[]
    for wl in w_liq_values
        p = remake(prob; p = [compiled.w_liq => wl])
        s = solve(p)
        push!(outflows, s[compiled.q_liq_out][end])
    end
    # Outflow should increase (or stay zero) with more liquid
    for i in 2:length(outflows)
        @test outflows[i] >= outflows[i - 1]
    end

    # More ice reduces the irreducible retention S_r*(1-θ_ice), so at the same θ_liq
    # the outflow actually increases with more ice. This is because higher ice content
    # reduces the pore space available for capillary retention (Eq. 5.18).
    prob_lo_ice = remake(prob; p = [compiled.w_liq => 20.0, compiled.w_ice => 1.0])
    sol_lo_ice = solve(prob_lo_ice)
    prob_hi_ice = remake(prob; p = [compiled.w_liq => 20.0, compiled.w_ice => 40.0])
    sol_hi_ice = solve(prob_hi_ice)
    @test sol_hi_ice[compiled.q_liq_out][end] >= sol_lo_ice[compiled.q_liq_out][end]
end

@testitem "SoilHydraulicProperties - Qualitative Behavior" setup = [HydrologySetup] tags = [:hydrology] begin
    sys = SoilHydraulicProperties()
    compiled = mtkcompile(sys)

    # Higher sand content should increase k_sat (more permeable)
    prob = ODEProblem(
        compiled,
        [compiled.pct_sand => 30.0, compiled.pct_clay => 20.0, compiled.θ_i => 0.3],
        (0.0, 1.0)
    )
    sand_values = [20.0, 40.0, 60.0, 80.0]
    ksat_values = Float64[]
    for sand in sand_values
        p = remake(prob; p = [compiled.pct_sand => sand])
        s = solve(p)
        push!(ksat_values, s[compiled.k_sat][end])
    end
    for i in 2:length(ksat_values)
        @test ksat_values[i] > ksat_values[i - 1]
    end

    # Higher sand content should decrease porosity
    θsat_values = Float64[]
    for sand in sand_values
        p = remake(prob; p = [compiled.pct_sand => sand])
        s = solve(p)
        push!(θsat_values, s[compiled.θ_sat][end])
    end
    for i in 2:length(θsat_values)
        @test θsat_values[i] < θsat_values[i - 1]
    end

    # Higher clay content should increase B
    clay_values = [10.0, 20.0, 30.0, 40.0]
    B_values = Float64[]
    for clay in clay_values
        p = remake(prob; p = [compiled.pct_clay => clay])
        s = solve(p)
        push!(B_values, s[compiled.B][end])
    end
    for i in 2:length(B_values)
        @test B_values[i] > B_values[i - 1]
    end

    # ψ_sat should be negative for all sand contents
    for sand in sand_values
        p = remake(prob; p = [compiled.pct_sand => sand])
        s = solve(p)
        @test s[compiled.ψ_sat][end] < 0.0
    end
end

# ============================================================================
# RichardsEquation Tests
# ============================================================================

@testitem "RichardsEquation - Structural Verification" setup = [HydrologySetup] tags = [:hydrology] begin
    using MethodOfLines, DomainSets

    θ_val = 0.3
    result = RichardsEquation(;
        N_layers = 5, Δz_total = 2.0,
        pct_sand = 50.0, pct_clay = 20.0,
        θ_top_val = θ_val, θ_bottom_val = θ_val,
        θ_init_val = θ_val,
    )

    # Verify the return type
    @test haskey(result, :prob)
    @test haskey(result, :θ)
    @test haskey(result, :t_pde)
    @test haskey(result, :z_var)
    @test haskey(result, :k_sat_val)
    @test haskey(result, :θ_sat_val)
    @test haskey(result, :B_val)
    @test haskey(result, :ψ_sat_val)

    # Check computed soil properties
    @test result.k_sat_val ≈ 0.0070556e-3 * 10.0^(-0.884 + 0.0153 * 50.0) rtol = 1.0e-6
    @test result.θ_sat_val ≈ 0.489 - 0.00126 * 50.0 rtol = 1.0e-6
    @test result.B_val ≈ 2.91 + 0.159 * 20.0 rtol = 1.0e-6
    @test result.ψ_sat_val ≈ -10.0e-3 * 10.0^(1.88 - 0.0131 * 50.0) rtol = 1.0e-6

    # Verify the problem can be solved
    sol = solve(result.prob, saveat = 0.1)
    @test sol.retcode == :Success || sol.retcode == ReturnCode.Success
end

@testitem "RichardsEquation - Steady State Uniform" setup = [HydrologySetup] tags = [:hydrology] begin
    using MethodOfLines, DomainSets

    θ_val = 0.3
    # With uniform BCs and uniform initial condition, the profile should remain
    # nearly uniform (the gravity term redistributes slightly, but with identical
    # BCs at top and bottom the steady state should be close to uniform).
    result = RichardsEquation(;
        N_layers = 10, Δz_total = 2.0,
        pct_sand = 50.0, pct_clay = 20.0,
        θ_top_val = θ_val, θ_bottom_val = θ_val,
        θ_init_val = θ_val,
    )

    sol = solve(result.prob, saveat = [0.0, 1.0])
    θ_mat = sol[result.θ(result.t_pde, result.z_var)]

    θ_final = θ_mat[end, :]
    # Profile should remain close to uniform θ_val
    for v in θ_final
        @test v ≈ θ_val rtol = 0.1
    end
end

@testitem "RichardsEquation - Wetting from Top" setup = [HydrologySetup] tags = [:hydrology] begin
    using MethodOfLines, DomainSets

    Δz_total = 2.0
    θ_sat_val = 0.489 - 0.00126 * 50.0
    θ_dry = 0.1 * θ_sat_val
    θ_wet = 0.8 * θ_sat_val

    # Wet top, dry bottom: water should flow from top into the column
    result = RichardsEquation(;
        N_layers = 10, Δz_total = Δz_total,
        pct_sand = 50.0, pct_clay = 20.0,
        θ_top_val = θ_wet,
        θ_bottom_val = θ_dry,
        θ_init_val = θ_dry,
    )

    tspan = (0.0, 50000.0)
    prob2 = remake(result.prob; tspan = tspan)
    sol = solve(prob2, saveat = [0.0, 50000.0])
    θ_mat = sol[result.θ(result.t_pde, result.z_var)]

    θ_initial = θ_mat[1, :]
    θ_final = θ_mat[end, :]

    # Upper portion of profile should be wetter than initial dry condition
    # (water entered from the wet top boundary; wetting front may not have
    # reached the deepest layers within the simulation time)
    n = length(θ_final)
    n_upper = div(n, 3)  # Check top third of the profile
    for i in 2:n_upper
        @test θ_final[i] > θ_dry
    end

    # All interior points should be at least as wet as initial condition
    for i in 2:(n - 1)
        @test θ_final[i] >= θ_dry - 1.0e-10
    end

    # Profile should be monotonically decreasing from top to bottom
    # (wetting front propagates down)
    for i in 1:(n - 1)
        @test θ_final[i] >= θ_final[i + 1] - 1.0e-6
    end
end

@testitem "RichardsEquation - Gravity Drainage" setup = [HydrologySetup] tags = [:hydrology] begin
    using MethodOfLines, DomainSets

    Δz_total = 2.0
    θ_sat_val = 0.489 - 0.00126 * 50.0
    θ_mid = 0.4 * θ_sat_val

    # Uniform BCs with gravity: the gravity term ∂K/∂z drives downward flow
    # At steady state with same BCs, the profile develops a slight bulge
    # because gravity pushes water down while diffusion resists.
    result = RichardsEquation(;
        N_layers = 10, Δz_total = Δz_total,
        pct_sand = 50.0, pct_clay = 20.0,
        θ_top_val = θ_mid,
        θ_bottom_val = θ_mid,
        θ_init_val = θ_mid,
    )

    tspan = (0.0, 100000.0)
    prob2 = remake(result.prob; tspan = tspan)
    sol = solve(prob2, saveat = [0.0, 100000.0])
    θ_mat = sol[result.θ(result.t_pde, result.z_var)]

    θ_final = θ_mat[end, :]
    n = length(θ_final)

    # With equal BCs and gravity, the system should be near the boundary values
    # (gravity redistributes but BCs constrain the endpoints)
    @test θ_final[1] ≈ θ_mid rtol = 0.05
    @test θ_final[end] ≈ θ_mid rtol = 0.05
end

@testitem "RichardsEquation - Different Soil Types" setup = [HydrologySetup] tags = [:hydrology] begin
    using MethodOfLines, DomainSets

    θ_val = 0.25

    # Sandy soil should have higher k_sat and faster dynamics
    result_sand = RichardsEquation(;
        N_layers = 5, Δz_total = 2.0,
        pct_sand = 80.0, pct_clay = 5.0,
        θ_top_val = θ_val, θ_bottom_val = θ_val,
        θ_init_val = θ_val,
    )

    # Clay soil should have lower k_sat
    result_clay = RichardsEquation(;
        N_layers = 5, Δz_total = 2.0,
        pct_sand = 10.0, pct_clay = 50.0,
        θ_top_val = θ_val, θ_bottom_val = θ_val,
        θ_init_val = θ_val,
    )

    # Sandy soil has higher k_sat
    @test result_sand.k_sat_val > result_clay.k_sat_val

    # Sandy soil has lower porosity
    @test result_sand.θ_sat_val < result_clay.θ_sat_val

    # Both should produce solvable problems
    sol_sand = solve(result_sand.prob, saveat = 0.1)
    sol_clay = solve(result_clay.prob, saveat = 0.1)
    @test sol_sand.retcode == :Success || sol_sand.retcode == ReturnCode.Success
    @test sol_clay.retcode == :Success || sol_clay.retcode == ReturnCode.Success
end
