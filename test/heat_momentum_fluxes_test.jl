@testsnippet HeatMomentumSetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy

    """
    Create default parameter values for a typical urban canyon scenario.
    Returns a compiled system and dictionary of parameter => value pairs.
    """
    function default_params(sys)
        compiled = mtkcompile(sys)
        # Typical mid-latitude urban canyon
        defaults = Dict(
            compiled.H_W => 1.0,             # H/W = 1 (symmetric canyon)
            compiled.H_canyon => 10.0,        # 10 m building height
            compiled.W_roof => 0.25,          # 25% roof fraction
            compiled.f_prvrd => 0.5,          # 50% of road is pervious
            compiled.H_w_est => 5.0,          # Wind estimated at 5 m height
            compiled.u_atm => 3.0,            # 3 m/s zonal wind
            compiled.v_atm => 2.0,            # 2 m/s meridional wind
            compiled.θ_atm => 300.0,          # 300 K potential temperature
            compiled.q_atm => 0.01,           # 10 g/kg specific humidity
            compiled.ρ_atm => 1.2,            # ~1.2 kg/m³ air density
            compiled.z_atm_m => 30.0,         # 30 m momentum reference height
            compiled.z_atm_h => 30.0,         # 30 m heat reference height
            compiled.z_atm_w => 30.0,         # 30 m moisture reference height
            compiled.T_g_roof => 310.0,       # Warm roof
            compiled.T_g_prvrd => 305.0,      # Warm pervious road
            compiled.T_g_imprvrd => 308.0,    # Warm impervious road
            compiled.T_g_sunwall => 312.0,    # Warm sunlit wall
            compiled.T_g_shdwall => 303.0,    # Cool shaded wall
            compiled.q_g_roof => 0.015,       # Saturated specific humidity at roof temp
            compiled.q_g_prvrd => 0.012,      # Pervious road humidity
            compiled.q_g_imprvrd => 0.014,    # Saturated at impervious road temp
            compiled.f_wet_roof => 0.0,       # Dry roof
            compiled.f_wet_imprvrd => 0.0,    # Dry impervious road
            compiled.ζ_in => 0.1,             # Mildly stable MO stability parameter
        )
        return compiled, defaults
    end

    """
    Solve the system with given parameters.
    Returns the solution.
    """
    function solve_system(compiled, params)
        prob = ODEProblem(compiled, params, (0.0, 1.0))
        sol = solve(prob)
        return sol
    end
end

@testitem "Structural Verification" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()

    # Verify equation count
    @test length(equations(sys)) == 48

    # Verify unknown count
    @test length(unknowns(sys)) == 48

    # Verify key unknown names
    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))

    # Roughness parameters
    @test "d_canopy" in unk_names
    @test "z_0m_canopy" in unk_names
    @test "λ_p" in unk_names
    @test "λ_F" in unk_names

    # MO similarity variables
    @test "u_star" in unk_names
    @test "L_MO" in unk_names
    @test "V_a" in unk_names
    @test "V_r" in unk_names

    # Stability corrections
    @test "ψ_m_atm" in unk_names
    @test "ψ_m_0" in unk_names
    @test "ψ_h_atm" in unk_names
    @test "ψ_h_0" in unk_names
    @test "ψ_w_atm" in unk_names
    @test "ψ_w_0" in unk_names
    @test "ζ_0m" in unk_names

    # Resistances
    @test "r_am" in unk_names
    @test "r_ah" in unk_names
    @test "r_aw" in unk_names
    @test "r_s_u" in unk_names

    # Canyon wind
    @test "U_can" in unk_names
    @test "U_ac" in unk_names

    # Heat fluxes
    @test "H_roof" in unk_names
    @test "H_prvrd" in unk_names
    @test "H_imprvrd" in unk_names
    @test "H_sunwall" in unk_names
    @test "H_shdwall" in unk_names
    @test "H_total" in unk_names

    # Water vapor fluxes
    @test "E_roof" in unk_names
    @test "E_prvrd" in unk_names
    @test "E_imprvrd" in unk_names
    @test "E_total" in unk_names

    # Momentum fluxes
    @test "τ_x" in unk_names
    @test "τ_y" in unk_names

    # UCL diagnostics
    @test "T_ac" in unk_names
    @test "q_ac" in unk_names

    # Conductances
    @test "c_a_h" in unk_names
    @test "c_a_w" in unk_names

    # Convective velocity
    @test "U_c" in unk_names
    @test "θ_v_atm" in unk_names

    # Partial derivatives
    @test "dH_roof_dT" in unk_names
    @test "dH_prvrd_dT" in unk_names
    @test "dH_imprvrd_dT" in unk_names
    @test "dH_sunwall_dT" in unk_names
    @test "dH_shdwall_dT" in unk_names

    # Verify key parameter names
    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    @test "H_W" in param_names
    @test "H_canyon" in param_names
    @test "W_roof" in param_names
    @test "f_prvrd" in param_names
    @test "u_atm" in param_names
    @test "v_atm" in param_names
    @test "θ_atm" in param_names
    @test "T_g_roof" in param_names
    @test "ρ_atm" in param_names
    @test "ζ_in" in param_names

    # Verify system compiles
    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "Unit Verification" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()

    expected_units = Dict(
        "d_canopy" => u"m",
        "z_0m_canopy" => u"m",
        "u_star" => u"m/s",
        "L_MO" => u"m",
        "V_a" => u"m/s",
        "V_r" => u"m/s",
        "r_am" => u"s/m",
        "r_ah" => u"s/m",
        "r_aw" => u"s/m",
        "r_s_u" => u"s/m",
        "U_can" => u"m/s",
        "U_ac" => u"m/s",
        "H_roof" => u"W/m^2",
        "H_prvrd" => u"W/m^2",
        "H_imprvrd" => u"W/m^2",
        "H_sunwall" => u"W/m^2",
        "H_shdwall" => u"W/m^2",
        "H_total" => u"W/m^2",
        "E_roof" => u"kg/(m^2*s)",
        "E_prvrd" => u"kg/(m^2*s)",
        "E_imprvrd" => u"kg/(m^2*s)",
        "E_total" => u"kg/(m^2*s)",
        "τ_x" => u"kg/(m*s^2)",
        "τ_y" => u"kg/(m*s^2)",
        "T_ac" => u"K",
        "q_ac" => u"kg/kg",
        "c_a_h" => u"m/s",
        "c_a_w" => u"m/s",
        "U_c" => u"m/s",
        "θ_v_atm" => u"K",
        "dH_roof_dT" => u"W/(m^2*K)",
        "dH_prvrd_dT" => u"W/(m^2*K)",
        "dH_imprvrd_dT" => u"W/(m^2*K)",
        "dH_sunwall_dT" => u"W/(m^2*K)",
        "dH_shdwall_dT" => u"W/(m^2*K)",
    )

    for v in unknowns(sys)
        name = string(Symbolics.tosymbol(v, escape = false))
        if haskey(expected_units, name)
            u = ModelingToolkit.get_unit(v)
            @test u == expected_units[name]
        end
    end
end

@testitem "Roughness Length and Displacement Height" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    # With H/W = 1, plan area index λ_p = 1/(1+1) = 0.5
    @test sol[compiled.λ_p][end] ≈ 0.5 rtol = 1.0e-10

    # Frontal area index λ_F = (1 - λ_p)(H/W) = 0.5 * 1.0 = 0.5
    @test sol[compiled.λ_F][end] ≈ 0.5 rtol = 1.0e-10

    # Displacement height should be between 0 and H
    d = sol[compiled.d_canopy][end]
    @test d > 0.0
    @test d < 10.0  # H_canyon = 10

    # Eq. 3.55: d = H[1 + α^{-λ_p}(λ_p - 1)]
    d_expected = 10.0 * (1 + 4.43^(-0.5) * (0.5 - 1))
    @test d ≈ d_expected rtol = 1.0e-6

    # Roughness length should be positive and less than H - d
    z0m = sol[compiled.z_0m_canopy][end]
    @test z0m > 0.0
    @test z0m < 10.0 - d  # Must be less than H - d

    # Check z_0m formula (Eq. 3.57)
    z0m_expected = 10.0 * (1 - d_expected / 10.0) *
        exp(-(0.5 * 1.0 * 1.2 / 0.4^2 * (1 - d_expected / 10.0) * 0.5)^(-0.5))
    @test z0m ≈ z0m_expected rtol = 1.0e-6
end

@testitem "Roughness - Limiting H/W" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()

    # Test with small H/W (open terrain limit)
    compiled, params = default_params(sys)
    params[compiled.H_W] = 0.1
    sol = solve_system(compiled, params)

    # λ_p = 0.1/1.1 ≈ 0.0909
    @test sol[compiled.λ_p][end] ≈ 0.1 / 1.1 rtol = 1.0e-10

    # d_canopy should be small relative to H
    d = sol[compiled.d_canopy][end]
    @test d < 5.0  # Much less than H_canyon = 10

    # Test with large H/W (deep canyon)
    params[compiled.H_W] = 5.0
    sol2 = solve_system(compiled, params)

    # λ_p = 5/6 ≈ 0.833
    @test sol2[compiled.λ_p][end] ≈ 5.0 / 6.0 rtol = 1.0e-10

    # d_canopy should be larger for deeper canyons
    d2 = sol2[compiled.d_canopy][end]
    @test d2 > d
end

@testitem "Stability Functions - Stable vs Unstable" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Stable conditions: ζ_in > 0
    params_stable = copy(params)
    params_stable[compiled.ζ_in] = 0.5
    sol_stable = solve_system(compiled, params_stable)

    # In stable conditions, ψ_m should be negative (reduces transfer)
    @test sol_stable[compiled.ψ_m_atm][end] < 0.0
    @test sol_stable[compiled.ψ_h_atm][end] < 0.0

    # Unstable conditions: ζ_in < 0
    params_unstable = copy(params)
    params_unstable[compiled.ζ_in] = -0.5
    sol_unstable = solve_system(compiled, params_unstable)

    # In unstable conditions, ψ_m should be positive (enhances transfer)
    @test sol_unstable[compiled.ψ_m_atm][end] > 0.0
    @test sol_unstable[compiled.ψ_h_atm][end] > 0.0

    # Friction velocity should be higher in unstable conditions
    # (positive ψ_m reduces the log-law denominator, increasing u_*)
    @test sol_unstable[compiled.u_star][end] > sol_stable[compiled.u_star][end]
end

@testitem "Stability Functions - Very Unstable" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Very unstable: ζ < -1.574 (momentum) and ζ < -0.465 (heat)
    params_vu = copy(params)
    params_vu[compiled.ζ_in] = -3.0
    sol_vu = solve_system(compiled, params_vu)

    # ψ_m should still be positive and larger than for moderate instability
    @test sol_vu[compiled.ψ_m_atm][end] > 0.0

    # Moderately unstable
    params_mu = copy(params)
    params_mu[compiled.ζ_in] = -0.3
    sol_mu = solve_system(compiled, params_mu)

    @test sol_mu[compiled.ψ_m_atm][end] > 0.0

    # Very unstable should have larger ψ_m correction
    @test sol_vu[compiled.ψ_m_atm][end] > sol_mu[compiled.ψ_m_atm][end]
end

@testitem "Stability Functions - Very Stable" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Very stable: ζ > 1
    params_vs = copy(params)
    params_vs[compiled.ζ_in] = 2.0
    sol_vs = solve_system(compiled, params_vs)

    # ψ_m should be strongly negative
    @test sol_vs[compiled.ψ_m_atm][end] < 0.0

    # Moderately stable
    params_ms = copy(params)
    params_ms[compiled.ζ_in] = 0.3
    sol_ms = solve_system(compiled, params_ms)

    # Very stable should have more negative ψ_m
    @test sol_vs[compiled.ψ_m_atm][end] < sol_ms[compiled.ψ_m_atm][end]
end

@testitem "Aerodynamic Resistances" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    # All resistances must be positive
    @test sol[compiled.r_am][end] > 0.0
    @test sol[compiled.r_ah][end] > 0.0
    @test sol[compiled.r_aw][end] > 0.0
    @test sol[compiled.r_s_u][end] > 0.0

    # r_am = V_a / u_*² (Eq. 3.65)
    V_a = sol[compiled.V_a][end]
    u_star = sol[compiled.u_star][end]
    @test sol[compiled.r_am][end] ≈ V_a / u_star^2 rtol = 1.0e-6

    # Surface resistance (Eq. 3.68)
    ρ = 1.2  # ρ_atm
    C_p = 1004.64
    U_ac = sol[compiled.U_ac][end]
    r_s_expected = ρ * C_p / (11.8 + 4.2 * U_ac)
    @test sol[compiled.r_s_u][end] ≈ r_s_expected rtol = 1.0e-6
end

@testitem "Canyon Wind Speed" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Test skimming flow (H/W ≥ 1)
    params[compiled.H_W] = 2.0
    sol = solve_system(compiled, params)
    U_can_skimming = sol[compiled.U_can][end]

    # Canyon wind should be positive and less than reference wind
    V_r = max(sqrt(3.0^2 + 2.0^2), 1.0)
    @test U_can_skimming > 0.0
    @test U_can_skimming < V_r

    # U_ac = √(U_can² + u_*²) (Eq. 3.59)
    u_star = sol[compiled.u_star][end]
    @test sol[compiled.U_ac][end] ≈ sqrt(U_can_skimming^2 + u_star^2) rtol = 1.0e-6

    # Test isolated flow (H/W < 0.5)
    params[compiled.H_W] = 0.3
    sol2 = solve_system(compiled, params)
    U_can_isolated = sol2[compiled.U_can][end]
    @test U_can_isolated > 0.0
end

@testitem "UCL Temperature" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    T_ac = sol[compiled.T_ac][end]

    # UCL temperature should be between minimum and maximum surface/atm temperatures
    T_min = min(300.0, 310.0, 305.0, 308.0, 312.0, 303.0) # θ_atm and all T_g values
    T_max = max(300.0, 310.0, 305.0, 308.0, 312.0, 303.0)
    @test T_ac ≥ T_min
    @test T_ac ≤ T_max

    # T_ac is a weighted average (Eq. 3.75), so it should be strictly between extremes
    @test T_ac > T_min
    @test T_ac < T_max
end

@testitem "UCL Humidity" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    q_ac = sol[compiled.q_ac][end]

    # UCL humidity should be positive
    @test q_ac > 0.0

    # With dry roof and impervious road (f_wet = 0), the humidity
    # should still be bounded by atmospheric and surface values
    q_min = min(0.01, 0.012)  # q_atm and q_g_prvrd (only active surface)
    q_max = max(0.01, 0.012)
    @test q_ac ≥ q_min - 1.0e-6
    @test q_ac ≤ q_max + 1.0e-6
end

@testitem "Sensible Heat Flux Signs" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    T_ac = sol[compiled.T_ac][end]

    # H = -ρ C_p (T_ac - T_g) / r_s: flux is from surface to air
    # If T_g > T_ac, then (T_ac - T_g) < 0, so H = -ρ C_p (negative) / r_s > 0
    # Positive H means heat flux from surface to atmosphere (upward)

    # Roof is warmer than T_ac typically
    if T_ac < 310.0  # T_g_roof
        @test sol[compiled.H_roof][end] > 0.0
    end

    # Sensible heat from warmer surfaces should be positive (upward)
    # H_total should be positive when surfaces are warmer than atmosphere
    @test sol[compiled.H_total][end] > 0.0
end

@testitem "Wall Evaporation Zero" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    # Eqs. 3.79, 3.80: E_sunwall = E_shdwall = 0
    # These are not explicit variables, but the total E should only include
    # roof, pervious road, and impervious road contributions
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    # Total E = W_roof * E_roof + (1-W_roof)[f_prvrd * E_prvrd + (1-f_prvrd) * E_imprvrd]
    E_total = sol[compiled.E_total][end]
    E_roof = sol[compiled.E_roof][end]
    E_prvrd = sol[compiled.E_prvrd][end]
    E_imprvrd = sol[compiled.E_imprvrd][end]

    expected = 0.25 * E_roof + 0.75 * (0.5 * E_prvrd + 0.5 * E_imprvrd)
    @test E_total ≈ expected rtol = 1.0e-6
end

@testitem "Momentum Flux Direction" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    # τ_x = -ρ u_atm / r_am: when u_atm > 0, τ_x < 0 (downward momentum transfer)
    @test sol[compiled.τ_x][end] < 0.0

    # τ_y = -ρ v_atm / r_am: when v_atm > 0, τ_y < 0
    @test sol[compiled.τ_y][end] < 0.0

    # Magnitude proportional to wind components
    τ_x = sol[compiled.τ_x][end]
    τ_y = sol[compiled.τ_y][end]
    # τ_x/τ_y = u_atm/v_atm = 3/2
    @test τ_x / τ_y ≈ 3.0 / 2.0 rtol = 1.0e-6
end

@testitem "Friction Velocity Positivity" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    @test sol[compiled.u_star][end] > 0.0

    # u_* should be a reasonable value (typically 0.01 to 2 m/s)
    @test sol[compiled.u_star][end] < 5.0
    @test sol[compiled.u_star][end] ≥ 0.01
end

@testitem "Sensitivity to Wind Speed" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Low wind
    params_low = copy(params)
    params_low[compiled.u_atm] = 1.0
    params_low[compiled.v_atm] = 0.0
    sol_low = solve_system(compiled, params_low)

    # High wind
    params_high = copy(params)
    params_high[compiled.u_atm] = 10.0
    params_high[compiled.v_atm] = 0.0
    sol_high = solve_system(compiled, params_high)

    # Higher wind should give higher friction velocity
    @test sol_high[compiled.u_star][end] > sol_low[compiled.u_star][end]

    # Higher wind should give lower aerodynamic resistance
    @test sol_high[compiled.r_am][end] < sol_low[compiled.r_am][end]
end

@testitem "Wet Fraction Effect" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Dry surfaces
    params_dry = copy(params)
    params_dry[compiled.f_wet_roof] = 0.0
    params_dry[compiled.f_wet_imprvrd] = 0.0
    sol_dry = solve_system(compiled, params_dry)

    # Wet surfaces
    params_wet = copy(params)
    params_wet[compiled.f_wet_roof] = 1.0
    params_wet[compiled.f_wet_imprvrd] = 1.0
    sol_wet = solve_system(compiled, params_wet)

    # Wet surfaces should have larger evaporation magnitude
    @test abs(sol_wet[compiled.E_total][end]) ≥ abs(sol_dry[compiled.E_total][end]) - 1.0e-10
end

@testitem "Monin-Obukhov Length Consistency" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # L_MO = (z_atm_m - d) / ζ_in (Eq. 3.50 rearranged)
    sol = solve_system(compiled, params)

    z_atm_m = 30.0
    d = sol[compiled.d_canopy][end]
    ζ_in = 0.1
    L_expected = (z_atm_m - d) / ζ_in
    @test sol[compiled.L_MO][end] ≈ L_expected rtol = 1.0e-6

    # L_MO should be positive for stable conditions (ζ > 0)
    @test sol[compiled.L_MO][end] > 0.0

    # L_MO should be negative for unstable conditions (ζ < 0)
    params_unstable = copy(params)
    params_unstable[compiled.ζ_in] = -0.5
    sol_u = solve_system(compiled, params_unstable)
    @test sol_u[compiled.L_MO][end] < 0.0
end

@testitem "Neutral Stability Limit" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Near-neutral conditions: ζ ≈ 0 (use small positive value)
    params_neutral = copy(params)
    params_neutral[compiled.ζ_in] = 0.001
    sol = solve_system(compiled, params_neutral)

    # ψ corrections should be near zero
    @test abs(sol[compiled.ψ_m_atm][end]) < 0.1
    @test abs(sol[compiled.ψ_h_atm][end]) < 0.1

    # u_* ≈ k * V_a / ln((z-d)/z0) for neutral conditions
    V_a = sol[compiled.V_a][end]
    d = sol[compiled.d_canopy][end]
    z0m = sol[compiled.z_0m_canopy][end]
    u_star_neutral = 0.4 * V_a / log((30.0 - d) / z0m)
    @test sol[compiled.u_star][end] ≈ u_star_neutral rtol = 0.05
end

@testitem "Convective Velocity" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)

    # Stable conditions (ζ > 0): U_c should be zero
    params_stable = copy(params)
    params_stable[compiled.ζ_in] = 0.5
    sol_stable = solve_system(compiled, params_stable)
    @test sol_stable[compiled.U_c][end] ≈ 0.0 atol = 1.0e-10

    # V_a should equal V_r for stable conditions (U_c = 0)
    @test sol_stable[compiled.V_a][end] ≈ sol_stable[compiled.V_r][end] rtol = 1.0e-6

    # Unstable conditions (ζ < 0): U_c should be positive
    params_unstable = copy(params)
    params_unstable[compiled.ζ_in] = -0.5
    sol_unstable = solve_system(compiled, params_unstable)
    @test sol_unstable[compiled.U_c][end] > 0.0

    # V_a should be larger than V_r for unstable conditions (U_c > 0)
    @test sol_unstable[compiled.V_a][end] > sol_unstable[compiled.V_r][end]

    # More unstable conditions should give larger U_c
    params_very_unstable = copy(params)
    params_very_unstable[compiled.ζ_in] = -2.0
    sol_vu = solve_system(compiled, params_very_unstable)
    @test sol_vu[compiled.U_c][end] > sol_unstable[compiled.U_c][end]

    # Virtual potential temperature (θ_v = θ(1 + 0.61*q))
    θ_v = sol_stable[compiled.θ_v_atm][end]
    @test θ_v ≈ 300.0 * (1 + 0.61 * 0.01) rtol = 1.0e-6
end

@testitem "Partial Derivatives of Sensible Heat" setup = [HeatMomentumSetup] tags = [:heat_momentum] begin
    sys = HeatMomentumFluxes()
    compiled, params = default_params(sys)
    sol = solve_system(compiled, params)

    # All partial derivatives should be positive (warming surface increases upward heat flux)
    @test sol[compiled.dH_roof_dT][end] > 0.0
    @test sol[compiled.dH_prvrd_dT][end] > 0.0
    @test sol[compiled.dH_imprvrd_dT][end] > 0.0
    @test sol[compiled.dH_sunwall_dT][end] > 0.0
    @test sol[compiled.dH_shdwall_dT][end] > 0.0

    # Partial derivatives should have reasonable magnitudes
    # Typical values are O(1-100) W/(m²·K)
    for var in [
            compiled.dH_roof_dT, compiled.dH_prvrd_dT, compiled.dH_imprvrd_dT,
            compiled.dH_sunwall_dT, compiled.dH_shdwall_dT,
        ]
        @test sol[var][end] > 0.1    # Not too small
        @test sol[var][end] < 1000.0  # Not too large
    end

    # Verify numerical consistency: perturb T_g_roof and check that
    # the finite difference matches dH_roof_dT
    ΔT = 0.01
    params_plus = copy(params)
    params_plus[compiled.T_g_roof] = 310.0 + ΔT
    sol_plus = solve_system(compiled, params_plus)

    params_minus = copy(params)
    params_minus[compiled.T_g_roof] = 310.0 - ΔT
    sol_minus = solve_system(compiled, params_minus)

    dH_fd = (sol_plus[compiled.H_roof][end] - sol_minus[compiled.H_roof][end]) / (2 * ΔT)
    @test sol[compiled.dH_roof_dT][end] ≈ dH_fd rtol = 0.01
end
