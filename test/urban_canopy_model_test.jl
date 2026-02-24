@testsnippet UrbanCanopyModelSetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy

    """
    Create default parameter values for the UrbanCanopyModel.
    Returns a compiled system and dictionary of parameter => value pairs.
    """
    function ucm_default_params()
        sys = UrbanCanopyModel()
        compiled = mtkcompile(sys)
        defaults = Dict(
            # Met forcing
            compiled.S_atm => 500.0,        # W/m^2
            compiled.T_atm => 300.0,         # K
            compiled.P_atm => 101325.0,      # Pa
            compiled.q_atm => 0.01,          # kg/kg
            compiled.W_atm => 4.0,           # m/s
            compiled.P_precip => 0.0,        # kg/(m^2*s)
            compiled.μ_zen => 0.7854,        # rad (~45 deg)
            # Canyon geometry
            compiled.H_W => 1.0,
            compiled.W_roof => 0.25,
            compiled.H_canyon => 10.0,        # m
            compiled.f_prvrd => 0.5,
            compiled.H_w_est => 5.0,          # m
            compiled.z_0 => 0.5,              # m
            compiled.z_d => 5.0,              # m
            # Surface albedos (typical urban values)
            compiled.α_roof_dir_vis => 0.2,
            compiled.α_roof_dir_nir => 0.2,
            compiled.α_roof_dif_vis => 0.2,
            compiled.α_roof_dif_nir => 0.2,
            compiled.α_road_dir_vis => 0.08,
            compiled.α_road_dir_nir => 0.08,
            compiled.α_road_dif_vis => 0.08,
            compiled.α_road_dif_nir => 0.08,
            compiled.α_wall_dir_vis => 0.25,
            compiled.α_wall_dir_nir => 0.25,
            compiled.α_wall_dif_vis => 0.25,
            compiled.α_wall_dif_nir => 0.25,
            # Surface emissivities
            compiled.ε_roof => 0.9,
            compiled.ε_road => 0.95,
            compiled.ε_wall => 0.9,
            # Surface temperatures (warmer than atmosphere = UHI)
            compiled.T_g_roof => 310.0,       # K
            compiled.T_g_prvrd => 305.0,      # K
            compiled.T_g_imprvrd => 308.0,    # K
            compiled.T_g_sunwall => 312.0,    # K
            compiled.T_g_shdwall => 303.0,    # K
            # Surface humidities
            compiled.q_g_roof => 0.015,       # kg/kg
            compiled.q_g_prvrd => 0.012,      # kg/kg
            compiled.q_g_imprvrd => 0.014,    # kg/kg
            # Wetted fractions
            compiled.f_wet_roof => 0.0,
            compiled.f_wet_imprvrd => 0.0,
            # Stability
            compiled.ζ_in => 0.1,
        )
        return compiled, defaults
    end

    """
    Solve the UrbanCanopyModel with given parameters.
    Returns the solution.
    """
    function ucm_solve(compiled, params)
        prob = ODEProblem(compiled, params, (0.0, 1.0))
        sol = solve(prob)
        return sol
    end
end

@testitem "Structural Verification" setup = [UrbanCanopyModelSetup] tags = [:urban_canopy_model] begin
    sys = UrbanCanopyModel()

    # Verify subsystems exist
    subsys_names = Set(string(ModelingToolkit.nameof(s)) for s in ModelingToolkit.get_systems(sys))
    @test "offline" in subsys_names
    @test "atmos" in subsys_names
    @test "rad" in subsys_names
    @test "flux" in subsys_names
    @test length(ModelingToolkit.get_systems(sys)) == 4

    # Verify the system compiles
    compiled, defaults = ucm_default_params()
    @test compiled isa ModelingToolkit.AbstractSystem

    # Verify it can create an ODEProblem
    prob = ODEProblem(compiled, defaults, (0.0, 1.0))
    @test prob isa ODEProblem
end

@testitem "Unit Verification" setup = [UrbanCanopyModelSetup] tags = [:urban_canopy_model] begin
    compiled, defaults = ucm_default_params()

    # Solve and check that key outputs are accessible and have reasonable values
    sol = ucm_solve(compiled, defaults)

    # T_ac should be in Kelvin range
    T_ac = sol[compiled.flux.T_ac][end]
    @test 250.0 < T_ac < 350.0

    # H_total should be finite and in W/m^2 range
    H_total = sol[compiled.flux.H_total][end]
    @test isfinite(H_total)

    # E_total should be finite
    E_total = sol[compiled.flux.E_total][end]
    @test isfinite(E_total)

    # Radiation outputs should be accessible
    S_net = sol[compiled.rad.S_net_total][end]
    @test isfinite(S_net)
    @test S_net > 0  # Sun is shining

    L_net = sol[compiled.rad.L_net_total][end]
    @test isfinite(L_net)
end

@testitem "Subsystem Connectivity" setup = [UrbanCanopyModelSetup] tags = [:urban_canopy_model] begin
    compiled, defaults = ucm_default_params()

    # Solve with baseline T_atm = 300 K
    sol_base = ucm_solve(compiled, defaults)
    T_ac_base = sol_base[compiled.flux.T_ac][end]

    # Increase T_atm and verify T_ac responds
    params_warm = copy(defaults)
    params_warm[compiled.T_atm] = 310.0
    sol_warm = ucm_solve(compiled, params_warm)
    T_ac_warm = sol_warm[compiled.flux.T_ac][end]

    # Warmer atmosphere should increase T_ac
    @test T_ac_warm > T_ac_base

    # Increase S_atm and verify radiation responds
    params_bright = copy(defaults)
    params_bright[compiled.S_atm] = 800.0
    sol_bright = ucm_solve(compiled, params_bright)
    S_net_bright = sol_bright[compiled.rad.S_net_total][end]

    sol_dim = ucm_solve(compiled, defaults)
    S_net_base = sol_dim[compiled.rad.S_net_total][end]

    # More solar radiation should increase net absorbed solar
    @test S_net_bright > S_net_base
end

@testitem "Physical Consistency" setup = [UrbanCanopyModelSetup] tags = [:urban_canopy_model] begin
    compiled, defaults = ucm_default_params()
    sol = ucm_solve(compiled, defaults)

    T_ac = sol[compiled.flux.T_ac][end]
    H_total = sol[compiled.flux.H_total][end]
    E_total = sol[compiled.flux.E_total][end]
    T_atm_val = defaults[compiled.T_atm]

    # T_ac should be between atmospheric and surface temperatures
    T_surfaces = [
        defaults[compiled.T_g_roof],
        defaults[compiled.T_g_prvrd],
        defaults[compiled.T_g_imprvrd],
        defaults[compiled.T_g_sunwall],
        defaults[compiled.T_g_shdwall],
    ]
    T_min = min(T_atm_val, minimum(T_surfaces))
    T_max = max(T_atm_val, maximum(T_surfaces))
    @test T_min ≤ T_ac ≤ T_max

    # H_total should be positive when surfaces are warmer than atmosphere
    # (upward heat flux from warm surfaces)
    @test H_total > 0

    # Wind speed sensitivity: increasing wind should increase heat flux magnitude
    params_wind = copy(defaults)
    params_wind[compiled.W_atm] = 8.0  # Double the wind speed
    sol_wind = ucm_solve(compiled, params_wind)
    H_total_wind = sol_wind[compiled.flux.H_total][end]
    @test abs(H_total_wind) > abs(H_total)
end

@testitem "Stability Sensitivity" setup = [UrbanCanopyModelSetup] tags = [:urban_canopy_model] begin
    compiled, defaults = ucm_default_params()

    # Stable conditions (ζ > 0)
    params_stable = copy(defaults)
    params_stable[compiled.ζ_in] = 0.5
    sol_stable = ucm_solve(compiled, params_stable)
    T_ac_stable = sol_stable[compiled.flux.T_ac][end]

    # Unstable conditions (ζ < 0) → more mixing → T_ac closer to surface temps
    params_unstable = copy(defaults)
    params_unstable[compiled.ζ_in] = -0.5
    sol_unstable = ucm_solve(compiled, params_unstable)
    T_ac_unstable = sol_unstable[compiled.flux.T_ac][end]

    # Under unstable conditions with warmer surfaces, T_ac should differ from stable
    @test T_ac_stable != T_ac_unstable

    # Check that H_total also responds
    H_stable = sol_stable[compiled.flux.H_total][end]
    H_unstable = sol_unstable[compiled.flux.H_total][end]
    @test H_stable != H_unstable
end

@testitem "Radiation Consistency" setup = [UrbanCanopyModelSetup] tags = [:urban_canopy_model] begin
    compiled, defaults = ucm_default_params()
    sol_base = ucm_solve(compiled, defaults)

    S_net_base = sol_base[compiled.rad.S_net_total][end]
    L_net_base = sol_base[compiled.rad.L_net_total][end]

    # Net absorbed solar should be positive when sun is shining
    @test S_net_base > 0

    # Increasing surface temperature should increase net longwave (more emission)
    params_hot = copy(defaults)
    params_hot[compiled.T_g_roof] = 330.0
    params_hot[compiled.T_g_imprvrd] = 328.0
    params_hot[compiled.T_g_sunwall] = 332.0
    params_hot[compiled.T_g_shdwall] = 320.0
    sol_hot = ucm_solve(compiled, params_hot)
    L_net_hot = sol_hot[compiled.rad.L_net_total][end]

    # Hotter surfaces emit more LW → net LW should increase (more positive = more emission)
    @test L_net_hot > L_net_base

    # Increasing S_atm should increase S_net
    params_bright = copy(defaults)
    params_bright[compiled.S_atm] = 900.0
    sol_bright = ucm_solve(compiled, params_bright)
    S_net_bright = sol_bright[compiled.rad.S_net_total][end]
    @test S_net_bright > S_net_base
end
