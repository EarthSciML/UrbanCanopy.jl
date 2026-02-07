@testsnippet AlbedoRadSetup begin
    using Test
    using ModelingToolkit
    using DynamicQuantities
    using OrdinaryDiffEqDefault
    using UrbanCanopy

    function default_params(compiled)
        Dict(
            compiled.H_W => 1.0,
            compiled.W_roof => 0.25,
            compiled.S_atm_dir_vis => 200.0,
            compiled.S_atm_dir_nir => 150.0,
            compiled.S_atm_dif_vis => 80.0,
            compiled.S_atm_dif_nir => 60.0,
            compiled.L_atm_down => 350.0,
            compiled.μ_zen => 0.5236,  # ~30 degrees
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
            compiled.ε_roof => 0.9,
            compiled.ε_road => 0.95,
            compiled.ε_wall => 0.85,
            compiled.T_roof => 305.0,
            compiled.T_road => 310.0,
            compiled.T_sunwall => 315.0,
            compiled.T_shdwall => 300.0,
        )
    end

    function solve_system(compiled, params)
        prob = ODEProblem(compiled, params, (0.0, 1.0))
        solve(prob)
    end
end

@testitem "Structural Verification" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()

    @test length(equations(sys)) == 62
    @test length(unknowns(sys)) == 62
    @test length(parameters(sys)) == 31

    unk_names = Set(string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys))
    for name in [
            "Ψ_wall_sky", "Ψ_road_sky", "Ψ_wall_wall", "Ψ_road_wall",
            "S_net_roof", "S_net_road", "S_net_sunwall", "S_net_shdwall",
            "L_net_roof", "L_net_road", "L_net_total", "α_uc_dir_vis",
        ]
        @test name in unk_names
    end

    param_names = Set(string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys))
    for name in [
            "H_W", "W_roof", "S_atm_dir_vis", "L_atm_down", "μ_zen",
            "α_roof_dir_vis", "α_road_dir_vis", "α_wall_dir_vis",
            "ε_roof", "ε_road", "ε_wall", "T_roof", "T_road",
        ]
        @test name in param_names
    end

    compiled = mtkcompile(sys)
    @test compiled isa ModelingToolkit.AbstractSystem
end

@testitem "View Factors - Known Values" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)

    # For H/W = 1 (Eqs. 2.24-2.28)
    hw = 1.0
    ψ_ws_expected = (hw + 1 - sqrt(1 + hw^2)) / (2 * hw)
    ψ_rs_expected = sqrt(1 + hw^2) - hw
    ψ_rw_expected = (1 - ψ_rs_expected) / 2
    ψ_ww_expected = 1 - 2 * ψ_ws_expected

    @test sol[compiled.Ψ_wall_sky][end] ≈ ψ_ws_expected rtol = 1.0e-10
    @test sol[compiled.Ψ_road_sky][end] ≈ ψ_rs_expected rtol = 1.0e-10
    @test sol[compiled.Ψ_road_wall][end] ≈ ψ_rw_expected rtol = 1.0e-10
    @test sol[compiled.Ψ_wall_wall][end] ≈ ψ_ww_expected rtol = 1.0e-10

    # View factors from wall sum to 1
    @test sol[compiled.Ψ_wall_sky][end] + sol[compiled.Ψ_wall_road][end] +
        sol[compiled.Ψ_wall_wall][end] ≈ 1.0 rtol = 1.0e-10

    # View factors from road sum to 1
    @test sol[compiled.Ψ_road_sky][end] + 2 * sol[compiled.Ψ_road_wall][end] ≈ 1.0 rtol = 1.0e-10

    # Symmetry: Ψ_wall_sky = Ψ_wall_road (Eq. 2.26)
    @test sol[compiled.Ψ_wall_sky][end] ≈ sol[compiled.Ψ_wall_road][end] rtol = 1.0e-10

    # Symmetry: Ψ_road_sky = Ψ_sky_road (Eq. 2.25)
    @test sol[compiled.Ψ_road_sky][end] ≈ sol[compiled.Ψ_sky_road][end] rtol = 1.0e-10
end

@testitem "View Factors - Limiting Behaviors" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)

    # Low H/W (flat canyon): road sees mostly sky
    params_low = default_params(compiled)
    params_low[compiled.H_W] = 0.01
    sol_low = solve_system(compiled, params_low)
    @test sol_low[compiled.Ψ_road_sky][end] > 0.99
    @test sol_low[compiled.Ψ_wall_wall][end] < 0.01
    @test sol_low[compiled.Ψ_wall_sky][end] > 0.45

    # High H/W (deep canyon): road sees little sky, walls see each other
    params_high = default_params(compiled)
    params_high[compiled.H_W] = 10.0
    sol_high = solve_system(compiled, params_high)
    @test sol_high[compiled.Ψ_road_sky][end] < 0.11
    @test sol_high[compiled.Ψ_wall_wall][end] > 0.8
end

@testitem "Direct Beam Solar" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)

    hw = 1.0
    s_road = sol[compiled.S_road_dir_vis][end]
    s_sunwall = sol[compiled.S_sunwall_dir_vis][end]

    # Road and sunwall factors should be positive
    @test s_road > 0
    @test s_sunwall > 0

    # Energy conservation: road_factor + H/W * sunwall_factor ≈ 1 (Eq. 2.18)
    @test sol[compiled.S_road_factor][end] + hw * sol[compiled.S_sunwall_factor][end] ≈ 1.0 rtol = 0.01

    # At high zenith angle (nearly overhead), most radiation hits road
    params_overhead = copy(params)
    params_overhead[compiled.μ_zen] = 0.1
    sol_overhead = solve_system(compiled, params_overhead)
    @test sol_overhead[compiled.S_road_factor][end] > sol_overhead[compiled.S_sunwall_factor][end]
end

@testitem "Roof Solar Radiation" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)

    # Roof absorbed = Σ S_atm*(1-α) (Eqs. 2.34-2.35)
    expected = 200.0 * (1 - 0.2) + 150.0 * (1 - 0.2) +
        80.0 * (1 - 0.2) + 60.0 * (1 - 0.2)
    @test sol[compiled.S_net_roof][end] ≈ expected rtol = 1.0e-10
end

@testitem "Energy Conservation - Solar" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)
    hw = params[compiled.H_W]

    # For each waveband: incident = absorbed + reflected
    for (net_road, net_sw, net_sh, refl, s_road, s_sw, s_sh) in [
            (
                compiled.S_net_road_dir_vis, compiled.S_net_sunwall_dir_vis,
                compiled.S_net_shdwall_dir_vis, compiled.S_refl_sky_dir_vis,
                compiled.S_road_dir_vis, compiled.S_sunwall_dir_vis, nothing,
            ),
            (
                compiled.S_net_road_dir_nir, compiled.S_net_sunwall_dir_nir,
                compiled.S_net_shdwall_dir_nir, compiled.S_refl_sky_dir_nir,
                compiled.S_road_dir_nir, compiled.S_sunwall_dir_nir, nothing,
            ),
            (
                compiled.S_net_road_dif_vis, compiled.S_net_sunwall_dif_vis,
                compiled.S_net_shdwall_dif_vis, compiled.S_refl_sky_dif_vis,
                compiled.S_road_dif_vis, compiled.S_wall_dif_vis, compiled.S_wall_dif_vis,
            ),
            (
                compiled.S_net_road_dif_nir, compiled.S_net_sunwall_dif_nir,
                compiled.S_net_shdwall_dif_nir, compiled.S_refl_sky_dif_nir,
                compiled.S_road_dif_nir, compiled.S_wall_dif_nir, compiled.S_wall_dif_nir,
            ),
        ]
        incident = sol[s_road][end] + sol[s_sw][end] * hw
        if s_sh !== nothing
            incident += sol[s_sh][end] * hw
        end
        absorbed = sol[net_road][end] + (sol[net_sw][end] + sol[net_sh][end]) * hw
        @test incident ≈ absorbed + sol[refl][end] rtol = 1.0e-8
    end
end

@testitem "Roof Longwave Radiation" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)

    σ = 5.67e-8
    L_up_expected = 0.9 * σ * 305.0^4 + 0.1 * 350.0
    @test sol[compiled.L_up_roof][end] ≈ L_up_expected rtol = 1.0e-8

    @test sol[compiled.L_net_roof][end] ≈ L_up_expected - 350.0 rtol = 1.0e-8
end

@testitem "Canyon Albedo Decreases with H/W" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)

    hw_values = [0.5, 1.0, 2.0, 4.0]
    albedos = Float64[]

    for hw in hw_values
        params = default_params(compiled)
        params[compiled.H_W] = hw
        sol = solve_system(compiled, params)
        push!(albedos, sol[compiled.α_uc_dif_vis][end])
    end

    # Monotonically decreasing (radiation trapping, Figure 2.6)
    for i in 1:(length(albedos) - 1)
        @test albedos[i] > albedos[i + 1]
    end

    for α in albedos
        @test 0.0 < α < 1.0
    end
end

@testitem "Multi-Reflection Increases Absorption" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)

    # Canyon albedo < surface albedos due to multi-reflection trapping
    @test sol[compiled.α_uc_dif_vis][end] < 0.25  # wall albedo
    @test sol[compiled.α_uc_dif_vis][end] < 0.08  # road albedo

    # For very shallow canyon, canyon albedo ≈ road albedo
    params_flat = copy(params)
    params_flat[compiled.H_W] = 0.01
    sol_flat = solve_system(compiled, params_flat)
    @test sol_flat[compiled.α_uc_dif_vis][end] ≈ 0.08 atol = 0.01
end

@testitem "Canyon Net LW Varies with H/W" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)

    hw_values = [0.5, 1.0, 2.0]
    net_lw = Float64[]

    for hw in hw_values
        params = default_params(compiled)
        params[compiled.H_W] = hw
        sol = solve_system(compiled, params)
        push!(net_lw, sol[compiled.L_net_uc][end])
    end

    # All values should be positive (surfaces warmer than atmosphere)
    for lw in net_lw
        @test lw > 0
    end

    # Net LW should be finite and bounded
    for lw in net_lw
        @test isfinite(lw)
        @test lw < 2000.0  # reasonable upper bound
    end
end

@testitem "Symmetry - Equal Wall Temperatures" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)

    params = default_params(compiled)
    params[compiled.T_sunwall] = 310.0
    params[compiled.T_shdwall] = 310.0
    sol = solve_system(compiled, params)

    @test sol[compiled.L_emit_sunwall][end] ≈ sol[compiled.L_emit_shdwall][end] rtol = 1.0e-10
    @test sol[compiled.L_net_sunwall][end] ≈ sol[compiled.L_net_shdwall][end] rtol = 1.0e-8

    # For diffuse solar, both walls receive equal radiation → equal absorption
    @test sol[compiled.S_net_sunwall_dif_vis][end] ≈ sol[compiled.S_net_shdwall_dif_vis][end] rtol = 1.0e-8
end

@testitem "Shaded Wall Gets No Direct Beam" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)

    # Shaded wall absorbs less direct beam than sunlit wall
    @test sol[compiled.S_net_shdwall_dir_vis][end] < sol[compiled.S_net_sunwall_dir_vis][end]
    @test sol[compiled.S_net_shdwall_dir_nir][end] < sol[compiled.S_net_sunwall_dir_nir][end]

    # With zero albedos (no reflection), shaded wall gets no direct beam at all
    params_black = copy(params)
    params_black[compiled.α_wall_dir_vis] = 0.0
    params_black[compiled.α_road_dir_vis] = 0.0
    sol_black = solve_system(compiled, params_black)
    @test sol_black[compiled.S_net_shdwall_dir_vis][end] ≈ 0.0 atol = 1.0e-10
end

@testitem "Total Net Radiation" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)
    sol = solve_system(compiled, params)

    w_roof = 0.25
    expected_solar = w_roof * sol[compiled.S_net_roof][end] +
        (1 - w_roof) * sol[compiled.S_net_uc][end]
    @test sol[compiled.S_net_total][end] ≈ expected_solar rtol = 1.0e-8

    expected_lw = w_roof * sol[compiled.L_net_roof][end] +
        (1 - w_roof) * sol[compiled.L_net_uc][end]
    @test sol[compiled.L_net_total][end] ≈ expected_lw rtol = 1.0e-8
end

@testitem "Zero Albedo - Full Absorption" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)
    params = default_params(compiled)

    for key in keys(params)
        name = string(Symbolics.tosymbol(key, escape = false))
        if startswith(name, "α_")
            params[key] = 0.0
        end
    end

    sol = solve_system(compiled, params)

    # All canyon reflected to sky should be zero
    @test sol[compiled.S_refl_sky_dir_vis][end] ≈ 0.0 atol = 1.0e-10
    @test sol[compiled.S_refl_sky_dif_vis][end] ≈ 0.0 atol = 1.0e-10

    # Canyon albedos should be zero
    @test sol[compiled.α_uc_dif_vis][end] ≈ 0.0 atol = 1.0e-10

    # Total absorbed = total incident for direct VIS
    hw = params[compiled.H_W]
    total_incident = sol[compiled.S_road_dir_vis][end] +
        sol[compiled.S_sunwall_dir_vis][end] * hw
    total_absorbed = sol[compiled.S_net_road_dir_vis][end] +
        (
        sol[compiled.S_net_sunwall_dir_vis][end] +
            sol[compiled.S_net_shdwall_dir_vis][end]
    ) * hw
    @test total_absorbed ≈ total_incident rtol = 1.0e-10
end

@testitem "Road Net LW Decreases with H/W" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)

    # When all surfaces are at the same temperature, road net LW should
    # decrease with H/W as more LW from walls replaces sky view (Figure 2.7)
    hw_values = [0.5, 1.0, 2.0, 5.0]
    net_lw_road = Float64[]

    for hw in hw_values
        params = default_params(compiled)
        params[compiled.H_W] = hw
        # Set all temperatures equal
        params[compiled.T_roof] = 292.16
        params[compiled.T_road] = 292.16
        params[compiled.T_sunwall] = 292.16
        params[compiled.T_shdwall] = 292.16
        # Set all emissivities equal
        params[compiled.ε_roof] = 0.95
        params[compiled.ε_road] = 0.95
        params[compiled.ε_wall] = 0.95
        sol = solve_system(compiled, params)
        push!(net_lw_road, sol[compiled.L_net_road][end])
    end

    # Net LW for road decreases with H/W (more LW from warm walls reduces net loss)
    for i in 1:(length(net_lw_road) - 1)
        @test net_lw_road[i + 1] < net_lw_road[i]
    end
end

@testitem "Roof LW Independent of H/W" setup = [AlbedoRadSetup] tags = [:albedo_rad] begin
    sys = UrbanRadiation()
    compiled = mtkcompile(sys)

    net_lw_roof = Float64[]
    for hw in [0.5, 1.0, 2.0, 5.0]
        params = default_params(compiled)
        params[compiled.H_W] = hw
        sol = solve_system(compiled, params)
        push!(net_lw_roof, sol[compiled.L_net_roof][end])
    end

    # Roof LW is independent of canyon geometry
    @test all(x -> isapprox(x, net_lw_roof[1], rtol = 1.0e-10), net_lw_roof)
end
