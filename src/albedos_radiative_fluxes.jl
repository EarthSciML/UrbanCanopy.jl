export UrbanRadiation

"""
    UrbanRadiation(; name=:UrbanRadiation)

Urban albedos and radiative fluxes for the Community Land Model Urban (CLMU)
parameterization, following Chapter 2 of Oleson et al. (2010).

Computes the view factors, incident direct beam and diffuse solar radiation,
absorbed and reflected solar radiation (with closed-form multi-reflection solution
for the urban canyon), incident longwave radiation, and absorbed/reflected/emitted
longwave radiation for urban surfaces: roof, sunlit wall, shaded wall, and road.

The urban canyon is represented as an infinitely long canyon of width W bordered
by buildings of height H. Radiation interactions within the canyon account for
multiple reflections between walls and road using view factors. The iterative
multi-reflection scheme (Eqs. 2.63-2.85, 2.138-2.170) is solved in closed form
using the matrix geometric series: `S_absorbed = (I - A) * (I - T*A)^{-1} * S_0`,
where A is the albedo matrix, T is the view-factor transport matrix, and S_0 is
the initial incident radiation vector.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp. Chapter 2: Albedos and Radiative Fluxes (pp. 26-60).
"""
@component function UrbanRadiation(; name = :UrbanRadiation)

    # Physical Constants
    @constants begin
        σ_SB = 5.67e-8, [description = "Stefan-Boltzmann constant (Table 1.4)", unit = u"W/(m^2*K^4)"]
        π_val = Float64(π), [description = "Pi (dimensionless)"]
        one_rad = 1.0, [description = "Unit angle for non-dimensionalization", unit = u"rad"]
        one_wm2 = 1.0, [description = "Unit radiative flux", unit = u"W/m^2"]
        zero_wm2 = 0.0, [description = "Zero radiative flux", unit = u"W/m^2"]
    end

    # Canyon geometry
    @parameters begin
        H_W, [description = "Canyon height to width ratio H/W (dimensionless)"]
        W_roof, [description = "Roof fraction of urban surface (dimensionless)"]
    end

    # Solar radiation inputs
    @parameters begin
        S_atm_dir_vis, [description = "Incident direct beam visible solar radiation", unit = u"W/m^2"]
        S_atm_dir_nir, [description = "Incident direct beam NIR solar radiation", unit = u"W/m^2"]
        S_atm_dif_vis, [description = "Incident diffuse visible solar radiation", unit = u"W/m^2"]
        S_atm_dif_nir, [description = "Incident diffuse NIR solar radiation", unit = u"W/m^2"]
    end

    # Longwave and solar angle
    @parameters begin
        L_atm_down, [description = "Incident longwave radiation from atmosphere", unit = u"W/m^2"]
        μ_zen, [description = "Solar zenith angle", unit = u"rad"]
    end

    # Surface albedos (effective road albedo = weighted by pervious/impervious fraction)
    @parameters begin
        α_roof_dir_vis, [description = "Roof albedo, direct beam, visible (dimensionless)"]
        α_roof_dir_nir, [description = "Roof albedo, direct beam, NIR (dimensionless)"]
        α_roof_dif_vis, [description = "Roof albedo, diffuse, visible (dimensionless)"]
        α_roof_dif_nir, [description = "Roof albedo, diffuse, NIR (dimensionless)"]
        α_road_dir_vis, [description = "Effective road albedo, direct beam, visible (dimensionless)"]
        α_road_dir_nir, [description = "Effective road albedo, direct beam, NIR (dimensionless)"]
        α_road_dif_vis, [description = "Effective road albedo, diffuse, visible (dimensionless)"]
        α_road_dif_nir, [description = "Effective road albedo, diffuse, NIR (dimensionless)"]
        α_wall_dir_vis, [description = "Wall albedo, direct beam, visible (dimensionless)"]
        α_wall_dir_nir, [description = "Wall albedo, direct beam, NIR (dimensionless)"]
        α_wall_dif_vis, [description = "Wall albedo, diffuse, visible (dimensionless)"]
        α_wall_dif_nir, [description = "Wall albedo, diffuse, NIR (dimensionless)"]
    end

    # Surface emissivities and temperatures
    @parameters begin
        ε_roof, [description = "Roof emissivity (dimensionless)"]
        ε_road, [description = "Effective road emissivity (dimensionless)"]
        ε_wall, [description = "Wall emissivity (dimensionless)"]
        T_roof, [description = "Roof surface temperature", unit = u"K"]
        T_road, [description = "Road surface temperature", unit = u"K"]
        T_sunwall, [description = "Sunlit wall surface temperature", unit = u"K"]
        T_shdwall, [description = "Shaded wall surface temperature", unit = u"K"]
    end

    # View factors (Section 2.3, Eqs. 2.22-2.28)
    @variables begin
        Ψ_wall_sky(t), [description = "View factor: wall to sky (Eq. 2.24) (dimensionless)"]
        Ψ_wall_road(t), [description = "View factor: wall to road (Eq. 2.26) (dimensionless)"]
        Ψ_wall_wall(t), [description = "View factor: wall to wall (Eq. 2.28) (dimensionless)"]
        Ψ_road_sky(t), [description = "View factor: road to sky (Eq. 2.25) (dimensionless)"]
        Ψ_road_wall(t), [description = "View factor: road to wall (Eq. 2.27) (dimensionless)"]
        Ψ_sky_road(t), [description = "View factor: sky to road (Eq. 2.25) (dimensionless)"]
        Ψ_sky_wall(t), [description = "View factor: sky to wall (Eq. 2.24) (dimensionless)"]
    end

    # Direct beam solar intermediates
    @variables begin
        θ_0(t), [description = "Critical canyon orientation angle (Eq. 2.11)", unit = u"rad"]
        S_sunwall_factor(t), [description = "Sunlit wall direct beam factor (Eq. 2.16) (dimensionless)"]
        S_road_factor(t), [description = "Road direct beam factor (Eq. 2.17) (dimensionless)"]
    end

    # Incident solar radiation per surface (before multi-reflection)
    @variables begin
        S_road_dir_vis(t), [description = "Direct beam VIS on road (Eq. 2.17)", unit = u"W/m^2"]
        S_road_dir_nir(t), [description = "Direct beam NIR on road", unit = u"W/m^2"]
        S_sunwall_dir_vis(t), [description = "Direct beam VIS on sunlit wall (Eq. 2.16)", unit = u"W/m^2"]
        S_sunwall_dir_nir(t), [description = "Direct beam NIR on sunlit wall", unit = u"W/m^2"]
        S_road_dif_vis(t), [description = "Diffuse VIS on road (Eq. 2.30)", unit = u"W/m^2"]
        S_road_dif_nir(t), [description = "Diffuse NIR on road", unit = u"W/m^2"]
        S_wall_dif_vis(t), [description = "Diffuse VIS on wall (Eq. 2.31-2.32)", unit = u"W/m^2"]
        S_wall_dif_nir(t), [description = "Diffuse NIR on wall", unit = u"W/m^2"]
    end

    # Transport matrix elements: T maps reflected flux between canyon surfaces
    # with H/W area conversions (Section 2.5, pp. 42-43)
    @variables begin
        T_rw(t), [description = "Transport: wall to road = H/W * Ψ_wall_road (dimensionless)"]
        T_wr(t), [description = "Transport: road to wall = Ψ_road_wall / (H/W) (dimensionless)"]
        T_ww(t), [description = "Transport: wall to wall = Ψ_wall_wall (dimensionless)"]
    end

    # Multi-reflection determinants for each waveband×type
    @variables begin
        det_dir_vis(t), [description = "Determinant of (I-TA), direct beam VIS (dimensionless)"]
        det_dir_nir(t), [description = "Determinant of (I-TA), direct beam NIR (dimensionless)"]
        det_dif_vis(t), [description = "Determinant of (I-TA), diffuse VIS (dimensionless)"]
        det_dif_nir(t), [description = "Determinant of (I-TA), diffuse NIR (dimensionless)"]
    end

    # Net absorbed solar radiation per surface and waveband
    @variables begin
        S_net_roof(t), [description = "Net absorbed solar by roof (Eqs. 2.34-2.35)", unit = u"W/m^2"]
        S_net_road_dir_vis(t), [description = "Absorbed solar: road, direct VIS", unit = u"W/m^2"]
        S_net_road_dir_nir(t), [description = "Absorbed solar: road, direct NIR", unit = u"W/m^2"]
        S_net_road_dif_vis(t), [description = "Absorbed solar: road, diffuse VIS", unit = u"W/m^2"]
        S_net_road_dif_nir(t), [description = "Absorbed solar: road, diffuse NIR", unit = u"W/m^2"]
        S_net_sunwall_dir_vis(t), [description = "Absorbed solar: sunwall, direct VIS", unit = u"W/m^2"]
        S_net_sunwall_dir_nir(t), [description = "Absorbed solar: sunwall, direct NIR", unit = u"W/m^2"]
        S_net_sunwall_dif_vis(t), [description = "Absorbed solar: sunwall, diffuse VIS", unit = u"W/m^2"]
        S_net_sunwall_dif_nir(t), [description = "Absorbed solar: sunwall, diffuse NIR", unit = u"W/m^2"]
        S_net_shdwall_dir_vis(t), [description = "Absorbed solar: shdwall, direct VIS", unit = u"W/m^2"]
        S_net_shdwall_dir_nir(t), [description = "Absorbed solar: shdwall, direct NIR", unit = u"W/m^2"]
        S_net_shdwall_dif_vis(t), [description = "Absorbed solar: shdwall, diffuse VIS", unit = u"W/m^2"]
        S_net_shdwall_dif_nir(t), [description = "Absorbed solar: shdwall, diffuse NIR", unit = u"W/m^2"]
        S_net_road(t), [description = "Total net absorbed solar by road", unit = u"W/m^2"]
        S_net_sunwall(t), [description = "Total net absorbed solar by sunlit wall per wall area", unit = u"W/m^2"]
        S_net_shdwall(t), [description = "Total net absorbed solar by shaded wall per wall area", unit = u"W/m^2"]
        S_net_uc(t), [description = "Net absorbed solar by urban canyon per ground area", unit = u"W/m^2"]
        S_net_total(t), [description = "Total net absorbed solar (Eq. 2.95)", unit = u"W/m^2"]
    end

    # Reflected solar to sky (for canyon albedo)
    @variables begin
        S_refl_sky_dir_vis(t), [description = "Reflected to sky: direct VIS", unit = u"W/m^2"]
        S_refl_sky_dir_nir(t), [description = "Reflected to sky: direct NIR", unit = u"W/m^2"]
        S_refl_sky_dif_vis(t), [description = "Reflected to sky: diffuse VIS", unit = u"W/m^2"]
        S_refl_sky_dif_nir(t), [description = "Reflected to sky: diffuse NIR", unit = u"W/m^2"]
    end

    # Canyon albedos
    @variables begin
        α_uc_dir_vis(t), [description = "Canyon direct beam albedo, VIS (Eq. 2.93) (dimensionless)"]
        α_uc_dir_nir(t), [description = "Canyon direct beam albedo, NIR (dimensionless)"]
        α_uc_dif_vis(t), [description = "Canyon diffuse albedo, VIS (Eq. 2.94) (dimensionless)"]
        α_uc_dif_nir(t), [description = "Canyon diffuse albedo, NIR (dimensionless)"]
    end

    # Longwave radiation variables
    @variables begin
        L_up_roof(t), [description = "Upward LW from roof (Eq. 2.102)", unit = u"W/m^2"]
        L_net_roof(t), [description = "Net LW for roof (Eq. 2.101)", unit = u"W/m^2"]
        L_emit_road(t), [description = "Emitted LW from road (Eq. 2.111)", unit = u"W/m^2"]
        L_emit_sunwall(t), [description = "Emitted LW from sunlit wall (Eq. 2.112)", unit = u"W/m^2"]
        L_emit_shdwall(t), [description = "Emitted LW from shaded wall (Eq. 2.113)", unit = u"W/m^2"]
        det_lw(t), [description = "Determinant of (I-TA) for longwave (dimensionless)"]
        L_net_road(t), [description = "Net LW for road", unit = u"W/m^2"]
        L_net_sunwall(t), [description = "Net LW for sunlit wall per wall area", unit = u"W/m^2"]
        L_net_shdwall(t), [description = "Net LW for shaded wall per wall area", unit = u"W/m^2"]
        L_net_uc(t), [description = "Net LW for urban canyon per ground area (Eq. 2.172)", unit = u"W/m^2"]
        L_net_total(t), [description = "Total net LW (Eq. 2.175)", unit = u"W/m^2"]
    end

    eqs = Equation[]

    # ===== View Factors (Section 2.3, Eqs. 2.22-2.28) =====
    push!(eqs, Ψ_sky_wall ~ (H_W + 1 - sqrt(1 + H_W^2)) / (2 * H_W))        # Eq. 2.24
    push!(eqs, Ψ_wall_sky ~ Ψ_sky_wall)                                        # Eq. 2.22/2.24
    push!(eqs, Ψ_road_sky ~ sqrt(1 + H_W^2) - H_W)                            # Eq. 2.25
    push!(eqs, Ψ_sky_road ~ Ψ_road_sky)                                        # Eq. 2.25
    push!(eqs, Ψ_wall_road ~ Ψ_wall_sky)                                       # Eq. 2.26
    push!(eqs, Ψ_road_wall ~ (1 - Ψ_road_sky) / 2)                            # Eq. 2.27
    push!(eqs, Ψ_wall_wall ~ 1 - Ψ_wall_sky - Ψ_wall_road)                    # Eq. 2.28

    # ===== Transport matrix elements (Section 2.5, pp. 42-43) =====
    # T[road,wall] = H/W * Ψ_wall_road  (wall→road area conversion)
    # T[wall,road] = Ψ_road_wall / (H/W)  (road→wall area conversion)
    # T[wall,wall] = Ψ_wall_wall
    push!(eqs, T_rw ~ H_W * Ψ_wall_road)
    push!(eqs, T_wr ~ Ψ_road_wall / H_W)
    push!(eqs, T_ww ~ Ψ_wall_wall)

    # ===== Direct Beam Solar (Section 2.2, Eqs. 2.11, 2.14-2.17) =====
    push!(eqs, θ_0 ~ asin(min(1 / (H_W * tan(μ_zen / one_rad)), 1)) * one_rad) # Eq. 2.11
    push!(
        eqs, S_sunwall_factor ~ 2 * (
            (1 / H_W) * (1 / 2 - θ_0 / (π_val * one_rad)) +  # Eq. 2.16
                (1 / π_val) * tan(μ_zen / one_rad) * (1 - cos(θ_0 / one_rad))
        )
    )
    push!(eqs, S_sunwall_dir_vis ~ S_sunwall_factor * S_atm_dir_vis)
    push!(eqs, S_sunwall_dir_nir ~ S_sunwall_factor * S_atm_dir_nir)
    push!(
        eqs, S_road_factor ~ 2 * θ_0 / (π_val * one_rad) -                             # Eq. 2.17
            (2 / π_val) * H_W * tan(μ_zen / one_rad) * (1 - cos(θ_0 / one_rad))
    )
    push!(eqs, S_road_dir_vis ~ S_road_factor * S_atm_dir_vis)
    push!(eqs, S_road_dir_nir ~ S_road_factor * S_atm_dir_nir)

    # ===== Diffuse Solar (Section 2.4, Eqs. 2.29-2.32) =====
    push!(eqs, S_road_dif_vis ~ S_atm_dif_vis * Ψ_sky_road)   # Eq. 2.30
    push!(eqs, S_road_dif_nir ~ S_atm_dif_nir * Ψ_sky_road)
    push!(eqs, S_wall_dif_vis ~ S_atm_dif_vis * Ψ_sky_wall)   # Eqs. 2.31-2.32
    push!(eqs, S_wall_dif_nir ~ S_atm_dif_nir * Ψ_sky_wall)

    # ===== Roof Absorbed Solar (Eqs. 2.34-2.35) =====
    push!(
        eqs, S_net_roof ~ S_atm_dir_vis * (1 - α_roof_dir_vis) +
            S_atm_dir_nir * (1 - α_roof_dir_nir) +
            S_atm_dif_vis * (1 - α_roof_dif_vis) +
            S_atm_dif_nir * (1 - α_roof_dif_nir)
    )

    # ===== Canyon Multi-Reflection Solar (Section 2.5, closed-form) =====
    # The iterative multi-reflection (Eqs. 2.63-2.85) is equivalent to the
    # matrix geometric series: S_abs = (I-A) * (I-T*A)^{-1} * S_0
    #
    # With M = I - T*A, surfaces [road, sunwall, shdwall], and
    # b = a_w*T_rw, c = a_r*T_wr, d = a_w*T_ww:
    #
    #   M = [[1, -b, -b], [-c, 1, -d], [-c, -d, 1]]
    #   det(M) = (1 - d^2) - 2bc(1+d)
    #   adj(M) = [[1-d^2,   b(1+d), b(1+d)],
    #             [c(1+d),  1-bc,   d+bc  ],
    #             [c(1+d),  d+bc,   1-bc  ]]
    #   M^{-1} = adj(M) / det(M)
    for (
            a_r, a_w, S_r, S_sw, S_sh, det_var,
            net_road, net_sw, net_sh, refl_sky,
        ) in [
            (
                α_road_dir_vis, α_wall_dir_vis, S_road_dir_vis, S_sunwall_dir_vis, zero_wm2,
                det_dir_vis,
                S_net_road_dir_vis, S_net_sunwall_dir_vis, S_net_shdwall_dir_vis, S_refl_sky_dir_vis,
            ),
            (
                α_road_dir_nir, α_wall_dir_nir, S_road_dir_nir, S_sunwall_dir_nir, zero_wm2,
                det_dir_nir,
                S_net_road_dir_nir, S_net_sunwall_dir_nir, S_net_shdwall_dir_nir, S_refl_sky_dir_nir,
            ),
            (
                α_road_dif_vis, α_wall_dif_vis, S_road_dif_vis, S_wall_dif_vis, S_wall_dif_vis,
                det_dif_vis,
                S_net_road_dif_vis, S_net_sunwall_dif_vis, S_net_shdwall_dif_vis, S_refl_sky_dif_vis,
            ),
            (
                α_road_dif_nir, α_wall_dif_nir, S_road_dif_nir, S_wall_dif_nir, S_wall_dif_nir,
                det_dif_nir,
                S_net_road_dif_nir, S_net_sunwall_dif_nir, S_net_shdwall_dif_nir, S_refl_sky_dif_nir,
            ),
        ]
        local b = a_w * T_rw
        local c = a_r * T_wr
        local d = a_w * T_ww

        push!(eqs, det_var ~ (1 - d^2) - 2 * b * c * (1 + d))

        # Q = M^{-1} * S_0 (total incoming after all reflections)
        local Q_road = ((1 - d^2) * S_r + b * (1 + d) * S_sw + b * (1 + d) * S_sh) / det_var
        local Q_sw = (c * (1 + d) * S_r + (1 - b * c) * S_sw + (d + b * c) * S_sh) / det_var
        local Q_sh = (c * (1 + d) * S_r + (d + b * c) * S_sw + (1 - b * c) * S_sh) / det_var

        # Absorbed = (1 - albedo) * total_incoming
        push!(eqs, net_road ~ (1 - a_r) * Q_road)
        push!(eqs, net_sw ~ (1 - a_w) * Q_sw)
        push!(eqs, net_sh ~ (1 - a_w) * Q_sh)

        # Reflected to sky = incident - absorbed (energy conservation)
        push!(
            eqs, refl_sky ~ (S_r + (S_sw + S_sh) * H_W) -
                (net_road + (net_sw + net_sh) * H_W)
        )
    end

    # Sum absorbed solar by surface
    push!(
        eqs, S_net_road ~ S_net_road_dir_vis + S_net_road_dir_nir +
            S_net_road_dif_vis + S_net_road_dif_nir
    )
    push!(
        eqs, S_net_sunwall ~ S_net_sunwall_dir_vis + S_net_sunwall_dir_nir +
            S_net_sunwall_dif_vis + S_net_sunwall_dif_nir
    )
    push!(
        eqs, S_net_shdwall ~ S_net_shdwall_dir_vis + S_net_shdwall_dir_nir +
            S_net_shdwall_dif_vis + S_net_shdwall_dif_nir
    )

    # Canyon net absorbed solar per ground area
    push!(eqs, S_net_uc ~ S_net_road + (S_net_sunwall + S_net_shdwall) * H_W)

    # Canyon albedos (Eqs. 2.93-2.94): α = reflected / incident per ground area
    push!(
        eqs, α_uc_dir_vis ~ ifelse(
            S_road_dir_vis + S_sunwall_dir_vis * H_W > zero_wm2,
            S_refl_sky_dir_vis / (S_road_dir_vis + S_sunwall_dir_vis * H_W), 0.0
        )
    )
    push!(
        eqs, α_uc_dir_nir ~ ifelse(
            S_road_dir_nir + S_sunwall_dir_nir * H_W > zero_wm2,
            S_refl_sky_dir_nir / (S_road_dir_nir + S_sunwall_dir_nir * H_W), 0.0
        )
    )
    push!(
        eqs, α_uc_dif_vis ~ ifelse(
            S_road_dif_vis + 2 * S_wall_dif_vis * H_W > zero_wm2,
            S_refl_sky_dif_vis / (S_road_dif_vis + 2 * S_wall_dif_vis * H_W), 0.0
        )
    )
    push!(
        eqs, α_uc_dif_nir ~ ifelse(
            S_road_dif_nir + 2 * S_wall_dif_nir * H_W > zero_wm2,
            S_refl_sky_dif_nir / (S_road_dif_nir + 2 * S_wall_dif_nir * H_W), 0.0
        )
    )

    # Total absorbed solar (Eq. 2.95)
    push!(eqs, S_net_total ~ W_roof * S_net_roof + (1 - W_roof) * S_net_uc)

    # ===== Longwave Radiation (Sections 2.6-2.7) =====

    # Roof LW (Eqs. 2.101-2.102)
    push!(eqs, L_up_roof ~ ε_roof * σ_SB * T_roof^4 + (1 - ε_roof) * L_atm_down)
    push!(eqs, L_net_roof ~ L_up_roof - L_atm_down)

    # Emitted LW (Eqs. 2.109-2.113)
    push!(eqs, L_emit_road ~ ε_road * σ_SB * T_road^4)
    push!(eqs, L_emit_sunwall ~ ε_wall * σ_SB * T_sunwall^4)
    push!(eqs, L_emit_shdwall ~ ε_wall * σ_SB * T_shdwall^4)

    # Canyon LW multi-reflection (closed-form, Section 2.7)
    # Same matrix structure as solar, with reflectivity (1-ε) as "albedo".
    # Source vector includes atmospheric LW and emitted LW from other surfaces.
    begin
        local a_r_lw = 1 - ε_road
        local a_w_lw = 1 - ε_wall
        local b_lw = a_w_lw * T_rw
        local c_lw = a_r_lw * T_wr
        local d_lw = a_w_lw * T_ww

        push!(eqs, det_lw ~ (1 - d_lw^2) - 2 * b_lw * c_lw * (1 + d_lw))

        # Initial incoming from atmosphere + emitted from other canyon surfaces
        local S0_road = L_atm_down * Ψ_sky_road +
            T_rw * (L_emit_sunwall + L_emit_shdwall)
        local S0_sw = L_atm_down * Ψ_sky_wall +
            T_wr * L_emit_road + T_ww * L_emit_shdwall
        local S0_sh = L_atm_down * Ψ_sky_wall +
            T_wr * L_emit_road + T_ww * L_emit_sunwall

        # Q = M^{-1} * S_0 (total incoming after all reflections)
        local Q_road_lw = (
            (1 - d_lw^2) * S0_road +
                b_lw * (1 + d_lw) * S0_sw +
                b_lw * (1 + d_lw) * S0_sh
        ) / det_lw
        local Q_sw_lw = (
            c_lw * (1 + d_lw) * S0_road +
                (1 - b_lw * c_lw) * S0_sw +
                (d_lw + b_lw * c_lw) * S0_sh
        ) / det_lw
        local Q_sh_lw = (
            c_lw * (1 + d_lw) * S0_road +
                (d_lw + b_lw * c_lw) * S0_sw +
                (1 - b_lw * c_lw) * S0_sh
        ) / det_lw

        # Net LW = emitted - absorbed = L_emit - ε * Q
        push!(eqs, L_net_road ~ L_emit_road - ε_road * Q_road_lw)
        push!(eqs, L_net_sunwall ~ L_emit_sunwall - ε_wall * Q_sw_lw)
        push!(eqs, L_net_shdwall ~ L_emit_shdwall - ε_wall * Q_sh_lw)
    end

    # Canyon net LW per ground area (Eq. 2.172)
    push!(eqs, L_net_uc ~ L_net_road + (L_net_sunwall + L_net_shdwall) * H_W)

    # Total net LW (Eq. 2.175)
    push!(eqs, L_net_total ~ W_roof * L_net_roof + (1 - W_roof) * L_net_uc)

    return System(eqs, t; name)
end
