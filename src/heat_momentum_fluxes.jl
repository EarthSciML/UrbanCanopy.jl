export HeatMomentumFluxes

"""
    HeatMomentumFluxes(; name=:HeatMomentumFluxes)

Urban heat and momentum fluxes for the Community Land Model Urban (CLMU)
parameterization, following Chapter 3 of Oleson et al. (2010).

Computes the Monin-Obukhov similarity theory parameters, roughness length and
displacement height, wind speed in the urban canyon, aerodynamic and surface
resistances, sensible and latent heat fluxes for the five urban surfaces (roof,
sunlit wall, shaded wall, pervious road, impervious road), total momentum fluxes,
and the UCL air temperature and specific humidity.

The Monin-Obukhov stability parameter ζ is taken as an input parameter, since
in the original CLM implementation it is computed through an iterative procedure
(Section 3.2.3). Given ζ, the system computes all stability corrections, friction
velocity, aerodynamic resistances, and surface fluxes as a directed algebraic system
without circular dependencies.

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp. Chapter 3: Heat and Momentum Fluxes (pp. 61-89).
"""
@component function HeatMomentumFluxes(; name = :HeatMomentumFluxes)

    # ===== Physical Constants (Table 1.4) =====
    @constants begin
        C_p = 1004.64, [description = "Specific heat capacity of dry air (Table 1.4)", unit = u"J/(kg*K)"]
        g_acc = 9.80616, [description = "Gravitational acceleration (Table 1.4)", unit = u"m/s^2"]
        k_vk = 0.4, [description = "von Karman constant (Table 1.4) (dimensionless)"]
    end

    # Unit reference constants for non-integer exponents and dimensional comparisons
    @constants begin
        one_ms = 1.0, [description = "Unit velocity", unit = u"m/s"]
        one_K = 1.0, [description = "Unit temperature", unit = u"K"]
        one_kgkg = 1.0, [description = "Unit specific humidity", unit = u"kg/kg"]
        zero_ms = 0.0, [description = "Zero velocity", unit = u"m/s"]
        π_val = Float64(π), [description = "Pi (dimensionless)"]
    end

    # Roughness length empirical constants
    @constants begin
        α_rough = 4.43, [description = "Empirical coefficient for displacement height (Eq. 3.55) (dimensionless)"]
        B_drag = 1.0, [description = "Drag correction coefficient (Eq. 3.57) (dimensionless)"]
        C_D = 1.2, [description = "Depth-integrated mean drag coefficient (Eq. 3.57) (dimensionless)"]
    end

    # Convective boundary layer constants
    @constants begin
        β_conv = 1.0, [description = "Convective velocity coefficient (Eq. 3.29) (dimensionless)"]
        z_i = 1000.0, [description = "Convective boundary layer height (Eq. 3.30)", unit = u"m"]
        zero_61 = 0.61, [description = "Virtual temperature coefficient (dimensionless)"]
        one_m = 1.0, [description = "Unit length for non-dimensionalization", unit = u"m"]
    end

    # Surface resistance coefficients (Eq. 3.68, Rowley et al. 1930)
    @constants begin
        h_c_free = 11.8, [description = "Free convection coefficient in surface resistance (Eq. 3.68)", unit = u"W/(m^2*K)"]
        h_c_forced = 4.2, [description = "Forced convection coefficient in surface resistance (Eq. 3.68)", unit = u"W*s/(m^3*K)"]
    end

    # ===== Parameters =====

    # Canyon geometry (Table 1.3)
    @parameters begin
        H_W, [description = "Canyon height to width ratio H/W (dimensionless)"]
        H_canyon, [description = "Canyon (building/roof) height", unit = u"m"]
        W_roof, [description = "Roof fraction of urban surface (dimensionless)"]
        f_prvrd, [description = "Fraction of road that is pervious (dimensionless)"]
        H_w_est, [description = "Height at which canyon wind speed is estimated (Table 1.3)", unit = u"m"]
    end

    # Atmospheric forcing (Table 1.1)
    @parameters begin
        u_atm, [description = "Zonal wind at reference height", unit = u"m/s"]
        v_atm, [description = "Meridional wind at reference height", unit = u"m/s"]
        θ_atm, [description = "Atmospheric potential temperature (Eq. 3.8)", unit = u"K"]
        q_atm, [description = "Atmospheric specific humidity", unit = u"kg/kg"]
        ρ_atm, [description = "Atmospheric air density (Eq. 3.9)", unit = u"kg/m^3"]
    end

    # Reference heights
    @parameters begin
        z_atm_m, [description = "Reference height for wind (momentum)", unit = u"m"]
        z_atm_h, [description = "Reference height for temperature (sensible heat)", unit = u"m"]
        z_atm_w, [description = "Reference height for humidity (water vapor)", unit = u"m"]
    end

    # Surface temperatures (from Chapter 4 solution, previous timestep)
    @parameters begin
        T_g_roof, [description = "Roof surface temperature", unit = u"K"]
        T_g_prvrd, [description = "Pervious road surface temperature", unit = u"K"]
        T_g_imprvrd, [description = "Impervious road surface temperature", unit = u"K"]
        T_g_sunwall, [description = "Sunlit wall surface temperature", unit = u"K"]
        T_g_shdwall, [description = "Shaded wall surface temperature", unit = u"K"]
    end

    # Surface specific humidities
    @parameters begin
        q_g_roof, [description = "Roof surface saturated specific humidity", unit = u"kg/kg"]
        q_g_prvrd, [description = "Pervious road surface specific humidity", unit = u"kg/kg"]
        q_g_imprvrd, [description = "Impervious road surface saturated specific humidity", unit = u"kg/kg"]
    end

    # Wetted fractions
    @parameters begin
        f_wet_roof, [description = "Wetted fraction of roof (Eqs. 3.82-3.83) (dimensionless)"]
        f_wet_imprvrd, [description = "Wetted fraction of impervious road (Eqs. 3.82-3.83) (dimensionless)"]
    end

    # Monin-Obukhov stability parameter (from iterative solution, Section 3.2.3)
    @parameters begin
        ζ_in, [description = "Monin-Obukhov stability parameter ζ = (z_atm_m - d)/L, from iterative solution (Eqs. 3.47-3.50) (dimensionless)"]
    end

    # ===== Roughness Length and Displacement Height (Section 3.2.1) =====

    @variables begin
        λ_p(t), [description = "Plan area index (Eq. 3.56) (dimensionless)"]
        λ_F(t), [description = "Frontal area index (Eq. 3.58) (dimensionless)"]
        d_canopy(t), [description = "Canopy displacement height (Eq. 3.55)", unit = u"m"]
        z_0m_canopy(t), [description = "Canopy roughness length for momentum (Eq. 3.57)", unit = u"m"]
    end

    # ===== MO-derived quantities =====

    @variables begin
        V_a(t), [description = "Effective wind speed (Eq. 3.25)", unit = u"m/s"]
        V_r(t), [description = "Reference level atmospheric wind (Eq. 3.63)", unit = u"m/s"]
        u_star(t), [description = "Friction velocity (Eq. 3.26)", unit = u"m/s"]
        L_MO(t), [description = "Monin-Obukhov length (Eq. 3.17)", unit = u"m"]
    end

    # ===== Stability corrections =====

    @variables begin
        ψ_m_atm(t), [description = "Momentum stability correction at ζ (dimensionless)"]
        ψ_m_0(t), [description = "Momentum stability correction at ζ_0m (dimensionless)"]
        ψ_h_atm(t), [description = "Heat stability correction at ζ_h (dimensionless)"]
        ψ_h_0(t), [description = "Heat stability correction at ζ_0h (dimensionless)"]
        ψ_w_atm(t), [description = "Moisture stability correction at ζ_w (dimensionless)"]
        ψ_w_0(t), [description = "Moisture stability correction at ζ_0w (dimensionless)"]
        ζ_0m(t), [description = "Stability parameter at z_0m (dimensionless)"]
        ζ_h(t), [description = "Stability parameter for heat at z_atm_h (dimensionless)"]
        ζ_0h(t), [description = "Stability parameter at z_0h (dimensionless)"]
        ζ_w(t), [description = "Stability parameter for moisture at z_atm_w (dimensionless)"]
        ζ_0w(t), [description = "Stability parameter at z_0w (dimensionless)"]
    end

    # ===== Resistances =====

    @variables begin
        r_am(t), [description = "Aerodynamic resistance for momentum (Eq. 3.65)", unit = u"s/m"]
        r_ah(t), [description = "Aerodynamic resistance for sensible heat (Eq. 3.66)", unit = u"s/m"]
        r_aw(t), [description = "Aerodynamic resistance for water vapor (Eq. 3.67)", unit = u"s/m"]
        r_s_u(t), [description = "Surface resistance for all urban surfaces (Eq. 3.68)", unit = u"s/m"]
    end

    # ===== Canyon Wind Speed =====

    @variables begin
        U_can(t), [description = "Canyon horizontal wind speed (Eqs. 3.60-3.62)", unit = u"m/s"]
        U_ac(t), [description = "Canyon wind speed (Eq. 3.59)", unit = u"m/s"]
    end

    # ===== UCL Temperature and Humidity =====

    @variables begin
        T_ac(t), [description = "UCL air temperature (Eq. 3.75)", unit = u"K"]
        q_ac(t), [description = "UCL air specific humidity (Eq. 3.93)", unit = u"kg/kg"]
    end

    # ===== Heat Fluxes =====

    @variables begin
        H_roof(t), [description = "Sensible heat flux from roof (Eq. 3.69)", unit = u"W/m^2"]
        H_prvrd(t), [description = "Sensible heat flux from pervious road (Eq. 3.70)", unit = u"W/m^2"]
        H_imprvrd(t), [description = "Sensible heat flux from impervious road (Eq. 3.71)", unit = u"W/m^2"]
        H_sunwall(t), [description = "Sensible heat flux from sunlit wall (Eq. 3.72)", unit = u"W/m^2"]
        H_shdwall(t), [description = "Sensible heat flux from shaded wall (Eq. 3.73)", unit = u"W/m^2"]
        H_total(t), [description = "Total sensible heat flux (Eq. 3.74)", unit = u"W/m^2"]
    end

    @variables begin
        E_roof(t), [description = "Water vapor flux from roof (Eq. 3.76)", unit = u"kg/(m^2*s)"]
        E_prvrd(t), [description = "Water vapor flux from pervious road (Eq. 3.77)", unit = u"kg/(m^2*s)"]
        E_imprvrd(t), [description = "Water vapor flux from impervious road (Eq. 3.78)", unit = u"kg/(m^2*s)"]
        E_total(t), [description = "Total water vapor flux (Eq. 3.81)", unit = u"kg/(m^2*s)"]
    end

    @variables begin
        τ_x(t), [description = "Zonal momentum flux (Eq. 3.6)", unit = u"kg/(m*s^2)"]
        τ_y(t), [description = "Meridional momentum flux (Eq. 3.7)", unit = u"kg/(m*s^2)"]
    end

    # Conductances for UCL solution
    @variables begin
        c_a_h(t), [description = "Sensible heat conductance UCL to atmosphere (1/r_ah)", unit = u"m/s"]
        c_a_w(t), [description = "Latent heat conductance UCL to atmosphere (1/r_aw)", unit = u"m/s"]
    end

    # Convective velocity (Eqs. 3.29-3.30)
    @variables begin
        U_c(t), [description = "Convective velocity (Eq. 3.29)", unit = u"m/s"]
        θ_v_atm(t), [description = "Virtual potential temperature at reference height (dimensionless)", unit = u"K"]
    end

    # ===== Partial Derivatives of Fluxes (Section 3.2.4, Eqs. 3.95-3.104) =====

    @variables begin
        dH_roof_dT(t), [description = "∂H_roof/∂T_g,roof (Eq. 3.95)", unit = u"W/(m^2*K)"]
        dH_prvrd_dT(t), [description = "∂H_prvrd/∂T_g,prvrd (Eq. 3.96)", unit = u"W/(m^2*K)"]
        dH_imprvrd_dT(t), [description = "∂H_imprvrd/∂T_g,imprvrd (Eq. 3.97)", unit = u"W/(m^2*K)"]
        dH_sunwall_dT(t), [description = "∂H_sunwall/∂T_g,sunwall (Eq. 3.98)", unit = u"W/(m^2*K)"]
        dH_shdwall_dT(t), [description = "∂H_shdwall/∂T_g,shdwall (Eq. 3.99)", unit = u"W/(m^2*K)"]
    end

    eqs = Equation[]

    # ===== Roughness Length and Displacement Height (Section 3.2.1, Eqs. 3.55-3.58) =====

    push!(eqs, λ_p ~ H_W / (H_W + 1))                                            # Eq. 3.56
    push!(eqs, λ_F ~ (1 - λ_p) * H_W)                                            # Eq. 3.58 (with B_S/B_L = λ_p)
    push!(eqs, d_canopy ~ H_canyon * (1 + (α_rough)^(-λ_p) * (λ_p - 1)))         # Eq. 3.55
    push!(
        eqs, z_0m_canopy ~ H_canyon * (1 - d_canopy / H_canyon) *               # Eq. 3.57
            exp(-((0.5 * B_drag * C_D / k_vk^2 * (1 - d_canopy / H_canyon) * λ_F))^(-0.5))
    )

    # ===== Wind Speed =====
    push!(eqs, V_r ~ max(sqrt(u_atm^2 + v_atm^2), one_ms))                       # Eq. 3.63

    # ===== Stability corrections from ζ_in =====
    # ζ_in = (z_atm_m - d_canopy) / L_MO, so L_MO = (z_atm_m - d_canopy) / ζ_in
    push!(eqs, L_MO ~ (z_atm_m - d_canopy) / ζ_in)                               # Eq. 3.50

    # ζ at different levels: scale by height ratios
    push!(eqs, ζ_0m ~ z_0m_canopy / L_MO)
    push!(eqs, ζ_h ~ (z_atm_h - d_canopy) / L_MO)
    push!(eqs, ζ_0h ~ z_0m_canopy / L_MO)  # z_0h = z_0m for canopy
    push!(eqs, ζ_w ~ (z_atm_w - d_canopy) / L_MO)
    push!(eqs, ζ_0w ~ z_0m_canopy / L_MO)  # z_0w = z_0m for canopy

    # ψ_m stability correction function (Eqs. 3.31, 3.33-3.37)
    function _ψ_m(ζ_val)
        x = (1 - 16 * ζ_val)^(1 / 4)
        ψ_unstable = 2 * log((1 + x) / 2) + log((1 + x^2) / 2) - 2 * atan(x) + π_val / 2  # Eq. 3.37
        ψ_stable = -5 * ζ_val

        # Very unstable transition at ζ = -1.574
        ζ_m_trans = -1.574
        x_trans = (1 - 16 * ζ_m_trans)^(1 / 4)
        ψ_m_at_trans = 2 * log((1 + x_trans) / 2) + log((1 + x_trans^2) / 2) - 2 * atan(x_trans) + π_val / 2
        ψ_very_unstable = ψ_m_at_trans + log(ζ_val / ζ_m_trans) + 1.14 * ((-ζ_val)^(1 / 3) - (-ζ_m_trans)^(1 / 3))  # From Eq. 3.33

        # Very stable: ψ_m(ζ>1) = -4ln(ζ) - ζ - 4 (derived from Eq. 3.36)
        ψ_very_stable = -4 * log(ζ_val) - ζ_val - 4

        return ifelse(
            ζ_val < ζ_m_trans, ψ_very_unstable,
            ifelse(
                ζ_val < 0, ψ_unstable,
                ifelse(ζ_val ≤ 1, ψ_stable, ψ_very_stable)
            )
        )
    end

    # ψ_h = ψ_w stability correction function (Eqs. 3.32, 3.38-3.46)
    function _ψ_h(ζ_val)
        x = (1 - 16 * ζ_val)^(1 / 4)
        ψ_unstable = 2 * log((1 + x^2) / 2)  # Eq. 3.46
        ψ_stable = -5 * ζ_val

        # Very unstable transition at ζ = -0.465
        ζ_h_trans = -0.465
        x_trans = (1 - 16 * ζ_h_trans)^(1 / 4)
        ψ_h_at_trans = 2 * log((1 + x_trans^2) / 2)
        ψ_very_unstable = ψ_h_at_trans + log(ζ_val / ζ_h_trans) + 0.8 * ((-ζ_h_trans)^(-1 / 3) - (-ζ_val)^(-1 / 3))  # From Eq. 3.38

        # Very stable: ψ_h(ζ>1) = -4ln(ζ) - ζ - 4 (derived from Eq. 3.41)
        ψ_very_stable = -4 * log(ζ_val) - ζ_val - 4

        return ifelse(
            ζ_val < ζ_h_trans, ψ_very_unstable,
            ifelse(
                ζ_val < 0, ψ_unstable,
                ifelse(ζ_val ≤ 1, ψ_stable, ψ_very_stable)
            )
        )
    end

    push!(eqs, ψ_m_atm ~ _ψ_m(ζ_in))
    push!(eqs, ψ_m_0 ~ _ψ_m(ζ_0m))
    push!(eqs, ψ_h_atm ~ _ψ_h(ζ_h))
    push!(eqs, ψ_h_0 ~ _ψ_h(ζ_0h))
    push!(eqs, ψ_w_atm ~ _ψ_h(ζ_w))   # φ_w = φ_h, so ψ_w = ψ_h (Eq. 3.32)
    push!(eqs, ψ_w_0 ~ _ψ_h(ζ_0w))

    # ===== Virtual Potential Temperature =====
    push!(eqs, θ_v_atm ~ θ_atm * (1 + zero_61 * q_atm))                          # θ_v,atm = θ̄_atm(1 + 0.61*q_atm)

    # ===== Convective Velocity (Eqs. 3.29-3.30) =====
    # U_c = 0 for stable (ζ ≥ 0), U_c = β*w* for unstable (ζ < 0)
    # w* = (-g*u*_star*θ_v*_star*z_i / θ_v,atm)^{1/3} (Eq. 3.30)
    # Since ζ_in is given as input, we can compute L and from it θ_v*:
    #   L = u*²θ_v,atm / (k*g*θ_v*) (Eq. 3.17 rearranged), so θ_v* = u*²θ_v,atm/(k*g*L)
    # Then w* = (-g * u* * θ_v* * z_i / θ_v,atm)^{1/3}
    #         = (-g * u* * u*² θ_v,atm / (k*g*L) * z_i / θ_v,atm)^{1/3}
    #         = (-u*³ * z_i / (k * L))^{1/3}
    # For unstable conditions L < 0, so -z_i/(k*L) > 0 and the cube root is real.
    # We use the simplified approach: compute U_c after u_star is known.
    # To break the V_a → u_star → U_c → V_a loop, we use V_r (without U_c) to get
    # the initial u_star, then compute U_c from that u_star and include it in V_a.
    # This is consistent with the iterative approach described in Section 3.2.3.
    begin
        # Compute u_star from V_r first (no convective velocity) for the U_c calculation
        local u_star_init = max(
            k_vk * V_r /
                (log((z_atm_m - d_canopy) / z_0m_canopy) - ψ_m_atm + ψ_m_0), 0.01 * one_ms
        )

        # w* = u* * (-z_i * ζ_in / (k * (z_atm_m - d)))^{1/3}  (Eq. 3.30 simplified)
        # The ratio X = z_i / (z_atm_m - d) * (-ζ_in) / k is dimensionless.
        # Non-dimensionalize z_i and (z_atm_m - d) with one_m to help the unit checker.
        local X_dimless = (z_i / one_m) / ((z_atm_m - d_canopy) / one_m) * (-ζ_in) / k_vk
        local X_pos = max(X_dimless, 0.0)

        # w* = u* * X^{1/3} for unstable (ζ < 0), 0 for stable
        local w_star = ifelse(
            ζ_in < 0,
            u_star_init * X_pos^(1 / 3),
            zero_ms
        )

        push!(eqs, U_c ~ β_conv * w_star)                                         # Eq. 3.29
    end

    # ===== Effective Wind Speed (Eq. 3.25) =====
    push!(eqs, V_a ~ max(sqrt(u_atm^2 + v_atm^2 + U_c^2), one_ms))              # Eq. 3.25

    # ===== Friction Velocity (Eq. 3.26) =====
    push!(
        eqs, u_star ~ max(
            k_vk * V_a /
                (log((z_atm_m - d_canopy) / z_0m_canopy) - ψ_m_atm + ψ_m_0), 0.01 * one_ms
        )
    )

    # ===== Aerodynamic Resistances (Eqs. 3.65-3.67) =====
    push!(eqs, r_am ~ V_a / u_star^2)                                             # Eq. 3.65

    push!(
        eqs, r_ah ~ (1 / (k_vk^2 * V_a)) *                                     # Eq. 3.66
            (log((z_atm_m - d_canopy) / z_0m_canopy) - ψ_m_atm + ψ_m_0) *
            (log((z_atm_h - d_canopy) / z_0m_canopy) - ψ_h_atm + ψ_h_0)
    )

    push!(
        eqs, r_aw ~ (1 / (k_vk^2 * V_a)) *                                     # Eq. 3.67
            (log((z_atm_m - d_canopy) / z_0m_canopy) - ψ_m_atm + ψ_m_0) *
            (log((z_atm_w - d_canopy) / z_0m_canopy) - ψ_w_atm + ψ_w_0)
    )

    # ===== Canyon Wind Speed (Section 3.2.2, Eqs. 3.59-3.62) =====
    begin
        local log_ratio = log(max((H_canyon - d_canopy) / z_0m_canopy, 1.0 + 1.0e-10)) /
            log(max((z_atm_m - d_canopy) / z_0m_canopy, 1.0 + 1.0e-10))
        local exp_decay = exp(-0.5 * H_W * (1 - H_w_est / H_canyon))

        local U_skim = V_r * (2 / π_val) * log_ratio * exp_decay                  # Eq. 3.60
        local U_isol = V_r * log_ratio * exp_decay                                 # Eq. 3.61
        local U_wake = V_r * (1 + 2 * (2 / π_val - 1) * (H_W - 0.5)) * log_ratio * exp_decay  # Eq. 3.62

        push!(
            eqs, U_can ~ ifelse(
                H_W ≥ 1, U_skim,
                ifelse(H_W ≥ 0.5, U_wake, U_isol)
            )
        )
    end

    # W_can = u_* (Masson 2000); U_ac = √(U_can² + W_can²) (Eq. 3.59)
    push!(eqs, U_ac ~ sqrt(U_can^2 + u_star^2))

    # ===== Surface Resistance (Eq. 3.68) =====
    push!(eqs, r_s_u ~ ρ_atm * C_p / (h_c_free + h_c_forced * U_ac))

    # ===== UCL Air Temperature (Eq. 3.75) =====
    # Conductances: c = weight / r_s_u for each surface
    push!(eqs, c_a_h ~ 1 / r_ah)
    push!(eqs, c_a_w ~ 1 / r_aw)

    begin
        local c_roof = W_roof / r_s_u
        local c_prvrd = (1 - W_roof) * f_prvrd / r_s_u
        local c_imprvrd = (1 - W_roof) * (1 - f_prvrd) / r_s_u
        local c_sunwall = (1 - W_roof) * H_W / r_s_u
        local c_shdwall = (1 - W_roof) * H_W / r_s_u

        push!(
            eqs, T_ac ~ (
                c_a_h * θ_atm + c_roof * T_g_roof + c_prvrd * T_g_prvrd +  # Eq. 3.75
                    c_imprvrd * T_g_imprvrd + c_sunwall * T_g_sunwall + c_shdwall * T_g_shdwall
            ) /
                (c_a_h + c_roof + c_prvrd + c_imprvrd + c_sunwall + c_shdwall)
        )
    end

    # ===== UCL Air Specific Humidity (Eq. 3.93) =====
    begin
        local cw_roof = W_roof / r_s_u
        local cw_prvrd = (1 - W_roof) * f_prvrd / r_s_u
        local cw_imprvrd = (1 - W_roof) * (1 - f_prvrd) / r_s_u

        push!(
            eqs, q_ac ~ (
                c_a_w * q_atm + cw_roof * f_wet_roof * q_g_roof +            # Eq. 3.93
                    cw_prvrd * q_g_prvrd + cw_imprvrd * f_wet_imprvrd * q_g_imprvrd
            ) /
                (c_a_w + f_wet_roof * cw_roof + cw_prvrd + f_wet_imprvrd * cw_imprvrd)
        )
    end

    # ===== Sensible Heat Fluxes (Eqs. 3.69-3.74) =====
    push!(eqs, H_roof ~ -ρ_atm * C_p * (T_ac - T_g_roof) / r_s_u)                # Eq. 3.69
    push!(eqs, H_prvrd ~ -ρ_atm * C_p * (T_ac - T_g_prvrd) / r_s_u)              # Eq. 3.70
    push!(eqs, H_imprvrd ~ -ρ_atm * C_p * (T_ac - T_g_imprvrd) / r_s_u)          # Eq. 3.71
    push!(eqs, H_sunwall ~ -ρ_atm * C_p * (T_ac - T_g_sunwall) / r_s_u)          # Eq. 3.72
    push!(eqs, H_shdwall ~ -ρ_atm * C_p * (T_ac - T_g_shdwall) / r_s_u)          # Eq. 3.73

    push!(
        eqs, H_total ~ W_roof * H_roof + (1 - W_roof) *                         # Eq. 3.74
            (
            f_prvrd * H_prvrd + (1 - f_prvrd) * H_imprvrd +
                H_W * H_sunwall + H_W * H_shdwall
        )
    )

    # ===== Water Vapor Fluxes (Eqs. 3.76-3.81) =====
    push!(eqs, E_roof ~ -ρ_atm * f_wet_roof * (q_ac - q_g_roof) / r_s_u)         # Eq. 3.76
    push!(eqs, E_prvrd ~ -ρ_atm * (q_ac - q_g_prvrd) / r_s_u)                    # Eq. 3.77
    push!(eqs, E_imprvrd ~ -ρ_atm * f_wet_imprvrd * (q_ac - q_g_imprvrd) / r_s_u)  # Eq. 3.78
    # E_sunwall = 0 (Eq. 3.79), E_shdwall = 0 (Eq. 3.80)

    push!(
        eqs, E_total ~ W_roof * E_roof + (1 - W_roof) *                         # Eq. 3.81
            (f_prvrd * E_prvrd + (1 - f_prvrd) * E_imprvrd)
    )

    # ===== Momentum Fluxes (Eqs. 3.6-3.7) =====
    push!(eqs, τ_x ~ -ρ_atm * u_atm / r_am)                                      # Eq. 3.6
    push!(eqs, τ_y ~ -ρ_atm * v_atm / r_am)                                      # Eq. 3.7

    # ===== Partial Derivatives of Sensible Heat Fluxes (Eqs. 3.95-3.99) =====
    # These are needed for the implicit soil temperature calculation (Chapter 4)
    # and for the flux update equations (Eqs. 3.106-3.107).
    #
    # The general form is (Eq. 3.95):
    # ∂H_g/∂T_g = ρ_atm * C_p * (c_a^h + Σ c_i/W_i [excluding surface g]) * c_g / W_g
    #             / (c_a^h + Σ c_i/W_i [all surfaces])
    # where c_i = conductance for surface i, W_i = area weight for surface i

    begin
        # Define conductances per unit area (c_i / W_i = 1 / r_s_u for all surfaces)
        # and area-weighted conductances c_i = W_i / r_s_u
        local cₕ_roof = W_roof / r_s_u
        local cₕ_prvrd = (1 - W_roof) * f_prvrd / r_s_u
        local cₕ_imprvrd = (1 - W_roof) * (1 - f_prvrd) / r_s_u
        local cₕ_sunwall = (1 - W_roof) * H_W / r_s_u
        local cₕ_shdwall = (1 - W_roof) * H_W / r_s_u

        # Sum of all conductances (denominator in Eq. 3.75)
        local cₕ_sum = c_a_h + cₕ_roof + cₕ_prvrd + cₕ_imprvrd + cₕ_sunwall + cₕ_shdwall

        # For each surface, the numerator of dH/dT includes all OTHER conductances + c_a_h
        push!(
            eqs, dH_roof_dT ~ ρ_atm * C_p *                                    # Eq. 3.95
                (c_a_h + cₕ_prvrd + cₕ_imprvrd + cₕ_sunwall + cₕ_shdwall) *
                cₕ_roof / (W_roof * cₕ_sum)
        )

        push!(
            eqs, dH_prvrd_dT ~ ρ_atm * C_p *                                   # Eq. 3.96
                (c_a_h + cₕ_roof + cₕ_imprvrd + cₕ_sunwall + cₕ_shdwall) *
                cₕ_prvrd / (((1 - W_roof) * f_prvrd) * cₕ_sum)
        )

        push!(
            eqs, dH_imprvrd_dT ~ ρ_atm * C_p *                                 # Eq. 3.97
                (c_a_h + cₕ_roof + cₕ_prvrd + cₕ_sunwall + cₕ_shdwall) *
                cₕ_imprvrd / (((1 - W_roof) * (1 - f_prvrd)) * cₕ_sum)
        )

        push!(
            eqs, dH_sunwall_dT ~ ρ_atm * C_p *                                 # Eq. 3.98
                (c_a_h + cₕ_roof + cₕ_prvrd + cₕ_imprvrd + cₕ_shdwall) *
                cₕ_sunwall / (((1 - W_roof) * H_W) * cₕ_sum)
        )

        push!(
            eqs, dH_shdwall_dT ~ ρ_atm * C_p *                                 # Eq. 3.99
                (c_a_h + cₕ_roof + cₕ_prvrd + cₕ_imprvrd + cₕ_sunwall) *
                cₕ_shdwall / (((1 - W_roof) * H_W) * cₕ_sum)
        )
    end

    return System(eqs, t; name)
end
