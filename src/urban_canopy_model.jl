export UrbanCanopyModel, UrbanCanopyModelCoupler

"""
    UrbanCanopyModelCoupler

Coupler type for the [`UrbanCanopyModel`](@ref) component, enabling coupling
with external data sources (e.g., GEOS-FP) via `EarthSciMLBase.couple2`.
"""
struct UrbanCanopyModelCoupler
    sys::Any
end

"""
    UrbanCanopyModel(; name=:UrbanCanopyModel)

Top-level composed diagnostic urban canopy model following the Community Land Model
Urban (CLMU) parameterization of Oleson et al. (2010).

This component wires together the four diagnostic CLMU subsystems:
- [`OfflineModeForcing`](@ref): Partitions raw meteorological forcing into derived quantities
- [`CLMUAtmosphere`](@ref): Computes atmospheric density, vapor pressure, and reference height
- [`UrbanRadiation`](@ref): Computes absorbed/reflected solar and longwave radiation
- [`HeatMomentumFluxes`](@ref): Computes sensible/latent heat fluxes and UCL air temperature

Surface temperatures are prescribed as parameters (diagnostic mode). The prognostic
thermal solvers (Chapter 4) and hydrology (Chapter 5) are left for future work.

Key outputs (accessed via subsystem dot notation):
- `flux.T_ac` [K]: Urban canopy layer air temperature
- `flux.q_ac` [kg/kg]: UCL specific humidity
- `flux.H_total` [W/m^2]: Total sensible heat flux
- `flux.E_total` [kg/(m^2*s)]: Total evaporation
- `rad.S_net_total` [W/m^2]: Net absorbed solar radiation
- `rad.L_net_total` [W/m^2]: Net longwave radiation

**Reference**: Oleson, K.W., G.B. Bonan, J.J. Feddema, M. Vertenstein, and E. Kluzek,
2010: Technical Description of an Urban Parameterization for the Community Land Model
(CLMU). NCAR Technical Note NCAR/TN-480+STR, National Center for Atmospheric Research,
Boulder, CO, 168 pp.
"""
@component function UrbanCanopyModel(; name = :UrbanCanopyModel)

    # ===== Top-Level Parameters (External Inputs) =====

    # Meteorological forcing (fed to OfflineModeForcing)
    @parameters begin
        S_atm, [description = "Total incident solar radiation (Eq. 6.1 output)", unit = u"W/m^2"]
        T_atm, [description = "Atmospheric temperature", unit = u"K"]
        P_atm, [description = "Atmospheric pressure", unit = u"Pa"]
        q_atm, [description = "Atmospheric specific humidity", unit = u"kg/kg"]
        W_atm, [description = "Atmospheric wind speed", unit = u"m/s"]
        P_precip, [description = "Total precipitation rate", unit = u"kg/(m^2*s)"]
    end

    # Solar geometry (fed to UrbanRadiation)
    @parameters begin
        μ_zen, [description = "Solar zenith angle", unit = u"rad"]
    end

    # Canyon geometry (shared by UrbanRadiation + HeatMomentumFluxes)
    @parameters begin
        H_W, [description = "Canyon height to width ratio H/W (dimensionless)"]
        W_roof, [description = "Roof fraction of urban surface (dimensionless)"]
        H_canyon, [description = "Canyon (building/roof) height", unit = u"m"]
        f_prvrd, [description = "Fraction of road that is pervious (dimensionless)"]
        H_w_est, [description = "Height at which canyon wind speed is estimated", unit = u"m"]
    end

    # Surface roughness (fed to CLMUAtmosphere)
    @parameters begin
        z_0, [description = "Roughness length", unit = u"m"]
        z_d, [description = "Displacement height", unit = u"m"]
    end

    # Surface albedos (fed to UrbanRadiation)
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

    # Surface emissivities (fed to UrbanRadiation)
    @parameters begin
        ε_roof, [description = "Roof emissivity (dimensionless)"]
        ε_road, [description = "Effective road emissivity (dimensionless)"]
        ε_wall, [description = "Wall emissivity (dimensionless)"]
    end

    # Surface temperatures (fed to UrbanRadiation and HeatMomentumFluxes)
    @parameters begin
        T_g_roof, [description = "Roof surface temperature", unit = u"K"]
        T_g_prvrd, [description = "Pervious road surface temperature", unit = u"K"]
        T_g_imprvrd, [description = "Impervious road surface temperature", unit = u"K"]
        T_g_sunwall, [description = "Sunlit wall surface temperature", unit = u"K"]
        T_g_shdwall, [description = "Shaded wall surface temperature", unit = u"K"]
    end

    # Surface humidities (fed to HeatMomentumFluxes)
    @parameters begin
        q_g_roof, [description = "Roof surface saturated specific humidity", unit = u"kg/kg"]
        q_g_prvrd, [description = "Pervious road surface specific humidity", unit = u"kg/kg"]
        q_g_imprvrd, [description = "Impervious road surface saturated specific humidity", unit = u"kg/kg"]
    end

    # Wetted fractions (fed to HeatMomentumFluxes)
    @parameters begin
        f_wet_roof, [description = "Wetted fraction of roof (dimensionless)"]
        f_wet_imprvrd, [description = "Wetted fraction of impervious road (dimensionless)"]
    end

    # Stability (fed to HeatMomentumFluxes)
    @parameters begin
        ζ_in, [description = "Monin-Obukhov stability parameter (dimensionless)"]
    end

    # ===== Step 1: Create Subsystems =====

    offline = OfflineModeForcing(; name = :offline)
    atmos = CLMUAtmosphere(; name = :atmos)
    rad = UrbanRadiation(; name = :rad)
    flux = HeatMomentumFluxes(; name = :flux)

    # ===== Step 2: Convert Subsystem Input Parameters to Variables =====

    # OfflineModeForcing: all 6 input params are driven by parent params
    offline = param_to_var(offline, :S_atm, :T_atm, :P_atm, :q_atm, :W_atm, :P_precip)

    # CLMUAtmosphere: only params actually used in its equations
    # (u_atm, v_atm, θ_atm, L_atm_down, q_rain, q_sno, S_atm_dir/dif_* are declared
    # but unused in CLMUAtmosphere equations, so MTK drops them from the system)
    atmos = param_to_var(
        atmos,
        :q_atm, :P_atm, :T_atm,
        :z_prime_atm, :z_0, :z_d
    )

    # UrbanRadiation: radiation inputs, geometry, surface temps, emissivities, albedos
    rad = param_to_var(
        rad,
        :S_atm_dir_vis, :S_atm_dir_nir, :S_atm_dif_vis, :S_atm_dif_nir,
        :L_atm_down, :μ_zen, :H_W, :W_roof,
        :α_roof_dir_vis, :α_roof_dir_nir, :α_roof_dif_vis, :α_roof_dif_nir,
        :α_road_dir_vis, :α_road_dir_nir, :α_road_dif_vis, :α_road_dif_nir,
        :α_wall_dir_vis, :α_wall_dir_nir, :α_wall_dif_vis, :α_wall_dif_nir,
        :ε_roof, :ε_road, :ε_wall,
        :T_roof, :T_road, :T_sunwall, :T_shdwall
    )

    # HeatMomentumFluxes: atmospheric inputs, geometry, surface conditions
    flux = param_to_var(
        flux,
        :u_atm, :v_atm, :θ_atm, :q_atm, :ρ_atm,
        :z_atm_m, :z_atm_h, :z_atm_w,
        :H_W, :H_canyon, :W_roof, :f_prvrd, :H_w_est,
        :T_g_roof, :T_g_prvrd, :T_g_imprvrd, :T_g_sunwall, :T_g_shdwall,
        :q_g_roof, :q_g_prvrd, :q_g_imprvrd,
        :f_wet_roof, :f_wet_imprvrd,
        :ζ_in
    )

    # ===== Step 3: Coupling Equations =====

    eqs = Equation[]

    # --- Parent params → OfflineModeForcing ---
    append!(
        eqs, [
            offline.S_atm ~ S_atm,
            offline.T_atm ~ T_atm,
            offline.P_atm ~ P_atm,
            offline.q_atm ~ q_atm,
            offline.W_atm ~ W_atm,
            offline.P_precip ~ P_precip,
        ]
    )

    # --- OfflineModeForcing outputs + parent params → CLMUAtmosphere ---
    # Only parameters that CLMUAtmosphere actually uses in its equations
    append!(
        eqs, [
            atmos.z_prime_atm ~ offline.z_prime_atm,
            atmos.T_atm ~ T_atm,
            atmos.P_atm ~ P_atm,
            atmos.q_atm ~ q_atm,
            atmos.z_0 ~ z_0,
            atmos.z_d ~ z_d,
        ]
    )

    # --- OfflineModeForcing outputs → UrbanRadiation ---
    append!(
        eqs, [
            rad.S_atm_dir_vis ~ offline.S_atm_dir_vis,
            rad.S_atm_dir_nir ~ offline.S_atm_dir_nir,
            rad.S_atm_dif_vis ~ offline.S_atm_dif_vis,
            rad.S_atm_dif_nir ~ offline.S_atm_dif_nir,
            rad.L_atm_down ~ offline.L_atm_down,
        ]
    )

    # --- Parent params → UrbanRadiation (geometry, surface props) ---
    append!(
        eqs, [
            rad.μ_zen ~ μ_zen,
            rad.H_W ~ H_W,
            rad.W_roof ~ W_roof,
            # Albedos
            rad.α_roof_dir_vis ~ α_roof_dir_vis,
            rad.α_roof_dir_nir ~ α_roof_dir_nir,
            rad.α_roof_dif_vis ~ α_roof_dif_vis,
            rad.α_roof_dif_nir ~ α_roof_dif_nir,
            rad.α_road_dir_vis ~ α_road_dir_vis,
            rad.α_road_dir_nir ~ α_road_dir_nir,
            rad.α_road_dif_vis ~ α_road_dif_vis,
            rad.α_road_dif_nir ~ α_road_dif_nir,
            rad.α_wall_dir_vis ~ α_wall_dir_vis,
            rad.α_wall_dir_nir ~ α_wall_dir_nir,
            rad.α_wall_dif_vis ~ α_wall_dif_vis,
            rad.α_wall_dif_nir ~ α_wall_dif_nir,
            # Emissivities
            rad.ε_roof ~ ε_roof,
            rad.ε_road ~ ε_road,
            rad.ε_wall ~ ε_wall,
            # Surface temperatures (UrbanRadiation uses T_roof, T_road, T_sunwall, T_shdwall)
            rad.T_roof ~ T_g_roof,
            rad.T_road ~ T_g_imprvrd,  # Road temp uses impervious as representative
            rad.T_sunwall ~ T_g_sunwall,
            rad.T_shdwall ~ T_g_shdwall,
        ]
    )

    # --- OfflineModeForcing + CLMUAtmosphere → HeatMomentumFluxes ---
    append!(
        eqs, [
            flux.u_atm ~ offline.u_atm,
            flux.v_atm ~ offline.v_atm,
            flux.θ_atm ~ offline.θ_atm,
            flux.q_atm ~ q_atm,
            flux.ρ_atm ~ atmos.ρ_atm,
            flux.z_atm_m ~ atmos.z_atm,
            flux.z_atm_h ~ atmos.z_atm,
            flux.z_atm_w ~ atmos.z_atm,
        ]
    )

    # --- Parent params → HeatMomentumFluxes (geometry, surface conditions) ---
    append!(
        eqs, [
            flux.H_W ~ H_W,
            flux.H_canyon ~ H_canyon,
            flux.W_roof ~ W_roof,
            flux.f_prvrd ~ f_prvrd,
            flux.H_w_est ~ H_w_est,
            flux.T_g_roof ~ T_g_roof,
            flux.T_g_prvrd ~ T_g_prvrd,
            flux.T_g_imprvrd ~ T_g_imprvrd,
            flux.T_g_sunwall ~ T_g_sunwall,
            flux.T_g_shdwall ~ T_g_shdwall,
            flux.q_g_roof ~ q_g_roof,
            flux.q_g_prvrd ~ q_g_prvrd,
            flux.q_g_imprvrd ~ q_g_imprvrd,
            flux.f_wet_roof ~ f_wet_roof,
            flux.f_wet_imprvrd ~ f_wet_imprvrd,
            flux.ζ_in ~ ζ_in,
        ]
    )

    # ===== Step 4: Return Composed System =====

    return System(
        eqs, t; systems = [offline, atmos, rad, flux], name,
        metadata = Dict(CoupleType => UrbanCanopyModelCoupler)
    )
end
