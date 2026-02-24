@testitem "GEOS-FP Coupling Structure" tags = [:geosfp_coupling] begin
    using UrbanCanopy, EarthSciData
    using Dates, ModelingToolkit, EarthSciMLBase

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15,
    )

    coupled = couple(
        UrbanCanopyModel(),
        GEOSFP("4x5", domain),
    )

    sys = convert(System, coupled)

    # The coupled system should compile without errors
    @test sys isa ModelingToolkit.AbstractSystem
end

@testitem "GEOS-FP Coupling Equations" tags = [:geosfp_coupling] begin
    using UrbanCanopy, EarthSciData
    using Dates, ModelingToolkit, EarthSciMLBase

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15,
    )

    coupled = couple(
        UrbanCanopyModel(),
        GEOSFP("4x5", domain),
    )

    sys = convert(System, coupled)
    obs = observed(sys)

    # Helper: find observed equation where the LHS symbol name contains `lhs_term`,
    # then check that its RHS string contains `rhs_term`.
    function has_coupling(obs, lhs_term, rhs_term)
        for eq in obs
            lhs_str = string(eq.lhs)
            rhs_str = string(eq.rhs)
            # Match top-level UCM equation (e.g. "UrbanCanopyModel₊S_atm(t)")
            # but not subsystem propagations (e.g. "UrbanCanopyModel₊offline₊S_atm(t)")
            if occursin("UrbanCanopyModel", lhs_str) && occursin(lhs_term, lhs_str) &&
                    count("₊", lhs_str) == 1 && occursin(rhs_term, rhs_str)
                return true
            end
        end
        return false
    end

    # Check that GEOS-FP solar radiation is coupled to UCM S_atm
    @test has_coupling(obs, "S_atm", "SWGDN")

    # Check that GEOS-FP 2m temperature is coupled to UCM T_atm
    @test has_coupling(obs, "T_atm", "T2M")

    # Check that GEOS-FP surface pressure is coupled to UCM P_atm
    @test has_coupling(obs, "P_atm", "PS")

    # Check that GEOS-FP 2m specific humidity is coupled to UCM q_atm
    @test has_coupling(obs, "q_atm", "QV2M")

    # Check that GEOS-FP wind components are coupled to UCM W_atm (via sqrt(U10M^2 + V10M^2))
    @test has_coupling(obs, "W_atm", "U10M") || has_coupling(obs, "W_atm", "V10M")

    # Check that GEOS-FP precipitation is coupled to UCM P_precip
    @test has_coupling(obs, "P_precip", "PRECTOT")

    # Check that Monin-Obukhov stability (ζ_in) coupling exists (uses USTAR and HFLUX)
    @test has_coupling(obs, "ζ_in", "USTAR") || has_coupling(obs, "ζ_in", "HFLUX")

    # Check that solar zenith angle (μ_zen) coupling exists
    @test has_coupling(obs, "μ_zen", "acos")
end
