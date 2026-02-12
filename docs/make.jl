using UrbanCanopy
using Documenter

makedocs(;
    modules = [UrbanCanopy],
    authors = "EarthSciML authors and contributors",
    repo = "https://github.com/EarthSciML/UrbanCanopy.jl/blob/{commit}{path}#{line}",
    sitename = "UrbanCanopy.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://urbancanopy.earthsci.dev",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "CLMU Introduction" => "clmu_introduction.md",
        "Offline Mode" => "offline_mode.md",
        "Albedos and Radiative Fluxes" => "albedos_radiative_fluxes.md",
        "Heat and Momentum Fluxes" => "heat_momentum_fluxes.md",
        "Roof, Wall, Road, Snow Temperatures" => "roof_wall_road_snow_temperatures.md",
        "Hydrology" => "hydrology.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo = "github.com/EarthSciML/UrbanCanopy.jl")
