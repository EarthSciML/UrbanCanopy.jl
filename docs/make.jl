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
        "Albedos and Radiative Fluxes" => "albedos_radiative_fluxes.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo = "github.com/EarthSciML/UrbanCanopy.jl")
