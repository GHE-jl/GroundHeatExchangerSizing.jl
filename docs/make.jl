using GroundHeatExchangerSizing
using Documenter

# Make `using GroundHeatExchangerSizing` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    GroundHeatExchangerSizing,
    :DocTestSetup,
    :(using GroundHeatExchangerSizing);
    recursive = true,
)

makedocs(;
    modules = [GroundHeatExchangerSizing],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "GroundHeatExchangerSizing.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/GroundHeatExchangerSizing.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Sizing theory" => [
            "Overview" => "theory/overview.md",
            "Alternative ASHRAE equation" => "theory/alternative_equation.md",
            "Outlet transfer function" => "theory/outlet_transfer_function.md",
            "Optimisation" => "theory/optimization.md",
        ],
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/GroundHeatExchangerSizing.jl",
    devbranch = "main",
)
