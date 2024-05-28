using Documenter
using MvNormalCalibration
using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "assets", "refs.bib");
    style=:authoryear
)
makedocs(
    sitename = "MvNormalCalibration",
    format = Documenter.HTML(
        assets = [
            # I just like justified text better :)
            "./assets/justified-text.css",
        ]
    ),
    clean = true,
    checkdocs = :exports,
    modules = [MvNormalCalibration],
    repo = Remotes.GitHub("RomeoV", "MvNormalCalibration.jl");
    pages = [
        "index.md",
        "mathematical_background.md",
        "api_reference.md",
    ],
    plugins=[bib]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/RomeoV/MvNomalCalibration.jl.git";
)
