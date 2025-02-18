using Documenter, DocumenterTools, MaterialPointSolver

makedocs(
    modules = [MaterialPointSolver],
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        prettyurls = true
    ),
    clean = false,
    sitename = "MaterialPointSolver.jl",
    authors = "Zenan Huo",
    pages = [
        "Home" => "index.md",
        "Interface" => Any[
            "interface/structs.md",
            "interface/kernelfunctions.md",
            "interface/mpmprocess.md",
            "interface/backendagnostic.md",
            "interface/workflow.md",
            "interface/import.md",
            "interface/export.md"
        ],
        "Tutorials" => Any[
            "tutorials/example1.md",
            "tutorials/example2.md",
            "tutorials/example3.md"
        ],
        "Advanced Topics" => Any[
            "advancedtopics/customkernel.md",
            "advancedtopics/customconstitutive.md",
            "advancedtopics/custommpmprocess.md",
            "advancedtopics/customextension.md",
            "advancedtopics/customfriction.md"
        ],
        "utils.md"
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/LandslideSIM/MaterialPointSolver.jl.git",
    target = "build",
    branch = "gh-pages",
    versions=["stable" => "v^", "dev" => "dev"]
)