using Documenter, DocumenterVitepress, MaterialPointSolver, WGLMakie

makedocs(
    modules = [MaterialPointSolver],
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/LandslideSIM/MaterialPointSolver.jl",
        devbranch = "main",
        devurl = "dev"
    ),
    sitename = "MaterialPointSolver.jl",
    authors = "Zenan Huo",
    pages = [
        "Home" => "index.md",
        "Introduction" => Any[
            "introduction/getstarted.md",
            "introduction/why.md",
            "introduction/ecosystem.md",
        ],
        "DataType System" => Any[
            "datatype/backgroundgrid.md",
        ],
        "Advanced Topics" => Any[
            "advancedtopics/plugin.md",
        ],
        "Examples" => Any[
            "examples/example1.md",
        ],
        "Others" => Any[
            "others/gpu.md",
        ]
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/LandslideSIM/MaterialPointSolver.jl",
    target = "build",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)