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
            "introduction/tips.md",
        ],
        "Interface" => Any[
            "interface/structs.md",
            "interface/backendagnostic.md",
            "interface/kernelfunctions.md",
            "interface/mpmprocess.md",
            "interface/workflow.md",
            "interface/imexport.md",
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
        "Useful Tools" => Any[
            "utils/debug.md",
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