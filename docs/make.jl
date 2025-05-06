using Atlans
using Documenter
using DocumenterMarkdown


makedocs(;
    sitename = "Atlans.jl",
    repo = "https://gitlab.com/deltares/subsidence/atlans.jl.git",
    format = Markdown(),
    authors = "Deltares and contributors",
    modules = [Atlans],
)
