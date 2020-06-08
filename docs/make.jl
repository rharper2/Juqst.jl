import Pkg; Pkg.add("Documenter")
using Documenter,Juqst

makedocs(sitename="Juqst Documentation",pages = [
        "index.md",
        "chp.md",
        "open-systems.md",
        "channels.md",
        "probability-algs.md"
    ])
