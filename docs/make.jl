import Pkg; Pkg.add("Documenter")
using Documenter,Juqst

push!(LOAD_PATH,"/Users/robin/Dropbox/Juqst.jl")

makedocs(sitename="Juqst Documentation",pages = [
        "index.md",
        "probability-algs.md",
        "chp.md",
        "open-systems.md",
        "channels.md"
    ])

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/rharper2/Juqst.jl.git",
    #target = "build",
)


