push!(LOAD_PATH,"../src/")
using Documenter, Zurcher

makedocs(modules = [Zurcher], sitename = "Zurcher.jl")

deploydocs(repo = "github.com/floswald/Zurcher.jl.git", devbranch = "main")
