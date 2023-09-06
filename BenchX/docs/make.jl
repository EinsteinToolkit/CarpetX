# Generate documentation with this command:
# (cd docs && julia make.jl)

push!(LOAD_PATH, "..")

using Documenter
using BenchX

makedocs(; sitename="BenchX", format=Documenter.HTML(), modules=[BenchX])

# deploydocs(; repo="github.com/eschnett/BenchX.jl.git", devbranch="main", push_preview=true)
