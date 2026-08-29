using CSV
using CairoMakie
using Printf
using SixelTerm



orderstring(params) = Dict(2 => "", 4 => "-n4", 8 => "-n8")[params.order]
resolutionstring(params) = Dict(32 => "", 64 => "-64", 128 => "-128")[params.res]
simname(params) = "standing$(orderstring(params))$(resolutionstring(params))"

simdir(params) = "/Users/eschnett/simulations/$(simname(params))/output-0000/$(simname(params))"
iterstring(params) = @sprintf "it%06d" params.iter

paramset = [
    (; iter, res, order)
    for order in [2, 4, 8]
        for res in [32, 64, 128]
            for iter in (res÷32) * (32:32:32)
                ]
println("Parameters:")
for params in paramset
    println("    $params")
end

coords = Dict(
    params => CSV.File("$(simdir(params))/dgcoordinatesx-coords.$(iterstring(params)).x.tsv")
    for params in paramset)
u = Dict(
    params => CSV.File("$(simdir(params))/dgwavetoyx-u.$(iterstring(params)).x.tsv")
    for params in paramset)
f = Dict(
    params => CSV.File("$(simdir(params))/dgwavetoyx-f.$(iterstring(params)).x.tsv")
    for params in paramset)
u_rhs = Dict(
    params => CSV.File("$(simdir(params))/dgwavetoyx-u_rhs.$(iterstring(params)).x.tsv")
    for params in paramset)
u_error = Dict(
    params => CSV.File("$(simdir(params))/dgwavetoyx-u_error.$(iterstring(params)).x.tsv")
    for params in paramset)

fig = Figure(; size=(800, 450))
ax = Axis(fig[1, 1]; title="Solution", xlabel="x", ylabel="u")
lines = Dict(
    params => scatterlines!(ax, coords[params].coordx, u[params].u)
    for params in paramset)
Legend(fig[1, 2],
       [lines[params] for params in paramset],
       ["o$(params.order)/i$(params.iter)/n$(params.res)" for params in paramset])
display(fig)

foreach([2, 4, 8]) do order
    paramset′ = [params for params in paramset if params.order == order]
    fig = Figure(; size=(800, 450))
    ax = Axis(fig[1, 1]; title="Convergence test (order=$order)", xlabel="x", ylabel="scaled error")
    lines = Dict(
        let
            α = float(params.res) ^ (params.order-1)
            params => scatterlines!(ax, coords[params].coordx, α * u_error[params].u_error)
        end
        for params in paramset′)
    Legend(fig[1, 2],
           [lines[params] for params in paramset′],
           ["o$(params.order)/i$(params.iter)/n$(params.res)" for params in paramset′])
    display(fig)
end

nothing
