using CSV
using CairoMakie
using DataFrames
using SixelTerm

@info "Evaluating Crusher benchmarks"

frames = CSV.read("data/benchmark-z4c+grhydrox-crusher-hip-sim-hip.csv", DataFrame)

set_theme!(Theme(; fontsize=30))
fig = Figure(; resolution=(1280, 960))

ax = Axis(
    fig[1, 1];
    title="Performance on Crusher (GRaM-X, Z4c, CarpetX)",
    xlabel="# nodes",
    xscale=log10,
    xticks=2 .^ (0:7),
    ylabel="cell updates / sec / node",
)

ncellss = frames[!, "\$ncells_x"]
nodess = frames.nodes

legend = []

# weak scaling
for ncells_per_node in [maximum(ncellss .^ 3 ./ nodess)]
    @info "ncells_per_node=$ncells_per_node"
    frames1 = filter(row -> row["\$ncells_x"]^3 / row.nodes == ncells_per_node, frames)
    frames1 = sort(frames1, :nodes)
    nodes = frames1.nodes
    cell_updates_per_second_per_node = frames1.cell_updates_per_second ./ frames1.nodes
    obj = scatterlines!(ax, nodes, cell_updates_per_second_per_node; color=:red, linestyle=:dash, linewidth=5, markersize=20)
    push!(legend, obj => "ncells/node=$(round(Int, cbrt(ncells_per_node)))³ (weak scaling)")
    ylims!(; low=0)
end

# strong scaling
for ncells in sort(unique(ncellss))
    @info "ncells=$ncells"
    frames1 = filter(row -> row["\$ncells_x"] == ncells, frames)
    frames1 = sort(frames1, :nodes)
    nodes = frames1.nodes
    cell_updates_per_second_per_node = frames1.cell_updates_per_second ./ frames1.nodes
    obj = scatterlines!(ax, nodes, cell_updates_per_second_per_node; linewidth=5, markersize=20)
    push!(legend, obj => "ncells=$(ncells)³ (strong scaling)")
    ylims!(; low=0)
end

Legend(fig[2, 1], [obj for (obj, leg) in legend], [leg for (obj, leg) in legend]; tellheight=true, tellwidth=false)

display(fig)

@info "Done."
