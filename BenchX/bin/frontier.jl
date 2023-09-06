using CSV
using CairoMakie
using DataFrames
using SixelTerm

@info "Evaluating Frontier benchmarks"

# Frontier is described on <https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html#>.
# Each node has 8 GPUs.
# Each GPU has a peak performance of 26.5 TFlop/s and a peak memory bandwidth of of 1.6 TByte/s.

# 38M cell updates/sec correspond to 700 kFlop and 42 kByte per cell update.

# frames = CSV.read("data/benchmark-z4c+asterx-sim-hip-benchmark-z4c+asterx-nlevels1-9ee1f994.csv", DataFrame)
# frames = CSV.read("data/benchmark-z4c+asterx-sim-hip-benchmark-z4c+asterx-nlevels4-0071c23b.csv", DataFrame)
# frames = CSV.read("data/benchmark-z4c+asterx-i0001-sim-hip-benchmark-z4c+asterx-nlevels1-660dbcca.csv", DataFrame)
# frames = CSV.read("data/benchmark-z4c+asterx-i0001-sim-hip-benchmark-z4c+asterx-nlevels4-431c3f43.csv", DataFrame)
frames = CSV.read("data/benchmark-z4c+asterx-i0002-sim-hip-benchmark-z4c+asterx-nlevels1-79ab1644.csv", DataFrame)

set_theme!(Theme(; fontsize=30))
fig = Figure(; resolution=(1536, 1152))

nlevels = 1
# TODO: Filter by `nlevels`

ax = Axis(
    fig[1, 1];
    title="Performance on OLCF Frontier (AsterX, Z4c, CarpetX, $nlevels levels AMR)",
    xlabel="# nodes",
    xscale=log10,
    xticks=2 .^ (0:8),
    ylabel="cell updates / sec / node",
)

ncells3s = frames.nlevels .* frames.ncells_x .* frames.ncells_y .* frames.ncells_z
nodess = frames.nodes

legend = []

# weak scaling
for ncells_per_node in [maximum(ncells3s ./ nodess)]
    @info "ncells_per_node=$ncells_per_node"
    frames1 = filter(row -> row.nlevels * row.ncells_x * row.ncells_y * row.ncells_z / row.nodes == ncells_per_node, frames)
    frames1 = sort(frames1, :nodes)
    nodes = frames1.nodes
    cell_updates_per_second_per_node = frames1.cell_updates_per_second ./ frames1.nodes
    obj = scatterlines!(ax, nodes, cell_updates_per_second_per_node; color=:red, linestyle=:dash, linewidth=5, markersize=20)
    push!(legend, obj => "# cells/node=$(round(Int, cbrt(ncells_per_node)))³ (weak scaling)")
    ylims!(; low=0)
end

# strong scaling
for ncells3 in sort(unique(ncells3s))
    @info "ncells=$(round(Int, cbrt(ncells3)))"
    frames1 = filter(row -> row.nlevels * row.ncells_x * row.ncells_y * row.ncells_z == ncells3, frames)
    frames1 = sort(frames1, :nodes)
    nodes = frames1.nodes
    cell_updates_per_second_per_node = frames1.cell_updates_per_second ./ frames1.nodes
    obj = scatterlines!(ax, nodes, cell_updates_per_second_per_node; linewidth=5, markersize=20)
    push!(legend, obj => "# cells=$nlevels⋅$(round(Int, cbrt(ncells3 / nlevels)))³ (strong scaling)")
    ylims!(; low=0)
end

Legend(fig[1, 2], [obj for (obj, leg) in legend], [leg for (obj, leg) in legend]; tellheight=true, tellwidth=false)

display(fig)

@info "Done."
