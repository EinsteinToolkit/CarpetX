using BenchX
using CSV
using DataFrames

# TODO:
# 1. create and store set of benchmarks, then run them independently, and analyze them independently
# 2. put benchmark description CSV into simulation directory
# 3. distinguish between "physics" and "discretization" in benchmarks; allow varying discretizations

physics = PhysicsParams(;
    configuration="sim-amd",
    name="carpetx",
    parfile="repos/CarpetX/BenchX/par/benchmark-carpetx.par",
    extra_physics_params=Dict("\$nlevels" => 1),
)

run = RunParams(;
    machine="symmetry-amd",
    queue="amdq",
    nodes=1,
    cores_per_node=64,
    pus_per_core=1,
    processes=1,
    threads_per_process=64,
    smts_per_thread=1,
    walltime_seconds=3600,
    extra_run_params=Dict(
        "\$ncells_x" => 128,
        "\$ncells_y" => 128,
        "\$ncells_z" => 128,
        "\$max_grid_size_x" => 32,
        "\$max_grid_size_y" => 32,
        "\$max_grid_size_z" => 32,
        "\$max_tile_size_x" => 1024^3,
        "\$max_tile_size_y" => 16,
        "\$max_tile_size_z" => 16,
    ),
)

# benchmark = Benchmark(physics, run)
# 
# status = get_run_status(benchmark)
# if status ≡ rs_unknown
#     submit_run(benchmark)
# end
# if status ≢ rs_finished
#     wait_for_run(benchmark)
# end
# 
# timing = read_run_timing(benchmark)

benchmarks = Benchmark[]
for nlevels in [1]
    for log_ncells in 21:24     # 21 => 128³
        log_ncells_z = log_ncells ÷ 3
        log_ncells_y = (log_ncells - log_ncells_z) ÷ 2
        log_ncells_x = log_ncells - log_ncells_z - log_ncells_y
        ncells_x = 2^log_ncells_x
        ncells_y = 2^log_ncells_y
        ncells_z = 2^log_ncells_z
        total_ncells = nlevels * ncells_x * ncells_y * ncells_z
        for nodes in [1]
            if true       # nodes * 2^23 ≤ total_ncells ≤ nodes * 2^26
                for max_grid_size_z in [32, 64, 128], max_grid_size_y in [32, 64, 128], max_grid_size_x in [32, 64, 128]
                    if max_grid_size_x ≥ max_grid_size_y ≥ max_grid_size_z
                        for max_tile_size_z in [2, 4, 8, 16, 32], max_tile_size_y in [2, 4, 8, 16, 32]
                            if max_tile_size_y ≥ max_tile_size_z
                                for threads_per_process in [1, 2, 4, 8, 16, 32, 64]
                                    @assert nodes % run.nodes == 0
                                    @assert run.threads_per_process % threads_per_process == 0
                                    processes =
                                        (nodes ÷ run.nodes) * (run.threads_per_process ÷ threads_per_process) * run.processes

                                    physics′ = PhysicsParams(;
                                        name=physics.name,
                                        configuration=physics.configuration,
                                        config_id=physics.config_id,
                                        build_id=physics.build_id,
                                        parfile=physics.parfile,
                                        extra_physics_params=copy(physics.extra_physics_params),
                                    )
                                    physics′.extra_physics_params["\$nlevels"] => nlevels

                                    run′ = RunParams(;
                                        machine=run.machine,
                                        queue=run.queue,
                                        walltime_seconds=run.walltime_seconds,
                                        # nodes=run.nodes,
                                        nodes=nodes,
                                        cores_per_node=run.cores_per_node,
                                        pus_per_core=run.pus_per_core,
                                        # processes=run.processes,
                                        # threads_per_process=run.threads_per_process,
                                        smts_per_thread=run.smts_per_thread,
                                        processes=processes,
                                        threads_per_process=threads_per_process,
                                        extra_run_params=copy(run.extra_run_params),
                                    )
                                    run′.extra_run_params["\$ncells_x"] = ncells_x
                                    run′.extra_run_params["\$ncells_y"] = ncells_y
                                    run′.extra_run_params["\$ncells_z"] = ncells_z
                                    run′.extra_run_params["\$max_grid_size_x"] = max_grid_size_x
                                    run′.extra_run_params["\$max_grid_size_y"] = max_grid_size_y
                                    run′.extra_run_params["\$max_grid_size_z"] = max_grid_size_z
                                    run′.extra_run_params["\$max_tile_size_y"] = max_tile_size_y
                                    run′.extra_run_params["\$max_tile_size_z"] = max_tile_size_z

                                    benchmark = Benchmark(physics′, run′)
                                    push!(benchmarks, benchmark)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

benchmarks′ = Benchmark[]
for benchmark in benchmarks
    status = get_run_status(benchmark)
    status ∈ (rs_missing, rs_broken) && push!(benchmarks′, benchmark)
end
submit_runs(benchmarks′)

wait_for_runs(benchmarks)

results = BenchmarkResult[]
for benchmark in benchmarks
    timing = read_run_timing(benchmark)
    result = BenchmarkResult(benchmark, timing)
    println(result)
    push!(results, result)
end

sort!(results; by=result -> -result.timing.cell_updates_per_second)

undollar(str::AbstractString) = replace(str, r"\$" => s"")
undollar(sym::Symbol) = Symbol(undollar(String(sym)))

dataframe = DataFrame(;
    # PhysicsParams
    name=AbstractString[],
    configuration=AbstractString[],
    config_id=AbstractString[],
    build_id=AbstractString[],
    parfile=AbstractString[],
    [Symbol(undollar(key)) => typeof(val)[] for (key, val) in sort!(collect(physics.extra_physics_params); by=first)]...,
    # RunParams
    machine=AbstractString[],
    queue=AbstractString[],
    walltime_seconds=Float64[],
    nodes=Int[],
    cores_per_node=Int[],
    pus_per_core=Int[],
    processes=Int[],
    threads_per_process=Int[],
    smts_per_thread=Int[],
    [Symbol(undollar(key)) => typeof(val)[] for (key, val) in sort!(collect(run.extra_run_params); by=first)]...,
    # Timing
    submitted=Bool[],
    success=Bool[],
    evolution_seconds=Float64[],
    evolution_compute_seconds=Float64[],
    evolution_output_seconds=Float64[],
    evolution_cell_updates=Float64[],
    evolution_iterations=Int[],
    cells=Float64[],
    cell_updates_per_second=Float64[],
)
for result in results
    dict = merge(
        convert(Dict{Symbol,Any}, result.benchmark.physics),
        convert(Dict{Symbol,Any}, result.benchmark.run),
        convert(Dict{Symbol,Any}, result.timing),
    )
    dict = Dict(undollar(k) => v for (k, v) in dict)
    push!(dataframe, dict)
end

physics_name = make_physics_name(physics)
CSV.write("$physics_name.csv", dataframe)
