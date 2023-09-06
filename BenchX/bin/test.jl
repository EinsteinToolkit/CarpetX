using BenchX
using CSV
using DataFrames

# TODO:
# 1. create and store set of benchmarks, then run them independently, and analyze them independently
# 2. put benchmark description CSV into simulation directory
# 3. distinguish between "physics" and "discretization" in benchmarks; allow varying discretizations

physics = PhysicsParams(;
    # configuration="sim",
    # configuration="sim-llvm",
    # configuration="sim-llvm-32",
    configuration="sim-hip",

    # name="z4c",
    # parfile="repos/cactusamrex/BenchX/par/benchmark-z4c.par",
    # name="z4c+grhydrox-i0002",
    # parfile="repos/cactusamrex/BenchX/par/benchmark-z4c+grhydrox.par",
    name="z4c+asterx-i0002",
    parfile="repos/cactusamrex/BenchX/par/benchmark-z4c+asterx.par",
    extra_physics_params=Dict(
        "\$nlevels" => 1,
        # "\$nlevels" => 4,
        # "\$nlevels" => 8,
    ),
)

run = RunParams(;
    # machine="symmetry",
    # machine="symmetry-llvm",
    # machine="symmetry-llvm-32",
    # queue="debugq",
    # nodes=1,
    # cores_per_node=40,
    # pus_per_core=1,
    # processes=2,
    # threads_per_process=20,
    # smts_per_thread=1,

    # machine="crusher-hip",
    machine="frontier-hip",
    queue="batch",
    nodes=1,
    cores_per_node=64,
    pus_per_core=1,
    processes=8,
    threads_per_process=1,
    smts_per_thread=1,
    walltime_seconds=3600,
    extra_run_params=Dict(
        "\$ncells_x" => 256,
        "\$ncells_y" => 256,
        "\$ncells_z" => 256,
        # "\$ncells_x" => 512,
        # "\$ncells_y" => 512,
        # "\$ncells_z" => 512,
        # "\$ncells_x" => 1024,
        # "\$ncells_y" => 1024,
        # "\$ncells_z" => 1024,
        # "\$ncells_x" => 2048,
        # "\$ncells_y" => 2048,
        # "\$ncells_z" => 2048,
        # "\$max_grid_size" => 32,
        "\$max_grid_size" => 128,
        # "\$max_grid_size" => 256,
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
for nlevels in [1] # [1, 4]    
    for log_ncells in 24:33     # 24 => 256³, 33 => 2048³
        log_ncells_z = log_ncells ÷ 3
        log_ncells_y = (log_ncells - log_ncells_z) ÷ 2
        log_ncells_x = log_ncells - log_ncells_z - log_ncells_y
        ncells_x = 2^log_ncells_x
        ncells_y = 2^log_ncells_y
        ncells_z = 2^log_ncells_z
        total_ncells = nlevels * ncells_x * ncells_y * ncells_z
        for nodes in [1, 2, 4, 8, 16, 32, 64, 128, 256] # [1, 2, 4, 8, 16, 32, 64, 128, 256]
            if nodes * 2^23 ≤ total_ncells ≤ nodes * 2^26
                for max_tile_size_y in [128] # [2, 4, 8, 16, 32]
                    for max_tile_size_z in [128] # [2, 4, 8, 16, 32]
                        if max_tile_size_y ≥ max_tile_size_z
                            for threads_per_process in [1] # [1, 2, 4, 5, 8, 10, 20, 40]
                                @assert nodes % run.nodes == 0
                                @assert run.threads_per_process % threads_per_process == 0
                                processes = (nodes ÷ run.nodes) * (run.threads_per_process ÷ threads_per_process) * run.processes

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
