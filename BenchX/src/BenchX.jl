module BenchX

using SHA
using YAML

################################################################################

export PhysicsParams
"""
    struct PhysicsParams

"Physics" parameters define a setup that is to be benchmarks. These
parameters completely describe, in detail, the complete setup that
make a simulation reproducible. For example, the numerical resolution
would be such a physics parameter, while the number of threads used to
run a simulation would not be.

See [`RunParams`](@ref).
"""
struct PhysicsParams
    name::AbstractString
    configuration::AbstractString
    config_id::AbstractString
    build_id::AbstractString
    parfile::AbstractString
    extra_physics_params::Dict{AbstractString,Union{AbstractString,Float64,Int}}
    function PhysicsParams(;
        name::AbstractString,
        configuration::AbstractString,
        config_id::Union{Nothing,AbstractString}=nothing,
        build_id::Union{Nothing,AbstractString}=nothing,
        parfile::AbstractString,
        extra_physics_params::Dict,
    )
        if config_id ≡ nothing
            cactusdir = find_cactusdir()
            config_id = chomp(read(joinpath(cactusdir, "configs", configuration, "CONFIG-ID"), String))
        end
        if build_id ≡ nothing
            cactusdir = find_cactusdir()
            build_id = chomp(read(joinpath(cactusdir, "configs", configuration, "BUILD-ID"), String))
        end
        return new(name, configuration, config_id, build_id, parfile, extra_physics_params)
    end
end
function Base.show(io::IO, physics::PhysicsParams)
    println(io, "PhysicsParams:")
    println(io, "    name: ", physics.name)
    println(io, "    configuration: ", physics.configuration)
    println(io, "    config_id: ", physics.config_id)
    println(io, "    build_id: ", physics.build_id)
    println(io, "    parfile: ", physics.parfile)
    println(io, "    extra_physics_params:")
    for (key, val) in sort!(collect(physics.extra_physics_params))
        println(io, "        ", key, ": ", val)
    end
    return nothing
end
function Base.convert(::Type{DICT}, physics::PhysicsParams) where {DICT<:AbstractDict}
    dict = DICT(
        :name => physics.name,
        :configuration => physics.configuration,
        :config_id => physics.config_id,
        :build_id => physics.build_id,
        :parfile => physics.parfile,
    )
    for (key, val) in physics.extra_physics_params
        dict[Symbol(key)] = val
    end
    return dict
end
function Base.hash(physics::PhysicsParams, h::UInt)
    return hash(
        physics.name,
        hash(
            physics.configuration,
            hash(physics.config_id, hash(physics.build_id, hash(physics.parfile, hash(physics.extra_physics_params, h)))),
        ),
    )
end

export RunParams
"""
    struct RunParams

"Run" parameters define how a physics setup is to be run. These
parameters completely describe, in detail, how a simulation is to be
executed on a particular machine. This includes e.g. the name of the
machine, the number of nodes, processes, threads, etc.

See also: [`PhysicsParams`](@ref).
"""
struct RunParams
    machine::AbstractString
    queue::AbstractString
    walltime_seconds::Float64
    nodes::Int
    cores_per_node::Int
    pus_per_core::Int
    processes::Int
    threads_per_process::Int
    smts_per_thread::Int
    extra_run_params::Dict{AbstractString,Union{AbstractString,Float64,Int}}
    function RunParams(;
        machine::AbstractString,
        queue::AbstractString,
        walltime_seconds::Union{AbstractFloat,Integer},
        nodes::Int,
        cores_per_node::Int,
        pus_per_core::Int,
        processes::Int,
        threads_per_process::Int,
        smts_per_thread::Int,
        extra_run_params::Dict,
    )
        return new(
            machine,
            queue,
            walltime_seconds,
            nodes,
            cores_per_node,
            pus_per_core,
            processes,
            threads_per_process,
            smts_per_thread,
            extra_run_params,
        )
    end
end
function Base.convert(::Type{DICT}, run::RunParams) where {DICT<:AbstractDict}
    dict = DICT(
        :machine => run.machine,
        :queue => run.queue,
        :walltime_seconds => run.walltime_seconds,
        :nodes => run.nodes,
        :cores_per_node => run.cores_per_node,
        :pus_per_core => run.pus_per_core,
        :processes => run.processes,
        :threads_per_process => run.threads_per_process,
        :smts_per_thread => run.smts_per_thread,
    )
    for (key, val) in run.extra_run_params
        dict[Symbol(key)] = val
    end
    return dict
end
function Base.show(io::IO, run::RunParams)
    println(io, "RunParams:")
    println(io, "    machine: ", run.machine)
    println(io, "    queue: ", run.queue)
    println(io, "    walltime_seconds: ", run.walltime_seconds)
    println(io, "    nodes: ", run.nodes)
    println(io, "    cores_per_node: ", run.cores_per_node)
    println(io, "    pus_per_core: ", run.pus_per_core)
    println(io, "    processes: ", run.processes)
    println(io, "    threads_per_process: ", run.threads_per_process)
    println(io, "    smts_per_thread: ", run.smts_per_thread)
    println(io, "    extra_run_params:")
    for (key, val) in sort!(collect(run.extra_run_params))
        println(io, "        ", key, ": ", val)
    end
    return nothing
end
function Base.hash(run::RunParams, h::UInt)
    return hash(
        run.machine,
        hash(
            run.queue,
            hash(
                run.walltime_seconds,
                hash(
                    run.nodes,
                    hash(
                        run.cores_per_node,
                        hash(
                            run.pus_per_core,
                            hash(
                                run.processes,
                                hash(run.threads_per_process, hash(run.smts_per_thread, hash(run.extra_run_params, h))),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
end

export Benchmark
"""
    struct Benchmark

A benchmark setup, consisting of [`PhysicsParams`](@ref) and
[`RunParams`](@ref).
"""
struct Benchmark
    physics::PhysicsParams
    run::RunParams
end
function Base.show(io::IO, benchmark::Benchmark)
    println(io, "Benchmark:")
    print(io, benchmark.physics)
    print(io, benchmark.run)
    return nothing
end
Base.hash(benchmark::Benchmark, h::UInt) = hash(benchmark.physics, hash(benchmark.run, h))

export Timing
"""
    struct Timing

The timing result of running a [`Benchmark`](@ref). This includes the
total run time, but also some other interesting quantities.
"""
struct Timing
    submitted::Bool
    success::Bool
    evolution_seconds::Float64
    evolution_compute_seconds::Float64
    evolution_output_seconds::Float64
    evolution_cell_updates::Float64
    evolution_iterations::Int
    cells::Float64
    cell_updates_per_second::Float64
    function Timing(;
        submitted::Bool,
        success::Bool,
        evolution_seconds::Real=NaN,
        evolution_compute_seconds::Real=NaN,
        evolution_output_seconds::Real=NaN,
        evolution_cell_updates::Real=NaN,
        evolution_iterations::Integer=0,
    )
        evolution_seconds = Float64(evolution_seconds)
        evolution_compute_seconds = Float64(evolution_compute_seconds)
        evolution_output_seconds = Float64(evolution_output_seconds)
        evolution_cell_updates = Float64(evolution_cell_updates)
        evolution_iterations = Int(evolution_iterations)
        if isnan(evolution_seconds) ||
            isnan(evolution_compute_seconds) ||
            isnan(evolution_output_seconds) ||
            isnan(evolution_cell_updates) ||
            evolution_iterations == 0
            @assert isnan(evolution_seconds) &&
                isnan(evolution_compute_seconds) &&
                isnan(evolution_output_seconds) &&
                isnan(evolution_cell_updates) &&
                evolution_iterations == 0
        end
        cells = evolution_cell_updates / evolution_iterations
        cell_updates_per_second = evolution_cell_updates / evolution_compute_seconds
        return new(
            submitted,
            success,
            evolution_seconds,
            evolution_compute_seconds,
            evolution_output_seconds,
            evolution_cell_updates,
            evolution_iterations,
            cells,
            cell_updates_per_second,
        )
    end
end

function Base.show(io::IO, timing::Timing)
    println(io, "Timing:")
    println(io, "    submitted: ", timing.submitted)
    println(io, "    success: ", timing.success)
    println(io, "    evolution_seconds: ", timing.evolution_seconds)
    println(io, "    evolution_compute_seconds: ", timing.evolution_compute_seconds)
    println(io, "    evolution_output_seconds: ", timing.evolution_output_seconds)
    println(io, "    evolution_cell_updates: ", timing.evolution_cell_updates)
    println(io, "    evolution_iterations: ", timing.evolution_iterations)
    println(io, "    number of cells (problem size): ", timing.cells)
    println(io, "    cell updates per second (strong performance): ", timing.cell_updates_per_second)
    return nothing
end
function Base.convert(::Type{DICT}, timing::Timing) where {DICT<:AbstractDict}
    return DICT(
        :submitted => timing.submitted,
        :success => timing.success,
        :evolution_seconds => timing.evolution_seconds,
        :evolution_compute_seconds => timing.evolution_compute_seconds,
        :evolution_output_seconds => timing.evolution_output_seconds,
        :evolution_cell_updates => timing.evolution_cell_updates,
        :evolution_iterations => timing.evolution_iterations,
        :cells => timing.cells,
        :cell_updates_per_second => timing.cell_updates_per_second,
    )
end

export BenchmarkResult
"""
    struct BenchmarkResult

A benchmark result, consisting of a [`Benchmark`](@ref) and a
[`Timing`](@ref).
"""
struct BenchmarkResult
    benchmark::Benchmark
    timing::Timing
end
function Base.show(io::IO, result::BenchmarkResult)
    println(io, "BenchmarkResult:")
    print(io, result.benchmark)
    print(io, result.timing)
    return nothing
end

################################################################################

export find_cactusdir
"""
    dir = find_cactusdir()
    dir::AbstractString

Find the main Cactus directory.
"""
function find_cactusdir()
    dir = pwd()
    while true
        isdir(joinpath(dir, "simfactory")) && return dir
        @assert dir ≠ "/"
        dir = splitdir(dir)[1]
    end
    @assert false
end

export make_physics_name
"""
    name = make_physics_name(physics::PhysicsParams)
    name::AbstractString

Create a unique name for the physics configuration. This name can be
used as a file name.
"""
function make_physics_name(physics::PhysicsParams; add_tag::Bool=true)
    parfile = physics.parfile
    parfile = replace(parfile, r"^.*/" => s"")
    parfile = replace(parfile, r"(\.par)?$" => s"")

    settings = [replace("$key$val", r"\$" => s"") for (key, val) in sort!(collect(physics.extra_physics_params); by=first)]
    settings::AbstractVector

    name = join(["benchmark", physics.name, physics.configuration, parfile, settings...], "-")

    # Limit length of simulation name
    maxlen = 100
    if length(name) > maxlen
        name = name[1:maxlen]
    end

    if add_tag
        tag = bytes2hex(sha256(string(physics)))[1:8]
        name = "$name-$tag"
    end

    return name
end

export make_simulation_name
"""
    name = make_simulation_name(benchmark::Benchmark)
    name::AbstractString

Create a unique name for the benchmark. This name can be used as
SimFactory simulation name.
"""
function make_simulation_name(benchmark::Benchmark; add_tag::Bool=true)
    physics_name = make_physics_name(benchmark.physics; add_tag=false)

    settings = [replace("$key$val", r"\$" => s"") for (key, val) in sort!(collect(benchmark.run.extra_run_params); by=first)]
    settings::AbstractVector

    name = join(
        [
            physics_name,
            benchmark.run.machine,
            benchmark.run.queue,
            "n$(benchmark.run.nodes)",
            "c$(benchmark.run.cores_per_node)",
            "v$(benchmark.run.pus_per_core)",
            "p$(benchmark.run.processes)",
            "t$(benchmark.run.threads_per_process)",
            "s$(benchmark.run.smts_per_thread)",
            settings...,
        ],
        "-",
    )

    # Limit length of simulation name
    maxlen = 200
    if length(name) > maxlen
        name = name[1:maxlen]
    end

    if add_tag
        tag = bytes2hex(sha256(string(benchmark)))[1:8]
        name = "$name-$tag"
    end

    return name
end

export RunStatus, rs_error, rs_missing, rs_queued, rs_running, rs_finished, rs_broken
"""
    Enum RunStatus:
        rs_error
        rs_missing
        rs_queued
        rs_running
        rs_finished
        rs_broken

The state of a simulation that is being run.
"""
@enum RunStatus rs_error rs_missing rs_queued rs_running rs_finished rs_broken
function Base.show(io::IO, status::RunStatus)
    status ≡ rs_error && return print(io, "error")
    status ≡ rs_missing && return print(io, "missing")
    status ≡ rs_queued && return print(io, "queued")
    status ≡ rs_running && return print(io, "running")
    status ≡ rs_finished && return print(io, "finished")
    status ≡ rs_broken && return print(io, "broken")
    @assert false
    return nothing
end
Base.show(io::IO, ::MIME"text/plain", status::RunStatus) = show(io, status)

export get_run_status
"""
    status = get_run_status(benchmark::Benchmark)
    status::RunStatus

Determine the status of a benchmark. This queries SimFactory, and it
take about a second to do so.

See [`RunStatus`](@ref).
"""
function get_run_status(benchmark::Benchmark)
    cactusdir = find_cactusdir()
    name = make_simulation_name(benchmark)

    @info "Querying simulation $name..."

    status = rs_error

    cmd = Cmd(
        Cmd([joinpath("simfactory", "bin", "sim"), "--machine=$(benchmark.run.machine)", "get-output-dir", name]); dir=cactusdir
    )
    output_dir = chomp(read(cmd, String))
    if !ispath(output_dir)
        status = rs_missing
    end

    if status ≡ rs_error
        for iter in 1:3
            iter ≠ 1 && sleep(1)
            output = read(
                Cmd(
                    Cmd([joinpath("simfactory", "bin", "sim"), "--machine=$(benchmark.run.machine)", "list-simulations", name]);
                    dir=cactusdir,
                ),
                String,
            )
            for line in split(output, "\n")
                if match(r"PRESUBMITTED|QUEUED", line) ≢ nothing
                    status = rs_queued
                    break
                elseif match(r"RUNNING", line) ≢ nothing
                    status = rs_running
                    break
                elseif match(r"FINISHED|INACTIVE", line) ≢ nothing
                    status = rs_finished
                    break
                end
            end
            status ≢ rs_error && break
        end
    end

    if status ≡ rs_error
        # The simulation exists, but Simfactory doesn't tell us anything about it. Assume the simulation is broken.
        status = rs_broken
    end

    @info "    $status"
    @assert status ≢ rs_error

    return status
end

"""
    struct NoCmd <: Base.AbstractCmd

Similar to a `Cmd`, but doesn't actually do anything.
"""
struct NoCmd <: Base.AbstractCmd end
Base.show(io::IO, ::NoCmd) = print(io, "NoCmd()")
Base.:(==)(::NoCmd, ::NoCmd) = true
Base.hash(::NoCmd, h::UInt) = hash(0xdf3a71e1, h)
Base.ignorestatus(::NoCmd) = NoCmd()
Base.wait(::NoCmd) = nothing

"""
    cmd = function submit_run_nowait(benchmark::Benchmark)
    cmd::AbstractCmd

Submit a benchmark. Do not wait for the submission to complete (nor
for the run to complete). Use `wait(cmd)` to wait for the submission
to complete.

Submitting a job takes about a second. This is the time to perform the
job submission itself, not the time to run the benchmark, which is
much longer. When multiple jobs need to be submitted at the same time,
then it is convenient to submit them asynchronously. This speeds up
job submission.

This is an internal function. Usually one would use
[`submit_run`](@ref) or [`submit_runs`](@ref) instead.
"""
function submit_run_nowait(benchmark::Benchmark)
    # Don't double-submit
    status = get_run_status(benchmark)

    if status ≡ rs_broken
        # Delete broken simulations manually
        @info "Manually deleting $name..."
        cactusdir = find_cactusdir()
        name = make_simulation_name(benchmark)
        cmd = Cmd(
            Cmd([joinpath("simfactory", "bin", "sim"), "--machine=$(benchmark.run.machine)", "get-output-dir", name]); dir=cactusdir
        )
        output_dir = chomp(read(cmd, String))
        rm(output_dir; force=true, recursive=true)
        status = rs_missing
    end

    status ≢ rs_missing && return NoCmd()

    cactusdir = find_cactusdir()
    name = make_simulation_name(benchmark)

    replacements = [
        ["--replace=$key=$val" for (key, val) in sort!(collect(benchmark.physics.extra_physics_params); by=first)]
        ["--replace=$key=$val" for (key, val) in sort!(collect(benchmark.run.extra_run_params); by=first)]
    ]

    walltime_seconds = round(Int, benchmark.run.walltime_seconds)
    hours = walltime_seconds ÷ 3600
    minutes = walltime_seconds % 3600 ÷ 60
    seconds = walltime_seconds % 60
    walltime = "$hours:$minutes:$seconds"

    nodes = benchmark.run.nodes
    cores_per_node = benchmark.run.cores_per_node
    pus_per_core = benchmark.run.pus_per_core
    processes = benchmark.run.processes
    threads_per_process = benchmark.run.threads_per_process
    smts_per_thread = benchmark.run.smts_per_thread

    threads = processes * threads_per_process * smts_per_thread
    @assert threads % nodes == 0
    threads_per_node = threads ÷ nodes
    cores = nodes * cores_per_node
    threads_per_core = cld(threads, cores)

    @info "Submitting $name..."
    # @info "    --configuration=$(benchmark.physics.configuration)"
    # @info "    --parfile=$(benchmark.physics.parfile)"
    # for replacement in replacements
    #     @info "    $replacement"
    # end
    # @info "    --queue=$(benchmark.run.queue)"
    # @info "    --walltime=$walltime"
    # @info "    --ppn=$cores_per_node"
    # @info "    --procs=$threads"
    # @info "    --ppn-used=$threads_per_node"
    # @info "    --num-threads=$threads_per_process"
    # @info "    --num-smt=$threads_per_core"

    cmd = Cmd([
        joinpath("simfactory", "bin", "sim"),
        "--machine=$(benchmark.run.machine)",
        "submit",
        name,
        "--configuration=$(benchmark.physics.configuration)",
        "--parfile=$(benchmark.physics.parfile)",
        replacements...,
        "--queue=$(benchmark.run.queue)",
        "--walltime=$walltime",
        "--ppn=$cores_per_node",
        "--procs=$threads",
        "--ppn-used=$threads_per_node",
        "--num-threads=$threads_per_process",
        "--num-smt=$threads_per_core",
    ])
    # @info cmd

    cmd = run(Cmd(cmd; dir=cactusdir); wait=false)

    return cmd
end

export submit_run
"""
    submit_run(benchmark::Benchmark)

Submit a benchmark. Use [`submit_runs`](@ref) to submit multiple
benchmarks simultaneously. This speeds up the submission process.

Use [`get_run_status`](@ref) to see the current state of a submitted
benchmark, and [`wait_for_run`](@ref) to wait for a benchmark to
complete.
"""
function submit_run(benchmark::Benchmark)
    cmd = submit_run_nowait(benchmark::Benchmark)
    @info "Waiting for submission to finish..."
    wait(cmd)
    @assert cmd.success()
    @info "Done."
    return nothing
end

export submit_runs
"""
    submit_runs(benchmarks::AbstractVector{Benchmark})

Submit several benchmarks at the same time. You can also use
[`submit_run`](@ref) to submit a single benchmark. Submitting multiple
benchmarks simultaneously speeds up the submission process.

Use [`get_run_status`](@ref) to see the current state of a submitted
benchmark, and [`wait_for_runs`](@ref) to wait for the benchmarks to
complete.
"""
function submit_runs(benchmarks::AbstractVector{Benchmark})
    cmds = submit_run_nowait.(benchmarks)
    @info "Waiting for all submissions to finish..."
    wait.(cmds)
    sleep(10)
    @info "Done."
    return nothing
end

export wait_for_run
"""
    status = wait_for_run(benchmark::Bnchmark)
    status::RunStatus

Wait for a benchmark to complete. This might wait for a long time
(minutes or hours) depending on the state of the queuing system.

You can use [`wait_for_runs`](@ref) to wait for multiple benchmarks.

See [`RunStatus`](@ref).
"""
function wait_for_run(benchmark::Benchmark)
    name = make_simulation_name(benchmark)
    @info "Waiting for $name..."
    delay_seconds = 1
    max_delay_seconds = 60
    while true
        status = get_run_status(benchmark)
        status ≡ rs_missing && return status
        status ≡ rs_finished && return status
        status ≡ rs_broken && return status
        @info "    Waiting $delay_seconds seconds..."
        sleep(delay_seconds)
        # Back off exponentially
        delay_seconds = min(max_delay_seconds, 2 * delay_seconds)
    end
end

export wait_for_runs
"""
    wait_for_runs(benchmarks::AbstractVector{Benchmark})

Wait for a benchmark to complete. This might wait for a long time
(minutes or hours) depending on the state of the queuing system.

You can use [`wait_for_run`](@ref) to wait for a single benchmarks.
"""
function wait_for_runs(benchmarks::AbstractVector{Benchmark})
    for benchmark in benchmarks
        status = wait_for_run(benchmark)
        @assert status ∈ (rs_missing, rs_finished, rs_broken)
    end
    return nothing
end

export read_run_timing
"""
    timing = read_run_timing(benchmark::Benchmark)
    timing::Timing

Read and interpret the timing result from a completed benchmark. The
benchmark must be in the `rs_finished` [`RunStatus`](@ref). See
[`get_run_status`](@ref), [`wait_for_run`](@ref).

See [`Timing`](@ref).
"""
function read_run_timing(benchmark::Benchmark)
    status = wait_for_run(benchmark)
    if status ∈ (rs_missing, rs_broken)
        # Run does not exist
        @info "   Run does not exist"
        return Timing(; submitted=false, success=false)
    end

    cactusdir = find_cactusdir()
    name = make_simulation_name(benchmark)

    cmd = Cmd(
        Cmd([joinpath("simfactory", "bin", "sim"), "--machine=$(benchmark.run.machine)", "get-output-dir", name]); dir=cactusdir
    )
    output_dir = chomp(read(cmd, String))

    parfile = benchmark.physics.parfile
    parfile = replace(parfile, r"^.*/" => s"")
    parfile = replace(parfile, r"(\.par)?$" => s"")

    timing = try
        @show filename = joinpath(output_dir, parfile, "performance.yaml")
        @show yaml = YAML.load_file(filename)
        @show performance = yaml["performance"]::Dict

        # Find last two iterations
        @show iters = sort!(collect(keys(performance)))
        @show iter1 = iters[end - 1]
        @show iter2 = iters[end - 0]
        @show perf1 = performance[iter1]
        @show perf2 = performance[iter2]
        perf(key) = perf2[key] - perf1[key]
        Timing(;
            submitted=true,
            success=true,
            evolution_seconds=perf("evolution-seconds"),
            evolution_compute_seconds=perf("evolution-compute-seconds"),
            evolution_output_seconds=perf("evolution-output-seconds"),
            evolution_cell_updates=perf("evolution-cell-updates"),
            evolution_iterations=perf("evolution-iterations"),
        )
    catch
        # Could not read output; run probably failed
        @info "   Could not read output; run probably failed"
        return Timing(; submitted=true, success=false)
    end

    return timing
end

end
