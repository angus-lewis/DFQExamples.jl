if length(ARGS)∉(7,8)
    println("usage: (assuming you are running from DFQExamples.jl directory)")
    println("  julia -p n_procs empirical_cdfs.jl ARGS[1] ... ARGS[7] [ARGS[8]]")
    println("  where:")
    println("    n_procs: the number of workers to run (and hence number of output files as each worker writes its own data)")
    println("             n_procs+1 cpus will be used in total")
    println("    ARGS[1]: path to include model definition")
    println("    ARGS[2]: path to include simulation definitions (initial_condition_generator, stopping time)")
    println("    ARGS[3]: n_sims, number of sims to do")
    println("    ARGS[4]: n_boot, number of bootstrap resamples to do")
    println("    ARGS[5]: seed for rng")
    println("    ARGS[6]: filepath to write CDF to")
    println("    ARGS[7]: filepath to write bootstrap sample of CDFs to")
    println("    ARGS[8]: (opional) true/false, to specify whether to record time or position of fluid at the stopping time of the sim, default false")
    throw(ArgumentError("Wrong number of arguments"))
end

using Distributed
const n_procs = nprocs()-1
@show ps = procs()

@show pwd()
@everywhere begin 
    using Pkg; Pkg.activate(".")
    using DiscretisedFluidQueues
    using FluidFluidQueues
    using Random
    using StatsBase
    using JSON
    using TickTock
end 

@everywhere include($(ARGS[1])) # define model 
@everywhere include($(ARGS[2])) # define initial_condition_generator, stopping time, and points to evaluate the CDF, 
const n_sims = parse(Int, ARGS[3]) # define number of sims to do 
const n_boot = parse(Int, ARGS[4]) # number of bootstrap sample to do
const seed = parse(Int, ARGS[5]) # set rng seed
const rng = MersenneTwister(seed)
const record_time = length(ARGS)>7 ? parse(Bool, ARGS[8]) : false # specify whether to record time or position of sim

@everywhere begin
    function _simulate(model::FluidFluidQueue, tau, ic) 
        rest, ~ = simulate(model, tau, ic)
        X = rest.X
        φ = Int.(rest.φ)
        t = rest.t
        return (X=X, φ=φ, t=t)
    end
    _simulate(model, tau, ic) = simulate(model, tau, ic)
end
@everywhere function empirical_cdf(rng, ic_generator, model, n_sims, tau, x::AbstractRange, record_time=false)
    println("Simulating...")
    cdf = zeros(Int, length(x)+1, typeof(model)===BoundedFluidQueue ? 2 : 4)
    f, h = first(x), step(x)
    n_laps = 0
    total_laps = n_sims÷100_000_000
    mean_lap_time = 0.0
    tick()
    t0 = peektimer()
    for n in 1:n_sims
        (n%500_000_000 == 0) && begin 
            n_laps += 1
            t1 = peektimer()
            @info "msg @ time "*string(TickTock.now())
            println("  Simulating... progress on this worker: ", n/n_sims * 100, "%, elapsed: ", t1, " secs, last lap: ", t1-t0, " secs.")
            mean_lap_time = ((n_laps-1)*mean_lap_time + t1-t0)/(n_laps)
            println("  Simulating... predicted time remaining on this worker: ", (total_laps-n_laps)*mean_lap_time, " seconds.")
            println()
            t0 = t1
        end
        # simulate from model 
        ic = ic_generator(rng,1,model)
        sim_tau = _simulate(model, tau, ic)
        Xτ, φτ, tτ = only(sim_tau.X), only(sim_tau.φ), only(sim_tau.t)
        quantity_of_interest = record_time ? tτ : Xτ
        # add to cdf
        idx = ceil(Int, (quantity_of_interest - f) / h)+1
        # idx = Xτ<x[idx] ? idx-1 : idx
        cdf[idx, φτ] += 1
    end
    @info "msg @ time "*string(TickTock.now())
    tock()
    return cumsum(cdf,dims=1)./n_sims
end

@everywhere begin 
    @inline function batch_boot!(cdf,rng,pdf_ish,n)
        idxs = wsample(rng, 1:length(pdf_ish), pdf_ish[:], n)
        for i in idxs
            cdf[i] += 1
        end
        return nothing
    end

    function boot1(rng,pdf_ish,n)
        cdf_boot = zeros(Int, size(pdf_ish))
        # process in batches of 1_000_000 samples
        for _ in 1:(n÷1_000_000)
            batch_boot!(cdf_boot, rng, pdf_ish, 1_000_000)
        end
        batch_boot!(cdf_boot, rng, pdf_ish, n%1_000_000)
        return cumsum(cdf_boot, dims=1)./n
    end

    function boot(filepath, rng, cdf_paths, n_boot, ssize)
        jsons = [JSON.parsefile(pth) for pth in cdf_paths]
        cdf = zeros(length(jsons[1][1]),length(jsons[1]))
        for data in jsons 
            cdf .+= hcat(data...)
        end
        cdf ./= length(jsons)
        println("Bootstrapping...")
        boot_samples = Array{Array{Array{Float64,1},1},1}(undef,n_boot)
        pdf_ish = [cdf[1,:]'; diff(cdf,dims=1)]
        n_laps = 0
        total_laps = n_boot÷50
        mean_lap_time = 0.0
        tick()
        t0 = peektimer()
        for n in 1:n_boot
            (n%5 == 0) && begin 
                n_laps += 1
                t1 = peektimer()
                @info "msg @ time "*string(TickTock.now())
                println("  Bootstrapping... progress on this worker: ", n/n_boot * 100, "%, elapsed: ", t1, " secs, last lap: ", t1-t0, " secs.")
                mean_lap_time = ((n_laps-1)*mean_lap_time + t1-t0)/(n_laps)
                println("  Bootstrapping... predicted time remaining on this worker: ", (total_laps-n_laps)*mean_lap_time, " seconds.")
                println()
                t0 = t1
            end
            samp = boot1(rng, pdf_ish, ssize)
            boot_samples[n] = [samp[:,i] for i in 1:size(cdf,2)]
        end
        tock()
        @info "msg @ time "*string(TickTock.now())
        println("Writing bootstrap on this worker...")
        JSON.write(filepath,JSON.json(boot_samples))
        println("...write complete.")
        return nothing 
    end
end

@everywhere function sim_and_write(write_path, rng, initial_condition_generator, model, n_sims, tau, x, record_time)
    ecdf = empirical_cdf(rng, initial_condition_generator, model, n_sims, tau, x, record_time)
    JSON.write(write_path, JSON.json(ecdf))
    return nothing 
end

@time @sync begin 
    for n in n_procs+1:-1:3
        local_rng = MersenneTwister(rand(rng, UInt32)) # make sure rngs on each thread are distinct
        local_write_file = "proc_"*string(ps[n])*"_"*ARGS[6]
        @spawnat ps[n] sim_and_write(local_write_file, rng, initial_condition_generator, model, n_sims÷n_procs, tau, x, record_time)
    end
    local_rng = MersenneTwister(rand(rng, UInt32)) # make sure rngs on each thread are distinct
    local_write_file = "proc_"*string(ps[2])*"_"*ARGS[6]
    @spawnat ps[2] sim_and_write(local_write_file, local_rng, initial_condition_generator, model, n_sims-(n_procs-1)*(n_sims÷n_procs), tau, x, record_time)
end

local_cdf_paths = "proc_".*string.(ps[2:end]).*"_".*ARGS[6]
println("Bootstrapping on workers...")
@time @sync begin 
    for n in n_procs+1:-1:3
        local_rng = MersenneTwister(rand(rng, UInt32)) # make sure rngs on each thread are distinct
        local_write_file = "proc_"*string(ps[n])*"_"*ARGS[7]
        @spawnat ps[n] boot(local_write_file, local_rng, local_cdf_paths, n_boot÷n_procs, n_sims)
    end
    local_rng = MersenneTwister(rand(rng, UInt32)) # make sure rngs on each thread are distinct
    local_write_file = "proc_"*string(ps[2])*"_"*ARGS[7]
    @spawnat ps[2] boot(local_write_file, local_rng, local_cdf_paths, n_boot-(n_procs-1)*(n_boot÷n_procs), n_sims)
end
println("...workers complete.")
