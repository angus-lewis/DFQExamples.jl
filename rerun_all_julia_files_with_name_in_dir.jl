using CSV, DataFrames, Distributed

function rerun_all(dir,filename_to_run)
    # @sync begin
    for (root, dirs, files) in walkdir(dir)
        for file in files
            if (file==filename_to_run)
                fpath = joinpath(root,file)
                try 
                    #p = addprocs(1)
                    #@spawnat p[1] :(using CSV, DataFrames; 
                    include(fpath); println("Success")
                    #)
                catch
                    @warn "Could not rerun file at "*fpath
                end
            end
        end
    end
    # end
    rmprocs(procs()[2:end])
end

rerun_all(".","phd_ch_5_figs.jl")