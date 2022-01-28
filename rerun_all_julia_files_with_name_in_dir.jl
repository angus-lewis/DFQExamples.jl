using CSV, DataFrames

function rerun_all(dir,filename_to_run)
    for (root, dirs, files) in walkdir(dir)
        for file in files
            if (file==filename_to_run)
                fpath = joinpath(root,file)
                try 
                    include(fpath)
                    println("Success")
                catch 
                    @warn "Could not rerun file at "*fpath
                end
            end
        end
    end
end

rerun_all(".","phd_ch_5_figs.jl")