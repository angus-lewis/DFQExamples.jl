using CSV, DataFrames

function rename_csv_col_Order_1_to_Unif(pth)
    for (root, dirs, files) in walkdir(pth)
        for file in files
            fname, extension = splitext(file)
            if ((extension==".csv")&&(!occursin("sim",fname))&&!occursin("sample",fname))
                fpath = joinpath(root,file)
                df = CSV.read(fpath,DataFrame)
                try 
                    rename!(df,Dict(:Order_1 => :Unif))
                    CSV.write(fpath,df)
                catch
                    @warn "file at "*fpath*" has col names "*join(names(df),", ")*". Did not overwrite"
                    @show names(df)
                end
            end
        end
    end
end

rename_csv_col_Order_1_to_Unif(".")