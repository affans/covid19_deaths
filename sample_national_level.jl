## Affan Shoukat, 2020
## this file goes through each state file `UST_{st}_posterior_y/z.csv` 
## `UST` files are created by Seyed from model output)
## and samples the empirical distribution to get the national level results. 
## the results are stored in files "national_y, national_z.csv" 

using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using DelimitedFiles
using TimerOutputs
using Random

const to = TimerOutput()

function read_df()
    #vector of states 
    validstates = ("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
    df_of_states = Array{DataFrame, 1}(undef, 50)
    for (i, vs) in enumerate(validstates)
        fn = "/data/actualdeaths_covid19/UST_$(vs)_posterior_z.csv"
        df_of_states[i] = CSV.File(fn, header=false) |> DataFrame
    end 
    return df_of_states
end

# fast preallocated
function _sample_for_national(mydat)
    cols = size(mydat[1])[2] ## get the number of columns for any matrix - this represents the time span of the analysis.
    full_df = zeros(Float64, cols, 12000)  # time x samples
    for simi = 1:12000        
        for st = 1:50 
            sl_df = mydat[st]
            for t = 1:cols
                sl_df_t = @view sl_df[:, t] 
                r_val = rand(sl_df_t)
                full_df[t, simi] += r_val
            end
        end   
    end
    return full_df
end

function process_national() 
    # make sure to change the filename "national_y, national_z" 
    # make sure to read the proper dataframes in the read_df function
    reset_timer!(to)
    dfs = @timeit to "read df" convert(Array{Matrix{Float64}, 1}, read_df())
    pdat = @timeit to "write file" _sample_for_national(dfs)
    writedlm("/data/actualdeaths_covid19/national_z.dat", pdat)
    return 1
end



