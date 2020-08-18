using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using DelimitedFiles
using TimerOutputs
using Random


function read_df()
    #vector of states 
    validstates = ("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
    df_of_states = Array{DataFrame, 1}(undef, 50)
    for (i, vs) in enumerate(validstates)
        fn = "/data/Seyed/UST_$(vs)_posterior_z.csv"
        df_of_states[i] = CSV.File(fn, header=false) |> DataFrame
    end 
    return df_of_states
end

const to = TimerOutput()

## very slow when used with map/pmap, requires df to be global
function main()
    national_vec = zeros(Float64, 141)        
    for st = 1:50 
        sl_df = dfs[st]
        for t = 1:141
            sl_df_t = @view sl_df[:, t] 
            r_val = rand(sl_df_t)
            national_vec[t] += r_val
        end
    end   
    return national_vec
end

# same as above, very slow, but argument is passed 
@everywhere function main_with_pass(mydat)
    national_vec = zeros(Float64, 141)        
    for st = 1:50 
        sl_df = mydat[st]
        for t = 1:141
            sl_df_t = @view sl_df[:, t] 
            r_val = rand(sl_df_t)
            national_vec[t] += r_val
        end
    end   
    return national_vec
end

# fast preallocated
function main_two(mydat)
    cols = size(mydat[1])[2] ## get the number of columns for any matrix
    full_df = zeros(Float64, cols, 20000)
    for simi = 1:20000        
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

function a() 
    #pdat = @timeit to "no arg" map(x -> main(), 1:100)
    #pdat = @timeit to "with arg" map(x -> main(), 1:100)
    reset_timer!(to)
    dfs = @timeit to "read df" convert(Array{Matrix{Float64}, 1}, read_df())
    pdat = @timeit to "write file" main_two(dfs)
    writedlm("national_z.dat", pdat)
    #pdat = pmap(x -> main(), 1:100)
    #pdat_mat = convert(Matrix, hcat(pdat...)')
    #return pdat_mat
    
    return 1
end

