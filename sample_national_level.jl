## Affan Shoukat, 2020
## processes the individual state posteriors at the national level 
## the order of functions are basically how they should be executed
## 1) run create_all_USTs() to generate state-specific UST files. 
##    - these are simply the posterior files, but padded to have uniform time 
## 2) run sample_national_level() (depends on read_UST(), and _sample_for_national())
##    - this function then uses the UST files to generate national_z and national_y files 
## 3) run national_statistics() to get means/cred intervals of these national_z / national_y files
##    - creates a single file called "antional_processed_zy.csv"
##    - and samples the empirical distribution to get the national level results. 
##    - the results are stored in files "national_y, national_z.csv" 
## to do: create function: pipeline() that runs the three functions at the same time

using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using DelimitedFiles
using TimerOutputs
using Random
using RCall
using Distributions
using Gnuplot

@rlibrary bayestestR

const to = TimerOutput()
const list_of_states = ["AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY"]

function create_all_USTs() 
    timepoints = 304  # starting at january 22 to november 20

    map(list_of_states) do vs
        println("processing $vs")
        fn_y = "/data/actualdeaths_covid19/st_$(vs)_00_posterior_y.dat"
        fn_z = "/data/actualdeaths_covid19/st_$(vs)_00_posterior_z.dat"
        df_y = CSV.File(fn_y, header=true) |> DataFrame!
        df_z = CSV.File(fn_z, header=true) |> DataFrame!
        n_zero_cols = timepoints - size(df_y)[2]
        new_df = DataFrame([Int64 for i = 1:n_zero_cols], [Symbol("c$i") for i = 1:n_zero_cols], 12000)
        new_df .= 0    
        newdf_y = hcat(new_df, df_y)
        newdf_z = hcat(new_df, df_z)
        # In Seyeds file the headers are off
        # but this script sets the headers by default.         
        CSV.write("/data/actualdeaths_covid19/UST_af_$(vs)_posterior_y.csv", newdf_y)
        CSV.write("/data/actualdeaths_covid19/UST_af_$(vs)_posterior_z.csv", newdf_z)        
    end
end

function read_UST(y_or_z) 
    df_of_states = Array{DataFrame, 1}(undef, length(list_of_states))
    for (i, vs) in enumerate(list_of_states)
        fn = "/data/actualdeaths_covid19/UST_af_$(vs)_posterior_$(y_or_z).csv"        
        df_of_states[i] = CSV.File(fn, header=true) |> DataFrame
    end 
    return df_of_states
end

# fast preallocated
# function _sample_for_national(mydat)
#     cols = size(mydat[1])[2] ## get the number of columns for any matrix - this represents the time span of the analysis.
#     full_df = zeros(Float64, cols, 12000)  # time x samples
#     for simi = 1:12000        
#         for st = 1:50 
#             sl_df = mydat[st]
#             for t = 1:cols
#                 sl_df_t = @view sl_df[:, t] 
#                 r_val = rand(sl_df_t)
#                 full_df[t, simi] += r_val
#             end
#         end   
#     end
#     return full_df
# end

function _sample_for_national(mydat)
    cols = size(mydat[1])[2] ## get the number of columns for any matrix - this represents the time span of the analysis.
    full_df = zeros(Float64, cols, 12000)  # time x samples

    for st = 1:length(mydat)
        sl_df = mydat[st]
        for t = 1:cols             
            reps = @view sl_df[:, t]
            full_df[t, :] += reps
        end
    end
    return full_df
end

function sample_national_level() 
    # make sure to change the filename "national_y, national_z" 
    # make sure to read the proper dataframes in the read_df function
    reset_timer!(to)
    df_y = @timeit to "read df z" convert(Array{Matrix{Float64}, 1}, read_UST("y"))
    df_z = @timeit to "read df y" convert(Array{Matrix{Float64}, 1}, read_UST("z"))
    pdat_y = @timeit to "write file y" _sample_for_national(df_y)
    pdat_z = @timeit to "write file z" _sample_for_national(df_z)
    writedlm("/data/actualdeaths_covid19/national_y.dat", pdat_y)    
    writedlm("/data/actualdeaths_covid19/national_z.dat", pdat_z)
end

function r_hdi(vdat)
    data = rand(Normal(0, 1), 1000)
    hdi = R"""
        library(bayestestR)
        time_ci = ci($vdat, method = "HDI", ci = 0.95)
        list(lo = time_ci[[2]], hi = time_ci[[3]])
    """
    hd_df = rcopy(hdi)
end

function get_national_statistics(df) 
    #df = CSV.File("/data/actualdeaths_covid19/national_y.dat", delim='\t', header=false) |> DataFrame!
    time_means = mean.(eachrow(df))
    lo_means = Float64[]
    hi_means = Float64[]
    map(eachrow(df)) do row 
        hdi = r_hdi(convert(Array, row))
        push!(lo_means, hdi[:lo])
        push!(hi_means, hdi[:hi])
    end
    DataFrame(time = 1:size(df)[1], means = time_means, lows = lo_means, highs = hi_means)
end 

function national_statistics() 
    #all_death_data = CSV.File("Death.csv", header=false) |> DataFrame! 
    #rename!(all_death_data, Symbol.(list_of_states))
    # we can just use our file now instead of Seyeds
    # but now we have total, confirmed, and probable
    all_death_data =  CSV.File("/data/actualdeaths_covid19/downloaded_data/incidence_deaths_lancetid.csv") |> DataFrame! 
    _deathcols =  ["death_$i" for i in list_of_states]
    all_death_data = select(all_death_data, _deathcols)

    df_y = CSV.File("/data/actualdeaths_covid19/national_y.dat", delim='\t', header=false) |> DataFrame!
    df_z = CSV.File("/data/actualdeaths_covid19/national_z.dat", delim='\t', header=false) |> DataFrame!
    proc_y = get_national_statistics(df_y) 
    proc_z = get_national_statistics(df_z)
    proc = innerjoin(proc_y, proc_z, on = :time, makeunique=true)
    proc[!, :raw] = sum.(eachrow(all_death_data))
    rename!(proc, [:time, :means_y, :lows_y, :highs_y, :means_z, :lows_z, :highs_z, :raw])
    
    # this code was used to check whether the julia code produces same as R code
    # checkdf = CSV.File("/data/actualdeaths_covid19/national_processed_zy.csv") |> DataFrame!
    # sum.(eachcol(checkdf)) .≈ sum.(eachcol(proc))
    # save the processed results in home directory
    CSV.write("/data/actualdeaths_covid19/national_processed_zy.csv", proc)
end

function create_national_plot()     
    est = CSV.File("/data/actualdeaths_covid19/national_processed_zy.csv", header=true) |> DataFrame!    
    est = est[1:end, :] 
    xvals = est.time    
    @gp "set grid" "set key left"
    #@gp "set term svg enhanced standalone mouse size 600,400"    
    #@gp :- "set size 1000, 1000"  ## dosn't work in the terminal, but maybe for saving?  
    @gp :- "unset colorbox"
    #@gp :- "set xtics 1" "set xlabel 'Week'" "set ylabel 'Difference"
    @gp :- "set xlabel 'Time (2020)'" "set ylabel 'Death Count'"
    @gp :- "set key title 'Legend'"
    @gp :- "set key top left Left reverse samplen 1" ## ???
    @gp :- "set label 'reported: $(round(sum(est.raw)))' at  graph 0.7, 0.95"
    @gp :- "set label 'cumulative mean sum: $(round(sum(est.means_y)))' at  graph 0.7, 0.9"
    @gp :- "set label 'cumulative lo sum: $(round(sum(est.lows_y)))' at  graph 0.7, 0.85"
    @gp :- "set label 'cumulative hi sum: $(round(sum(est.highs_y)))' at  graph 0.7, 0.80"
    #@gp :- "set xtics 1"
    @gp :- raw"""set xtics ('Mar 1' 39, 'Apr 1' 70, 'May 1' 100, 'Jun 1' 131, 'Jul 1' 161, 'Aug 1' 192, 'Sep 1' 223)"""
    @gp :- xvals est.lows_y est.highs_y "with filledcu notitle lw 2 lc rgb '#a1d99b' "
    @gp :- xvals est.means_y "with lines title 'Estimated mean' lw 2 lc rgb '#8031a354' "
    
    # add the z curve on the plot 
    @gp :- xvals est.means_z "with lines title 'Fitted curve to data' lc 'black'"
    @gp :- xvals est.raw "with linespoints title 'Reported Deaths' lc 'black' pt 8"
    @gp :- xvals est.raw "with lines notitle lc 'black'"
    save(term="svg enhanced standalone mouse size 1600,400", output="Figure_national.svg")
    #save(term="pngcairo", output="ex1.png")
    #save(term="pngcairo size 550,350 fontscale 0.8", output="assets/output.png")        
end

function create_trace_plot(st)
    ## creates a trace plot for each state. 
    #println("am i not revising")
    est = CSV.File("/data/actualdeaths_covid19/st_$(st)_00_traceplots.dat", header=true) |> DataFrame!    
    xvals = 1:nrow(est)
    @gp "reset session"
    @gp :- "set grid" "set key left"    
    @gp :- "unset colorbox"
    @gp :- "set xlabel 'iteration'" "set ylabel 'parameter value'"
    @gp :- "set title '$st'"
    #@gp :- "set size 1600, 1480"
    @gp :- "set multiplot layout 5,1 rowsfirst margins 0.10,0.95,0.05,0.9 spacing 0.05,0.0" #l r b t

    @gp :- "unset xlabel" 
    @gp :- "unset ylabel"
    @gp :- "unset xtics" 
    @gp :- 1 xvals est[!, :a0] "with lines title '$(string(:a0))' lw 0.1 lc 'black'" :-         
    @gp :- 2 xvals est[!, :a1] "with lines title '$(string(:a1))' lw 0.1 lc 'black'" :-    
    @gp :- 3 xvals est[!, :b0] "with lines title '$(string(:b0))' lw 0.1 lc 'black'" :-
    @gp :- "set ylabel 'parameter value'"
    @gp :- 4 xvals est[!, :b1] "with lines title '$(string(:b1))' lw 0.1 lc 'black'" :-    
    @gp :- "unset ylabel"
    @gp :- 5 xvals est[!, :b2] "with lines title '$(string(:b2))' lw 0.1 lc 'black'"     
    @gp :- "set xtics"    
    #save(term="svg enhanced standalone mouse size 1600,800", output="trace_plot.svg")    
    fname = "/data/actualdeaths_covid19/trace_plots_$st.ps"
    save(term="postscript eps enhanced color size 6,6", output=fname) # size in inches for postscript
    println("updated size 3")
    #return est
    return 1
end

function get_all_traceplots() 
    map(list_of_states) do x 
        println("working $x")
        create_trace_plot(x)
    end
end
