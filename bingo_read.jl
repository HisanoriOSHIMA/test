using LinearAlgebra
using Statistics
using Plots
using Pkg 
using StatsBase

function read_File(N, Nsample, free)

    dirname = "/Users/ooshimahisanori/Desktop/free_study/bingo/"
    fname = dirname*string(N)*"_sample"*string(Nsample)*"_"*free*".txt"

    data_out = zeros(Int64, Nsample)
    stats_info = zeros(Float64, 7)

    out = open(fname, "r")
        dt = readlines(out)
    close(out)

    data_out .= parse.(Int64, dt)

    stats_info[1] = N
    stats_info[2] = mean(data_out)
    stats_info[3] = std(data_out)
    stats_info[4] = mode(data_out)
    stats_info[5] = median(data_out)
    stats_info[6] = skewness(data_out)
    stats_info[7] = kurtosis(data_out)
    println("N = "*string(stats_info[1])*", mean = "*string(stats_info[2])*", std = "*string(stats_info[3])*", mode = "*string(stats_info[4])*", median = "*string(stats_info[5]))



    return data_out, stats_info
end



function hist(free, Nsample)
    N = 3
    results = read_File(N, Nsample, free)[1]
    histogram(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,26, step=1))

    N = 5
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 7
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 9
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 11
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 13
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 15
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 17
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 19
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 21
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 23
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 25
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 27
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 29
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 31
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 33
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 35
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))

    N = 37
    results = read_File(N, Nsample, free)[1]
    histogram!(results, ylabel = "Frequency", xlabel = "BINGO time", legend = false, bin=range(0,4000, step=1))



    Nvec = Vector(3:2:37)
    info = zeros(Float64, length(Nvec), 7)
    for iN in 1:length(Nvec)
        info[iN,:] .= read_File(Nvec[iN], Nsample, free)[2] 
    end

    Plots.plot!(info[:,4], 10^4*4.55 ./ info[:,3], linewidth = 3, marker = :star)

    dirname = "/Users/ooshimahisanori/Desktop/free_study/bingo/"
    #Plots.savefig(dirname*"BINGO"*"_sample"*string(Nsample)*"_"*free*".png")
end


function ratio(quant1, quant2, free, Nsample, Nvec)

    v = zeros(Float64, length(Nvec), 2)
    for iN in 1:length(Nvec)
        println(iN)
        results = read_File(Nvec[iN], Nsample, free)[2]
        if quant1 == "mean"
            q1 = results[2]
        elseif quant1 == "std"
            q1 = results[3]
        elseif quant1 == "var"
            q1 = results[3]^2
        elseif quant1 == "mode"
            q1 = results[4]
        elseif quant1 == "median"
            q1 = results[5]
        elseif quant1 == "skewness"
            q1 = results[6]
        elseif quant1 == "kurtosis"
            q1 = results[7]
        end
        if quant2 == "mean"
            q2 = results[2]
        elseif quant2 == "std"
            q2 = results[3]
        elseif quant2 == "var"
            q2 = results[3]^2
        elseif quant2 == "mode"
            q2 = results[4]
        elseif quant2 == "median"
            q2 = results[5]
        elseif quant2 == "skewness"
            q2 = results[6]
        elseif quant2 == "kurtosis"
            q2 = results[7]
        end

        v[iN, 1] = Nvec[iN]
        v[iN, 2] = (q1-q2)/results[3]
    end

    Plots.plot(v[:,1], v[:,2], linewidth=3, marker=:circ, legend=false, xlabel=raw"$N$", ylabel=quant1*"/"*quant2)

    Plots.plot!(x->1, linewidth=2)

    dirname = "/Users/ooshimahisanori/Desktop/free_study/bingo/"
    #Plots.savefig(dirname*"BINGO"*"_sample"*string(Nsample)*"_"*free*"_"*quant1*quant2*".png")
end


function quant(free, Nsample, quantity, Nvec)
    v = zeros(Float64, length(Nvec), 2)
    for iN in 1:length(Nvec)
        println(iN)
        results = read_File(Nvec[iN], Nsample, free)[2]

        if quantity == "mean"
            q = results[2]
        elseif quantity == "std"
            q = results[3]
        elseif quantity == "var"
            q = results[3].^2
        elseif quantity == "mode"
            q = results[4]
        elseif quantity == "median"
            q = results[5]
        elseif quantity == "skewness"
            q = results[6]
        elseif quantity == "kurtosis"
            q = results[7]
        end

        v[iN, 1] = Nvec[iN]
        v[iN, 2] = q
    end

    if quantity == "mean"
        ylb = quantity
    elseif quantity == "std"
        ylb = "standard deviation"
    elseif quantity == "var"
        ylb = "variance"
    elseif quantity == "mode"
        ylb = quantity
    elseif quantity == "median"
        ylb = quantity
    elseif quantity == "skewness"
        ylb = quantity
    elseif quantity == "kurtosis"
        ylb = quantity
    end


    Plots.plot(v[:,1], v[:,2], linewidth=3, marker=:circ, legend=false, xlabel=raw"$N$", ylabel=ylb)
end


free = "free"
       #"no_free"
Nsample = 10^5

hist(free, Nsample)
#ratio("mode", "mean", free, Nsample, Vector(3:2:37))
#quant(free, Nsample, "kurtosis", Vector(3:2:37))