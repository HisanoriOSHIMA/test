using LinearAlgebra
using Statistics
using Plots
using Pkg 
using StatsBase



function create(N, j, step)
    column = zeros(Int64, N)
    for i in 1:N
        r = rand(step*(j-1)+1:step*j)
        column[i] = r
    end

    return column
end

function judge(card_column, M)
    n = 0
    for i in 1:M
        for k in i+1:M
            if card_column[i] == card_column[k]
                n = 1
            end
        end
    end

    return n
end


function bingo_card(N, step, free)
    card = zeros(Int64, N,N)
    for j in 1:N
        card_column = create(N, j, step)
        n = judge(card_column, N)
        while n == 1
            card_column = create(N, j, step)
            n = judge(card_column, N)
        end
        card[:,j] .= card_column
    end

    if free == "free"
        card[div(N+1,2), div(N+1,2)] = 0
    elseif free == "no_free"
    end

    return card
end



function suji(M)

    v = sample(1:M, M, replace=false)

    return v
end

 #j回目のナンバーコール
function bingo(card, v, j, M)
    num = v[j]

    for i in 1:M
        for j in 1:M
            if card[i,j] == num
                card[i,j] = 0
            end
        end
    end


    for i in 1:M #行
        if sum(card[i,:]) < 0.5
            card = zeros(Int64, M,M)
            break
        end
    end

    for i in 1:M #列
        if sum(card[:,i]) < 0.5
            card = zeros(Int64, M,M)
            break
        end
    end

    #斜め
    S1 = 0
    S2 = 0
    for i in 1:M
        S1 += card[i,i]
        S2 += card[i,M-i+1]
    end

    if S1 < 0.5 || S2 < 0.5
        card = zeros(Int64, M,M)
    else 
        return card
    end
    
    
end


function exe(N, Nsample, free)
    step = N * 3
    K = step * N

    results = zeros(Int64, Nsample)
    for l in 1:Nsample
        card = bingo_card(N, step, free)
        v = suji(K)
        for j in 1:K
            if sum(card[1,:]) > 0
                card = bingo(card, v, j, N)
            else
                #println("j = ", j)
                results[l] = j
                break
            end
        end
    end
    println("Average time = ", mean(results))

    return results
end




N = 9
Nsample = 10^7
free = "free"
       #"no_free"

results = exe(N, Nsample, free)

dirname = "/Users/ooshimahisanori/Desktop/free_study/bingo/"
#dirname = "/data2/"

fname = dirname*string(N)*"_sample"*string(Nsample)*"_"*free*".txt"

out = open(fname, "w")
Base.print_array(out, results)
close(out)



#=
results = exe(3, Nsample, free)
histogram(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(5, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(7, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(9, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(11, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(13, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(15, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(17, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
=#
#=
results = exe(19, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(21, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(23, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(25, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(27, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
results = exe(29, Nsample, free)
histogram!(results, ylabel = "Frequency", xlabel = "BINGO", legend = false)
=#



#Plots.savefig("/Users/ooshimahisanori/Desktop/free_study/bingo/BINGO"*string(N)*"_sample"*string(Nsample)*"_"*free*".png")
#Plots.savefig("/Users/ooshimahisanori/Desktop/free_study/bingo/BINGO"*"_sample"*string(Nsample)*"_"*free*".png")

            

