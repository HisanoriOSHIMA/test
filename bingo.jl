using LinearAlgebra
using Statistics
using Plots
using StatsBase


#Assume N×N size bingo-card.
#Row j of the card contains numbers from step*(j-1) to step*j.
#The numbers are chosen from a uniform distribution so that there are no duplicates.
#This function creates row j of the card.
function create(N, j, step)
    column = zeros(Int64, N)
    for i in 1:N
        r = sample(step*(j-1):step*j, N, replace=false)
        column[i] = r[i]
    end

    return column
end

#This function creates N×N size bingo-card.
#free : decides whether to make a hole in the middle of the matrix (free space).
function bingo_card(N, step, free)
    card = zeros(Int64, N,N)
    for j in 1:N
        card_column = create(N, j, step)
        card[:,j] .= card_column
    end

    if free == "free"
        card[div(N+1,2), div(N+1,2)] = 0
    elseif free == "no_free"
    end

    return card
end

#N numbers without duplicates are drawn from a uniform distribution.
function call_numbers(N)

    v = sample(1:N, N, replace=false)

    return v
end

#j-th number call
function bingo(card, v, j, N)
    num = v[j]

    for i in 1:N
        for j in 1:N
            if card[i,j] == num
                card[i,j] = 0
            end
        end
    end


    for i in 1:N #Bingo in a column
        if sum(card[i,:]) < 0.5
            card = zeros(Int64, N,N)
            break
        end
    end

    for i in 1:N #Bingo in a row
        if sum(card[:,i]) < 0.5
            card = zeros(Int64, N,N)
            break
        end
    end

    #Bingo in a diagonal row
    S1 = 0
    S2 = 0
    for i in 1:N
        S1 += card[i,i]
        S2 += card[i,N-i+1]
    end

    if S1 < 0.5 || S2 < 0.5
        card = zeros(Int64, N,N)
    else 
        return card
    end
    
end

#Nsample : Number of game attempts
#Record the number of calls until one gets a bingo.
function exe(N, Nsample, free)
    step = N * 3
    K = step * N

    results = zeros(Int64, Nsample)
    for l in 1:Nsample
        card = bingo_card(N, step, free)
        v = call_numbers(K)
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




N = 3
Nsample = 10^5
free = "free"
       #"no_free"

results = exe(N, Nsample, free)

dirname = "/Users/ooshimahisanori/Desktop/free_study/bingo/"
fname = dirname*string(N)*"_sample"*string(Nsample)*"_"*free*".txt"

out = open(fname, "w")
Base.print_array(out, results)
close(out)