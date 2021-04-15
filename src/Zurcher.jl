
module Zurcher

using LinearAlgebra  # for norm()
using Plots
using DataFrames
using CSV
using StatsPlots
using StatsBase  # counts
using Optim
using NLopt
using PrettyTables

using JuMP 
using Ipopt

"""
    make_trans(θ, n)

Build a transition matrix of mileage states for `n` mileage bins and a vector θ of probabilities 
    to move to each possible bin (or to stay in the same bin).

"""
function make_trans(θ, n)
    transition = zeros(n, n);
    p = [θ ; 1 - sum(θ)]
    if any(p .< 0)
        @warn "negative probability in make_trans" maxlog=1
    end
    np = length(p)

    # standard block
    for i = 0:n - np
        transition[i + 1,(i + 1):(i + np)] = p
    end

    for k in 1:(np-1)
        transition[n - k,(n - k):n] = [p[1:k]...; 1 - sum(p[1:k])]
    end
    transition[n,n] = 1.0
    return transition
end

"""
# Harold Zurcher DataType

* This type represents an instance of Harold Zurcher, the decision maker in Rust (1987)
* The default constructor builds `Harold` for the values in table X of the paper.

"""
mutable struct Harold
    # parameters
    n::Int
    maxmiles::Int
    RC::Float64 
    c::Float64 
    θ::Vector{Float64}
    β::Float64 

    # numerical settings
    tol::Float64

    # state space
    mileage::Vector{Float64}
    transition::Matrix{Float64}

    function Harold(;n = 175,maxmiles = 450, RC = 11.7257,c = 2.45569,
                      θ = [0.107, 0.5152 ,0.3622, 0.0143,0.0009],
                      β = 0.999, tol =  1e-12)

        this = new()   # create new empty Harold instance
        this.n = n
        this.maxmiles = maxmiles
        this.RC = RC
        this.c = c
        this.θ = θ
        this.β = β
        this.tol = tol

        # build state space
        this.mileage = collect(0.0:n - 1)
        this.transition = make_trans(θ, n)
        return this
    end
end


"""
    bellman(h::Harold, ev0::Vector)

Bellman Operator `T(EV)`. Notice the centering of the log sum function by `M`.
"""
function bellman(h::Harold, ev0::Vector)
    maintainance = -0.001 .* h.mileage .* h.c  # StepRange of length n
    v0 = maintainance .+ h.β .* ev0  # vector of length n
    v1 = -h.RC + maintainance[1] + h.β * ev0[1]  #  a scalar. if you replace, you start at zero miles
    M  = maximum(vcat(v1, v0))  #  largest value in both vs
    logsum = M .+ log.(exp.(v0 .- M) .+ exp(v1 - M))
    ev1 = h.transition * logsum #  matrix multiplication
    ev1
end

"""
    ccp(h::Harold, ev0::Vector)

Conditional Choice Probability function returning the proability of replacing the engine at a given mileage state `x`.
"""
function ccp(h::Harold, ev::Vector)
    maintainance = -0.001 .* h.mileage .* h.c  # StepRange of length n
    v0 = maintainance .+ h.β .* ev  # vector of length n
    v1 = -h.RC + maintainance[1] + h.β * ev[1] 
    1.0 ./ (exp.(v0 .- v1) .+ 1)
end

"""
    vfi(h::Harold)

Value Function Iterator. Will iterate on [`bellman`](@ref) until the absolute norm is less than value `h.tol`.
"""
function vfi(h::Harold)
    if h.β >= 1.0
        throw(ArgumentError("value function iteration will not converge with β >= $(h.β)"))
    end
    ev0 = zeros(h.n)  # starting value
    err = 100.0
    iters = 0

    while err > h.tol
        iters += 1
        ev1 = bellman(h, ev0)
        err = norm(abs.(ev0 .- ev1))
        ev0[:] .= ev1   #  [:] do not reallocate a new object
    end
    return (ev0, iters)
end

"""
    runit(;n=90, β=0.9999,c=2.4, θ=[0.3, 0.68, 0.02])

Test run.
"""
function runit(;n=90, β=0.9999,c=2.4, θ=[0.3, 0.68, 0.02])
    z = Harold(n=n, β=β, c=c, θ=θ)
    sol, iters = vfi(z)
    pr  = ccp(z, sol)
    sol, pr, iters, z
end

"""
    plotit(;n=90, β=0.9)

Test plot.
"""
function plotit(;n=90, β=0.9)
    sol, pr, iters, z = runit(n=n, β=β)
    plot(plot(z.mileage, sol, title="Value Function"), 
            plot(z.mileage, pr,  title="Probability of Replacement"),
            xlab="Miles", leg=false)
end

"""
    simit(; T=500,n=500, θ=[0.3,0.6,0.1])

Simulate one bus.
"""
function simit(; T=500,n=500, θ=[0.3,0.6,0.1])
    z = Harold(n=n, θ=θ)  # need go higher with miles
    sol, iters = vfi(z)
    pr  = ccp(z, sol)

    P = cumsum(z.θ ./ sum(z.θ))

    miles = Int[]
    push!(miles, 1)  # start on first grid point of mileage
    a = Int[]  # 0/1 keep or replace
    push!(a, 0)  #  keep in first period

    for it in 2:T
        action = rand() < pr[miles[end]] ? 1 : 0
        push!(a, action)
            # update miles
        if action == 1
            push!(miles, 1)  # go back to first state
        else
            next_miles = findfirst(rand() .< P)  # index of first `true`
            push!(miles, miles[end] + next_miles)
        end
    end

    plot(1:T, miles, xlab="Period", ylab="miles",
                title="Simulating Harold Zurcher's Decisions",
                leg=false, lw=2)
end

"""
    dataplots()

return a dict of plots and tables from the data.
"""
function dataplots()
    d = CSV.read(joinpath(@__DIR__,"buses.csv"), DataFrame, header = false)
    select!(d, 1,2,5,6,7)
    rename!(d, [:busid, :bustype, :replaced, :miles_lag, :miles])

    # group by bus group (variable `:groupid`)
    gd = groupby(d, :bustype)

    # on the gd compute num of obs and compute mean of `:replace`
    cg = combine(gd, nrow => :nobs, 
                :busid => (x -> length(unique(x))) => :nbuses, 
                :replaced => sum => :n_replaced,
                :replaced => mean => :mean_replaced,
                :miles => mean => :mean_miles,
                :miles => maximum => :max_miles,
                )

    # output dict 
    di = Dict()

    di[:tab1] = pretty_table(String,cg, formatters = ft_printf(["%1.4f","%6.1f"], 5:6))

        # total buses?
    sum(cg.nbuses)         

        # subset to all but groups 1 and 2
    d2 = d[in.(d.bustype, Ref([3,4])), :]

    # group g2 and compute summary stats
    g2 = groupby(d2, :bustype)
    cg2 = combine(g2, nrow => :nobs, 
    :busid => (x -> length(unique(x))) => :nbuses, 
    :replaced => sum => :n_replaced,
    :replaced => mean => :mean_replaced,
    :miles => mean => :mean_miles,
    :miles => maximum => :max_miles,
    )

    stats_replacement = combine(filter(r -> r.replaced .== 1,d2), nrow => :nobs, 
        :miles_lag => mean => :mean_miles,
        :miles_lag => maximum => :max_miles,
    )
    di[:tab2] = pretty_table(String,stats_replacement, formatters = ft_printf(["%6.1f","%6.1f"], 2:3))

    di[:miles] = @df d2 scatter(:miles_lag, :replaced, leg = false)

    # discretized data
    dd = busdata(Harold())

    # plot average replacement indicator by mileage state
    b = groupby(dd, :x)
    bh = combine(b, :d => mean => :mean_replaced)
    di[:disc_miles] = @df bh bar(:x, :mean_replaced, xlab = "mileage state", leg = false)

    # plot replacement probs by group and mileage state
    b2 = groupby(dd, [:x, :bustype])
    bh = combine(b2, :d => mean => :mean_replaced, nrow)
    sort!(bh, [:x, :bustype])
    di[:disc_miles_repl] = @df bh bar(:x, :mean_replaced, group=:bustype, title="avg replacement indicator", bar_position=:dodge, ylab="share replaced", xlab="miles", alpha=0.9, legend = :topleft)
    di[:disc_miles_repl34] = @df filter(:mean_replaced => x -> x .> 0, bh) plot(:x, :mean_replaced, group=:bustype, title="avg replacement indicator", ylab="share replaced", xlab="milage state", alpha=0.9, legend = :topleft)
    @df bh scatter(:x, :mean_replaced, group=:bustype,xlab="mileage state", ylab="share replaced")

    di[:disc_miles_states] = @df bh groupedbar(:x, :nrow, group=:bustype, bar_position=:stack, xlab="mileage state", ylab="number of buses", title="mileage states by bus groups", bar_width=1)
    di[:disc_miles_states_hist] = @df bh groupedhist(:x, group=:bustype, bar_position=:dodge, xlab="miles",bar_width = 7,bins = 15)

    di
end

@doc raw"""
log likelihood function

```math
l_n(\theta,EV_\theta) = \log \mathcal{L}_n(\theta,EV_\theta) = \sum_{i=1}^{162}\sum_{t=2}^{T_i} \left( \log P(d_{i,t}|x_{i,t}) + \log \pi(x_{i,t}|x_{i,t-1},d_{i,t-1}) \right)
```
"""
function loglikelihood(z,p_model,
                       replace_data::BitArray,
                       miles_data::Vector{Int64},
                       miles_increase_data::Vector{Int64};do_θ = false)

    # check data arrays are consistent
    nt = length(replace_data)
    @assert(nt == length(miles_data))
    @assert(nt == length(miles_increase_data))

    # 1. get model-implied Probability of replace at miles_data
    prob_replace = p_model[miles_data]  
    prob_keep    = 1 .- prob_replace 

    # 2. for each observed discrete choice (replace or not)
    # compute the model implied probability of that happening
    logL = log.(prob_keep .* (.!(replace_data)) .+ (prob_replace .* replace_data))

    # 3. compute likelihood of mileage transition
    # adjust \theta for last element
    # this is tricky because implicitly need to enforce constraint that all p sum to 1
    # can get negative numbers here for the last element.
    if do_θ
        p_miles = clamp!([z.θ ; 1 .- sum(z.θ)], 0.0, 1.0)  # the clamp!() constrains the values to [0,1] mechanically
        logL = logL .+ log.(p_miles[miles_increase_data .+ 1])
    end
        
    mean((-1) * logL)
end

"""
    nested_likelihood(x::Vector{Float64}, h::Harold, d::DataFrame)

Outer loop of NFXP proceedure. takes candidate vector `x` of parameters, solves model, and evaluates the log likelihood.
"""
function nested_likelihood(x::Vector{Float64}, h::Harold, d::DataFrame)
    # update Harold
    h.RC = x[1]
    h.c  = x[2]
    if length(x) > 2
        h.θ  = x[3:end]
        h.transition = make_trans(h.θ,h.n)
    end

    # compute structural model 
    sol, iters = vfi(h)
    pr  = ccp(h, sol)

    # evaluate likelihood function 
    loglikelihood(h, pr, d.d, d.x, d.dx1, do_θ = length(x) > 2)
end


"""
    busdata(z::Harold; bustype = 4)

Data Loader.
"""
function busdata(z::Harold; bustype = 4) 
    d = CSV.read(joinpath(@__DIR__,"buses.csv"), DataFrame, header = false)
    select!(d, 1,2,5,7)
    rename!(d, [:id, :bustype, :d1, :odometer])

    d = filter(x -> x.bustype .<= bustype, d)

    # discretize odometer
    transform!(d, :odometer => (x -> Int.(ceil.(x .* z.n ./ (z.maxmiles * 1000)))) => :x)

    # replacement indicator
    dd = [d.d1[2:end] ; 0]

    # mileage increases
    dx1 = d.x .- [0;d.x[1:(end-1)]]
    dx1 = dx1 .* (1 .- d.d1) .+ d.x .* d.d1

    # make new dataframe
    df = [select(d, :id, :x, :bustype) DataFrame(dx1 = dx1, d = BitArray(dd))]

    # get rid of first observation for each bus
    idx = df.id .== [0; df.id[1:end-1]]
    df = df[idx,:]
end

function partialMLE(dx1::Vector) 
    c = counts(dx1) ./ length(dx1)
    c[1:(end-1)]  # forget about smallest category
end

"""
    nfxp(; β = 0.9, is_silent = false, doθ = false,n = 175, θ = [0.107, 0.5152 ,0.3622, 0.0143,0.0009])

Nested Fixed Point estimation. Can choose to do partial MLE (don't estimate the θs for mileage transition).
"""
function nfxp(; β = 0.9, is_silent = false, doθ = false,n = 175, θ = [0.107, 0.5152 ,0.3622, 0.0143,0.0009])
    z = Harold(RC = 0.0, c = 0.0, β = β, n = n, θ = θ)

    # get data
    d = busdata(z)
    # get initial transition probabilities
    if doθ
        p_0 = partialMLE(d.dx1)
        x_0 = [0.0, 0.0, p_0...]  # starting value with thetas
    else
        x_0 = [0.0, 0.0]  # starting value with thetas
    end

    # # optimize likelihood which calls the structural model

    r = Optim.optimize( x -> nested_likelihood(x, z, d), x_0 ,BFGS(),  Optim.Options(show_trace  = !is_silent ))
    o = Optim.minimizer(r)
    if doθ
        (RC = o[1], θc = o[2], θp = o[3:end])
    else
        (RC = o[1], θc = o[2])
    end

end


"""
    mpec(; β = 0.9, is_silent = false, doθ = false,n = 175, θ = [0.107, 0.5152 ,0.3622, 0.0143,0.0009])

MPEC estimation. Can choose to do partial MLE (don't estimate the θs for mileage transition).
"""
function mpec(; β = 0.9, is_silent = false, doθ = false,n = 175, θ = [0.107, 0.5152 ,0.3622, 0.0143,0.0009])

    h = Harold(β = β, n = n, θ = θ)
    N = h.n # number of mileage states

    d = busdata(h)
    p_0 = partialMLE(d.dx1)  # get a starting value for mileage
    J = length(p_0) + 1 # the partialMLE function discards the smallest bin, but later on sums up all probs

    # add time period index to each bus 
    gd = groupby(d, :id)
    dd = transform(gd, :id => (x -> 1:length(x)) => :it)
    icounter = combine(gd, nrow)
    I = nrow(icounter)

    # Jump model
    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, MOI.Silent(), is_silent)

    # variables
    @variable(m, θc >= 0 , start = 0.0)
    @variable(m, RC >= 0 , start = 0.0)
    if doθ
        p_0 = partialMLE(d.dx1)  # get a starting value for mileage
        J = length(p_0) + 1 # the partialMLE function discards the smallest bin, but later on sums up all probs 
        @variable(m, θlast[1:(J)] >= 0)
        set_start_value.(θlast[1:(J-1)], p_0)
        set_start_value(θlast[J], 0.0)
        @constraint(m, sum(θlast) == 1)

    else
        J = length(h.θ) + 1
        θlast = [h.θ; 1 - sum(h.θ)]
    end
    @variable(m, -50.0 <= EV[1:N] <= 50.0)  # need to help solver here!

    # expressions
    @NLexpression(m, opcost[i = 1:N], -0.001 * h.mileage[i] * θc)
    @NLexpression(m, VK[i = 1:N], opcost[i] + h.β * EV[i])  # value of keep
    @NLexpression(m, VR         , -RC + opcost[1] + h.β * EV[1])  # value of replace
    @NLexpression(m, diffV[i = 1:N], VR - VK[i] )  # payoff difference
    @NLexpression(m, pkeep[i = 1:N], 1 / (1 + exp(diffV[i])))
   
    # objective function
    # This is the likelihood function 
    if doθ
        @NLobjective(m, Max, 
        sum( log( (gd[i][it,:d]==false) * pkeep[ gd[i][it,:x] ] + (gd[i][it,:d]==true) * (1 - pkeep[ gd[i][it,:x] ]) )
             + log(θlast[gd[i][it,:dx1] + 1]) 
        for i in 1:I, it in 1:icounter[i,:nrow] )  )

    else
        @NLobjective(m, Max, 
        sum( log( (gd[i][it,:d]==false) * pkeep[ gd[i][it,:x] ] + (gd[i][it,:d]==true) * (1 - pkeep[ gd[i][it,:x] ]) )
             for i in 1:I, it in 1:icounter[i,:nrow] )  )
    end

    # 1. when state can move up any number of slots: easy
    @NLconstraint(m, evcon[i = 1:(N-J+1)], 
    EV[i] == sum( log(exp(VK[i + j]) + exp(VR)) * θlast[j+1] for j in 0:(J-1)  ))

    # 2. not all state progressions are possible: need       
    @NLconstraint(m, evconJ[i = (N-J+2):(N-1)], 
        EV[i] == sum(log( exp(VK[i + j]) + exp(VR)) * θlast[j+1] for j in 0:(N-i-1)) + 
                     log( exp(VK[N])     + exp(VR)) * ( 1 - sum( θlast[k+1] for k in 0:(N-i-1)) ) 
                    # these are all equivalent:
                    # log( exp(VK[N])     + exp(VR)) * sum( θlast[k+1] for k in (N-i):(J-1) ) 
                    #  sum(log( exp(VK[N])     + exp(VR)) *  θlast[k+1] for k in (N-i):(J-1) ) 
        )
    
    # 3. bellman equation at the final state
    @NLconstraint(m, evconN,
        EV[N] == log(exp(VK[N]) + exp(VR))
    )

    JuMP.optimize!(m)

    if doθ
        (RC = value(RC), θc = value(θc), θp = value.(θlast))
    else
        (RC = value(RC), θc = value(θc))
    end
    

end

function runall()
    d = Dict()
    for ns in [(90, [0.3489, 0.6394]), (175, [0.107, 0.5152 ,0.3622, 0.0143,0.0009])]
        d["n=$(ns[1])"] = Dict()
        for βs in [0.99, 0.995, 0.999]
            @info "running for n=$(ns[1]), β=$βs"
            d["n=$(ns[1])"]["β=$βs"] = Dict()
            t1 = @elapsed r1 = mpec(β = βs, doθ = false, is_silent = true,n = ns[1], θ = ns[2])
            t2 = @elapsed r2 = nfxp(β = βs, doθ = false, is_silent = true,n = ns[1], θ = ns[2])
            t3 = @elapsed r3 = mpec(β = βs, doθ = true, is_silent = true,n = ns[1])
            t4 = @elapsed r4 = nfxp(β = βs, doθ = true, is_silent = true,n = ns[1])

            d["n=$(ns[1])"]["β=$βs"][:mpec]      = Dict(:time => t1, :res => r1)
            d["n=$(ns[1])"]["β=$βs"][:nfxp]      = Dict(:time => t2, :res => r2)
            d["n=$(ns[1])"]["β=$βs"][:mpec_full] = Dict(:time => t3, :res => r3)
            d["n=$(ns[1])"]["β=$βs"][:nfxp_full] = Dict(:time => t4, :res => r4)
        end
    end
    
    d
end

end



