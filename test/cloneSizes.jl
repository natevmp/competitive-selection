include("../src/competitiveSelection.jl")
using .CompetitiveSelection
using Statistics
using CairoMakie

##

paramsPop = Dict(
    :N => Int(1E5),
    :α => 1,
    # :s => 0.12,
    :s => 0.2,
    :T => 90,
    :μ => 10,
    # :μ => 0,
    :σ => 0.1,
    :sType => "fixed",
    # :sType => "free",
)
runs = 1000

##
@time solEns = CompetitiveSelection.evolvePop(paramsPop, runs=runs)
##
# visualize individual clones

simID = 3
# xTot_t = [sum(solEns[simID][t][solEns[simID][t] .> 0.01]) for t in 1:length(solEns[simID])]
fig = Figure(resolution=(.8*700, .8*600))
Axis(
    fig[1,1],
    xlabel="age",
    ylabel="variant size",
    title="selection type = "*paramsPop[:sType]*"; s = "*string(paramsPop[:s])*"; μ = "*string(paramsPop[:μ]),
)
for i in 1:length(solEns[simID][end])
    lines!([solEns[simID][t][i] for t in 1:length(solEns[simID])])
end
# lines!(xTot_t, color=:black, linestyle=:dash)
display(fig)
##
save("Figures/cloneSizesTime.png", fig, px_per_unit = 2)


##
_f, nV_t_f = CompetitiveSelection.variantsSizeDistribution(solEns, xMin=0.01, xMax=nothing, bins=100, normalized=true)
sum(nV_t_f[paramsPop[:T], :])

##

fig1 = Figure(resolution=(.8*700, .8*600))
fig1[1,1] = Axis(
    fig1[1,1],
    xlabel="vaf",
    ylabel="density of clones",
    # yscale=log10,
    # xscale=log10,
)
xlims!(0.01,0.15)
# ylims!(1E-3,1)
for t in [20,25,30,35,40,45]
    lines!(_f, nV_t_f[t,:], label="t = "*string(t))
end
axislegend()
display(fig1)

save("Figures/cloneSizeDistributions.png", fig1, px_per_unit = 2)

##

detectThreshold = 0.01

nDetect_sim_t = Array{Float64}(undef, runs, length(solEns[1]))
for simId in 1:runs
    for tInd in 1:paramsPop[:T]
        detectableVariants_vid = solEns[simId][tInd] .> detectThreshold # is a boolean vec
        nDetect_sim_t[simId,tInd] = sum(detectableVariants_vid)
    end
end
nDetectAv_t = vec(mean(nDetect_sim_t,dims=1))
nDetectMax_t = vec(maximum(nDetect_sim_t,dims=1))

##

fig2 = Figure(resolution=(.8*700, .8*600))
Axis(
    fig2[1,1],
    # yscale=Makie.pseudolog10,
    xlabel="age",
    ylabel="av number of variants > $detectThreshold",
    title="selection type = "*paramsPop[:sType]*"; s = "*string(paramsPop[:s])*"; μ = "*string(paramsPop[:μ]),
    # xscale=identity,
)
lines!(nDetectAv_t)
lines!(nDetectMax_t)
display(fig2)

save("Figures/cloneNumberTime.png", fig2, px_per_unit = 2)

##
# =========== Clonal fraction of expanded clones ==============

vafTotDetect_sim_t = Array{Float64}(undef, runs, length(solEns[simID]))
vafDetect_sim_t_Vid = Array{Vector{Float64}, 2}(undef, length(solEns), length(solEns[simID].t))
for simId in 1:runs
    for tInd in 1:paramsPop[:T]
        isDetectable_vid = solEns[simId][tInd] .> detectThreshold #boolean vec
        vafTotDetect_sim_t[simId, tInd] = sum(solEns[simId][tInd][isDetectable_vid])
        vafDetect_sim_t_Vid[simId,tInd] = solEns[simId][tInd][isDetectable_vid]
    end
end

##
# lines(vec(mean(vafTotDetect_sim_t,dims=1))[1:end-1])

##
simId = 3

_t = 0:params[:T]
_tBInd = 1:5:(params[:T])

t_T_Vid = Vector{Vector{Int}}(undef, length(_tBInd))
vaf_T_Vid = Vector{Vector{Float64}}(undef, length(_tBInd))
grp_T_Vid = similar(t_T_Vid)
for (i, tBInd) in enumerate(_tBInd)
    vaf_id = vafDetect_sim_t_Vid[simId, tBInd]
    t = _t[tBInd]
    if length(vaf_id)==0
        t_T_Vid[i] = [t,]
        vaf_T_Vid[i] = [0.,]
        grp_T_Vid[i] = [1,]
    else
        t_T_Vid[i] = fill(t, length(vaf_id))
        vaf_T_Vid[i] = sort(vaf_id, rev=true)
        grp_T_Vid[i] = Vector(1:length(vaf_id))
    end
end
t_barInd = vcat(t_T_Vid...)
vaf_barInd = vcat(vaf_T_Vid...)
grp_barInd = vcat(grp_T_Vid...)
#

fig4 = Figure(resolution=(.8*700, .8*600))
Axis(
    fig4[1,1],
    xlabel="age",
    ylabel="clonal fraction of clones > $detectThreshold",
)
barplot!(t_barInd, vaf_barInd,
    stack=grp_barInd,
    color=grp_barInd,
    colormap=:tab10,
    label="single simulation"
)
lines!(_t[1:end-1], vec(mean(vafTotDetect_sim_t,dims=1))[1:end-1],
    color=:black,
    label="average over simulations"
)
axislegend( position=:lt )
display(fig4)

save("Figures/clonalFractions_age.png", fig4, px_per_unit = 2)
##
# =============Find space of variant numbers ==================

detectThreshold = 0.01
runs = 1000
_μ = range(1,50,length=10)
_α = range(1,15,length=10)

params = Dict(
    :N => Int(1E5),
    :s => 0.3,
    :T => 90,
    :sType => "fixed",
)
##
nMax_μ_α = Array{Float64,2}(undef, length(_μ), length(_α))
@showprogress for (i,α) in enumerate(_α)
    for (j,μ) in enumerate(_μ)
        # @time sol = CompetitiveSelection.evolvePop(paramsPop)
        params[:α] = α
        params[:μ] = μ
        solEns = CompetitiveSelection.evolvePop(params, runs=runs)
        nDetect_sim_t = Array{Float64}(undef, runs, length(solEns[1]))
        _t = solEns[1].t
        for simId in 1:runs
            for tInd in 1:length(solEns[1].t)
                detectableVariants_vid = solEns[simId][tInd] .> detectThreshold
                nDetect_sim_t[simId,tInd] = sum(detectableVariants_vid)
            end
        end
        nDetectAv_t = vec(mean(nDetect_sim_t,dims=1))
        nMax_μ_α[j,i] = nDetectAv_t[end-1]
    end
end

##
using JLD2
@load "Data/cloneNumbersAlphaMu.jld2"


##
using Dierckx, ProgressMeter

nMaxSplF_α_μ = Spline2D(_α, _μ, transpose(nMax_μ_α))

_αi = range(_α[1],_α[end], length=100)
_μi = range(_μ[1], _μ[end],length=100)

nMaxInt_α_μ = evalgrid(nMaxSplF_α_μ, _αi, _μi)

##

fig3 = Figure(resolution = (.8*800, .8*500))
Axis(
    fig3[1,1],
    xlabel="μ",
    ylabel="α",
)
hm = heatmap!(_μi, _αi, transpose(nMaxInt_α_μ))
Colorbar(fig3[1,2], hm, label="expected number of variants")
display(fig3)

save("Figures/maxCloneNumber_alphaMu.png", fig3, px_per_unit = 2)

##
fig4 = Figure(resolution=(.8*700, .8*400))
Axis(
    fig4[1,1],
    # yscale=Makie.pseudolog10,
    xlabel="μ (population mutation rate)",
    ylabel="number of variants at CH",
)
lines!(_μi, nMaxInt_α_μ[1,:])
display(fig4)
save("Figures/cloneNumber_Mu.png", fig4, px_per_unit = 2)

## ==================== plot theme ================
using ColorSchemes
myTheme = Theme(
    resolution = (.8*800, .8*600),
)
set_theme!(myTheme)