include("../src/competitiveSelection.jl")
using .CompetitiveSelection
using Statistics, JLD2

# Find space of variant numbers

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

nMax_μ_α = Array{Float64,2}(undef, length(_μ), length(_α))
@time for (i,α) in enumerate(_α)
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
        nMax_μ_α[j,i] = nDetectAv_t[80]
    end
end

jldsave("cloneNumbersAlphaMu.jld2"; nMax_μ_α, _μ, _α)