using JLD2
using GpABC
using CairoMakie

##

# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_particles500.jld2", "abcRes", "ctrlParams")
abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages6575_particles500.jld2", "abcRes", "ctrlParams")
# abcRes, ctrlParams = load("./Data/SimResults/abcResult_kolmogorov-smirnov_ages7080_cumulative_particles500.jld2", "abcRes", "ctrlParams")

##

ctrlParams[:metric]
ctrlParams[:tBounds]

##
myTheme = Theme(fontsize=20)
set_theme!(myTheme)

sets = length(abcRes.n_tries)
set = 3

fig = Figure(resolution=(1300,800))
Axis(fig[1,1],
    xlabel="s",
    ylabel="number of accepted particles"
)
hist!(abcRes.population[set][:,1])
Axis(fig[1,2],
    xlabel="σ",
)
hist!(abcRes.population[set][:,2])
Axis(fig[1,3],
    xlabel="μ",
)
hist!(abcRes.population[set][:,3])
Label(fig[0,:], 
    string(ctrlParams[:metric])*" distance of size distributions; ages "*
    string(ctrlParams[:tBounds][1])*"-"*string(ctrlParams[:tBounds][2]), 
    fontsize=30
)

Axis(fig[2,1],
    title="s",
    xlabel="σ",
    ylabel="μ"
)
scatter!(abcRes.population[set][:,2], abcRes.population[set][:,3])

Axis(fig[2,2],
    title="σ",
    xlabel="s",
    ylabel="μ"
)
scatter!(abcRes.population[set][:,1], abcRes.population[set][:,3])

Axis(fig[2,3],
    title="μ",
    xlabel="s",
    ylabel="σ"
)
scatter!(abcRes.population[set][:,1], abcRes.population[set][:,2])

display(fig)

##
abcRes.population[set]