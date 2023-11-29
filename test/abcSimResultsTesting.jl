using JLD2
using GpABC
using CairoMakie

##

# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages6070_particles500.jld2", "abcRes", "ctrlParams")
# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages6575_particles500.jld2", "abcRes", "ctrlParams")
# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages7080_particles500.jld2", "abcRes", "ctrlParams")

# abcRes, ctrlParams = load("./Data/SimResults/abcResult_kolmogorov-smirnov_ages6070_cumulative_particles500.jld2", "abcRes", "ctrlParams")
# abcRes, ctrlParams = load("./Data/SimResults/abcResult_kolmogorov-smirnov_ages6575_cumulative_particles500.jld2", "abcRes", "ctrlParams")
# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages7076_particles500.jld2", "abcRes", "ctrlParams")

abcRes, ctrlParams = load("./Data/SimResults/abcResult_Sampled_chisquare_ages7076_particles500.jld2", "abcRes", "ctrlParams")
# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages7076_particles500.jld2", "abcRes", "ctrlParams")


##

propertynames(abcRes)
abcRes.threshold_schedule
abcRes.n_tries
abcRes.n_accepted
abcRes.n_params
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
    string(ctrlParams[:metric])*" distance of size distributions, threshold="*string(abcRes.threshold_schedule[set])*"; ages "*
    string(ctrlParams[:tBounds][1])*"-"*string(ctrlParams[:tBounds][2]), 
    fontsize=30
)

Axis(fig[2,1],
    # title="s",
    xlabel="σ",
    ylabel="μ"
)
scatter!(abcRes.population[set][:,2], abcRes.population[set][:,3])

Axis(fig[2,2],
    # title="σ",
    xlabel="s",
    ylabel="μ"
)
scatter!(abcRes.population[set][:,1], abcRes.population[set][:,3])

Axis(fig[2,3],
    # title="μ",
    xlabel="s",
    ylabel="σ"
)
scatter!(abcRes.population[set][:,1], abcRes.population[set][:,2])

display(fig)

##
abcRes.population[set]


## ======================== Fixed fitness ===================================

# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages7080_SFixed_particles400.jld2", "abcRes", "ctrlParams")
abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages7076_SigmaFixed_particles500.jld2", "abcRes", "ctrlParams")

myTheme = Theme(fontsize=20)
set_theme!(myTheme)

abcRes.n_params
println("n tries ", abcRes.n_tries)
println("n accpeted ", abcRes.n_accepted)
println("thresholds ", ctrlParams[:threshold])
# ctrlParams
##
set = 3

fig = Figure(resolution=(1300,700))
Axis(fig[1,1],
    xlabel="σ",
    ylabel="number of accepted particles"
)
hist!(abcRes.population[set][:,1])
Axis(fig[1,2],
    xlabel="μ",
)
hist!(abcRes.population[set][:,2])
Label(fig[0,:], 
    string(ctrlParams[:metric])*" distance of size distributions, threshold="*string(abcRes.threshold_schedule[set])*"; ages "*
    string(ctrlParams[:tBounds][1])*"-"*string(ctrlParams[:tBounds][2]), 
    fontsize=30
)
Axis(fig[1,3],
    # title="s",
    xlabel="σ",
    ylabel="μ"
)
scatter!(abcRes.population[set][:,1], abcRes.population[set][:,2])

display(fig)

## ======================== Fixed fitness variance ===================================

# abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages7080_SFixed_particles400.jld2", "abcRes", "ctrlParams")
abcRes, ctrlParams = load("./Data/SimResults/abcResult_chisquare_ages7076_SigmaFixed_particles500.jld2", "abcRes", "ctrlParams")

myTheme = Theme(fontsize=20)
set_theme!(myTheme)

abcRes.n_params
println("n tries ", abcRes.n_tries)
println("n accpeted ", abcRes.n_accepted)
println("thresholds ", ctrlParams[:threshold])
# ctrlParams
##
set = 1

fig = Figure(resolution=(1300,700))
Axis(fig[1,1],
    xlabel="s",
    ylabel="number of accepted particles"
)
hist!(abcRes.population[set][:,1])
Axis(fig[1,2],
    xlabel="μ",
)
hist!(abcRes.population[set][:,2])
Label(fig[0,:], 
    string(ctrlParams[:metric])*" distance of size distributions, threshold="*string(abcRes.threshold_schedule[set])*"; ages "*
    string(ctrlParams[:tBounds][1])*"-"*string(ctrlParams[:tBounds][2]), 
    fontsize=30
)
Axis(fig[1,3],
    # title="s",
    xlabel="s",
    ylabel="μ"
)
scatter!(abcRes.population[set][:,1], abcRes.population[set][:,2])

display(fig)