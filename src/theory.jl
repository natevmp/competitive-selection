module Theory

export fLogistic, fIndFitnessWithAvCompetition, fCompetition

# fLogistic(t, t0, r) = 1 / ( 1 + (1 - 1/Nf)/(1/Nf) * exp(-r*(t-t0)) )
# fAntiLog(t, t0, r, b0) = 1 / ( 1 + (1-b0)/b0 * exp(-r*(t-t0)) )

fLogistic(t, t0, r, x0) = 1 / ( 1 + (1-x0)/x0 * exp(-r*(t-t0)) )

fLogistic(t, t0, r, x0, xF) = xF / ( 1 + (xF-x0)/x0 * exp(-r*(t-t0)) )

fLogistic(t, t0, r, x0, xF, z) = (xF-z/r) / ( 1 + (xF-z/r-x0)/x0 * exp(-(r-z)*(t-t0)) )

compCoverage(t, α, s, μ, N) = 1 - (1 + μ/(N*α*s))/( 1 + (μ/(N*α*s))*(exp(t*(α*s + μ/N))) )

fCompetition(t, t0, α, s, μ, N) = (1/N + exp(-(α*s + μ/N)*t0)*α*s/μ) / (1 + α*s/μ*N*exp(-(α*s + μ/N)*t))

fIndFitnessWithAvCompetition(t, t0, α, s, a, μ, N) = exp( -α*(a-s)*(t-t0) )*( α*a*exp(-(α*a+μ/N)*t0) +  μ/N ) /
    ( N * (α*a*exp(-(α*a+μ/N)*t) + μ/N) )


end
