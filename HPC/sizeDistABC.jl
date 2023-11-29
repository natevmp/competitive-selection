include("../src/dataStructuring.jl")
using .DataStructuring
include("../src/analysisTools.jl")
using .AnalysisTools
include("../src/abcSims.jl")
using .ABCSims
using JLD2, CSV, DataFrames


function main()
    
    include(ARGS[1])
    
    println("Available threads: ", Threads.nthreads())

    # ----- Data --------
    df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
    DataStructuring.structureData!(df)
    dfVid = DataStructuring.sortDataPerVariant(df)
    x_sid = [x for x in Iterators.flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
    t_sid = [t for t in Iterators.flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
    xMask_sid = x_sid.>0
    _fEdges, nData_f = AnalysisTools.sizeDistDens(
        x_sid[xMask_sid], t_sid[xMask_sid], tBounds; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax]
    )
    nCumData_f = [sum(nData_f[1:i]) for i in eachindex(nData_f)]
    data_f = ctrlParams[:cumulativeDist] ? nCumData_f : nData_f

    # ----- ABC --------
    @time abcRes = ABCSims.runABC(data_f, ctrlParams)

    # --------- save data -------------
    println("saving data")
    filename = begin
        "abcResult_Sampled_"*string(ctrlParams[:metric])*
        "_ages"*string(tBounds[1])*string(tBounds[2])*
        # "_cumulative"*string(ctrlParams[:cumulativeDist])*
        (ctrlParams[:cumulativeDist] ? "_cumulative" : "")*
        (!haskey(ctrlParams, :sBounds) ? "_SFixed" : "")*
        (!haskey(ctrlParams, :σBounds) ? "_SFixed" : "")*
        "_particles"*string(ctrlParams[:nParticles])*
        ".jld2"
    end
    jldsave(filename; abcRes, ctrlParams)

end

main()