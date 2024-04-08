module MullerPlot

export mullerBounds, buildFamilyArray, variantBirthTimes

function buildFamilyArray(parentVid_vid::Vector{T} where T<:Integer)
    childVid_child_Vid0 = [Int64[] for _ in 1:(length(parentVid_vid)+1)]
    for (vid, parId) in enumerate(parentVid_vid)
        push!(childVid_child_Vid0[1+parId], vid)
    end
    return childVid_child_Vid0
end

function sizeVariantRec(vid, n_vid, _child_Vid)
    childSizes = 0
    for child in _child_Vid[vid]
        childSizes += sizeVariantRec(child, n_vid, _child_Vid)
    end
    n_vid[vid] + childSizes
end

function getLowerCoordinateFromParent(parent, xL_vid, xU_vid, n_vid, _child_Vid)
    Δparent = xU_vid[parent]-xL_vid[parent]
    Δchild = sizeVariantRec(parent, n_vid, _child_Vid) - n_vid[parent]
    xL_vid[parent] + (Δparent-Δchild)/2
end

function setCoordinates!(parent, xL_vid, xU_vid, n_vid, _child_Vid)
    x0 = getLowerCoordinateFromParent(parent, xL_vid, xU_vid, n_vid, _child_Vid)
    for child in _child_Vid[parent]
        xL_vid[child] = x0
        xU_vid[child] = x0 + sizeVariantRec(child, n_vid, _child_Vid)
        x0 = xU_vid[child]
    end
    for child in _child_Vid[parent]
        setCoordinates!(child, xL_vid, xU_vid, n_vid, _child_Vid)
    end
end

"""
    mullerBounds(n_t_vid, _child_Vid)

Return lower and upper bounds `xL_t_vid` and `xU_t_vid` of variant bands for creating a Muller plot. Takes two arguments:
    `n_t_vid`: two dimensional Array whose elements are variant clone sizes. The first dimension `_t` is the timepoints at which the measurements occured; the second dimension `_vid` identifies the variant (i.e. the variant id). The first index refers to the wild type.
    `_child_Vid`: A nested array indexed by the variant id's. Its elements are 1D Vectors which contain the variant id's of each variant's children.
"""
function mullerBounds(n_t_vid::Array{T,2} where T<:Real, vid_child_Vid::Vector{Vector{V}} where V<:Integer)
    _t = 1:size(n_t_vid,1)
    xL_t_vid = Array{Float64}(undef, length(_t), size(vid_child_Vid,1))
    xU_t_vid = Array{Float64}(undef, length(_t), size(vid_child_Vid,1))
    xL_t_vid[:,1] .= 0
    xU_t_vid[:,1] .= 1
    for t in eachindex(_t)
        setCoordinates!(1, (@view xL_t_vid[t,:]), (@view xU_t_vid[t,:]), n_t_vid[t,:], vid_child_Vid)
    end
    return xL_t_vid, xU_t_vid
end

function mullerBounds(n_t_vid::Array{T,2} where T<:Real, vid_child_Vid::Vector{Vector{V}} where V<:Integer, visible_vid::AbstractVector{Bool})
    vidT_child_VidT = reduceCloneSpace(vid_child_Vid, visible_vid)
    mullerBounds(n_t_vid[:, visible_vid], vidT_child_VidT)
end

function mullerBounds(n_t_vid::Array{T,2} where T<:Real, parentVid_vid::Vector{V} where V<:Integer)
    childVid_child_Vid = buildFamilyArray(parentVid_vid)
    mullerBounds(n_t_vid, childVid_child_Vid)
end

function mullerBounds(n_t_vid::Array{T,2} where T<:Real, parentVid_vid::Vector{V} where V<:Integer, visible_vid::AbstractVector{Bool})
    childVid_child_Vid = buildFamilyArray(parentVid_vid)
    mullerBounds(n_t_vid, childVid_child_Vid, visible_vid)
end

function variantBirthTimes(n_t_vid)
    tBirth_vid = Vector{Int64}(undef, size(n_t_vid,2))
    for (vid, n_t) in enumerate(eachcol(n_t_vid))
        tBirth_vid[vid] = findfirst(n_t .> 0)
    end
    return tBirth_vid
end

function transformVidCoordinate(vidIn::Integer, vid1_vid2)
    for (vid2, vid1) in enumerate(vid1_vid2)
        if vid1==vidIn
            return vid2
        end
    end
    return nothing
end

function reduceCloneSpace(vid_child_Vid, visible_vid::AbstractVector{Bool})
    nClones2 = sum(visible_vid)
    _vid = range(0,length(visible_vid)-1)
    childVid2_child_Vid2 = Vector{Vector{Int}}(undef, nClones2)
    for (parentVid2, vid_child) in enumerate(vid_child_Vid[visible_vid])
        childVid2_child = Int[]
        for childVid in vid_child
            childVid2 = transformVidCoordinate(childVid, _vid[visible_vid])
            isnothing(childVid2) && continue
            push!(childVid2_child, childVid2)
        end
        childVid2_child_Vid2[parentVid2] = childVid2_child
    end
    return childVid2_child_Vid2
end

struct Variant
    vid::Integer
end

end

