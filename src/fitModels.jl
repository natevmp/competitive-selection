
# -------------------------- Model definitions -----------------------
abstract type ModelShape end
mutable struct LogisticModel <: ModelShape 
    γ::Float64
    t0::Float64
    x0::Float64
    xF::Float64
    LogisticModel() = new(0.1, 1.,0.1, 0.5)
end
mutable struct MaxAdjustedLogisticModel <: ModelShape
    γ::Float64
    t0::Float64
    x0::Float64
    xF::Float64
    MaxAdjustedLogisticModel() = new(0.1, 1., 0.1, 0.25)
end

abstract type ExpansionType end
struct PositiveExpansion <: ExpansionType end
struct NegativeExpansion <: ExpansionType end

struct ModelFit{S<:ModelShape, E<:ExpansionType}
    shape::S
    expansion::E
end
ModelFit(shape, exp) = ModelFit{typeof(shape), typeof(exp)}(shape, exp)

# --------------------------- Initial values ---------------------------
# b0[1]=t0, b0[2]=γ
function setInitialParamsMLE!(model::ModelFit{LogisticModel, PositiveExpansion}, _t, vaf_t, Nf)
    model.shape.t0 = 20.
    model.shape.γ = 0.1
    model.shape.x0 = 1/(2Nf)
    return [model.shape.t0, model.shape.γ]
end
# b0[1]=t0, b0[2]=γ, b0[3]=x0
function setInitialParamsMLE!(model::ModelFit{LogisticModel, NegativeExpansion}, _t, vaf_t, Nf)
    model.shape.t0 = 1.
    model.shape.γ = -1.
    model.shape.x0 = vaf_t[1]
    return [model.shape.t0, model.shape.γ, model.shape.x0]
end
# b0[1]=t0, b0[2]=γ, b0[3]=xF
function setInitialParamsMLE!(model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion}, _t, vaf_t, Nf)
    model.shape.t0 = 20.
    model.shape.γ = 0.5
    model.shape.x0 = 1/(2Nf)
    model.shape.xF = 0.25
    return [model.shape.t0, model.shape.γ, model.shape.xF]
end
# b0[1]=t0, b0[2]=γ, b0[3]=x0, b0[4]=xF
function setInitialParamsMLE!(model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}, _t, vaf_t, Nf)
    model.shape.t0 = 20.
    model.shape.γ = -0.5
    model.shape.x0 = vaf_t[1]
    model.shape.xF = 0.4
    return [model.shape.t0, model.shape.γ, model.shape.x0, model.shape.xF]
end

# ----------------------- Box constraints for MLE ------------------------
function boundsMLE(model::ModelFit{LogisticModel, PositiveExpansion}, _t, vaf_t)
    # b0[1]=t0, b0[2]=γ
    ([-10, 0], 
    [_t[1], 5])
end

function boundsMLE(model::ModelFit{LogisticModel, NegativeExpansion}, _t, vaf_t)
    # b0[1]=t0, b0[2]=γ, b0[3]=x0
    ([-10, -10, 0],
    [_t[1], 0, 0.49])
end

function boundsMLE(model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion}, _t, vaf_t)
    # b0[1]=t0, b0[2]=γ, b0[3]=xF
    ([0, 0, 1E-3],
    [_t[1], 2, 0.49])
end

# function boundsMLE(model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}, _t, vaf_t)
#     ([-1E6, -Inf, 0, 1E-4],
#     [_t[1], 0, 0.45, 0.49])
# end

# ------------------ linear constraints for MLE -----------------
function constraintsMLE(model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}, _t)
    # b0[1]=t0, b0[2]=γ, b0[3]=x0, b0[4]=xF
    constraintF(res, x, p) = (res .= [x[1], x[2], x[3], x[4], x[3] - x[4]]) # constraint function
    return (
        constraintF,
        [0, -10, 0, 0, -0.5], # lower constraint
        [_t[1], 0, 0.48, 0.49, -0.01] # upper constraint
    )
end

# ------------------------- growth functions -------------------------
function fLogModel(t,β,model::ModelFit{LogisticModel, PositiveExpansion})
    fLogistic(t, β[1], β[2], model.shape.x0, 0.5)
end

function fLogModel(t,β,model::ModelFit{LogisticModel, NegativeExpansion})
    fLogistic(t, β[1], β[2], β[3], 0.5)
end

function fLogModel(t,β,model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion})
    fLogistic(t, β[1], β[2], model.shape.x0, β[3])
end

function fLogModel(t,β,model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion})
    fLogistic(t, β[1], β[2], β[3], β[4])
end 

# --------------------- parameter adjustment methods ---------------------------
function addFitParams!(model::ModelFit{LogisticModel, PositiveExpansion}, _β)
    model.shape.t0 = _β[1]
    model.shape.γ = _β[2]
end
function addFitParams!(model::ModelFit{LogisticModel, NegativeExpansion}, _β)
    model.shape.t0 = _β[1]
    model.shape.γ = _β[2]
    model.shape.x0 = _β[3]
end
function addFitParams!(model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion}, _β)
    model.shape.t0 = _β[1]
    model.shape.γ = _β[2]
    model.shape.xF = _β[3]
end
function addFitParams!(model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}, _β)
    model.shape.t0 = _β[1]
    model.shape.γ = _β[2]
    model.shape.x0 = _β[3]
    model.shape.xF = _β[4]
end

