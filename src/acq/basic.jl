# assumes we want to maximize f, with the global solution above zero.
# This is because it is easy to model f(blah) = 0 with zero-mean GPs prior.


function stdCDF(x::T)::T where T
    return SpecialFunctions.erfc(-x/sqrttwo(T))/2
end

function stdPDF(x::T)::T where T
    return 1/sqrt(twopi(T)) * exp(-(x^2)/2)
end

# noiseless expected improvement.

abstract type AcquisitionStrategy end
abstract type ElementaryAcq <: AcquisitionStrategy end
abstract type CompoundAcq <: AcquisitionStrategy end

struct ExpectedImprovement{T <: AbstractFloat} <: ElementaryAcq
    f_best::T
end

function ExpectedImprovement(y::Vector{T}) where T <: AbstractFloat
    return ExpectedImprovement(maximum(y))
end

function evalacq(
    info::ExpectedImprovement{T},
    μ::T,
    σ::T,
    )::T where T <: AbstractFloat

    #μ, σ = getposterior(x)

    z = (μ-info.f_best)/σ

    c_z = stdCDF(z)
    p_z = stdPDF(z)

    return σ*(z*c_z + p_z)
end

struct ProbabilityOfImprovement{T <: AbstractFloat} <: ElementaryAcq
    f_best::T
end

function ProbabilityOfImprovement(y::Vector{T}) where T <: AbstractFloat
    return ProbabilityOfImprovement(maximum(y))
end

function evalacq(
    info::ProbabilityOfImprovement{T},
    μ::T,
    σ::T,
    )::T where T <: AbstractFloat

    #μ, σ = getposterior(x)

    z = (info.f_best - μ)/σ

    c_z = stdCDF(z)

    return c_z
end

struct UpperConfidenceBound{T} <: ElementaryAcq
    β::T
end

# Eq. 45 from (Gutmann 2016) D is dimension, M is number of points.
function UpperConfidenceBound(ϵ::T, D::Int, M::Int) where T # default to ϵ = 0.1.
    tmp = convert(T, M^((D/2)+2) * twopi(T) /(3*ϵ))
    β = 2*log(tmp)
    return UpperConfidenceBound(β)
end

function evalacq(info::UpperConfidenceBound{T}, μ_x::T, σ_x::T)::T where T <: AbstractFloat
    return μ_x + info.β*σ_x
end

struct MaskedAcquisition{MT <: Function, AT <: ElementaryAcq} <: CompoundAcq
    mask::MT # mask: (μ_x, σ_x) |-> a real number.
    acq::AT
end

# abstract type BatchEIMethod end
# struct SPSA{CT} <: BatchEIMethod
#     SPSAConfig::CT
#     seed::Int
# end

# function maxbatchEI(config::SPSA, x::Vector{T}, f_best::T)::T where T <: AbstractFloat
#     #
#     f = hh->EI(x, min(f_best,hh))
#     runSPSA(config, x, f_best)

# end

####################### for GPs.

function evalacq(
    info::CompoundAcq,
    η::RK.Problem,
    x::Vector{T},
    args...
    ) where T

    acq_eval = evalacq(info.acq, η, x, args...)
    mask_eval::T = info.mask(x)

    return acq_eval*mask_eval
end

function evalacq(
    info::ElementaryAcq,
    η::RK.Problem,
    x::Vector{T},
    ) where T

    m, v = RK.queryGP(x, η)
    return evalacq(info, m, sqrt(v))
end

function evalacq(
    info::ElementaryAcq,
    η::RK.Problem,
    x::Vector{T},
    s::T,
    ) where T

    m, v = RK.queryGP(x, η)
    m = s*m

    #return evalacq(info, m, v)
    return evalacq(info, m, sqrt(v))
end

##### maximize acq function.

function optimacq(
    info::AcquisitionStrategy,
    η::RK.Problem,
    C::BoxConstraint{T},
    runminimization,
    config,
    ) where T <: AbstractFloat

    lbs, ubs = C.lbs, C.ubs

    h = xx->(-evalacq(info, η, xx))
    return runminimization(h, lbs, ubs, config)
end

############ for use with an optimization or samplinga lgorithm that requns a population of the acqusition function.

# use this to get a variety of points, given their acquisition function evals, y.
function maxinquantile(y::Vector{T}, probs::Union{Vector{T}, LinRange}) where T
    
    inds = Vector{Int}(undef, length(probs))

    Q_ub = Statistics.quantile(y, one(T))
    for i in eachindex(probs)

        q = Statistics.quantile(y, probs[i])

        
        _, k = findmax( 
            (q <= y[n] <= Q_ub) ? y[n] : convert(T, -Inf)
            for n in eachindex(y) #if q <= y[n] <= Q_ub
        )
        inds[i] = k
        Q_ub = q
    end

    inds = unique(inds)
    return inds
end


# function getvarietyviaquantile(rewards::Vector{T}, probs) where T
    
#     #probs = LinRange(0.9, 0.6, 10)
#     inds = maxinquantile(rewards, probs)
    
#     return inds
# end