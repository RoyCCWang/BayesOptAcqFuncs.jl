module BayesOptAcqFuncs

import Statistics

using LinearAlgebra
#using Serialization

import RKHSRegularization
const RK = RKHSRegularization

# import SingleLinkagePartitions
# const SL = SingleLinkagePartitions

import SpecialFunctions


# constant values.
function twopi(::Type{Float32})::Float32
    return 6.2831855f0 #convert(T, 2*π)
end

function twopi(::Type{Float64})::Float64
    return 6.283185307179586 #convert(T, 2*π)
end

function twopi(::Type{T})::T where T <: AbstractFloat
    return convert(T, 2*π)
end

function onepi(::Type{Float32})::Float32
    return 3.1415927f0
end

function onepi(::Type{Float64})::Float64
    return 3.141592653589793
end

function onepi(::Type{T})::T where T <: AbstractFloat
    return convert(T, π)
end

function sqrttwo(::Type{Float32})::Float32
    return 1.4142135f0
end

function sqrttwo(::Type{Float64})::Float64
    return 1.4142135623730951
end

function sqrttwo(::Type{T})::T where T <: AbstractFloat
    return convert(T, sqrt(2))
end


include("types.jl")
include("preprocessing.jl")

#include("./constraints/box.jl")
include("./acq/basic.jl")
include("./acq/frontend.jl")

export optimacq

end
