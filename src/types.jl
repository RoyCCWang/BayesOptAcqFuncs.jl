

abstract type ConstraintContainer end

struct BoxConstraint{T <: Real} <: ConstraintContainer
    lbs::Vector{T}
    ubs::Vector{T}
end
