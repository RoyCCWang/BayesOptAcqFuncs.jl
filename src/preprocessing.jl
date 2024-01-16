# move to active learning repo.

abstract type ObjectiveInfo end

struct Minimization{T <: AbstractFloat} <: ObjectiveInfo
    lb::T
    ub::T
    new_ub::T
end


# construct an objective for a maximization problem, with values in [0,1].
function setupobjective(info::Minimization{T}, f) where T

    lb, ub, s = info.lb, info.ub, info.new_ub
    b = -s/(ub-lb)

    h = xx->min2max(f, xx, ub, b)::T

    return h
end
# h = xx->min2max(sin, xx, 1.0, -1.23/2)
# X = LinRange(-10, 10, 5000);
# @show maximum(h.(X)), minimum(h.(X))
# PythonPlot.plot(X, y, X, y2)

function min2max(f::Function, x, ub::T, b::T)::T where T <: AbstractFloat
    f_x::T = f(x)
    return (f_x - ub)*b
end

