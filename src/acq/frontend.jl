
######## 

abstract type PostProcessOption end

# maximize each quantile interval, which is diff([1; probs]).
struct MaxQuantileInfo{T} <: PostProcessOption
    probs::T # quantile probabilities.
end

function postacq(info::MaxQuantileInfo, agents::Vector{Vector{T}}, rewards::Vector{T}) where T
    inds = maxinquantile(rewards, info.probs)
    return agents[inds], rewards[inds]
end

function postacq(::Nothing, x::Vector{Vector{T}}, y::Vector{T}) where T
    return x, y
end

############

# so far it is just optimacq() for front end.

# # load a GP model from disk, find maximum.
# function maximizeacq(
#     acq::AcquisitionStrategy,
#     load_path::String, # no error checking.
#     save_path::String, # no error checking.
#     constraints::ConstraintContainer,
#     runminimization,
#     optim_config,
#     post_process::Union{PostProcessOption, Nothing};
#     verbose = false
#     )
    
#     mgp = deserialize(load_path)

#     x_star, agents, rewards = optimacq(
#         acq,
#         mgp,
#         constraints,
#         runminimization,
#         optim_config,
#     )
#     rewards = -rewards

#     # get the candidates, xs, and their acquisition function evals, as.
#     xs, as = postacq(post_process, agents, rewards)

#     # save to disk..
#     serialize(save_path, (xs,as))
    
#     if verbose
#         println("saved $save_path")
#     end

#     return nothing
# end