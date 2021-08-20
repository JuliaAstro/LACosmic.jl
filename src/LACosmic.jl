module LACosmic

using Photometry.Background: LocationEstimator, RMSEstimator
using StaticArrays
# using StaticKernels

export lacosmic


include("subsample.jl")
include("core.jl")


end
