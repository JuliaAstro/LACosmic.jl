module LACosmic

using Photometry.Background: LocationEstimator, RMSEstimator
using StaticArrays
using StaticKernels


include("subsample.jl")
include("core.jl")


end
