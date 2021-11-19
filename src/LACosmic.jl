module LACosmic

using Statistics
using StaticArrays
using StaticKernels
using StaticKernels: ExtensionConstant, ExtensionReplicate

export lacosmic


include("subsample.jl")
include("core.jl")


end
