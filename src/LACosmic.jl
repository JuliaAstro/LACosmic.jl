module LACosmic

using BlockArrays
using FillArrays
using Parameters
using Photometry.Background: LocationEstimator, RMSEstimator
using StaticArrays
using StaticKernels


function lacosmic(data::AbstractMatrix, contrast=2.0, sigma_clip=4.5, objlim=2, niter=4, gain=nothing, readnoise=nothing)
    block_size = 2
    subimg = extend(subsample(data, block_size), StaticKernels.ExtensionSymmetric())
    kernel = @kernel w -> max(0, w[0,-1] + w[-1,0] - 4*w[0,0] + w[1,0] + w[0,1])
    conv_img = similar(data, size(kernel, subimg))
    laplace_img = similar(data)

    for iteration in 1:maxiter
        map!(kernel, conv_img, subimg)
        block_reduce!(sum, laplace_img, conv_img)
    end

end


function block_reduce!(op, out, sub; block_size=2)
    block_dims = @SVector [Fill(block_size, size(sub, dim) ÷ block_size) for dim in 1:ndims(sub)]
    block_arr = PseudoBlockArray(sub, block_dims...)
    return map!(op, out, eachblock(block_arr))
end


end
