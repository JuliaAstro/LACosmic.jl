# subsample arrays without copying
using Base: @propagate_inbounds
import Base: getindex, setindex!

struct SubsampledArray{T,N,PT<:AbstractArray{T,N}} <: AbstractArray{T,N}
    parent::PT
    block_size::Int
end

"""
    LACosmic.subsample(array, block_size=2)

Subsample `array` by the given factor without copying or allocating. This effectively treats each pixel as `block_size` x `block_size` pixels. The value of each pixel is not normalized or averaged. This is a view into the parent array, and if the data was changed this view would change subsequently.

# Examples

```jldoctest
julia> arr = [1 2; 3 4]

julia> sub = LACosmic.subsample(arr)
4×4 LACosmic.SubsampledArray{Int64, 2, Matrix{Int64}}:
 1  1  2  2
 1  1  2  2
 3  3  4  4
 3  3  4  4

julia> all(sub[1:2, 1:2] .=== arr[1, 1])
true

julia> sub[3, 3] = 5;

julia> sub
4×4 LACosmic.SubsampledArray{Int64, 2, Matrix{Int64}}:
 1  1  2  2
 1  1  2  2
 3  3  5  5
 3  3  5  5
```
"""
subsample(parent, block_size=2) = SubsampledArray(parent, block_size)

Base.parent(arr::SubsampledArray) = arr.parent
Base.size(arr::SubsampledArray) = size(parent(arr)) .* arr.block_size

@inline block_idx(idx, size) = ceil(Int, idx / size)

@propagate_inbounds function getindex(arr::SubsampledArray{T,N}, idxs::Vararg{Int,N}) where {T,N}
    scaled_idxs = block_idx.(idxs, arr.block_size)
    return getindex(parent(arr), scaled_idxs...)
end

@propagate_inbounds function setindex!(arr::SubsampledArray{T,N}, value, idxs::Vararg{Int,N}) where {T,N}
    scaled_idxs = block_idx.(idxs, arr.block_size)
    return setindex!(parent(arr), value, scaled_idxs...)
end
