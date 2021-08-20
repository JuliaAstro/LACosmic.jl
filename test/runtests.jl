using LACosmic
using Test

@testset "subsample" for T in [Int, Float64, Float32, BigFloat, BigInt]
    data = T.(reshape(1:25, 5, 5))
    sub = @inferred LACosmic.subsample(data)

    @test size(sub) == (10, 10)
    @test parent(sub) === data

    @inbounds for idx in CartesianIndices(data)
        idx_11 = CartesianIndex(idx.I[1] * 2 - 1, idx.I[2] * 2 - 1)
        @test data[idx] === @inferred sub[idx_11]
        sub[idx_11] = 0
        @test data[idx] == 0

        idx_12 = CartesianIndex(idx.I[1] * 2 - 1, idx.I[2] * 2)
        @test data[idx] === @inferred sub[idx_12]
        sub[idx_12] = 0
        @test data[idx] == 0

        idx_21 = CartesianIndex(idx.I[1] * 2, idx.I[2] * 2 - 1)
        @test data[idx] === @inferred sub[idx_21]
        sub[idx_21] = 0
        @test data[idx] == 0

        idx_22 = CartesianIndex(idx.I[1] * 2, idx.I[2] * 2)
        @test data[idx] === @inferred sub[idx_22]
        sub[idx_22] = 0
        @test data[idx] == 0
    end

end
