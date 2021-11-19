using Distributions
using ImageFiltering
using LACosmic
using PSFModels: Gaussian
using StableRNGs
using Test

rng = StableRNG(125588)

function make_data(rng, T, N; N_sources=100, N_cosmics=100)
    imdata = fill(200.0, (N, N))

    # Add some fake sources
    for _ in 1:N_sources
        x = rand(rng, Uniform(1, N + 1))
        y = rand(rng, Uniform(1, N + 1))
        brightness = rand(rng, Uniform(1000, 30000)) / (2π * 3.5^2)
		model = Gaussian(T; x=x, y=y, fwhm=3.5, amp=brightness)
        imdata .+= model[axes(imdata)...]
	end

    # Add the poisson noise
    imdata .= rand.(rng, Poisson.(imdata))

    # Add readnoise
    imdata .+= rand(rng, Normal(0, 10), (N, N))

	clean_image = T.(imdata)
	
    # Add Nc fake cosmic rays
    crmask = falses((N, N))
	for i in 1:N_cosmics
    	cr_x = round(Int, rand(rng, Uniform(6, N - 5)))
    	cr_y = round(Int, rand(rng, Uniform(6, N - 5)))
    	cr_brightnesses = rand(rng, Uniform(1000, 30000))
    	imdata[cr_y, cr_x] += cr_brightnesses
    	crmask[cr_y, cr_x] = true
	end

    # Make a mask where the detected cosmic rays should be
    return (image=T.(imdata), clean_image, mask=crmask)
end

make_data(T, N; kwargs...) = make_data(rng, T, N; kwargs...)
make_data(N; kwargs...) = make_data(Float32, N; kwargs...)

@testset "LACosmic.jl" begin
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

    @testset "dilation - $k" for k in [3, 5]
        mask = randn(rng, 1001, 1001) .< 0.05
        kernel = ones(k, k)
        expected = mapwindow(any, mask, (k, k), border=Fill(false))
        dilated = @inferred LACosmic.dilate!(copy(mask), k)
        @test dilated == expected
    end

    @testset "rebin" for block_size in [2, 3]
        data = randn(rng, Float32, 512, 512) .+ 100
        sub = @inferred LACosmic.subsample(data, block_size)
        binned = @inferred LACosmic.rebin(sub, block_size)
        @test binned ≈ data

        binned_sum = @inferred LACosmic.rebin(sum, sub, block_size)
        @test binned_sum ≈ block_size^2 .* data
    end

    @testset "cleaning" begin
        input = make_data(1001)
        output, mask = lacosmic(input.image; sigma_clip=6, neighbor_thresh=1, readnoise=10)
        output2, mask2 = lacosmic(output; sigma_clip=6, neighbor_thresh=1, readnoise=10)
        @test iszero(count(mask2))
    end

    @testset "LACosmic algorithm" for T in [Float32, Float64]
        input = make_data(T, 1001)
        output, mask = @inferred lacosmic(input.image; sigma_clip=6, neighbor_thresh=1, readnoise=10)
        @test eltype(output) === T
        @test mask == input.mask
        @test all(output[input.mask] .!= input.image[input.mask])
    end
end