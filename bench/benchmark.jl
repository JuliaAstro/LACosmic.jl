using BenchmarkTools
using CSV
using Distributions
using LACosmic
using ProgressLogging
using PSFModels: gaussian
using PyCall
using Random
using Statistics

ENV["OMP_NUM_THREADS"] = 1
@info "OMP NUM THREADS=$(ENV["OMP_NUM_THREADS"])"

scrap = pyimport("astroscrappy")


# generate data
rng = Random.MersenneTwister(444219)

function make_data(rng, T, N; N_sources=N ÷ 10, N_cosmics=N ÷ 10)
    imdata = fill(T(200), (N, N))

    # Add some fake sources
    for _ in 1:N_sources
        x = rand(rng, Uniform(1, N + 1))
        y = rand(rng, Uniform(1, N + 1))
        brightness = rand(rng, Uniform(1000, 30000))
		model = gaussian(T; x=x, y=y, fwhm=3.5, amp=brightness)
        imdata .+= model[axes(imdata)...]
	end

    # Add the poisson noise
    imdata .= rand.(rng, Poisson.(imdata))

    # Add readnoise
    imdata .+= rand(rng, Normal(0, 10), (N, N))

	clean_image = T.(imdata)
	
    # Add Nc fake cosmic rays
    crmask = falses((N, N))
	for _ in 1:N_cosmics
    	cr_x = round(Int, rand(rng, Uniform(6, N - 5)))
    	cr_y = round(Int, rand(rng, Uniform(6, N - 5)))
    	cr_brightnesses = rand(rng, Uniform(1000, 30000))
    	imdata[cr_y, cr_x] += T(cr_brightnesses)
    	crmask[cr_y, cr_x] = true
	end

    # Make a mask where the detected cosmic rays should be
    return (image=imdata, clean_image, mask=crmask)
end

data_sizes = 101:100:1001

results = []
@progress for N_data in data_sizes
    data = make_data(rng, Float64, N_data)
    
    bench_jl = @benchmark lacosmic($(data.image);
        sigma_clip=6,
        contrast=5,
        neighbor_thresh=1)
    bench_py = @benchmark scrap.detect_cosmics($(data.image);
        sigclip=6,
        objlim=5,
        sigfrac=1,
        readnoise=0,
        sepmed=false)
    
    res = (N_data, t_python=median(bench_py).time/1e9, t_julia=median(bench_jl).time/1e9)
    @info "" res...
    push!(results, res)
end

CSV.write(joinpath(@__DIR__, "benchmark_results.csv"), results)