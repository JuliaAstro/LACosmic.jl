# Examples

## Setup

You will need the following packages installed to replicate this tutorial

```julia
julia> ]add Distributions LACosmic Plots PSFModels
```

## Removing bad pixels with LACosmic.jl

First, let's create some fake data with Gaussian sources

```@example clean
using Distributions
using PSFModels: gaussian
using Random

function make_data(rng, N; N_sources=20, N_cosmics=20)
    imdata = fill(200.0, (N, N))

    # Add some fake sources
    for _ in 1:N_sources
        x = rand(rng, Uniform(1, N + 1))
        y = rand(rng, Uniform(1, N + 1))
        brightness = rand(rng, Uniform(1000, 30000)) / (2π * 3.5^2)
		model = gaussian(;x, y, fwhm=3.5, amp=brightness)
        imdata .+= model[axes(imdata)...]
	end

    # Add the poisson noise
    imdata .= rand.(rng, Poisson.(imdata))

    # Add readnoise
    imdata .+= rand(rng, Normal(0, 10), (N, N))

    clean_image = copy(imdata)
	
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
    return (image=imdata, clean_image, mask=crmask)
end

rng = MersenneTwister(808)
data = make_data(rng, 201)
```

let's inspect it

```@example clean
using Plots

function imshow(image; kwargs...)
	axy, axx = axes(image)
	heatmap(axy, axx, image; 
        aspect_ratio=1,
        ticks=false,
        xlim=extrema(axx),
        ylim=extrema(axy),
        kwargs...)
end

plot(
    imshow(log10.(data.clean_image), title="original image"),
    imshow(log10.(data.image), title="image w/cosmics"),
    size=(775, 350)
)
```

now we can clean it using [`lacosmic`](@ref)

```@example clean
using LACosmic

clean_image, mask = lacosmic(data.image, sigma_clip=6, contrast=5, neighbor_thresh=1)

plot(
    imshow(log10.(data.clean_image), title="original image"),
    imshow(log10.(clean_image), title="cleaned image"),
    size=(775, 350)
)
```

```@example clean
plot(
    imshow(data.mask, title="true cosmics", cbar=false),
    imshow(mask, title="detected cosmics", cbar=false),
    size=(700, 400)
)
```

```@example clean
data.mask == mask
```
