# static kernels does laplace filter faster than imfiltering
const LAPLACE_KERNEL = @kernel w -> -w[0,-1] - w[-1,0] + 4 * w[0,0] - w[1,0] - w[0,1]
const MEDFILT3_KERNEL = Kernel{(-1:1,-1:1)}(@inline w -> median(Tuple(w)))
const MEDFILT5_KERNEL = Kernel{(-2:2,-2:2)}(@inline w -> median(Tuple(w)))
const MEDFILT7_KERNEL = Kernel{(-3:3,-3:3)}(@inline w -> median(Tuple(w)))

"""
    lacosmic(data::AbstractMatrix; 
        noise=nothing,
        gain=1,
        background=0,
        readnoise=0,
        mask=falses(size(data)),
        sigma_clip=4.5,
        contrast=5,
        neighbor_thresh=0.3,
        maxiter=4,
        saturation_level=2^16,
        block_size=2)

Laplacian cosmic ray detection (LACosmic). This algorithm implements the algorithm presented in [lacosmicx](https://github.com/cmccully/lacosmicx). The return values are the cleaned image and the bad pixel mask. The image cleaning is done via median interpolation.

# Parameters
- `noise` is the pre-determined estimate of the data noise (square root of variance), if any
- `gain` is the image gain in electrons per data number
- `background` is pre-determined image background, if any
- `readnoise` is the read noise of the image in electrons
- `mask` is an input bad pixel mask, where `true` represents a bad pixel
- `sigma_clip` is the Laplacian signal-to-noise ratio for flagging bad pixels
- `contrast` is the minimum contrast required to flag a bad pixel in the ratio of the Laplacian image to the fine-structure image
- `neighbor_thresh` is the fractional detection limit for cosmic rays surrounding other cosmic rays. Should be a number between 0 and 1.
- `maxiter` is the maximum number of iterations used for detecting bad pixels
- `saturation_level` is the saturation value in electrons
- `block_size` is the subsampling factor for the Laplacian filter image.

# Examples
```jldoctest
julia> image = 100 .* randn(1001, 1001) .+ 1000;

julia> clean_image, mask = lacosmic(image, gain=4);
```

# References
> [van Dokkum, P.G. (2001)](https://ui.adsabs.harvard.edu/abs/2001PASP..113.1420V/abstract) - "Cosmic-Ray Rejection by Laplacian Edge Detection"
"""
function lacosmic(
        data::AbstractMatrix{T};
        noise = nothing,
        gain = one(T),
        background = zero(T),
        readnoise = zero(T),
        mask = falses(size(data)),
        sigma_clip = 4.5,
        contrast = 5,
        neighbor_thresh = 0.3,
        maxiter = 4,
        saturation_level = T(2^16),
        block_size = 2) where T
    # convert image to electrons using gain and background if previously subtracted
    if gain != one(T)
        gain_corrected = data .* gain
    else
        gain_corrected = data
    end
    clean_image = @. gain_corrected - background * gain
    if !isnothing(noise) # noise provided, convert to electrons
        clean_var = @. (noise * gain)^2 - background * gain
    else # noise inferred
        clean_var = @. clean_image + readnoise^2
    end
    

    # update mask with saturated stars
    background_img = map(MEDFILT5_KERNEL, extend(gain_corrected, ExtensionReplicate()))
    sat_mask = @. (gain_corrected ≥ saturation_level) & (background_img > (saturation_level * 0.1))
    mask = dilate!(copy(mask), 3) .| dilate!(dilate!(sat_mask, 5), 5)

    good_data = @view clean_image[.!mask]

    # get the default background level for large cosmic rays
    background_level = median(good_data)

    # define cosmic ray mask
    ray_mask = falses(size(data))


    # loop for `maxiter` or if no more cosmic rays are found
    for iteration in 1:maxiter
        # subsample, convolve, clip, and rebin
        sub_img = subsample(clean_image, block_size)
        conv_img = map(LAPLACE_KERNEL, extend(sub_img, ExtensionConstant(zero(T))))
        snr = rebin(conv_img, block_size)

        # build Laplacian S/N map
        var_img = map(MEDFILT5_KERNEL, extend(clean_var, ExtensionReplicate()))
        @. snr /= block_size * sqrt(var_img + background * gain)
        # remove large structions
        snr_prime = snr .- map(MEDFILT5_KERNEL, extend(snr, ExtensionReplicate()))

        # build fine structure image
        med3 = map(MEDFILT3_KERNEL, extend(clean_image, ExtensionReplicate()))
        med7 = map(MEDFILT7_KERNEL, extend(med3, ExtensionReplicate()))
        # clip fine structure image since it is a divisor (similar to IRAF)
        f = @. max((med3 - med7) / sqrt(var_img + background * gain), T(0.01))
        
        # find candidate cosmic rays
        cosmics = @. !mask & (snr_prime > sigma_clip) & ((snr_prime / f) > contrast)

        # determine neighborhood
        cosmics = dilate!(cosmics, 3)
        @. cosmics &= !mask & (snr_prime > sigma_clip)
        cosmics = dilate!(cosmics, 3)
        @. cosmics &= !mask & (snr_prime > sigma_clip * neighbor_thresh)

        # update cosmic ray mask
        ray_mask .|= cosmics
        num_cosmics = count(cosmics)
        @debug "iteration $iteration: found $num_cosmics bad pixels this iteration"
        
        # if no more cosmics, break
        if iszero(num_cosmics)
            break
        # otherwise, clean the image
        else
            clean_image = clean(clean_image, ray_mask, mask, 5, background_level)
            clean_var = clean(clean_var, ray_mask, mask, 5, background_level)
        end
    end
    out_image = @. clean_image / gain + background
    return out_image, ray_mask
end

_padded_inds(idx, ax, half_width) = max(first(ax), (idx - half_width)):min(last(ax), (idx + half_width))

function dilate!(mask, size)
    half_width = Int((size - 1) / 2)
    @inbounds for idx in findall(mask)
        rows = _padded_inds(idx.I[1], axes(mask, 1), half_width)
        cols = _padded_inds(idx.I[2], axes(mask, 1), half_width)
        mask[rows, cols] .= true
    end
    return mask
end

function clean(reduction, image, crmask, mask, size, background_level)
    half_width = Int((size - 1) / 2)
    out = copy(image)
    @inbounds for idx in findall(crmask)
        rows = _padded_inds(idx.I[1], axes(mask, 1), half_width)
        cols = _padded_inds(idx.I[2], axes(mask, 1), half_width)
        window = @view image[rows, cols]
        m = @views @. !(crmask[rows, cols] | mask[rows, cols])
        masked_stat = reduction(view(window, m))
        out[idx] = max(masked_stat, background_level)
    end
    return out
end

clean(image, crmask, mask, size, background_level) = clean(mean, image, crmask, mask, size, background_level)