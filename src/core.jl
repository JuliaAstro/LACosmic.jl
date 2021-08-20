using ImageFiltering
using Statistics
using StaticArrays

function lacosmic(
        data::AbstractMatrix;
        mask=falses(size(data)),
        contrast=2.0,
        sigma_clip=4.5,
        objlim=2,
        maxiter=4,
        gain=1,
        background=0,
        saturation_level=2^16,
        readnoise=0,
        block_size=2)
    # convert image to electrons using gain and background if previously subtracted
    clean_image = muladd.(data, gain, background)

    # update mask with saturated stars
    mask = @. mask | (clean_image > saturation_level)
    good_data = @view clean_image[.!mask]

    # get the default background level for large cosmic rays
    background_level = median(good_data)

    # define cosmic ray mask
    ray_mask = falses(size(data))


    # loop for `maxiter` or if no more cosmic rays are found
    for iteration in 1:maxiter
        # subsample, convolve, clip, and rebin
        sub_img = subsample(clean_image, block_size)
        conv_img = imfilter(sub_img, Kernel.Laplacian((false, false)), "symmetric")
        Lplus = rebin(conv_img)

        # build Laplacian S/N map
        medfilt = mapwindow(median, clean_image, (5, 5))
        snr =  @. Lplus / (block_size * sqrt(medfilt + readnoise^2))
        # remove large structions
        snr_prime = snr .- mapwindow(median, snr, (5, 5))

        # build fine structure image
        f = mapwindow(median, clean_image, (3, 3))
        back = mapwindow(median, clean_image, (7, 7))
        f = @. (f - medfilt) / (block_size * sqrt(medfilt + readnoise^2))
        
        # find candidate cosmic rays
        cosmics = @. !mask & (snr_prime > sigma_clip) & ((snr_prime / f) > objlim)

        # determine neighborhood
        cosmics = mapwindow(any, cosmics, (3, 3))
        @. cosmics &= !mask & (snr_prime > sigma_clip)
        cosmics = mapwindow(any, cosmics, (3, 3))
        @. cosmics &= !mask & (snr_prime > (sigma_clip * 0.3))
        ray_mask .|= cosmics

        num_cosmics = count(cosmics)
        @debug "iteration $iteration: found $num_cosmics bad pixels this iteration"
        
        # if no more cosmics, break
        if iszero(num_cosmics)
            break
        # otherwise, clean the image
        else
            clean_image = mapwindow(w -> max(mean(w), background_level), clean_image, (5, 5), indices=findall(ray_mask))
        end
    end
    return clean_image, ray_mask
end

function rebin(arr)
    dims = Int.(size(arr) ./ 2)
    out = similar(arr, dims)
    @inbounds for idx in CartesianIndices(out)
        val = max(0, arr[idx.I[1] * 2 - 1, idx.I[2] * 2 - 1]) +
              max(0, arr[idx.I[1] * 2 - 1, idx.I[2] * 2]) +
              max(0, arr[idx.I[1] * 2, idx.I[2] * 2 - 1]) +
              max(0, arr[idx.I[1] * 2, idx.I[2] * 2])
        out[idx] = val
    end
    return out
end