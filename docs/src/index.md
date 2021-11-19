```@meta
CurrentModule = LACosmic
```

# LACosmic.jl

[![Code](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/JuliaAstro/LACosmic.jl)
[![Build Status](https://github.com/JuliaAstro/LACosmic.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaAstro/LACosmic.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/L/LACosmic.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/JuliaAstro/LACosmic.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/LACosmic.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Laplacian cosmic-ray detection (L.A.Cosmic) in pure Julia.

## Installation

To use the LACosmic library, first install it using `Pkg`

```julia
julia> ]add LACosmic
```

## Usage

To import the library

```julia
using LACosmic
```

there is one exported function: [`lacosmic`](@ref)

```julia
clean_image, mask = lacosmic(image)
```

## Benchmarks

This code has been benchmarked against the Cython implementation in [Astro-SCRAPPY](https://github.com/astropy/astroscrappy). This benchmark simply computes the time it takes to run the LACosmic algorithm with different image sizes. The size is the length of one dimension of the image, so the expected scaling should be ``\propto N^2``. The code can be found in `bench/benchmark.jl`. Here is the information for my system-

```julia
Julia Version 1.6.0
Commit f9720dc2eb* (2021-03-24 12:55 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin20.3.0)
  CPU: Intel(R) Core(TM) i5-8259U CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, skylake)
Environment:
  OMP_NUM_THREADS = 1
  JULIA_NUM_THREADS = 1
```

```@example
using CSV, DataFrames, Plots # hide
using LACosmic # hide
benchdir = joinpath(dirname(pathof(LACosmic)), "..", "bench") # hide
results = DataFrame(CSV.File(joinpath(benchdir, "benchmark_results.csv"))) # hide

plot(results.N_data, [results.t_python results.t_julia]; # hide
    label=["Astro-SCRAPPY" "LACosmic.jl"], leg=:topleft, # hide
    xlabel="image size", ylabel="time (s)", # hide
    yscale=:log10, shape=:o) # hide
points = range(extrema(results.N_data)..., length=51) # hide
factor = results.t_julia[begin] / (results.N_data[begin])^2 # hide
plot!(points, t -> t^2 * factor, ls=:dash, c=:black, alpha=0.3, lab="N²") # hide

```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/LACosmic.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/LACosmic.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/JuliaAstro/LACosmic.jl/issues).
