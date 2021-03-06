# LACosmic.jl

[![Build Status](https://github.com/JuliaAstro/LACosmic.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaAstro/LACosmic.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/L/LACosmic.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/JuliaAstro/LACosmic.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/LACosmic.jl)
[![License](https://img.shields.io/badge/License-BSD-yellow.svg)](https://opensource.org/licenses/BSD-3-Clause)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/LACosmic.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/LACosmic.jl/dev)

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

there is one exported function: `lacosmic`

```julia
clean_image, mask = lacosmic(image)
```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/LACosmic.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/LACosmic.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/JuliaAstro/LACosmic.jl/issues).
