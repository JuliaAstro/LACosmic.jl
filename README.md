# LACosmic.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.org/LACosmic/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.org/LACosmic.jl/dev/)

[![CI](https://github.com/JuliaAstro/LACosmic.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaAstro/LACosmic.jl/actions/workflows/CI.yml)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/L/LACosmic.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/JuliaAstro/LACosmic.jl/graph/badge.svg?token=2lJWwb9zoO)](https://codecov.io/gh/JuliaAstro/LACosmic.jl)
[![License](https://img.shields.io/badge/License-BSD-yellow.svg)](https://opensource.org/licenses/BSD-3-Clause)

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
