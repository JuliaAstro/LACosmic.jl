```@meta
CurrentModule = LACosmic
```

# LACosmic.jl

[![Code](https://img.shields.io/badge/Code-GitHub-black.svg))](https://github.com/JuliaAstro/LACosmic.jl)
[![Build Status](https://github.com/JuliaAstro/LACosmic.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaAstro/LACosmic.jl/actions)
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

there is one exported function: `lacosmic`

```julia
clean_image, mask = lacosmic(image)
```

## API/Reference

```@index
```

```@autodocs
Modules = [LACosmic]
```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/LACosmic.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/LACosmic.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/JuliaAstro/LACosmic.jl/issues).
