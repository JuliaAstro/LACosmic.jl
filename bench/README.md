# benchmarks

Benchmarks against [astroscrappy](https://github.com/astropy/astroscrappy)

## setup

1. install python dependencies using poetry
   ```sh
   $ poetry install
   $ poetry shell
   ```
2. instantiate julia environment and develop local project
   ```sh
   $ julia --project=bench -e 'using Pkg; Pkg.instantiate(); Pkg.develop(path=pwd())'
   ```
   or, from within the Julia REPL
   ```julia
   pkg> activate bench
   pkg> instantiate
   pkg> dev .
   ```
3. ensure PyCall.jl is linked to appropriate python
   ```sh
   $ PYTHON=$(which python) julia --project=bench -e 'using Pkg; Pkg.build("PyCall")`
   ```
   or, from within the Julia REPL
   ```julia
   shell> which python
   julia> ENV["PYTHON"] = "/path/to/python"
   pkg> build PyCall
   ```

## usage

To run the benchmark, first make sure you have appropriately linked the version of python you are using with PyCall, as shown in the setup above. Once PyCall is linked, simply activating the Julia environment and executing `benchmark.jl` will run the benchmark. For fair comparisons, make sure the number of threads used is the same for each language-

```sh
$ export OMP_NUM_THREADS=1
$ export JULIA_NUM_THREADS=1
$ julia --project=bench bench/benchmark.jl
```
or, from within the Julia REPL
```julia
pkg> activate bench
julia> include("bench/benchmark.jl")
```
