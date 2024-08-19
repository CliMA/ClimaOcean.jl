using ClimaOcean
using Oceananigans
using Test

using ClimaOcean.DataWrangling
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.JRA55: JRA55_field_time_series

using Oceananigans.Architectures: architecture
using Oceananigans.OutputReaders: interpolate!

using ClimaOcean

test_architectures = [CPU()]
