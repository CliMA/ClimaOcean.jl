#####
##### Reading the output of an OMIP simulation
#####

"""
    OMIPOutput(directory = "."; prefix, label = prefix,
               averaging_window = (0, Inf),
               start_date = DateTime(1958, 1, 1),
               forcing_years = (1958, 2017))

The output that [`add_omip_diagnostics!`](@ref) wrote to `directory` with `prefix`, shown as `label` in the
figures. Time means and climatologies average the records of the `averaging_window` [s]. The simulation started
at `start_date` and cycles through the `forcing_years`, which places its records in the calendar of the
observations.
"""
struct OMIPOutput
    label :: String
    directory :: String
    prefix :: String
    averaging_window :: Tuple{Float64, Float64}
    start_date :: DateTime
    forcing_years :: Tuple{Int, Int}
    cache :: Dict{Symbol, Any}
end

OMIPOutput(directory = "."; prefix, label = prefix, averaging_window = (0, Inf),
           start_date = DateTime(1958, 1, 1), forcing_years = (1958, 2017)) =
    OMIPOutput(label, directory, prefix, Float64.(averaging_window), start_date, forcing_years, Dict{Symbol, Any}())

output_path(output::OMIPOutput, group) = joinpath(output.directory, output.prefix * "_" * group * ".jld2")

Base.summary(output::OMIPOutput) = string("OMIPOutput ", output.label, " from ", output_path(output, "*"))
Base.show(io::IO, output::OMIPOutput) = print(io, summary(output))
