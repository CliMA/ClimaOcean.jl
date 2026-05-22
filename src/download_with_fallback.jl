using Downloads: Downloads
using NumericalEarth.DataWrangling: DataWrangling

const ARTIFACTS_BASE_URL = "https://github.com/NumericalEarth/NumericalEarthArtifacts/releases/download/data-v1/"

function download_from_artifacts(filepath::AbstractString)
    filename = basename(filepath)
    fallback_url = ARTIFACTS_BASE_URL * filename
    @info "Downloading $filename from NumericalEarthArtifacts fallback…"
    mktemp(dirname(filepath)) do tmppath, tmpio
        close(tmpio)
        Downloads.download(fallback_url, tmppath)
        mv(tmppath, filepath; force=true)
    end
end

download_from_artifacts(filepaths::AbstractVector) =
    foreach(download_from_artifacts, unique(filepaths))

"""
    download_with_fallback(metadata; dataset_name = string(metadata.name))

Download the data backing `metadata` from its primary source; if that fails,
fall back to the `NumericalEarthArtifacts` GitHub release mirror and retry.

Use to make scripts and docs examples robust to upstream data-server outages
(e.g. the ECCO JPL drive). After a successful call, `set!(field, metadata)`
finds the file locally and reads it normally.
"""
function download_with_fallback(metadata; dataset_name = string(metadata.name))
    filepaths = DataWrangling.metadata_path(metadata)
    try
        return DataWrangling.download_dataset(metadata)
    catch e
        @warn "Primary download failed for $dataset_name; trying NumericalEarthArtifacts fallback…" exception=(e, catch_backtrace())
        download_from_artifacts(filepaths)
        return DataWrangling.download_dataset(metadata)
    end
end
