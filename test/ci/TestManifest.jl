module TestManifest

using TOML

export Manifest
export all_test_files
export files_for_group
export files_for_shard
export load_manifest
export owner_shard
export shard_ids
export validate_manifest

const VALID_GROUPS = Set(["public", "internal", "expert"])

struct TestEntry
    path::String
    group::String
    shard::String
end

struct Manifest
    shard_order::Vector{String}
    tests::Vector{TestEntry}
    documentation_smoke::String
    full_run_paths::Set{String}
    full_run_prefixes::Vector{String}
    documentation_paths::Set{String}
    documentation_prefixes::Vector{String}
    source_impacts::Dict{String,Vector{String}}
    fixture_impacts::Dict{String,Vector{String}}
end

function string_vector(value, label::AbstractString)
    value isa Vector || throw(ArgumentError("$label must be an array"))
    all(item -> item isa AbstractString, value) ||
        throw(ArgumentError("$label must contain only strings"))
    return String[String(item) for item in value]
end

function string_map(value, label::AbstractString)
    value isa AbstractDict || throw(ArgumentError("$label must be a table"))
    result = Dict{String,Vector{String}}()
    for (key, entries) in value
        result[String(key)] = string_vector(entries, "$label.$key")
    end
    return result
end

function load_manifest(path::AbstractString)
    raw = TOML.parsefile(path)
    get(raw, "version", nothing) == 1 ||
        throw(ArgumentError("unsupported shard manifest version"))

    shard_order = string_vector(raw["shard_order"], "shard_order")
    tests = TestEntry[]
    for item in raw["tests"]
        push!(tests, TestEntry(
            String(item["path"]),
            String(item["group"]),
            String(item["shard"]),
        ))
    end

    return Manifest(
        shard_order,
        tests,
        String(raw["documentation_smoke"]),
        Set(string_vector(raw["full_run_paths"], "full_run_paths")),
        string_vector(raw["full_run_prefixes"], "full_run_prefixes"),
        Set(string_vector(raw["documentation_paths"], "documentation_paths")),
        string_vector(raw["documentation_prefixes"], "documentation_prefixes"),
        string_map(raw["source_impacts"], "source_impacts"),
        string_map(raw["fixture_impacts"], "fixture_impacts"),
    )
end

shard_ids(manifest::Manifest) = copy(manifest.shard_order)
all_test_files(manifest::Manifest) = [entry.path for entry in manifest.tests]
files_for_group(manifest::Manifest, group::AbstractString) =
    [entry.path for entry in manifest.tests if entry.group == group]
files_for_shard(manifest::Manifest, shard::AbstractString) =
    [entry.path for entry in manifest.tests if entry.shard == shard]

function owner_shard(manifest::Manifest, test_path::AbstractString)
    normalized = startswith(test_path, "test/") ? test_path[6:end] : String(test_path)
    index = findfirst(entry -> entry.path == normalized, manifest.tests)
    return isnothing(index) ? nothing : manifest.tests[index].shard
end

function validate_manifest(manifest::Manifest, test_root::AbstractString)
    length(unique(manifest.shard_order)) == length(manifest.shard_order) ||
        throw(ArgumentError("shard ids must be unique"))
    paths = all_test_files(manifest)
    length(unique(paths)) == length(paths) ||
        throw(ArgumentError("test paths must be unique"))

    shard_set = Set(manifest.shard_order)
    for entry in manifest.tests
        entry.group in VALID_GROUPS ||
            throw(ArgumentError("invalid test group: $(entry.group)"))
        entry.shard in shard_set ||
            throw(ArgumentError("unknown shard $(entry.shard) for $(entry.path)"))
        isfile(joinpath(test_root, entry.path)) ||
            throw(ArgumentError("missing test file: $(entry.path)"))
    end

    for shard in manifest.shard_order
        isempty(files_for_shard(manifest, shard)) &&
            throw(ArgumentError("empty shard: $shard"))
    end

    manifest.documentation_smoke in paths ||
        throw(ArgumentError("documentation smoke test must belong to the complete suite"))

    for impacts in values(manifest.source_impacts)
        all(shard -> shard in shard_set, impacts) ||
            throw(ArgumentError("source impact references an unknown shard"))
    end
    for impacts in values(manifest.fixture_impacts)
        all(shard -> shard in shard_set, impacts) ||
            throw(ArgumentError("fixture impact references an unknown shard"))
    end
    return nothing
end

end
