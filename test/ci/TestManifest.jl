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

function unique_string_vector(value, label::AbstractString)
    result = string_vector(value, label)
    length(unique(result)) == length(result) ||
        throw(ArgumentError("$label values must be unique"))
    return result
end

function required_string(table::AbstractDict, key::AbstractString, label::AbstractString)
    haskey(table, key) || throw(ArgumentError("$label is required"))
    value = table[key]
    value isa AbstractString || throw(ArgumentError("$label must be a string"))
    return String(value)
end

function string_map(value, label::AbstractString)
    value isa AbstractDict || throw(ArgumentError("$label must be a table"))
    result = Dict{String,Vector{String}}()
    for (key, entries) in value
        result[String(key)] = string_vector(entries, "$label.$key")
    end
    return result
end

function validate_relative_path(path::AbstractString, label::AbstractString)
    isempty(path) && throw(ArgumentError("$label must not be empty"))
    occursin('\\', path) && throw(ArgumentError("$label must use forward slashes"))
    components = split(path, '/'; keepempty = true)
    any(isempty, components) &&
        throw(ArgumentError("$label must not contain empty components"))
    any(component -> component == "." || component == "..", components) &&
        throw(ArgumentError("$label must not contain '.' or '..' components"))

    first_component = first(components)
    if ncodeunits(first_component) >= 2
        bytes = codeunits(first_component)
        drive_letter =
            (0x41 <= bytes[1] && bytes[1] <= 0x5a) ||
            (0x61 <= bytes[1] && bytes[1] <= 0x7a)
        drive_letter && bytes[2] == 0x3a &&
            throw(ArgumentError("$label must not use a Windows drive prefix"))
    end
    return nothing
end

repository_path(root::AbstractString, path::AbstractString) =
    joinpath(root, split(path, '/')...)

function validate_policy_values(values, label::AbstractString)
    isempty(values) && throw(ArgumentError("$label must not be empty"))
    all(value -> !isempty(value), values) ||
        throw(ArgumentError("$label must not contain empty values"))
    length(unique(values)) == length(values) ||
        throw(ArgumentError("$label values must be unique"))
    return nothing
end

function validate_policy_files(values, label::AbstractString, repository_root::AbstractString)
    validate_policy_values(values, label)
    for path in values
        validate_relative_path(path, "$label entry")
        isfile(repository_path(repository_root, path)) ||
            throw(ArgumentError("$label references a missing file: $path"))
    end
    return nothing
end

function validate_policy_prefixes(
    prefixes,
    label::AbstractString,
    repository_root::AbstractString,
)
    validate_policy_values(prefixes, label)
    for prefix in prefixes
        endswith(prefix, "/") ||
            throw(ArgumentError("$label entry must end with '/': $prefix"))
        directory = chop(prefix; tail = 1)
        validate_relative_path(directory, "$label entry")
        isdir(repository_path(repository_root, directory)) ||
            throw(ArgumentError("$label references a missing directory: $prefix"))
    end
    return nothing
end

function validate_impact_map(
    impacts::Dict{String,Vector{String}},
    label::AbstractString,
    required_prefix::AbstractString,
    repository_root::AbstractString,
    shard_set::Set{String},
)
    isempty(impacts) && throw(ArgumentError("$label must not be empty"))
    for (path, shards) in impacts
        validate_relative_path(path, "$label key")
        startswith(path, required_prefix) ||
            throw(ArgumentError("$label key must start with $required_prefix: $path"))
        isfile(repository_path(repository_root, path)) ||
            throw(ArgumentError("$label key references a missing file: $path"))
        validate_policy_values(shards, "$label.$path")
        all(shard -> shard in shard_set, shards) ||
            throw(ArgumentError("$label references an unknown shard"))
    end
    return nothing
end

function load_manifest(path::AbstractString)
    raw = TOML.parsefile(path)
    get(raw, "version", nothing) == 1 ||
        throw(ArgumentError("unsupported shard manifest version"))

    shard_order = string_vector(raw["shard_order"], "shard_order")
    raw_tests = get(raw, "tests", nothing)
    raw_tests isa Vector || throw(ArgumentError("tests must be an array"))
    tests = TestEntry[]
    for (index, item) in enumerate(raw_tests)
        item isa AbstractDict || throw(ArgumentError("tests[$index] must be a table"))
        push!(tests, TestEntry(
            required_string(item, "path", "tests[$index].path"),
            required_string(item, "group", "tests[$index].group"),
            required_string(item, "shard", "tests[$index].shard"),
        ))
    end

    return Manifest(
        shard_order,
        tests,
        String(raw["documentation_smoke"]),
        Set(unique_string_vector(raw["full_run_paths"], "full_run_paths")),
        string_vector(raw["full_run_prefixes"], "full_run_prefixes"),
        Set(unique_string_vector(raw["documentation_paths"], "documentation_paths")),
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
    isempty(manifest.shard_order) && throw(ArgumentError("shard ids must not be empty"))
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
        validate_relative_path(entry.path, "test path")
        startswith(entry.path, "$(entry.group)/") ||
            throw(ArgumentError("test path does not match group: $(entry.path)"))
        endswith(entry.path, ".jl") ||
            throw(ArgumentError("test path must name a Julia file: $(entry.path)"))
        if entry.group == "public"
            entry.shard == "public" ||
                throw(ArgumentError("public test must belong to the public shard"))
        else
            startswith(entry.shard, "$(entry.group)-") ||
                throw(ArgumentError("test group does not match shard: $(entry.path)"))
        end
        isfile(repository_path(test_root, entry.path)) ||
            throw(ArgumentError("missing test file: $(entry.path)"))
    end

    for shard in manifest.shard_order
        isempty(files_for_shard(manifest, shard)) &&
            throw(ArgumentError("empty shard: $shard"))
    end

    manifest.documentation_smoke in paths ||
        throw(ArgumentError("documentation smoke test must belong to the complete suite"))

    repository_root = normpath(joinpath(test_root, ".."))
    validate_policy_files(manifest.full_run_paths, "full_run_paths", repository_root)
    validate_policy_prefixes(
        manifest.full_run_prefixes,
        "full_run_prefixes",
        repository_root,
    )
    validate_policy_files(
        manifest.documentation_paths,
        "documentation_paths",
        repository_root,
    )
    validate_policy_prefixes(
        manifest.documentation_prefixes,
        "documentation_prefixes",
        repository_root,
    )
    validate_impact_map(
        manifest.source_impacts,
        "source_impacts",
        "src/",
        repository_root,
        shard_set,
    )
    validate_impact_map(
        manifest.fixture_impacts,
        "fixture_impacts",
        "test/fixtures/",
        repository_root,
        shard_set,
    )

    fixture_directory = joinpath(test_root, "fixtures")
    expected_fixtures = Set(
        "test/fixtures/$file" for file in readdir(fixture_directory) if
        endswith(file, ".jl") && isfile(joinpath(fixture_directory, file))
    )
    Set(keys(manifest.fixture_impacts)) == expected_fixtures ||
        throw(ArgumentError("fixture_impacts must cover every test fixture exactly"))
    return nothing
end

end
