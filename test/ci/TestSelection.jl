module TestSelection

using ..TestManifest

export Selection
export matrix_json
export select_targets

struct Selection
    targets::Vector{String}
    documentation_only::Bool
    reasons::Vector{String}
end

is_prefixed(path::AbstractString, prefixes) =
    any(prefix -> startswith(path, prefix), prefixes)

function ordered_targets(selected::Set{String}, manifest::Manifest)
    return [shard for shard in shard_ids(manifest) if shard in selected]
end

function select_targets(changed_paths::Vector{String}, manifest::Manifest)
    isempty(changed_paths) && return Selection(
        shard_ids(manifest), false, ["empty diff: full fallback"])

    if all(path -> path in manifest.documentation_paths ||
                   is_prefixed(path, manifest.documentation_prefixes), changed_paths)
        return Selection(["documentation-smoke"], true, ["documentation-only diff"])
    end

    selected = Set{String}()
    reasons = String[]
    for path in changed_paths
        if path in manifest.full_run_paths || is_prefixed(path, manifest.full_run_prefixes)
            return Selection(shard_ids(manifest), false, ["full-run trigger: $path"])
        elseif startswith(path, "src/")
            impacts = get(manifest.source_impacts, path, nothing)
            isnothing(impacts) && return Selection(
                shard_ids(manifest), false, ["unknown source path: $path"])
            union!(selected, impacts)
            push!(reasons, "$path => $(join(impacts, ','))")
        elseif startswith(path, "test/fixtures/")
            impacts = get(manifest.fixture_impacts, path, nothing)
            isnothing(impacts) && return Selection(
                shard_ids(manifest), false, ["unknown fixture path: $path"])
            union!(selected, impacts)
            push!(reasons, "$path => $(join(impacts, ','))")
        elseif startswith(path, "test/")
            owner = owner_shard(manifest, path)
            isnothing(owner) && return Selection(
                shard_ids(manifest), false, ["unknown test path: $path"])
            push!(selected, owner)
            push!(reasons, "$path => $owner")
        elseif path in manifest.documentation_paths ||
               is_prefixed(path, manifest.documentation_prefixes)
            push!(reasons, "$path => documentation companion")
        else
            return Selection(shard_ids(manifest), false, ["unknown path: $path"])
        end
    end

    isempty(selected) && return Selection(
        shard_ids(manifest), false, ["empty source selection: full fallback"])
    return Selection(ordered_targets(selected, manifest), false, reasons)
end

function escape_json_string(value::AbstractString)
    io = IOBuffer()
    for character in value
        if character == '"'
            print(io, "\\\"")
        elseif character == '\\'
            print(io, "\\\\")
        elseif character == '\b'
            print(io, "\\b")
        elseif character == '\f'
            print(io, "\\f")
        elseif character == '\n'
            print(io, "\\n")
        elseif character == '\r'
            print(io, "\\r")
        elseif character == '\t'
            print(io, "\\t")
        elseif UInt32(character) <= 0x1f
            print(io, "\\u", lpad(string(UInt32(character); base = 16), 4, '0'))
        else
            print(io, character)
        end
    end
    return String(take!(io))
end

function matrix_json(targets::Vector{String})
    escaped = escape_json_string.(targets)
    return "[" * join(["\"$target\"" for target in escaped], ",") * "]"
end

end
