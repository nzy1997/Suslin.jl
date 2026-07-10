module TestRunner

using ..TestManifest

export requested_targets

function requested_targets(args::Vector{String}, manifest::Manifest)
    names = isempty(args) ? ["public", "internal"] :
        reduce(vcat, [filter(!isempty, split(arg, ',')) for arg in args])
    "all" in names && (names = ["public", "internal", "expert"])

    targets = Pair{String,Vector{String}}[]
    seen = Set{String}()
    for name in names
        name in seen && continue
        push!(seen, name)
        if name in ("public", "internal", "expert")
            push!(targets, name => files_for_group(manifest, name))
        elseif name == "documentation-smoke"
            push!(targets, name => [manifest.documentation_smoke])
        elseif startswith(name, "shard:")
            shard = name[7:end]
            shard in shard_ids(manifest) ||
                throw(ArgumentError("unknown test shard: $shard"))
            push!(targets, name => files_for_shard(manifest, shard))
        else
            throw(ArgumentError("unknown test target: $name"))
        end
    end
    return targets
end

end
