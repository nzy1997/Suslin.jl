module TestRunner

using ..TestManifest

export requested_targets

function timed_test_file(operation, path::AbstractString; io::IO = stdout)
    started = time_ns()
    try
        return operation()
    finally
        elapsed = (time_ns() - started) / 1_000_000_000
        println(io, "TEST_FILE_TIME\t", path, "\t", round(elapsed; digits = 3))
        flush(io)
    end
end

function requested_targets(args::Vector{String}, manifest::Manifest)
    names = isempty(args) ? ["public", "internal"] :
        reduce(vcat, [filter(!isempty, split(arg, ',')) for arg in args])
    expanded_names = String[]
    for name in names
        if name == "all"
            append!(expanded_names, ["public", "internal", "expert"])
        else
            push!(expanded_names, name)
        end
    end
    isempty(expanded_names) && throw(ArgumentError("no test targets requested"))

    targets = Pair{String,Vector{String}}[]
    seen = Set{String}()
    target_for_file = Dict{String,String}()
    for name in expanded_names
        name in seen && continue
        push!(seen, name)
        files = if name in ("public", "internal", "expert")
            files_for_group(manifest, name)
        elseif name == "documentation-smoke"
            [manifest.documentation_smoke]
        elseif startswith(name, "shard:")
            shard = name[7:end]
            shard in shard_ids(manifest) ||
                throw(ArgumentError("unknown test shard: $shard"))
            files_for_shard(manifest, shard)
        else
            throw(ArgumentError("unknown test target: $name"))
        end
        isempty(files) && throw(ArgumentError("test target has no files: $name"))
        for file in files
            if haskey(target_for_file, file)
                previous = target_for_file[file]
                throw(ArgumentError(
                    "test file $file is selected by both $previous and $name",
                ))
            end
        end
        push!(targets, name => files)
        for file in files
            target_for_file[file] = name
        end
    end
    return targets
end

end
