include(joinpath(@__DIR__, "TestManifest.jl"))
using .TestManifest

include(joinpath(@__DIR__, "TestSelection.jl"))
using .TestSelection

const MANIFEST_PATH = joinpath(@__DIR__, "shards.toml")
const TEST_ROOT = normpath(joinpath(@__DIR__, ".."))

function changed_paths(base::AbstractString, head::AbstractString)
    range = "$base...$head"
    output = read(`git diff --name-only --diff-filter=ACMRT $range`, String)
    return filter(!isempty, split(chomp(output), '\n'))
end

function argument_value(argument::AbstractString, option::AbstractString)
    value = argument[length(option) + 2:end]
    isempty(value) && throw(ArgumentError("$option must not be empty"))
    return value
end

function parse_arguments(arguments::Vector{String})
    base = nothing
    head = "HEAD"
    format = nothing
    github_output = nothing

    for argument in arguments
        if startswith(argument, "--base=")
            base = argument_value(argument, "--base")
        elseif startswith(argument, "--head=")
            head = argument_value(argument, "--head")
        elseif startswith(argument, "--format=")
            format = argument_value(argument, "--format")
        elseif startswith(argument, "--github-output=")
            github_output = argument_value(argument, "--github-output")
        else
            throw(ArgumentError("unknown option: $argument"))
        end
    end

    isnothing(base) && throw(ArgumentError("--base is required"))
    !isnothing(format) && format != "lines" &&
        throw(ArgumentError("unsupported format: $format"))
    return (; base, head, format, github_output)
end

function write_github_output(path::AbstractString, selection::Selection)
    open(path, "a") do io
        println(io, "matrix=$(matrix_json(selection.targets))")
        println(io, "documentation_only=$(selection.documentation_only)")
        println(io, "reason<<EOF")
        foreach(reason -> println(io, reason), selection.reasons)
        println(io, "EOF")
    end
    return nothing
end

function main(arguments::Vector{String})
    options = parse_arguments(arguments)
    manifest = load_manifest(MANIFEST_PATH)
    validate_manifest(manifest, TEST_ROOT)
    paths = String.(changed_paths(options.base, options.head))
    selection = select_targets(paths, manifest)

    if options.format == "lines"
        foreach(println, selection.targets)
    end
    if !isnothing(options.github_output)
        write_github_output(options.github_output, selection)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main(ARGS)
    catch err
        showerror(stderr, err, catch_backtrace())
        println(stderr)
        exit(1)
    end
end
