include(joinpath(@__DIR__, "TestManifest.jl"))
using .TestManifest

include(joinpath(@__DIR__, "TestSelection.jl"))
using .TestSelection

const MANIFEST_PATH = joinpath(@__DIR__, "shards.toml")
const TEST_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function parse_name_status(output::AbstractString)
    fields = String.(split(output, '\0'; keepempty = false))
    paths = String[]
    index = 1
    while index <= length(fields)
        status = fields[index]
        index += 1
        path_count = first(status) in ('R', 'C') ? 2 : 1
        index + path_count - 1 <= length(fields) ||
            throw(ArgumentError("malformed git name-status output"))
        append!(paths, fields[index:index + path_count - 1])
        index += path_count
    end
    return unique(paths)
end

function changed_paths(
    base::AbstractString,
    head::AbstractString;
    repository_root::AbstractString = REPOSITORY_ROOT,
)
    range = "$base...$head"
    output = read(
        `git -C $repository_root diff --name-status -z --find-renames --diff-filter=ACMRTD $range`,
        String,
    )
    return parse_name_status(output)
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

function github_delimiter(reasons::Vector{String})
    payload = replace(join(reasons, "\n"), "\r\n" => "\n", "\r" => "\n")
    payload_lines = Set(split(payload, '\n'; keepempty = true))
    delimiter = "EOF"
    suffix = 0
    while delimiter in payload_lines
        suffix += 1
        delimiter = "EOF_$suffix"
    end
    return delimiter
end

function write_github_output(path::AbstractString, selection::Selection)
    delimiter = github_delimiter(selection.reasons)
    open(path, "a") do io
        println(io, "matrix=$(matrix_json(selection.targets))")
        println(io, "documentation_only=$(selection.documentation_only)")
        println(io, "reason<<$delimiter")
        foreach(reason -> println(io, reason), selection.reasons)
        println(io, delimiter)
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
