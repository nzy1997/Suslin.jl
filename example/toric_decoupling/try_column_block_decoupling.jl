import Pkg

const DEFAULT_CASES = ["case_001", "case_004"]
const DEFAULT_BLOCKS = ["column_Q", "pair_mix_2_1"]
const DEFAULT_TORICBUILDER_DIR = joinpath(
    homedir(),
    "jcode",
    "topological-code-decoupling",
    "julia_code",
    "ToricBuilder",
)
const SUSLIN_PKG_ID = Base.PkgId(Base.UUID("47b4a9d4-42d3-414b-bb10-8db3ed793992"), "Suslin")
const TORICBUILDER_PKG_ID = Base.PkgId(Base.UUID("603dc7be-2f45-49e8-9590-fa6d91b7b2aa"), "ToricBuilder")
const _SCRIPT_MODULE = @__MODULE__
const _CACHE_TRANSFER_FIELDS = (
    :mat_after_coarse_graining,
    :result_matrix,
    :row_transformation,
    :column_transformation,
)

default_toricbuilder_dir() = get(ENV, "TORICBUILDER_DIR", DEFAULT_TORICBUILDER_DIR)

default_cache_dir(toricbuilder_dir=default_toricbuilder_dir()) =
    get(ENV, "TORICBUILDER_CACHE_DIR", joinpath(toricbuilder_dir, "example", "decouple_bbcodes_cache"))

function local_toricbuilder_available(;
    toricbuilder_dir=default_toricbuilder_dir(),
    cache_dir=default_cache_dir(toricbuilder_dir),
)
    return isdir(toricbuilder_dir) &&
        isfile(joinpath(toricbuilder_dir, "Project.toml")) &&
        isfile(joinpath(cache_dir, "case_001.jls")) &&
        isfile(joinpath(cache_dir, "case_004.jls"))
end

function main(args=ARGS; io=stdout)
    options = _parse_args(args)
    return run_smoke(
        options.cases;
        toricbuilder_dir = options.toricbuilder_dir,
        cache_dir = options.cache_dir,
        io,
    )
end

function run_smoke(
    cases=DEFAULT_CASES;
    toricbuilder_dir=default_toricbuilder_dir(),
    cache_dir=default_cache_dir(toricbuilder_dir),
    io=stdout,
)
    normalized_cases = string.(collect(cases))
    modules = _load_optional_modules(toricbuilder_dir)
    if modules.status != "PASS"
        rows = [_cache_error_row(case_id; block = "cache") for case_id in normalized_cases]
        _emit_rows(io, rows)
        return rows
    end

    rows = NamedTuple[]
    for case_id in normalized_cases
        append!(rows, _run_case_smoke(case_id, modules, cache_dir))
    end

    _emit_rows(io, rows)
    return rows
end

function _parse_args(args)
    cases = DEFAULT_CASES
    toricbuilder_dir = default_toricbuilder_dir()
    cache_dir = nothing

    for arg in args
        if startswith(arg, "--case=")
            raw_cases = split(arg[length("--case=")+1:end], ",")
            cases = [strip(case_id) for case_id in raw_cases if !isempty(strip(case_id))]
        elseif startswith(arg, "--toricbuilder-dir=")
            toricbuilder_dir = arg[length("--toricbuilder-dir=")+1:end]
        elseif startswith(arg, "--cache-dir=")
            cache_dir = arg[length("--cache-dir=")+1:end]
        else
            throw(ArgumentError("unsupported argument: $arg"))
        end
    end

    isempty(cases) && throw(ArgumentError("--case must name at least one cache case"))
    return (;
        cases,
        toricbuilder_dir,
        cache_dir = isnothing(cache_dir) ? default_cache_dir(toricbuilder_dir) : cache_dir,
    )
end

function _load_optional_modules(toricbuilder_dir::AbstractString)
    if !isdir(toricbuilder_dir) || !isfile(joinpath(toricbuilder_dir, "Project.toml"))
        return (; status = "FAIL", failure = "CACHE_ERROR")
    end

    _prepend_load_path!(_suslin_repo_root(), toricbuilder_dir)

    try
        Pkg.activate(_suslin_repo_root(); io=devnull)
        suslin = Base.require(SUSLIN_PKG_ID)
        toricbuilder = Base.require(TORICBUILDER_PKG_ID)
        return (;
            status = "PASS",
            failure = "NONE",
            Oscar = getfield(toricbuilder, :Oscar),
            Suslin = suslin,
            ToricBuilder = toricbuilder,
        )
    catch err
        err isa InterruptException && rethrow()
        return (; status = "FAIL", failure = "CACHE_ERROR")
    end
end

function _prepend_load_path!(paths...)
    for path in reverse(abspath.(paths))
        filter!(entry -> abspath(entry) != path, LOAD_PATH)
        pushfirst!(LOAD_PATH, path)
    end
    return LOAD_PATH
end

_suslin_repo_root() = normpath(joinpath(@__DIR__, "..", ".."))

function _run_case_smoke(case_id::AbstractString, modules, cache_dir::AbstractString)
    cached_case = _load_cache_case(case_id, modules.ToricBuilder, cache_dir)
    if cached_case.status != "PASS"
        return [_cache_error_row(case_id; block = "cache")]
    end

    cache = cached_case.case
    if _cache_case_unusable(cache)
        return [_cache_error_row(case_id; block = "cache")]
    end

    rows = NamedTuple[]
    for block in DEFAULT_BLOCKS
        block_matrix = _extract_block(cache, block, modules)
        if block_matrix.status != "PASS"
            push!(rows, _cache_error_row(case_id; block))
        else
            push!(rows, _route_block(case_id, block, block_matrix.matrix, modules))
        end
    end
    return rows
end

function _load_cache_case(case_id::AbstractString, ToricBuilder, cache_dir::AbstractString)
    path = joinpath(cache_dir, string(case_id, ".jls"))
    isfile(path) || return (; status = "FAIL", failure = "CACHE_ERROR")

    try
        return (;
            status = "PASS",
            failure = "NONE",
            case = Base.invokelatest(getfield(ToricBuilder, :load_cached_toric_case), path),
        )
    catch err
        err isa InterruptException && rethrow()
        return (; status = "FAIL", failure = "CACHE_ERROR")
    end
end

function _cache_case_unusable(cached_case)
    if !hasproperty(cached_case, :status) || getproperty(cached_case, :status) != :ok
        return true
    end
    if !hasproperty(cached_case, :transfer_result) || isnothing(getproperty(cached_case, :transfer_result))
        return true
    end

    transfer = getproperty(cached_case, :transfer_result)
    for field in _CACHE_TRANSFER_FIELDS
        if !hasproperty(transfer, field) || isnothing(getproperty(transfer, field))
            return true
        end
    end
    return false
end

function _extract_block(cached_case, block::AbstractString, modules)
    try
        if block == "column_Q"
            return (; status = "PASS", failure = "NONE", matrix = _decoder_basis_change(cached_case, modules).Q)
        elseif block == "pair_mix_2_1"
            return (; status = "PASS", failure = "NONE", matrix = _toric_pair_mix_2_1(cached_case, modules))
        end
        return (; status = "FAIL", failure = "CACHE_ERROR")
    catch err
        err isa InterruptException && rethrow()
        return (; status = "FAIL", failure = "CACHE_ERROR")
    end
end

function _decoder_basis_change(cached_case, modules)
    transfer = cached_case.transfer_result
    hasproperty(transfer, :Q_size) || throw(ArgumentError("missing Q_size"))
    Qsize = transfer.Q_size[1] ÷ 2
    column_transformation = transfer.column_transformation
    Q = _matrix_slice(column_transformation, 1:Qsize, 1:Qsize, modules)
    laurent_conjugate = getfield(modules.ToricBuilder, :laurent_conjugate)
    lower = _matrix_slice(column_transformation, (Qsize + 1):(2 * Qsize), (Qsize + 1):(2 * Qsize), modules)
    R = _matrix_base_ring(column_transformation, modules)
    Qinv = _matrix_from_entries(
        R,
        [
            Base.invokelatest(laurent_conjugate, _matrix_entry(lower, col, row))
            for row in 1:Qsize, col in 1:Qsize
        ],
        modules,
    )
    return (; Qsize, Q, Qinv)
end

function _toric_pair_mix_2_1(cached_case, modules)
    transfer = cached_case.transfer_result
    hasproperty(transfer, :product_state_num) || throw(ArgumentError("missing product_state_num"))
    hasproperty(transfer, :toric_num) || throw(ArgumentError("missing toric_num"))
    transfer.toric_num >= 2 || throw(ArgumentError("pair_mix_2_1 requires at least two toric factors"))

    decoder = _decoder_basis_change(cached_case, modules)
    R = _matrix_base_ring(transfer.column_transformation, modules)
    e = Base.invokelatest(getfield(modules.Oscar, :identity_matrix), R, decoder.Qsize)
    i = 2
    j = 1
    target_odd = 2 * transfer.product_state_num + 2 * i - 1
    source_odd = 2 * transfer.product_state_num + 2 * j - 1
    target_even = 2 * transfer.product_state_num + 2 * i
    source_even = 2 * transfer.product_state_num + 2 * j
    for row in 1:decoder.Qsize
        _set_matrix_entry!(
            e,
            row,
            target_odd,
            _add_ring_entries(_matrix_entry(e, row, target_odd), _matrix_entry(e, row, source_odd)),
        )
        _set_matrix_entry!(
            e,
            row,
            target_even,
            _add_ring_entries(_matrix_entry(e, row, target_even), _matrix_entry(e, row, source_even)),
        )
    end
    return Base.invokelatest(*, Base.invokelatest(*, decoder.Q, e), decoder.Qinv)
end

function _matrix_base_ring(A, modules)
    return Base.invokelatest(getfield(modules.Oscar, :base_ring), A)
end

function _add_ring_entries(left, right)
    return Base.invokelatest(+, left, right)
end

function _matrix_entry(A, row::Int, col::Int)
    return Base.invokelatest(getindex, A, row, col)
end

function _set_matrix_entry!(A, row::Int, col::Int, value)
    return Base.invokelatest(setindex!, A, value, row, col)
end

function _matrix_from_entries(R, entries, modules)
    return Base.invokelatest(getfield(modules.Oscar, :matrix), R, entries)
end

function _matrix_slice(A, rows, cols, modules)
    R = _matrix_base_ring(A, modules)
    return _matrix_from_entries(
        R,
        [_matrix_entry(A, row, col) for row in rows, col in cols],
        modules,
    )
end

function _route_block(case_id::AbstractString, block::AbstractString, A, modules)
    try
        profile = Base.invokelatest(getfield(modules.Suslin, :classify_laurent_determinant), A)
        det_token = string(profile.classification)
        if profile.classification == :one
            return _route_sl_core(case_id, block, A, det_token, modules)
        elseif profile.classification == :laurent_monomial_unit
            return _route_gl_certificate(case_id, block, A, det_token, modules)
        end
        return _unsupported_staged_row(
            case_id,
            block,
            A,
            det_token;
            normalization = "NORMALIZATION_SKIP",
            sl_core = "SL_CORE_SKIP",
            gl_cert = "GL_CERT_SKIP",
        )
    catch err
        err isa InterruptException && rethrow()
        if _is_staged_argument_error(err)
            return _unsupported_staged_row(
                case_id,
                block,
                A,
                "unknown";
                normalization = "NORMALIZATION_SKIP",
                sl_core = "SL_CORE_SKIP",
                gl_cert = "GL_CERT_SKIP",
            )
        end
        return _route_error_row(case_id, block, A, "unknown")
    end
end

function _route_sl_core(case_id::AbstractString, block::AbstractString, A, det_token::AbstractString, modules)
    try
        certificate = Base.invokelatest(getfield(modules.Suslin, :_factor_laurent_sl_column_peel), A)
        verified = string(Bool(certificate.verification.overall_ok))
        return _row(
            case_id,
            block,
            A,
            det_token;
            normalization = "NORMALIZATION_SKIP",
            sl_core = verified == "true" ? "SL_CORE_PASS" : "SL_CORE_FAIL",
            gl_cert = "GL_CERT_SKIP",
            factors = string(length(certificate.factors)),
            verified,
            status = verified == "true" ? "PASS" : "FAIL",
            failure = verified == "true" ? "NONE" : "ROUTE_ERROR",
        )
    catch err
        err isa InterruptException && rethrow()
        if _is_staged_argument_error(err)
            return _unsupported_staged_row(
                case_id,
                block,
                A,
                det_token;
                normalization = "NORMALIZATION_SKIP",
                sl_core = "UNSUPPORTED_STAGED",
                gl_cert = "GL_CERT_SKIP",
            )
        end
        return _route_error_row(case_id, block, A, det_token)
    end
end

function _route_gl_certificate(case_id::AbstractString, block::AbstractString, A, det_token::AbstractString, modules)
    normalization_status = "NORMALIZATION_SKIP"
    try
        Base.invokelatest(getfield(modules.Suslin, :normalize_laurent_gl_matrix), A)
        normalization_status = "NORMALIZATION_PASS"
        certificate = Base.invokelatest(getfield(modules.Suslin, :laurent_gl_factorization_certificate), A)
        verified = string(Base.invokelatest(
            getfield(modules.Suslin, :verify_laurent_gl_factorization_certificate),
            certificate,
        ))
        return _row(
            case_id,
            block,
            A,
            det_token;
            normalization = normalization_status,
            sl_core = "SL_CORE_SKIP",
            gl_cert = verified == "true" ? "GL_CERT_PASS" : "GL_CERT_FAIL",
            factors = string(length(certificate.core_factors)),
            verified,
            status = verified == "true" ? "PASS" : "FAIL",
            failure = verified == "true" ? "NONE" : "ROUTE_ERROR",
        )
    catch err
        err isa InterruptException && rethrow()
        if _is_staged_argument_error(err)
            return _unsupported_staged_row(
                case_id,
                block,
                A,
                det_token;
                normalization = normalization_status,
                sl_core = "SL_CORE_SKIP",
                gl_cert = "UNSUPPORTED_STAGED",
            )
        end
        return _row(
            case_id,
            block,
            A,
            det_token;
            normalization = normalization_status == "NORMALIZATION_PASS" ?
                normalization_status :
                "NORMALIZATION_FAIL",
            sl_core = "SL_CORE_SKIP",
            gl_cert = "GL_CERT_FAIL",
            factors = "0",
            verified = "false",
            status = "FAIL",
            failure = "ROUTE_ERROR",
        )
    end
end

function _is_staged_argument_error(err)
    err isa ArgumentError || return false
    message = sprint(showerror, err)
    return occursin("staged", message) || occursin("unsupported Laurent GL_n determinant", message)
end

function _cache_error_row(case_id::AbstractString; block::AbstractString)
    return (;
        case = string(case_id),
        block = string(block),
        size = "NA",
        det = "unknown",
        normalization = "NORMALIZATION_SKIP",
        sl_core = "SL_CORE_SKIP",
        gl_cert = "GL_CERT_SKIP",
        factors = "0",
        verified = "false",
        status = "FAIL",
        failure = "CACHE_ERROR",
    )
end

function _unsupported_staged_row(
    case_id::AbstractString,
    block::AbstractString,
    A,
    det_token::AbstractString;
    normalization::AbstractString,
    sl_core::AbstractString,
    gl_cert::AbstractString,
)
    return _row(
        case_id,
        block,
        A,
        det_token;
        normalization,
        sl_core,
        gl_cert,
        factors = "0",
        verified = "false",
        status = "WARN",
        failure = "UNSUPPORTED_STAGED",
    )
end

function _route_error_row(case_id::AbstractString, block::AbstractString, A, det_token::AbstractString)
    return _row(
        case_id,
        block,
        A,
        det_token;
        normalization = "NORMALIZATION_FAIL",
        sl_core = "SL_CORE_SKIP",
        gl_cert = "GL_CERT_FAIL",
        factors = "0",
        verified = "false",
        status = "FAIL",
        failure = "ROUTE_ERROR",
    )
end

function _row(
    case_id::AbstractString,
    block::AbstractString,
    A,
    det_token::AbstractString;
    normalization::AbstractString,
    sl_core::AbstractString,
    gl_cert::AbstractString,
    factors::AbstractString,
    verified::AbstractString,
    status::AbstractString,
    failure::AbstractString,
)
    nrow, ncol = Base.invokelatest(size, A)
    return (;
        case = string(case_id),
        block = string(block),
        size = string(nrow, "x", ncol),
        det = string(det_token),
        normalization = string(normalization),
        sl_core = string(sl_core),
        gl_cert = string(gl_cert),
        factors = string(factors),
        verified = string(verified),
        status = string(status),
        failure = string(failure),
    )
end

function _emit_rows(io, rows)
    for row in rows
        println(
            io,
            "TORIC_SMOKE case=$(row.case) block=$(row.block) size=$(row.size) det=$(row.det) " *
            "normalization=$(row.normalization) sl_core=$(row.sl_core) gl_cert=$(row.gl_cert) " *
            "factors=$(row.factors) verified=$(row.verified) status=$(row.status) failure=$(row.failure)",
        )
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
