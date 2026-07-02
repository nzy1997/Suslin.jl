using Test
using Oscar
using Suslin

const PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sln_driver_cases.jl")
const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")
const PARK_WOODBURN_SLN_REDUCTION_REF =
    "refs/arXiv-alg-geom9405003v1 Section \"Reduction to SL_3(k[x_1,...,x_m])\""
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"

const REQUIRED_SLN_DRIVER_CASE_IDS = Set([
    "sln-driver-legacy-recursive-column-peel-qq",
    "sln-driver-sl4-gf2-ecp-mainline",
    "sln-driver-sl5-gf2-two-step",
    "sln-driver-sl4-final-sl3-evidence-qq",
    "sln-driver-staged-missing-final-sl3-qq",
])

const REQUIRED_SLN_DRIVER_NEGATIVE_IDS = Set([
    "sln-driver-negative-det-not-one",
    "sln-driver-negative-unsupported-coefficient-ring",
    "sln-driver-negative-corrupt-peel-expectation",
    "sln-driver-negative-unknown-staged-reason",
    "sln-driver-negative-false-mainline-support",
])

const REQUIRED_SLN_DRIVER_FIELDS = (
    :id,
    :support_role,
    :expected_status,
    :route_provenance,
    :ring_constructor,
    :ring,
    :matrix,
    :expected_peel_count,
    :descent_dimensions,
    :peel_steps,
    :final_route,
    :source_refs,
    :consumer_issue_ids,
)

const ALLOWED_SLN_DRIVER_SUPPORT_ROLES = Set([
    :issue186_mainline,
    :staged_issue186_candidate,
    :legacy_regression,
])

const ALLOWED_SLN_DRIVER_STATUS = Set([
    :supported,
    :staged,
])

const ALLOWED_SLN_ECP_STATUS = Set([
    :replayed,
    :missing,
    :absent,
])

const ALLOWED_SLN_FINAL_ROUTE_STATUS = Set([
    :replayed,
    :missing,
    :legacy_regression,
    :absent,
])

const ALLOWED_SLN_STAGED_REASON_CODES = Set([
    :missing_ecp_evidence,
    :missing_final_sl3_route,
    :unsupported_coefficient_ring,
    :legacy_regression_only,
])

function _sln_entry_label(entry)
    return hasproperty(entry, :id) ? string(getproperty(entry, :id)) : "<missing id>"
end

function _sln_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _sln_nested_field(entry, parent::Symbol, child::Symbol)
    parent_value = _sln_field(entry, parent)
    hasproperty(parent_value, child) ||
        throw(ArgumentError("fixture $(_sln_entry_label(entry)) $(parent) metadata missing field $(child)"))
    return getproperty(parent_value, child)
end

function _sln_collect_sequence(value, label)
    value isa Tuple || value isa AbstractVector ||
        throw(ArgumentError("$(label) must be a tuple or vector"))
    return collect(value)
end

function _sln_require_matrix_over(matrix_value, R, n::Int, label)
    matrix_value isa AbstractAlgebra.MatElem ||
        throw(ArgumentError("$(label) must be a matrix"))
    nrows(matrix_value) == n && ncols(matrix_value) == n ||
        throw(ArgumentError("$(label) must be a $(n) x $(n) matrix"))
    base_ring(matrix_value) == R ||
        throw(ArgumentError("$(label) must be defined over the fixture ring"))
    return true
end

function _sln_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _sln_factor_product(factors, R, n::Int, label)
    product = identity_matrix(R, n)
    for factor in _sln_collect_sequence(factors, label)
        _sln_require_matrix_over(factor, R, n, "$(label) entry")
        product *= factor
    end
    return product
end

function _sln_expected_right_factors(after_left, R, n::Int)
    factors = typeof(identity_matrix(R, n))[]
    for col in 1:(n - 1)
        coeff = -after_left[n, col]
        coeff == zero(R) || push!(factors, elementary_matrix(n, n, col, coeff, R))
    end
    return factors
end

function _sln_extract_first_field(data, candidate_fields::Tuple)
    for candidate in candidate_fields
        hasproperty(data, candidate) && return getproperty(data, candidate)
    end
    return nothing
end

function _sln_sl3_catalog_module()
    if isdefined(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        sln_catalog_module = getfield(Main, :ParkWoodburnSLnDriverFixtureCatalog)
        if isdefined(sln_catalog_module, :ParkWoodburnSL3DriverFixtureCatalog)
            return getfield(sln_catalog_module, :ParkWoodburnSL3DriverFixtureCatalog)
        end
    end

    if !isdefined(Main, :ParkWoodburnSL3DriverFixtureCatalog)
        isfile(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH) || throw(
            ArgumentError("validator requires the Park-Woodburn SL3 driver catalog at $(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)")
        )
        include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)
    end
    return getfield(Main, :ParkWoodburnSL3DriverFixtureCatalog)
end

function _sln_sl3_cases_by_id()
    catalog_module = _sln_sl3_catalog_module()
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function _sln_assert_ring_metadata(entry)
    ring_constructor = _sln_field(entry, :ring_constructor)
    _sln_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use ordinary polynomial_ring metadata"))
    _sln_field(ring_constructor, :coefficient) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring constructor coefficient must be a string"))
    _sln_field(ring_constructor, :variables) isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be a tuple"))

    ring = _sln_field(entry, :ring)
    _sln_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    generator_names = _sln_field(ring, :generator_names)
    generators = _sln_field(ring, :generators)
    generator_names isa Tuple && generators isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring metadata must include generator_names and generators tuples"))
    length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator metadata is inconsistent"))
    ring_constructor.variables == generator_names ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables do not match ring metadata"))
    all(name -> name isa AbstractString, generator_names) ||
        throw(ArgumentError("fixture $(entry.id) ring generator names must be strings"))

    R = _sln_field(ring, :object)
    R isa MPolyRing || R isa PolyRing ||
        throw(ArgumentError("fixture $(entry.id) must use an ordinary polynomial ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("fixture $(entry.id) must use an exact coefficient type"))
    coefficient_ring(R) isa Field ||
        throw(ArgumentError("fixture $(entry.id) must be field-backed"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator parent mismatch"))
    return R
end

function _sln_assert_route_provenance(entry)
    provenance = _sln_field(entry, :route_provenance)
    route = _sln_extract_first_field(provenance, (:route, :driver, :algorithm))
    route isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) route_provenance must include symbolic route metadata"))
    _sln_field(provenance, :reduction_ref) == PARK_WOODBURN_SLN_REDUCTION_REF ||
        throw(ArgumentError("fixture $(entry.id) route_provenance must record the Park-Woodburn SLn reduction ref"))
    source = _sln_extract_first_field(provenance, (:source, :source_case_id, :description))
    source isa AbstractString || source isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) route_provenance must include source metadata"))
    return true
end

function _sln_assert_issue_ids(entry)
    issue_ids = _sln_field(entry, :consumer_issue_ids)
    issue_ids isa Tuple && !isempty(issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) must include consumer issue ids"))
    all(id -> id isa AbstractString && occursin(r"^#[0-9]+$", id), issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must be #number references"))
    issue_ids == ("#186",) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must be exactly (\"#186\",)"))
    return true
end

function _sln_assert_source_refs(entry; ecp_replayed::Bool, final_route_replayed::Bool)
    source_refs = _sln_field(entry, :source_refs)
    source_refs isa Tuple && !isempty(source_refs) ||
        throw(ArgumentError("fixture $(entry.id) must include source refs"))
    all(ref -> ref isa AbstractString, source_refs) ||
        throw(ArgumentError("fixture $(entry.id) source refs must be strings"))
    PARK_WOODBURN_SLN_REDUCTION_REF in source_refs ||
        throw(ArgumentError("fixture $(entry.id) must include the Park-Woodburn SLn reduction source ref"))
    ecp_replayed && !(PARK_WOODBURN_SECTION_4_REF in source_refs) &&
        throw(ArgumentError("fixture $(entry.id) with replayed ECP steps must include the Park-Woodburn Section 4 ref"))
    final_route_replayed && !(PARK_WOODBURN_SECTION_5_REF in source_refs) &&
        throw(ArgumentError("fixture $(entry.id) with replayed final SL3 evidence must include the Park-Woodburn Section 5 ref"))
    (PARK_WOODBURN_SECTION_4_REF in source_refs || PARK_WOODBURN_SECTION_5_REF in source_refs) ||
        throw(ArgumentError("fixture $(entry.id) must include a relevant Park-Woodburn Section 4 or 5 ref"))
    return true
end

function _sln_assert_support_metadata(entry)
    support_role = _sln_field(entry, :support_role)
    expected_status = _sln_field(entry, :expected_status)
    support_role isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) support_role must be a symbol"))
    expected_status isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) expected_status must be a symbol"))
    support_role in ALLOWED_SLN_DRIVER_SUPPORT_ROLES ||
        throw(ArgumentError("fixture $(entry.id) support_role must be recognized"))
    expected_status in ALLOWED_SLN_DRIVER_STATUS ||
        throw(ArgumentError("fixture $(entry.id) expected_status must be :supported or :staged"))

    if support_role == :issue186_mainline
        expected_status == :supported ||
            throw(ArgumentError("fixture $(entry.id) issue186 mainline entries must be supported"))
    elseif support_role == :staged_issue186_candidate
        expected_status == :staged ||
            throw(ArgumentError("fixture $(entry.id) staged issue186 candidates must remain staged"))
    elseif support_role == :legacy_regression
        expected_status == :staged ||
            throw(ArgumentError("fixture $(entry.id) legacy regression entries must remain staged"))
    end
    return true
end

function _sln_assert_staged_reason_codes(entry)
    expected_status = _sln_field(entry, :expected_status)
    support_role = _sln_field(entry, :support_role)
    codes = hasproperty(entry, :staged_reason_codes) ? entry.staged_reason_codes : ()
    codes isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) staged_reason_codes must be a tuple"))
    all(code -> code isa Symbol, codes) ||
        throw(ArgumentError("fixture $(entry.id) staged_reason_codes must be symbols"))
    all(code -> code in ALLOWED_SLN_STAGED_REASON_CODES, codes) ||
        throw(ArgumentError("fixture $(entry.id) staged_reason_codes include an unknown reason code"))

    if expected_status == :staged
        isempty(codes) &&
            throw(ArgumentError("fixture $(entry.id) staged entries must include staged_reason_codes"))
    else
        isempty(codes) ||
            throw(ArgumentError("fixture $(entry.id) supported entries must not carry staged_reason_codes"))
    end

    if :legacy_regression_only in codes
        support_role == :legacy_regression ||
            throw(ArgumentError("fixture $(entry.id) legacy regression reason is not mainline support"))
    end
    support_role == :legacy_regression && !(:legacy_regression_only in codes) &&
        throw(ArgumentError("fixture $(entry.id) legacy regression entries must carry :legacy_regression_only"))
    return codes
end

function _sln_assert_last_column_ecp(entry, step, R, d::Int, current_matrix)
    last_column = _sln_collect_sequence(_sln_field(step, :last_column), "fixture $(entry.id) peel step last_column")
    length(last_column) == d ||
        throw(ArgumentError("fixture $(entry.id) peel step last_column length must match dimension"))
    all(value -> parent(value) == R, last_column) ||
        throw(ArgumentError("fixture $(entry.id) peel step last_column entries must lie in the fixture ring"))
    last_column == [current_matrix[row, d] for row in 1:d] ||
        throw(ArgumentError("fixture $(entry.id) peel step last_column does not match the input matrix"))

    ecp = _sln_field(step, :last_column_ecp)
    status = _sln_field(ecp, :status)
    status isa Symbol && status in ALLOWED_SLN_ECP_STATUS ||
        throw(ArgumentError("fixture $(entry.id) peel step last_column_ecp status must be recognized"))
    target_column = _sln_target_column(R, d)
    _sln_field(ecp, :target_column) == target_column ||
        throw(ArgumentError("fixture $(entry.id) peel step ECP target column metadata is wrong"))

    if hasproperty(ecp, :source_case_id)
        source_case_id = ecp.source_case_id
        source_case_id === nothing || source_case_id isa AbstractString ||
            throw(ArgumentError("fixture $(entry.id) peel step ECP source_case_id must be a string or nothing"))
    end

    if status != :replayed
        return (; status, left_product = nothing)
    end

    certificate = _sln_field(ecp, :certificate)
    certificate isa Suslin.ECPColumnReductionCertificate ||
        throw(ArgumentError("fixture $(entry.id) replayed ECP metadata must include an ECP certificate"))
    Suslin.verify_ecp_column_reduction(certificate) ||
        throw(ArgumentError("fixture $(entry.id) replayed ECP certificate does not verify"))
    certificate.ring == R ||
        throw(ArgumentError("fixture $(entry.id) replayed ECP certificate ring mismatch"))
    certificate.original_column == last_column ||
        throw(ArgumentError("fixture $(entry.id) replayed ECP certificate original column mismatch"))
    certificate.final_column == target_column ||
        throw(ArgumentError("fixture $(entry.id) replayed ECP certificate final column mismatch"))

    left_product = _sln_factor_product(certificate.factors, R, d, "fixture $(entry.id) replayed ECP factors")
    left_product * matrix(R, d, 1, last_column) == target_column ||
        throw(ArgumentError("fixture $(entry.id) replayed ECP factors do not send the last column to e_d"))
    return (; status, left_product)
end

function _sln_assert_right_clearing(entry, step, R, d::Int, current_matrix, left_product)
    right_clearing = _sln_field(step, :right_clearing)
    _sln_field(right_clearing, :status) == :replayed ||
        throw(ArgumentError("fixture $(entry.id) right-clearing metadata must be replayed"))

    after_left = _sln_field(right_clearing, :after_left_matrix)
    _sln_require_matrix_over(after_left, R, d, "fixture $(entry.id) right-clearing after_left_matrix")
    if left_product !== nothing
        left_product * current_matrix == after_left ||
            throw(ArgumentError("fixture $(entry.id) right-clearing after_left_matrix does not replay the ECP left product"))
    end
    after_left[:, d:d] == _sln_target_column(R, d) ||
        throw(ArgumentError("fixture $(entry.id) right-clearing after_left_matrix must have normalized final column"))

    right_factors = _sln_collect_sequence(
        _sln_field(right_clearing, :right_factors),
        "fixture $(entry.id) right-clearing factors",
    )
    right_factors == _sln_expected_right_factors(after_left, R, d) ||
        throw(ArgumentError("fixture $(entry.id) right-clearing factors do not match the recorded after-left matrix"))
    right_product = _sln_factor_product(right_factors, R, d, "fixture $(entry.id) right-clearing factors")

    peeled_matrix = _sln_field(right_clearing, :peeled_matrix)
    _sln_require_matrix_over(peeled_matrix, R, d, "fixture $(entry.id) right-clearing peeled_matrix")
    after_left * right_product == peeled_matrix ||
        throw(ArgumentError("fixture $(entry.id) right-clearing factors do not replay the peeled matrix"))

    next_block = _sln_field(step, :next_block)
    _sln_require_matrix_over(next_block, R, d - 1, "fixture $(entry.id) peel step next_block")
    _sln_field(right_clearing, :next_block) == next_block ||
        throw(ArgumentError("fixture $(entry.id) right-clearing next_block metadata is inconsistent"))
    peeled_matrix == block_embedding(next_block, d, collect(1:(d - 1))) ||
        throw(ArgumentError("fixture $(entry.id) right-clearing peeled_matrix must be the embedded next block"))
    det(next_block) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) peel step next_block determinant must be one"))
    return next_block
end

function _sln_assert_peel_step(entry, step, R, current_matrix)
    d = _sln_field(step, :dimension)
    d isa Int && d >= 4 ||
        throw(ArgumentError("fixture $(entry.id) peel step dimension must be an integer at least 4"))
    nrows(current_matrix) == d && ncols(current_matrix) == d ||
        throw(ArgumentError("fixture $(entry.id) peel step dimension must match the current matrix"))

    input_matrix = _sln_field(step, :input_matrix)
    _sln_require_matrix_over(input_matrix, R, d, "fixture $(entry.id) peel step input_matrix")
    input_matrix == current_matrix ||
        throw(ArgumentError("fixture $(entry.id) peel step input_matrix does not match the current descent matrix"))

    ecp = _sln_assert_last_column_ecp(entry, step, R, d, current_matrix)
    next_block = _sln_assert_right_clearing(entry, step, R, d, current_matrix, ecp.left_product)
    return (; next_block, ecp_status = ecp.status)
end

function _sln_assert_peel_metadata(entry, R)
    matrix_value = _sln_field(entry, :matrix)
    n = nrows(matrix_value)
    peel_steps = _sln_collect_sequence(_sln_field(entry, :peel_steps), "fixture $(entry.id) peel_steps")
    expected_peel_count = _sln_field(entry, :expected_peel_count)
    expected_peel_count isa Int && expected_peel_count >= 1 ||
        throw(ArgumentError("fixture $(entry.id) expected_peel_count must be a positive integer"))
    length(peel_steps) == expected_peel_count ||
        throw(ArgumentError("fixture $(entry.id) expected_peel_count does not match peel_steps"))

    current_matrix = matrix_value
    ecp_statuses = Symbol[]
    observed_dimensions = Int[n]
    for step in peel_steps
        result = _sln_assert_peel_step(entry, step, R, current_matrix)
        push!(ecp_statuses, result.ecp_status)
        current_matrix = result.next_block
        push!(observed_dimensions, nrows(current_matrix))
    end
    observed_dimensions[end] == 3 ||
        throw(ArgumentError("fixture $(entry.id) recursive descent must end at an SL3 block"))

    descent_dimensions = Tuple(_sln_collect_sequence(
        _sln_field(entry, :descent_dimensions),
        "fixture $(entry.id) descent_dimensions",
    ))
    descent_dimensions == Tuple(observed_dimensions) ||
        throw(ArgumentError("fixture $(entry.id) descent_dimensions do not match peel_steps"))
    return (;
        final_block = current_matrix,
        ecp_statuses,
        all_ecp_replayed = all(status -> status == :replayed, ecp_statuses),
        any_ecp_replayed = any(status -> status == :replayed, ecp_statuses),
    )
end

function _sln_assert_final_route(entry, R, final_block)
    final_route = _sln_field(entry, :final_route)
    status = _sln_field(final_route, :status)
    status isa Symbol && status in ALLOWED_SLN_FINAL_ROUTE_STATUS ||
        throw(ArgumentError("fixture $(entry.id) final_route status must be recognized"))

    if status == :replayed
        case_id = _sln_extract_first_field(final_route, (:case_id, :sl3_case_id, :source_case_id))
        case_id isa AbstractString && !isempty(case_id) ||
            throw(ArgumentError("fixture $(entry.id) replayed final_route must include an SL3 case id"))
        sl3_cases = _sln_sl3_cases_by_id()
        haskey(sl3_cases, case_id) ||
            throw(ArgumentError("fixture $(entry.id) final_route references unknown SL3 case $(case_id)"))
        sl3_entry = sl3_cases[case_id]
        hasproperty(sl3_entry, :expected_status) && sl3_entry.expected_status == :supported ||
            throw(ArgumentError("fixture $(entry.id) final_route must reference a supported SL3 driver case"))
        _sln_require_matrix_over(sl3_entry.matrix, R, 3, "fixture $(entry.id) referenced SL3 matrix")
        sl3_entry.matrix == final_block ||
            throw(ArgumentError("fixture $(entry.id) final_route SL3 matrix does not match the peeled final block"))
        if hasproperty(final_route, :matrix)
            _sln_require_matrix_over(final_route.matrix, R, 3, "fixture $(entry.id) final_route matrix")
            final_route.matrix == final_block ||
                throw(ArgumentError("fixture $(entry.id) final_route matrix does not match the peeled final block"))
        end
        return true
    end

    if status == :missing
        entry.expected_status == :staged ||
            throw(ArgumentError("fixture $(entry.id) cannot be supported with missing final SL3 route evidence"))
        :missing_final_sl3_route in _sln_assert_staged_reason_codes(entry) ||
            throw(ArgumentError("fixture $(entry.id) missing final_route requires :missing_final_sl3_route"))
    elseif status == :legacy_regression
        entry.support_role == :legacy_regression ||
            throw(ArgumentError("fixture $(entry.id) legacy final_route status is only for legacy regressions"))
    elseif status == :absent
        entry.expected_status == :staged ||
            throw(ArgumentError("fixture $(entry.id) cannot be supported with absent final SL3 route evidence"))
    end
    return false
end

function _sln_assert_metadata(entry)
    for field in REQUIRED_SLN_DRIVER_FIELDS
        _sln_field(entry, field)
    end

    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    _sln_assert_support_metadata(entry)
    staged_reason_codes = _sln_assert_staged_reason_codes(entry)
    _sln_assert_route_provenance(entry)
    _sln_assert_issue_ids(entry)

    R = _sln_assert_ring_metadata(entry)
    matrix_value = _sln_field(entry, :matrix)
    matrix_value isa AbstractAlgebra.MatElem ||
        throw(ArgumentError("fixture $(entry.id) matrix must be a matrix"))
    n = nrows(matrix_value)
    n >= 4 || throw(ArgumentError("fixture $(entry.id) matrix size must be at least 4"))
    _sln_require_matrix_over(matrix_value, R, n, "fixture $(entry.id) matrix")
    all(idx -> parent(matrix_value[idx, idx]) == R, 1:n) ||
        throw(ArgumentError("fixture $(entry.id) matrix entries must lie in the fixture ring"))
    det(matrix_value) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) matrix determinant must be one"))

    peel_result = _sln_assert_peel_metadata(entry, R)
    final_route_replayed = _sln_assert_final_route(entry, R, peel_result.final_block)
    _sln_assert_source_refs(
        entry;
        ecp_replayed = peel_result.any_ecp_replayed,
        final_route_replayed,
    )

    if entry.support_role == :issue186_mainline
        peel_result.all_ecp_replayed ||
            throw(ArgumentError("fixture $(entry.id) issue186 mainline support requires replayed ECP for every peel step"))
        final_route_replayed ||
            throw(ArgumentError("fixture $(entry.id) issue186 mainline support requires replayed final SL3 evidence"))
    end
    if :missing_ecp_evidence in staged_reason_codes && peel_result.all_ecp_replayed
        throw(ArgumentError("fixture $(entry.id) claims missing ECP evidence but every peel step replayed"))
    end
    return true
end

function validate_park_woodburn_sln_driver_fixture(entry)
    _sln_assert_metadata(entry)
    return true
end

function validate_park_woodburn_sln_driver_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain valid cases"))
    isempty(catalog.negative_controls) && throw(ArgumentError("catalog must contain negative controls"))

    case_ids = [entry.id for entry in catalog.cases]
    control_ids = [entry.id for entry in catalog.negative_controls]
    all_ids = vcat(case_ids, control_ids)
    length(all_ids) == length(unique(all_ids)) ||
        throw(ArgumentError("catalog case and negative control ids must be unique"))
    issubset(REQUIRED_SLN_DRIVER_CASE_IDS, Set(case_ids)) ||
        throw(ArgumentError("catalog missing required SLn driver case ids"))
    issubset(REQUIRED_SLN_DRIVER_NEGATIVE_IDS, Set(control_ids)) ||
        throw(ArgumentError("catalog missing required SLn driver negative-control ids"))

    valid_case_ids = Set(case_ids)
    for entry in catalog.cases
        validate_park_woodburn_sln_driver_fixture(entry)
    end

    for entry in catalog.negative_controls
        hasproperty(entry, :base_case_id) &&
            entry.base_case_id isa AbstractString &&
            entry.base_case_id in valid_case_ids ||
            throw(ArgumentError("negative control $(entry.id) must record a valid base_case_id"))
        hasproperty(entry, :reason) &&
            entry.reason isa AbstractString &&
            !isempty(entry.reason) ||
            throw(ArgumentError("negative control $(entry.id) must record reason"))
        try
            validate_park_woodburn_sln_driver_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "Park-Woodburn SLn driver fixture catalog" begin
    include(PARK_WOODBURN_SLN_DRIVER_CATALOG_PATH)

    catalog = ParkWoodburnSLnDriverFixtureCatalog.catalog()
    @test validate_park_woodburn_sln_driver_fixture_catalog(catalog)

    entries = ParkWoodburnSLnDriverFixtureCatalog.cases_by_id()
    @test Set(keys(entries)) == Set(entry.id for entry in catalog.cases)
    @test issubset(REQUIRED_SLN_DRIVER_CASE_IDS, Set(keys(entries)))
    @test issubset(REQUIRED_SLN_DRIVER_NEGATIVE_IDS, Set(entry.id for entry in catalog.negative_controls))

    for entry in values(entries)
        @test validate_park_woodburn_sln_driver_fixture(entry)
        @test entry.consumer_issue_ids == ("#186",)
    end
    for entry in catalog.negative_controls
        @test_throws ArgumentError validate_park_woodburn_sln_driver_fixture(entry)
    end
end
