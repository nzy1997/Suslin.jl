using Test
using Oscar
using Suslin

const PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_sl3_driver_cases.jl")
const QUILLEN_MAINLINE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl")
const QUILLEN_PATCH_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"
const PARK_WOODBURN_SECTION_5_REF = "refs/arXiv-alg-geom9405003v1 Section 5"

const REQUIRED_SL3_DRIVER_CASE_IDS = Set([
    "sl3-driver-univariate-fast-local-qq",
    "sl3-driver-multivariate-monic-special-form-qq",
    "sl3-driver-quillen-mainline-evidence-gf2",
    "sl3-driver-legacy-quillen-patched-substitution-qq",
    "sl3-driver-det-one-no-witness-staged-qq",
])

const REQUIRED_SL3_DRIVER_NEGATIVE_IDS = Set([
    "sl3-driver-negative-det-not-one",
    "sl3-driver-negative-unsupported-coefficient-ring",
    "sl3-driver-negative-selected-variable-not-generator",
    "sl3-driver-negative-claimed-local-evidence-missing",
    "sl3-driver-negative-supported-without-witness",
])

const REQUIRED_DRIVER_FIELDS = (
    :id,
    :role,
    :expected_status,
    :ring_constructor,
    :ring,
    :matrix,
    :selected_variable,
    :local_form_status,
    :selected_variable_status,
    :supplied_witness_status,
    :upstream_evidence_status,
    :source_refs,
    :consumer_issue_ids,
)

const ALLOWED_DRIVER_STATUS = Set([
    :supported,
    :staged,
])

const ALLOWED_EVIDENCE_STATUS = Set([
    :passes,
    :replayed,
    :missing,
    :absent,
    :inapplicable,
])

function _sl3d_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _sl3d_require_matrix_over(matrix_value, R, n::Int, label)
    nrows(matrix_value) == n && ncols(matrix_value) == n ||
        throw(ArgumentError("$(label) must be a $(n) x $(n) matrix"))
    base_ring(matrix_value) == R ||
        throw(ArgumentError("$(label) must be defined over the fixture ring"))
    return true
end

function _sl3d_extract_first_symbolic_field(data, candidate_fields::Tuple)
    for candidate in candidate_fields
        hasproperty(data, candidate) && return getproperty(data, candidate)
    end
    return nothing
end

function _sl3d_monic_in_variable(p, selected_variable, selected_index::Int)
    R = parent(p)
    target_degree = degree(p, selected_index)
    if target_degree < 0
        return false
    end

    vars = collect(gens(R))
    leading = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        exponents[selected_index] == target_degree || continue
        term = R(coeff)
        for idx in eachindex(vars)
            idx == selected_index && continue
            exponent = exponents[idx]
            exponent == 0 || (term *= vars[idx]^exponent)
        end
        leading += term
    end
    return leading == one(R)
end

function _sl3d_quillen_patch_cases_by_id()
    isfile(QUILLEN_PATCH_CATALOG_PATH) ||
        throw(ArgumentError("validator requires the quillen patch catalog at $(QUILLEN_PATCH_CATALOG_PATH)"))
    isdefined(Main, :QuillenPatchFixtureCatalog) || include(QUILLEN_PATCH_CATALOG_PATH)
    catalog_module = getfield(Main, :QuillenPatchFixtureCatalog)
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function _sl3d_quillen_mainline_cases_by_id()
    isfile(QUILLEN_MAINLINE_CATALOG_PATH) ||
        throw(ArgumentError("validator requires the quillen mainline catalog at $(QUILLEN_MAINLINE_CATALOG_PATH)"))
    isdefined(Main, :QuillenMainlineFixtureCatalog) || include(QUILLEN_MAINLINE_CATALOG_PATH)
    catalog_module = getfield(Main, :QuillenMainlineFixtureCatalog)
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function _sl3d_selected_variable(entry, R)
    selected = _sl3d_field(entry, :selected_variable)
    for field in (:name, :generator, :index, :status)
        hasproperty(selected, field) || throw(
            ArgumentError("fixture $(entry.id) selected variable metadata missing field $(field)")
        )
    end
    selected.name isa AbstractString || selected.name isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) selected variable name must be a string-like symbol"))
    selected.index isa Int || throw(ArgumentError("fixture $(entry.id) selected variable index must be an integer"))
    selected.status isa Symbol || throw(ArgumentError("fixture $(entry.id) selected variable status must be a symbol"))

    generators = _sl3d_field(_sl3d_field(entry, :ring), :generators)
    generator_names = _sl3d_field(_sl3d_field(entry, :ring), :generator_names)
    all(idx -> idx > 0, [selected.index]) || throw(ArgumentError("fixture $(entry.id) selected variable index must be positive"))
    1 <= selected.index <= length(generators) ||
        throw(ArgumentError("fixture $(entry.id) selected variable index is out of bounds"))
    selected.generator == generators[selected.index] ||
        throw(ArgumentError("fixture $(entry.id) selected variable generator and index are inconsistent"))
    parent(selected.generator) == R ||
        throw(ArgumentError("fixture $(entry.id) selected variable must lie in the fixture ring"))
    String(selected.name) == String(generator_names[selected.index]) ||
        throw(ArgumentError("fixture $(entry.id) selected variable name must match ring generator metadata"))
    return selected
end

function _sl3d_assert_local_form_witness(entry, R, selected_variable, selected_index)
    status = _sl3d_field(entry, :local_form_status)
    status in ALLOWED_EVIDENCE_STATUS || throw(
        ArgumentError("fixture $(entry.id) local_form_status must be an evidence status symbol")
    )

    if status ∉ (:passes, :replayed)
        return false
    end

    if !hasproperty(entry, :local_form_witness)
        throw(ArgumentError("fixture $(entry.id) claims local-form support but is missing local_form_witness"))
    end
    witness = _sl3d_field(entry, :local_form_witness)
    monic_value = _sl3d_extract_first_symbolic_field(
        witness,
        (:entry, :polynomial, :p, :local_entry),
    )
    monic_value === nothing && throw(
        ArgumentError("fixture $(entry.id) local_form_witness must include a polynomial replay")
    )
    parent(monic_value) == R ||
        throw(ArgumentError("fixture $(entry.id) local-form replay polynomial must lie in the fixture ring"))
    _sl3d_monic_in_variable(monic_value, selected_variable, selected_index) ||
        throw(ArgumentError("fixture $(entry.id) local-form replay polynomial is not monic in selected variable"))
    if hasproperty(witness, :monic_entry_position)
        witness.monic_entry_position isa Int &&
            1 <= witness.monic_entry_position <= nrows(_sl3d_field(entry, :matrix)) ||
            throw(ArgumentError("fixture $(entry.id) local-form witness position is out of range"))
    end
    return true
end

function _sl3d_assert_supplied_witness(entry, R)
    status = _sl3d_field(entry, :supplied_witness_status)
    status in ALLOWED_EVIDENCE_STATUS || throw(
        ArgumentError("fixture $(entry.id) supplied_witness_status must be an evidence status symbol")
    )

    if status ∉ (:passes, :replayed)
        return false
    end

    hasproperty(entry, :supplied_witness) ||
        throw(ArgumentError("fixture $(entry.id) claims supplied-witness support but is missing supplied_witness"))
    witness = _sl3d_field(entry, :supplied_witness)
    if hasproperty(witness, :expected_matrix)
        _sl3d_require_matrix_over(_sl3d_field(witness, :expected_matrix), R, 3, "supplied_witness.expected_matrix")
        _sl3d_field(witness, :expected_matrix) == _sl3d_field(entry, :matrix) ||
            throw(ArgumentError("fixture $(entry.id) supplied witness replay matrix does not match the fixture matrix"))
    end
    if hasproperty(witness, :variable) && hasproperty(witness, :entry)
        _sl3d_field(witness, :variable) == _sl3d_field(entry, :selected_variable).generator ||
            throw(ArgumentError("fixture $(entry.id) supplied witness variable must match selected variable"))
    end
    return true
end

function _sl3d_replay_case_matrix(mainline_case)
    if hasproperty(mainline_case, :expected_global_product)
        return mainline_case.expected_global_product
    end
    if hasproperty(mainline_case, :expected) && hasproperty(mainline_case.expected, :global_correction)
        return mainline_case.expected.global_correction
    end
    if hasproperty(mainline_case, :patch_case) && hasproperty(mainline_case.patch_case, :target_matrix)
        return mainline_case.patch_case.target_matrix
    end
    return nothing
end

function _sl3d_assert_upstream_evidence(entry, R)
    status = _sl3d_field(entry, :upstream_evidence_status)
    status in ALLOWED_EVIDENCE_STATUS || throw(
        ArgumentError("fixture $(entry.id) upstream_evidence_status must be an evidence status symbol")
    )
    if status ∉ (:passes, :replayed)
        return false
    end

    hasproperty(entry, :upstream_evidence) ||
        throw(ArgumentError("fixture $(entry.id) claims upstream support but is missing upstream_evidence"))
    evidence = _sl3d_field(entry, :upstream_evidence)
    evidence_case_id = _sl3d_extract_first_symbolic_field(
        evidence,
        (:case_id, :mainline_case_id, :upstream_case_id),
    )
    evidence_case_id === nothing && throw(
        ArgumentError("fixture $(entry.id) upstream_evidence must include a case identifier")
    )
    evidence_case_id isa AbstractString && !isempty(evidence_case_id) ||
        throw(ArgumentError("fixture $(entry.id) upstream case id must be a non-empty string"))

    mainline_cases = _sl3d_quillen_mainline_cases_by_id()
    haskey(mainline_cases, evidence_case_id) ||
        throw(ArgumentError("fixture $(entry.id) upstream case id $(evidence_case_id) does not exist"))
    mainline_case = mainline_cases[evidence_case_id]

    replay_matrix = _sl3d_replay_case_matrix(mainline_case)
    replay_matrix === nothing &&
        throw(ArgumentError("fixture $(entry.id) upstream case $(evidence_case_id) does not expose a replayable target matrix"))
    _sl3d_require_matrix_over(replay_matrix, R, 3, "upstream replay matrix")
    replay_matrix == _sl3d_field(entry, :matrix) ||
        throw(ArgumentError("fixture $(entry.id) upstream replay matrix does not match the fixture matrix"))
    if hasproperty(evidence, :patch_case_id)
        patch_cases = _sl3d_quillen_patch_cases_by_id()
        patch_id = _sl3d_field(evidence, :patch_case_id)
        haskey(patch_cases, patch_id) ||
            throw(ArgumentError("fixture $(entry.id) upstream patch_case_id must exist in the Quillen patch catalog"))
    end
    return true
end

function _sl3d_assert_metadata(entry)
    for field in REQUIRED_DRIVER_FIELDS
        _sl3d_field(entry, field)
    end

    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    entry.role isa Symbol || throw(ArgumentError("fixture $(entry.id) role must be a symbol"))
    entry.expected_status isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) expected_status must be a symbol"))
    entry.expected_status in ALLOWED_DRIVER_STATUS ||
        throw(ArgumentError("fixture $(entry.id) expected_status must be :supported or :staged"))

    ring_constructor = _sl3d_field(entry, :ring_constructor)
    _sl3d_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use ordinary polynomial_ring metadata"))
    _sl3d_field(ring_constructor, :coefficient) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring constructor coefficient must be a string"))
    _sl3d_field(ring_constructor, :variables) isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be a tuple"))

    ring = _sl3d_field(entry, :ring)
    _sl3d_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    _sl3d_field(ring, :generator_names) isa Tuple || throw(
        ArgumentError("fixture $(entry.id) ring metadata must include generator_names tuple")
    )
    _sl3d_field(ring, :generators) isa Tuple || throw(
        ArgumentError("fixture $(entry.id) ring metadata must include generators tuple")
    )

    R = _sl3d_field(ring, :object)
    R isa MPolyRing || R isa PolyRing ||
        throw(ArgumentError("fixture $(entry.id) must use an ordinary polynomial ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("fixture $(entry.id) must use an exact coefficient type"))
    coefficient_ring(R) isa Field ||
        throw(ArgumentError("fixture $(entry.id) must be field-backed"))

    generator_names = _sl3d_field(ring, :generator_names)
    generators = _sl3d_field(ring, :generators)
    length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator metadata is inconsistent"))
    ring_constructor.variables == generator_names ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables do not match ring metadata"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator parent mismatch"))

    matrix = _sl3d_field(entry, :matrix)
    _sl3d_require_matrix_over(matrix, R, 3, "matrix")
    all(idx -> parent(matrix[idx, idx]) == R, 1:3) ||
        throw(ArgumentError("fixture $(entry.id) matrix entries must be in the fixture ring"))
    det(matrix) == one(R) || throw(ArgumentError("fixture $(entry.id) matrix determinant must be one"))

    selected_variable = _sl3d_selected_variable(entry, R)
    local_form_pass = _sl3d_assert_local_form_witness(entry, R, selected_variable.generator, selected_variable.index)
    supplied_pass = _sl3d_assert_supplied_witness(entry, R)
    upstream_pass = _sl3d_assert_upstream_evidence(entry, R)
    if entry.expected_status == :supported
        local_form_pass || upstream_pass ||
            throw(ArgumentError("fixture $(entry.id) supported status requires at least one replayable evidence path"))
    else
        (local_form_pass || supplied_pass || upstream_pass) || hasproperty(entry, :staged_reason) ||
            hasproperty(entry, :staged_status_reason) ||
            hasproperty(entry, :reason) ||
            throw(ArgumentError("fixture $(entry.id) staged status requires replay boundary metadata"))
    end

    selected_variable_status = _sl3d_field(entry, :selected_variable_status)
    selected_variable_status in ALLOWED_EVIDENCE_STATUS ||
        throw(ArgumentError("fixture $(entry.id) selected_variable_status must be an evidence status symbol"))
    selected_variable.status in ALLOWED_EVIDENCE_STATUS ||
        throw(ArgumentError("fixture $(entry.id) selected variable status must be an evidence status symbol"))
    selected_variable.status == selected_variable_status ||
        throw(ArgumentError("fixture $(entry.id) selected variable status metadata is inconsistent"))

    entry.source_refs isa Tuple && !isempty(entry.source_refs) ||
        throw(ArgumentError("fixture $(entry.id) must include source refs"))
    for ref in entry.source_refs
        ref isa AbstractString ||
            throw(ArgumentError("fixture $(entry.id) source ref entries must be strings"))
    end
    any(ref -> ref == PARK_WOODBURN_SECTION_3_REF || ref == PARK_WOODBURN_SECTION_5_REF, entry.source_refs) ||
        throw(ArgumentError("fixture $(entry.id) source refs must include Park-Woodburn section 3 or 5"))

    entry.consumer_issue_ids isa Tuple && !isempty(entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) must include consumer issue ids"))
    "#184" in entry.consumer_issue_ids ||
        throw(ArgumentError("fixture $(entry.id) must include consumer issue id #184"))
    all(id -> id isa AbstractString && startswith(id, "#"), entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must look like issue references"))
    return true
end

function validate_park_woodburn_sl3_driver_fixture(entry)
    _sl3d_assert_metadata(entry)
    return true
end

function validate_park_woodburn_sl3_driver_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain valid cases"))

    case_ids = [entry.id for entry in catalog.cases]
    control_ids = [entry.id for entry in catalog.negative_controls]
    all_ids = vcat(case_ids, control_ids)
    length(all_ids) == length(unique(all_ids)) ||
        throw(ArgumentError("catalog case and negative control ids must be unique"))

    for entry in catalog.cases
        validate_park_woodburn_sl3_driver_fixture(entry)
    end
    for entry in catalog.negative_controls
        hasproperty(entry, :base_case_id) &&
            entry.base_case_id isa AbstractString &&
            !isempty(entry.base_case_id) ||
            throw(ArgumentError("negative control $(entry.id) must record base_case_id"))
        hasproperty(entry, :reason) &&
            entry.reason isa AbstractString &&
            !isempty(entry.reason) ||
            throw(ArgumentError("negative control $(entry.id) must record reason"))
        try
            validate_park_woodburn_sl3_driver_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "Park-Woodburn SL3 driver fixture catalog" begin
    if !isdefined(Main, :QuillenMainlineFixtureCatalog)
        include(QUILLEN_MAINLINE_CATALOG_PATH)
    end
    if !isdefined(Main, :QuillenPatchFixtureCatalog)
        include(QUILLEN_PATCH_CATALOG_PATH)
    end
    include(PARK_WOODBURN_SL3_DRIVER_CATALOG_PATH)

    catalog = ParkWoodburnSL3DriverFixtureCatalog.catalog()
    validate_park_woodburn_sl3_driver_fixture_catalog(catalog)
    entries = ParkWoodburnSL3DriverFixtureCatalog.cases_by_id()
    negatives = catalog.negative_controls

    @test REQUIRED_SL3_DRIVER_CASE_IDS ⊆ Set(keys(entries))
    @test REQUIRED_SL3_DRIVER_NEGATIVE_IDS ⊆ Set(entry.id for entry in negatives)
    for entry in values(entries)
        @test "#184" in entry.consumer_issue_ids
        @test det(entry.matrix) == one(_sl3d_field(_sl3d_field(entry, :ring), :object))
    end
    for entry in negatives
        @test_throws ArgumentError validate_park_woodburn_sl3_driver_fixture(entry)
    end
end
