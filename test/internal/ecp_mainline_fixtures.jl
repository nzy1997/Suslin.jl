using Test
using Oscar
using Suslin

const ECP_MAINLINE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_mainline_cases.jl")
const PARK_WOODBURN_SECTION_4_REF = "refs/arXiv-alg-geom9405003v1 Section 4"

const REQUIRED_ECP_MAINLINE_CASE_IDS = Set([
    "ecp-mainline-gf2-hard-slice",
    "ecp-mainline-qq-link-bezout",
    "ecp-mainline-length4-coupled-qq",
    "ecp-mainline-monicity-change-gf2",
    "ecp-mainline-sl3-route-qq",
])

const REQUIRED_ECP_MAINLINE_NEGATIVE_IDS = Set([
    "ecp-mainline-negative-non-unimodular",
    "ecp-mainline-negative-corrupt-link-witness",
    "ecp-mainline-negative-selected-variable-not-generator",
    "ecp-mainline-negative-supported-without-evidence",
])

const REQUIRED_ECP_MAINLINE_FIELDS = (
    :id,
    :role,
    :expected_status,
    :ring_constructor,
    :ring,
    :column_entries,
    :column_order,
    :selected_variable,
    :monicity,
    :unimodularity,
    :support_evidence,
    :source_refs,
    :consumer_issue_ids,
)

const ALLOWED_ECP_MAINLINE_STATUS = Set([:supported, :staged])
const ALLOWED_ECP_MAINLINE_EVIDENCE_STATUS = Set([:passes, :replayed, :missing, :absent, :inapplicable])
const REPLAYABLE_ECP_MAINLINE_EVIDENCE_STATUS = Set([:passes, :replayed])
const ECP_MAINLINE_SUPPORT_STATUS_FIELDS = (
    :link_witness_status,
    :link_step_status,
    :lower_variable_status,
    :normality_status,
    :sl3_status,
)
const ECP_MAINLINE_EVIDENCE_FIELD_CANDIDATES = Dict(
    :link_witness_status => (:link_witness, :link_witness_evidence, :link_witness_record),
    :link_step_status => (:link_step, :link_step_certificate, :link_step_evidence),
    :lower_variable_status => (
        :lower_variable,
        :lower_variable_reduction,
        :lower_variable_certificate,
        :lower_reduction,
    ),
    :normality_status => (:normality, :normality_witness, :normality_certificate),
    :sl3_status => (:sl3, :sl3_route, :sl3_certificate, :sl3_evidence),
)

function _ecp_mainline_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _ecp_mainline_column(entry)
    column_entries = _ecp_mainline_field(entry, :column_entries)
    column_order = _ecp_mainline_field(entry, :column_order)
    return tuple((getproperty(column_entries, name) for name in column_order)...)
end

function _ecp_mainline_column_matrix(column, R)
    return matrix(R, length(column), 1, collect(column))
end

function _ecp_mainline_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _ecp_mainline_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _ecp_mainline_apply_factors(factors, column, R)
    return _ecp_mainline_factor_product(factors, R, length(column)) * _ecp_mainline_column_matrix(column, R)
end

function _ecp_mainline_monic_in_variable(p, R, variable_name::Symbol)
    gens_R = collect(gens(R))
    idx = findfirst(g -> Symbol(string(g)) == variable_name, gens_R)
    idx === nothing && throw(ArgumentError("fixture variable must be a generator of the ambient ring"))

    target_degree = degree(p, idx)
    target_degree < 0 && return false

    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        exponents[idx] == target_degree || continue
        term = R(coeff)
        for j in eachindex(gens_R)
            j == idx && continue
            exponent = exponents[j]
            exponent == 0 || (term *= gens_R[j]^exponent)
        end
        total += term
    end
    return total == one(R)
end

function _ecp_mainline_support_evidence_fields(field::Symbol)
    return get(
        ECP_MAINLINE_EVIDENCE_FIELD_CANDIDATES,
        field,
        (Symbol(replace(String(field), r"_status$" => "")),),
    )
end

function _ecp_mainline_support_evidence_status(evidence, field::Symbol, entry_id::AbstractString)
    hasproperty(evidence, field) ||
        throw(ArgumentError("fixture $(entry_id) support_evidence missing field $(field)"))
    status = getproperty(evidence, field)
    status isa Symbol ||
        throw(ArgumentError("fixture $(entry_id) support_evidence field $(field) must be a symbol"))
    status in ALLOWED_ECP_MAINLINE_EVIDENCE_STATUS ||
        throw(ArgumentError("fixture $(entry_id) support_evidence field $(field) must be an evidence status symbol"))
    return status
end

function _ecp_mainline_support_evidence_value(evidence, field::Symbol, entry_id::AbstractString)
    for candidate in _ecp_mainline_support_evidence_fields(field)
        hasproperty(evidence, candidate) || continue
        value = getproperty(evidence, candidate)
        value === nothing && continue
        return value
    end
    return nothing
end

function _ecp_mainline_is_replayable_status(status)
    return status in REPLAYABLE_ECP_MAINLINE_EVIDENCE_STATUS
end

function _ecp_mainline_replay_link_witness(entry, witness)
    if hasproperty(witness, :record)
        witness = witness.record
    elseif hasproperty(witness, :certificate)
        witness = witness.certificate
    end
    witness isa Suslin.ECPLinkWitnessRecord ||
        throw(ArgumentError("fixture $(entry.id) link witness evidence must be a replayable link witness record"))
    Suslin.verify_ecp_link_witness(witness) ||
        throw(ArgumentError("fixture $(entry.id) link witness evidence does not replay"))
    return true
end

function _ecp_mainline_replay_link_step(entry, link_step)
    if hasproperty(link_step, :record)
        link_step = link_step.record
    elseif hasproperty(link_step, :certificate)
        link_step = link_step.certificate
    end
    link_step isa Suslin.ECPLinkStepCertificate ||
        throw(ArgumentError("fixture $(entry.id) link step evidence must be a replayable link-step certificate"))
    Suslin.verify_ecp_link_step_certificate(link_step) ||
        throw(ArgumentError("fixture $(entry.id) link step evidence does not replay"))
    return true
end

function _ecp_mainline_replay_lower_variable(entry, lower_variable, column, R)
    lower_candidate = lower_variable
    if hasproperty(lower_candidate, :record)
        lower_candidate = lower_candidate.record
    elseif hasproperty(lower_candidate, :certificate)
        lower_candidate = lower_candidate.certificate
    end
    if lower_candidate isa Suslin.ECPColumnReductionCertificate
        Suslin.verify_ecp_column_reduction(lower_candidate) ||
            throw(ArgumentError("fixture $(entry.id) lower-variable evidence does not replay"))
        tuple(lower_candidate.original_column...) == column ||
            throw(ArgumentError("fixture $(entry.id) lower-variable certificate does not match the fixture column"))
        lower_candidate.ring == R ||
            throw(ArgumentError("fixture $(entry.id) lower-variable certificate must use the fixture ring"))
        return true
    end
    if hasproperty(lower_candidate, :factors)
        lower_candidate = lower_candidate.factors
    end
    Suslin._ecp_verified_lower_reduction(lower_candidate, column, R)
    return true
end

function _ecp_mainline_validate_supported_evidence(entry, support_evidence, column, R)
    link_witness_status = _ecp_mainline_support_evidence_status(support_evidence, :link_witness_status, entry.id)
    link_step_status = _ecp_mainline_support_evidence_status(support_evidence, :link_step_status, entry.id)
    lower_variable_status = _ecp_mainline_support_evidence_status(support_evidence, :lower_variable_status, entry.id)
    normality_status = _ecp_mainline_support_evidence_status(support_evidence, :normality_status, entry.id)
    sl3_status = _ecp_mainline_support_evidence_status(support_evidence, :sl3_status, entry.id)

    supported = entry.expected_status == :supported
    if supported
        link_step_status in REPLAYABLE_ECP_MAINLINE_EVIDENCE_STATUS ||
            throw(ArgumentError("fixture $(entry.id) supported entry must replay link step evidence"))
        lower_variable_status in REPLAYABLE_ECP_MAINLINE_EVIDENCE_STATUS ||
            throw(ArgumentError("fixture $(entry.id) supported entry must replay lower-variable evidence"))
    end

    if _ecp_mainline_is_replayable_status(link_witness_status)
        witness = _ecp_mainline_support_evidence_value(support_evidence, :link_witness_status, entry.id)
        witness === nothing &&
            throw(ArgumentError("fixture $(entry.id) link witness status requires replayable witness metadata"))
        _ecp_mainline_replay_link_witness(entry, witness)
    end

    if _ecp_mainline_is_replayable_status(link_step_status)
        link_step = _ecp_mainline_support_evidence_value(support_evidence, :link_step_status, entry.id)
        link_step === nothing &&
            throw(ArgumentError("fixture $(entry.id) link step status requires replayable link-step metadata"))
        _ecp_mainline_replay_link_step(entry, link_step)
    end

    if _ecp_mainline_is_replayable_status(lower_variable_status)
        lower_variable = _ecp_mainline_support_evidence_value(support_evidence, :lower_variable_status, entry.id)
        lower_variable === nothing &&
            throw(ArgumentError("fixture $(entry.id) lower-variable status requires replayable lower-variable metadata"))
        _ecp_mainline_replay_lower_variable(entry, lower_variable, column, R)
    end

    if _ecp_mainline_is_replayable_status(normality_status)
        normality = _ecp_mainline_support_evidence_value(support_evidence, :normality_status, entry.id)
        normality === nothing &&
            throw(ArgumentError("fixture $(entry.id) normality status requires replayable normality metadata"))
        if hasproperty(normality, :verification) && hasproperty(normality, :overall_ok)
            normality.overall_ok || throw(ArgumentError("fixture $(entry.id) normality evidence does not replay"))
        end
    end

    if _ecp_mainline_is_replayable_status(sl3_status)
        sl3 = _ecp_mainline_support_evidence_value(support_evidence, :sl3_status, entry.id)
        sl3 === nothing &&
            throw(ArgumentError("fixture $(entry.id) sl3 status requires replayable SL_3 metadata"))
        if hasproperty(sl3, :verification) && hasproperty(sl3, :overall_ok)
            sl3.overall_ok || throw(ArgumentError("fixture $(entry.id) SL_3 evidence does not replay"))
        end
    end

    return true
end

function _ecp_mainline_selected_variable(entry)
    selected = _ecp_mainline_field(entry, :selected_variable)
    for field in (:name, :generator, :index, :status)
        hasproperty(selected, field) || throw(
            ArgumentError("fixture $(entry.id) selected variable metadata missing field $(field)")
        )
    end
    selected.name isa AbstractString || selected.name isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) selected variable name must be string-like"))
    selected.index isa Int || throw(ArgumentError("fixture $(entry.id) selected variable index must be an integer"))
    selected.status isa Symbol || throw(ArgumentError("fixture $(entry.id) selected variable status must be a symbol"))
    return selected
end

function _ecp_mainline_check_ring(entry)
    ring_constructor = _ecp_mainline_field(entry, :ring_constructor)
    _ecp_mainline_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use polynomial_ring constructor metadata"))
    _ecp_mainline_field(ring_constructor, :coefficient) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring constructor coefficient must be a string"))
    _ecp_mainline_field(ring_constructor, :variables) isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be a tuple"))

    ring = _ecp_mainline_field(entry, :ring)
    _ecp_mainline_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    R = _ecp_mainline_field(ring, :object)
    R isa MPolyRing || R isa PolyRing ||
        throw(ArgumentError("fixture $(entry.id) must use an ordinary polynomial ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("fixture $(entry.id) must use an exact coefficient type"))
    coefficient_ring(R) isa Field ||
        throw(ArgumentError("fixture $(entry.id) must be field-backed"))
    generator_names = _ecp_mainline_field(ring, :generator_names)
    generators = _ecp_mainline_field(ring, :generators)
    generator_names isa Tuple && generators isa Tuple && length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generators metadata is inconsistent"))
    ring_constructor.variables == generator_names ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables do not match ring metadata"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator parent mismatch"))
    return R
end

function _ecp_mainline_assert_metadata(entry)
    for field in REQUIRED_ECP_MAINLINE_FIELDS
        _ecp_mainline_field(entry, field)
    end
    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    entry.role isa Symbol || throw(ArgumentError("fixture $(entry.id) role must be a symbol"))
    entry.expected_status in ALLOWED_ECP_MAINLINE_STATUS ||
        throw(ArgumentError("fixture $(entry.id) expected_status must be :supported or :staged"))
    entry.column_order isa Tuple && length(entry.column_order) >= 3 ||
        throw(ArgumentError("fixture $(entry.id) must contain at least three column entries"))
    entry.column_entries isa NamedTuple ||
        throw(ArgumentError("fixture $(entry.id) column_entries must be a NamedTuple"))
    entry.source_refs isa Tuple && PARK_WOODBURN_SECTION_4_REF in entry.source_refs ||
        throw(ArgumentError("fixture $(entry.id) must include the Park-Woodburn Section 4 source ref"))
    entry.consumer_issue_ids isa Tuple && "#185" in entry.consumer_issue_ids ||
        throw(ArgumentError("fixture $(entry.id) must include #185 as a consumer"))
    all(id -> id isa AbstractString && startswith(id, "#"), entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must look like issue references"))
    entry.support_evidence isa NamedTuple ||
        throw(ArgumentError("fixture $(entry.id) support_evidence must be a NamedTuple"))

    selected = _ecp_mainline_selected_variable(entry)
    R = _ecp_mainline_check_ring(entry)
    generators = _ecp_mainline_field(_ecp_mainline_field(entry, :ring), :generators)
    selected.index isa Int && 1 <= selected.index <= length(generators) ||
        throw(ArgumentError("fixture $(entry.id) selected variable index must be in range"))
    selected.generator == generators[selected.index] ||
        throw(ArgumentError("fixture $(entry.id) selected variable generator does not match its index"))
    parent(selected.generator) == R ||
        throw(ArgumentError("fixture $(entry.id) selected variable must lie in the fixture ring"))
    String(selected.name) == String(_ecp_mainline_field(_ecp_mainline_field(entry, :ring), :generator_names)[selected.index]) ||
        throw(ArgumentError("fixture $(entry.id) selected variable name does not match ring metadata"))

    column = _ecp_mainline_column(entry)
    length(column) == length(entry.column_order) ||
        throw(ArgumentError("fixture $(entry.id) column entries do not match column_order"))
    all(entry_name -> parent(getproperty(entry.column_entries, entry_name)) == R, entry.column_order) ||
        throw(ArgumentError("fixture $(entry.id) column entry has the wrong parent ring"))

    monicity = _ecp_mainline_field(entry, :monicity)
    hasproperty(monicity, :status) || throw(ArgumentError("fixture $(entry.id) monicity missing status"))
    hasproperty(monicity, :selected_entry) || throw(ArgumentError("fixture $(entry.id) monicity missing selected_entry"))
    hasproperty(monicity, :transformed_entry) || throw(ArgumentError("fixture $(entry.id) monicity missing transformed_entry"))
    hasproperty(monicity, :selected_variable) || throw(ArgumentError("fixture $(entry.id) monicity missing selected_variable"))

    unimodularity = _ecp_mainline_field(entry, :unimodularity)
    hasproperty(unimodularity, :status) || throw(ArgumentError("fixture $(entry.id) unimodularity missing status"))
    hasproperty(unimodularity, :witness) || throw(ArgumentError("fixture $(entry.id) unimodularity missing witness"))
    hasproperty(unimodularity, :coefficients) || throw(ArgumentError("fixture $(entry.id) unimodularity missing coefficients"))
    _ecp_mainline_support_evidence_status(_ecp_mainline_field(entry, :support_evidence), :link_step_status, entry.id)

    if entry.expected_status == :staged
        hasproperty(entry, :missing_evidence) ||
            throw(ArgumentError("fixture $(entry.id) staged entries must record missing_evidence"))
        entry.missing_evidence isa Tuple ||
            throw(ArgumentError("fixture $(entry.id) missing_evidence must be a tuple"))
    else
        hasproperty(entry, :missing_evidence) &&
            throw(ArgumentError("fixture $(entry.id) supported entries must not claim missing_evidence"))
    end

    return true
end

function _ecp_mainline_assert_unimodularity(entry)
    column = _ecp_mainline_column(entry)
    R = parent(column[1])
    unimodularity = _ecp_mainline_field(entry, :unimodularity)
    status = _ecp_mainline_field(unimodularity, :status)
    status in ALLOWED_ECP_MAINLINE_EVIDENCE_STATUS ||
        throw(ArgumentError("fixture $(entry.id) unimodularity status must be an evidence status symbol"))
    total = sum(_ecp_mainline_field(unimodularity, :coefficients)[idx] * column[idx] for idx in eachindex(column); init = zero(R))
    total == one(R) ||
        throw(ArgumentError("fixture $(entry.id) unimodularity witness does not reconstruct one"))
    if entry.expected_status == :supported
        status in (:passes, :replayed) ||
            throw(ArgumentError("fixture $(entry.id) supported entry must replay unimodularity"))
    end
    return true
end

function _ecp_mainline_assert_monicity(entry)
    monicity = _ecp_mainline_field(entry, :monicity)
    status = _ecp_mainline_field(monicity, :status)
    status in ALLOWED_ECP_MAINLINE_EVIDENCE_STATUS ||
        throw(ArgumentError("fixture $(entry.id) monicity status must be an evidence status symbol"))
    if entry.expected_status == :supported
        status in (:passes, :replayed) ||
            throw(ArgumentError("fixture $(entry.id) supported entry must replay monicity"))
        selected_variable = _ecp_mainline_field(entry, :selected_variable)
        selected_name = Symbol(String(_ecp_mainline_field(selected_variable, :name)))
        transformed_entry = _ecp_mainline_field(monicity, :transformed_entry)
        parent(transformed_entry) == _ecp_mainline_field(_ecp_mainline_field(entry, :ring), :object) ||
            throw(ArgumentError("fixture $(entry.id) monicity transformed entry must lie in the fixture ring"))
        _ecp_mainline_monic_in_variable(transformed_entry, parent(transformed_entry), selected_name) ||
            throw(ArgumentError("fixture $(entry.id) monicity transformed entry is not monic in the selected variable"))
    end
    return true
end

function _ecp_mainline_assert_stage_evidence(entry)
    support_evidence = _ecp_mainline_field(entry, :support_evidence)
    column = _ecp_mainline_column(entry)
    R = parent(column[1])
    for field in ECP_MAINLINE_SUPPORT_STATUS_FIELDS
        _ecp_mainline_support_evidence_status(support_evidence, field, entry.id)
    end
    _ecp_mainline_validate_supported_evidence(entry, support_evidence, column, R)
    if entry.expected_status == :staged
        hasproperty(entry, :missing_evidence) ||
            throw(ArgumentError("fixture $(entry.id) staged entry must declare missing_evidence"))
        entry.missing_evidence isa Tuple ||
            throw(ArgumentError("fixture $(entry.id) staged missing_evidence must be a tuple"))
        isempty(entry.missing_evidence) &&
            throw(ArgumentError("fixture $(entry.id) staged missing_evidence must not be empty"))
    end
    return true
end

function validate_ecp_mainline_fixture(entry)
    _ecp_mainline_assert_metadata(entry)
    _ecp_mainline_assert_unimodularity(entry)
    _ecp_mainline_assert_monicity(entry)
    _ecp_mainline_assert_stage_evidence(entry)
    return true
end

function validate_ecp_mainline_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain valid cases"))
    case_ids = [entry.id for entry in catalog.cases]
    control_ids = [entry.id for entry in catalog.negative_controls]
    all_ids = vcat(case_ids, control_ids)
    length(all_ids) == length(unique(all_ids)) ||
        throw(ArgumentError("catalog case and negative control ids must be unique"))
    for entry in catalog.cases
        validate_ecp_mainline_fixture(entry)
    end
    isempty(catalog.negative_controls) && throw(ArgumentError("catalog must contain negative controls"))
    for entry in catalog.negative_controls
        try
            validate_ecp_mainline_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "ECP mainline fixture catalog" begin
    include(ECP_MAINLINE_CATALOG_PATH)
    catalog = ECPMainlineFixtureCatalog.catalog()
    validate_ecp_mainline_fixture_catalog(catalog)
    entries = ECPMainlineFixtureCatalog.cases_by_id()
    negatives = Dict(entry.id => entry for entry in catalog.negative_controls)
    @test REQUIRED_ECP_MAINLINE_CASE_IDS ⊆ Set(keys(entries))
    @test REQUIRED_ECP_MAINLINE_NEGATIVE_IDS ⊆ Set(keys(negatives))
end
