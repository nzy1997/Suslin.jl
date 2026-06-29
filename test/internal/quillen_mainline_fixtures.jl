using Test
using Oscar
using Suslin

const QUILLEN_MAINLINE_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_mainline_cases.jl")
const REQUIRED_QUILLEN_MAINLINE_IDS = Set([
    "quillen-two-open-cover-qq",
    "quillen-nontrivial-multipliers-qq",
    "quillen-patched-substitution-witness-qq",
    "quillen-constructive-acceptance-gf2",
])
const REQUIRED_QUILLEN_MAINLINE_NEGATIVE_IDS = Set([
    "quillen-mainline-uncovered-denominator-control",
    "quillen-mainline-tampered-local-evidence-control",
])
const PARK_WOODBURN_SECTION_3_REF = "refs/arXiv-alg-geom9405003v1 Section 3"

function _qml_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("mainline fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _qml_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _qml_require_matrix_over(matrix_value, R, n::Int, label)
    nrows(matrix_value) == n && ncols(matrix_value) == n ||
        throw(ArgumentError("$(label) must be a square matrix of fixture size"))
    base_ring(matrix_value) == R ||
        throw(ArgumentError("$(label) must be defined over the fixture ring"))
    return true
end

function _qml_assert_patch_case(entry)
    patch_case = _qml_field(entry, :patch_case)
    for field in (
        :id,
        :kind,
        :stage_coverage,
        :ring_constructor,
        :ring,
        :size,
        :target_matrix,
        :base_matrix,
        :denominator_data,
        :local_factors,
        :expected,
        :patched_substitution_witness,
        :source_refs,
        :consumer_issue_ids,
    )
        _qml_field(patch_case, field)
    end

    patch_case.id == entry.id ||
        throw(ArgumentError("mainline fixture $(entry.id) patch_case id must match the wrapper id"))
    return patch_case
end

function _qml_assert_ring_metadata(entry, patch_case)
    ring_constructor = _qml_field(patch_case, :ring_constructor)
    _qml_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("mainline fixture $(entry.id) must use ordinary polynomial_ring metadata"))
    _qml_field(ring_constructor, :coefficient) isa AbstractString ||
        throw(ArgumentError("mainline fixture $(entry.id) ring constructor coefficient must be a string"))
    _qml_field(ring_constructor, :variables) isa Tuple ||
        throw(ArgumentError("mainline fixture $(entry.id) ring constructor variables must be a tuple"))

    ring = _qml_field(patch_case, :ring)
    _qml_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("mainline fixture $(entry.id) ring description must be a string"))
    R = _qml_field(ring, :object)
    R isa MPolyRing || R isa PolyRing ||
        throw(ArgumentError("mainline fixture $(entry.id) must use an ordinary polynomial ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("mainline fixture $(entry.id) must use an exact coefficient ring"))
    generator_names = _qml_field(ring, :generator_names)
    generators = _qml_field(ring, :generators)
    generator_names isa Tuple && generators isa Tuple && length(generator_names) == length(generators) ||
        throw(ArgumentError("mainline fixture $(entry.id) ring generator metadata is inconsistent"))
    ring_constructor.variables == generator_names ||
        throw(ArgumentError("mainline fixture $(entry.id) ring constructor variables do not match ring metadata"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("mainline fixture $(entry.id) ring generator parent mismatch"))

    size = _qml_field(patch_case, :size)
    size isa Int && size >= 2 ||
        throw(ArgumentError("mainline fixture $(entry.id) size must be at least 2"))
    _qml_require_matrix_over(_qml_field(patch_case, :target_matrix), R, size, "target matrix")
    _qml_require_matrix_over(_qml_field(entry, :expected_global_product), R, size, "expected global product")
    patch_case.target_matrix == entry.expected_global_product ||
        throw(ArgumentError("mainline fixture $(entry.id) expected global product must match the reused patch target matrix"))

    return R, size
end

function _qml_assert_denominator_cover(entry, R)
    denominator_cover = _qml_field(entry, :denominator_cover)
    raw_denominator_provenance = _qml_field(entry, :raw_denominator_provenance)
    for field in (:denominators, :multipliers, :coverage_terms, :coverage_sum)
        _qml_field(denominator_cover, field)
    end

    denominators = denominator_cover.denominators
    multipliers = denominator_cover.multipliers
    coverage_terms = denominator_cover.coverage_terms
    denominators isa Tuple && multipliers isa Tuple && coverage_terms isa Tuple ||
        throw(ArgumentError("mainline fixture $(entry.id) denominator cover collections must be tuples"))
    length(denominators) == length(multipliers) == length(coverage_terms) ||
        throw(ArgumentError("mainline fixture $(entry.id) denominator cover collections must align"))
    !isempty(denominators) ||
        throw(ArgumentError("mainline fixture $(entry.id) denominator cover must not be empty"))

    for denominator in denominators
        parent(denominator) == R ||
            throw(ArgumentError("mainline fixture $(entry.id) denominator cover denominator has wrong parent ring"))
    end
    for multiplier in multipliers
        parent(multiplier) == R ||
            throw(ArgumentError("mainline fixture $(entry.id) denominator cover multiplier has wrong parent ring"))
    end

    for idx in eachindex(coverage_terms)
        term = coverage_terms[idx]
        _qml_field(term, :denominator)
        _qml_field(term, :multiplier)
        term.denominator == denominators[idx] ||
            throw(ArgumentError("mainline fixture $(entry.id) denominator cover denominator does not match the explicit tuple entry"))
        term.multiplier == multipliers[idx] ||
            throw(ArgumentError("mainline fixture $(entry.id) denominator cover multiplier does not match the explicit tuple entry"))
        if hasproperty(term, :term)
            term.term == multipliers[idx] * denominators[idx] ||
                throw(ArgumentError("mainline fixture $(entry.id) denominator cover term does not replay"))
        end
        parent(term.denominator) == R ||
            throw(ArgumentError("mainline fixture $(entry.id) denominator cover denominator has wrong parent ring"))
        parent(term.multiplier) == R ||
            throw(ArgumentError("mainline fixture $(entry.id) denominator cover multiplier has wrong parent ring"))
    end

    cover_sum = zero(R)
    for term in coverage_terms
        cover_sum += term.multiplier * term.denominator
    end
    cover_sum == denominator_cover.coverage_sum ||
        throw(ArgumentError("mainline fixture $(entry.id) denominator cover sum does not replay"))
    denominator_cover.coverage_sum == one(R) ||
        throw(ArgumentError("mainline fixture $(entry.id) denominator cover must sum to one"))
    any(
        idx -> denominators[idx] == raw_denominator_provenance.denominator &&
            multipliers[idx] == raw_denominator_provenance.coverage_multiplier,
        eachindex(denominators),
    ) ||
        throw(ArgumentError("mainline fixture $(entry.id) raw denominator provenance must match a denominator cover pair"))
    return true
end

function _qml_assert_raw_denominator_provenance(entry, R)
    provenance = _qml_field(entry, :raw_denominator_provenance)
    for field in (:denominator, :coverage_multiplier, :exponent_l, :source_ref)
        _qml_field(provenance, field)
    end
    parent(provenance.denominator) == R ||
        throw(ArgumentError("mainline fixture $(entry.id) raw denominator has wrong parent ring"))
    parent(provenance.coverage_multiplier) == R ||
        throw(ArgumentError("mainline fixture $(entry.id) raw denominator coverage multiplier has wrong parent ring"))
    provenance.exponent_l isa Integer ||
        throw(ArgumentError("mainline fixture $(entry.id) raw denominator exponent_l must be an integer"))
    provenance.source_ref == PARK_WOODBURN_SECTION_3_REF ||
        throw(ArgumentError("mainline fixture $(entry.id) raw denominator provenance must cite Park-Woodburn Section 3"))
    return true
end

function _qml_assert_local_evidence(entry, patch_case, R, n::Int)
    local_evidence = _qml_field(entry, :local_evidence)
    for field in (:factors, :expected_product)
        _qml_field(local_evidence, field)
    end
    local_evidence.factors isa Tuple && !isempty(local_evidence.factors) ||
        throw(ArgumentError("mainline fixture $(entry.id) local evidence must include factors"))
    _qml_require_matrix_over(local_evidence.expected_product, R, n, "local evidence expected product")
    all(factor -> _qml_require_matrix_over(factor, R, n, "local evidence factor"), local_evidence.factors) ||
        throw(ArgumentError("mainline fixture $(entry.id) local evidence factors are malformed"))
    _qml_product(local_evidence.factors, R, n) == local_evidence.expected_product ||
        throw(ArgumentError("mainline fixture $(entry.id) local evidence factors do not replay"))
    local_evidence.expected_product == patch_case.expected.global_correction ||
        throw(ArgumentError("mainline fixture $(entry.id) local evidence must match the patch-case global correction"))
    return true
end

function _qml_assert_patched_chain(entry, patch_case, R, n::Int)
    patched_chain = _qml_field(entry, :patched_substitution_chain)
    for field in (:variable, :denominator, :exponent, :shift, :expected_matrix)
        _qml_field(patched_chain, field)
    end
    patched_chain.exponent isa Integer ||
        throw(ArgumentError("mainline fixture $(entry.id) patched substitution exponent must be an integer"))
    _qml_require_matrix_over(patched_chain.expected_matrix, R, n, "patched substitution expected matrix")
    patched_chain.variable == patch_case.substitution_variable ||
        throw(ArgumentError("mainline fixture $(entry.id) patched substitution variable must match the patch case"))
    if hasproperty(patch_case, :patched_substitution_witness) && patch_case.patched_substitution_witness !== nothing
        witness = patch_case.patched_substitution_witness
        actual = Suslin.patched_substitution(
            witness.matrix,
            witness.variable,
            witness.denominator,
            witness.exponent,
            witness.shift,
        )
        actual == witness.expected_matrix ||
            throw(ArgumentError("mainline fixture $(entry.id) patched substitution witness does not replay"))
    end
    return true
end

function _qml_assert_base_term_evidence(entry)
    base_term_evidence = _qml_field(entry, :base_term_evidence)
    for field in (:status, :source_ref)
        _qml_field(base_term_evidence, field)
    end
    base_term_evidence.status in (:assumes_identity, :supplied_factors, :staged) ||
        throw(ArgumentError("mainline fixture $(entry.id) base term evidence has an unsupported status"))
    base_term_evidence.source_ref == PARK_WOODBURN_SECTION_3_REF ||
        throw(ArgumentError("mainline fixture $(entry.id) base term evidence must cite Park-Woodburn Section 3"))
    return true
end

function _qml_assert_source_refs(entry)
    source_refs = _qml_field(entry, :source_refs)
    source_refs isa Tuple && !isempty(source_refs) ||
        throw(ArgumentError("mainline fixture $(entry.id) source refs must be a non-empty tuple"))
    PARK_WOODBURN_SECTION_3_REF in source_refs ||
        throw(ArgumentError("mainline fixture $(entry.id) source refs must include Park-Woodburn Section 3"))
    return true
end

function _qml_assert_consumer_ids(entry)
    consumer_issue_ids = _qml_field(entry, :consumer_issue_ids)
    consumer_issue_ids isa Tuple && !isempty(consumer_issue_ids) ||
        throw(ArgumentError("mainline fixture $(entry.id) consumer issue ids must be a non-empty tuple"))
    all(id -> id isa AbstractString && startswith(id, "#"), consumer_issue_ids) ||
        throw(ArgumentError("mainline fixture $(entry.id) consumer issue ids must look like issue references"))
    "#183" in consumer_issue_ids ||
        throw(ArgumentError("mainline fixture $(entry.id) must be consumable by #183"))
    return true
end

function validate_quillen_mainline_fixture(entry)
    patch_case = _qml_assert_patch_case(entry)
    R, n = _qml_assert_ring_metadata(entry, patch_case)
    _qml_assert_denominator_cover(entry, R)
    _qml_assert_raw_denominator_provenance(entry, R)
    _qml_assert_local_evidence(entry, patch_case, R, n)
    _qml_assert_patched_chain(entry, patch_case, R, n)
    _qml_assert_base_term_evidence(entry)
    _qml_assert_source_refs(entry)
    _qml_assert_consumer_ids(entry)
    return true
end

function validate_quillen_mainline_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain valid cases"))

    case_ids = [entry.id for entry in catalog.cases]
    length(case_ids) == length(unique(case_ids)) ||
        throw(ArgumentError("catalog valid case ids must be unique"))
    control_ids = [entry.id for entry in catalog.negative_controls]
    length(control_ids) == length(unique(control_ids)) ||
        throw(ArgumentError("catalog negative control ids must be unique"))
    combined_ids = vcat(case_ids, control_ids)
    length(combined_ids) == length(unique(combined_ids)) ||
        throw(ArgumentError("catalog case and negative control ids must be globally unique"))

    for entry in catalog.cases
        validate_quillen_mainline_fixture(entry)
    end
    isempty(catalog.negative_controls) &&
        throw(ArgumentError("catalog must contain negative controls"))
    for entry in catalog.negative_controls
        try
            validate_quillen_mainline_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "Quillen mainline fixture catalog" begin
    include(QUILLEN_MAINLINE_CATALOG_PATH)
    catalog = QuillenMainlineFixtureCatalog.catalog()
    validate_quillen_mainline_fixture_catalog(catalog)
    entries = QuillenMainlineFixtureCatalog.cases_by_id()
    negatives = Dict(entry.id => entry for entry in catalog.negative_controls)
    @test REQUIRED_QUILLEN_MAINLINE_IDS ⊆ Set(keys(entries))
    @test REQUIRED_QUILLEN_MAINLINE_NEGATIVE_IDS ⊆ Set(keys(negatives))
    @test length(entries) >= 4
    @test length(negatives) >= 2
    for entry in values(negatives)
        @test_throws ArgumentError validate_quillen_mainline_fixture(entry)
    end
end
