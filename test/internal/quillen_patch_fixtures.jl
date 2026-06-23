using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
const REQUIRED_QUILLEN_PATCH_FIELDS = (
    :id,
    :kind,
    :stage_coverage,
    :ring_constructor,
    :ring,
    :size,
    :substitution_variable,
    :target_matrix,
    :base_matrix,
    :denominator_data,
    :local_factors,
    :expected,
    :patched_substitution_witness,
    :source_refs,
    :consumer_issue_ids,
)

const REQUIRED_QUILLEN_PATCH_IDS = Set([
    "quillen-two-open-cover-qq",
    "quillen-nontrivial-multipliers-qq",
    "quillen-supplied-local-certificate-gf2",
    "quillen-patched-substitution-witness-qq",
    "quillen-constructive-acceptance-gf2",
])

const REQUIRED_QUILLEN_PATCH_NEGATIVE_IDS = Set([
    "quillen-uncovered-denominator-control",
    "quillen-tampered-local-factor-control",
])

function _quillen_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _quillen_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _quillen_require_matrix_over(matrix_value, R, n::Int, label)
    nrows(matrix_value) == n && ncols(matrix_value) == n ||
        throw(ArgumentError("$(label) must be a square matrix of fixture size"))
    base_ring(matrix_value) == R ||
        throw(ArgumentError("$(label) must be defined over the fixture ring"))
    return true
end

function _quillen_assert_metadata(entry)
    for field in REQUIRED_QUILLEN_PATCH_FIELDS
        _quillen_field(entry, field)
    end

    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    entry.kind isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) kind must be a symbol"))
    entry.stage_coverage isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) stage_coverage must be a symbol"))

    ring_constructor = _quillen_field(entry, :ring_constructor)
    _quillen_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use ordinary polynomial_ring metadata"))
    _quillen_field(ring_constructor, :coefficient) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring constructor coefficient must be a string"))
    _quillen_field(ring_constructor, :variables) isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be a tuple"))

    ring = _quillen_field(entry, :ring)
    _quillen_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    R = _quillen_field(ring, :object)
    R isa MPolyRing || R isa PolyRing ||
        throw(ArgumentError("fixture $(entry.id) must use an ordinary polynomial ring"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("fixture $(entry.id) must use an exact coefficient ring"))
    generator_names = _quillen_field(ring, :generator_names)
    generators = _quillen_field(ring, :generators)
    generator_names isa Tuple && generators isa Tuple && length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator metadata is inconsistent"))
    ring_constructor.variables == generator_names ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables do not match ring metadata"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator parent mismatch"))

    entry.size isa Int && entry.size >= 2 ||
        throw(ArgumentError("fixture $(entry.id) size must be at least 2"))
    _quillen_require_matrix_over(entry.target_matrix, R, entry.size, "target matrix")
    entry.base_matrix === nothing ||
        _quillen_require_matrix_over(entry.base_matrix, R, entry.size, "base matrix")
    entry.substitution_variable in generators ||
        throw(ArgumentError("fixture $(entry.id) substitution variable must be a ring generator"))

    entry.denominator_data isa Tuple && !isempty(entry.denominator_data) ||
        throw(ArgumentError("fixture $(entry.id) denominator_data must be a non-empty tuple"))
    for data in entry.denominator_data
        parent(_quillen_field(data, :denominator)) == R ||
            throw(ArgumentError("fixture $(entry.id) denominator has wrong parent ring"))
        parent(_quillen_field(data, :coverage_multiplier)) == R ||
            throw(ArgumentError("fixture $(entry.id) coverage multiplier has wrong parent ring"))
    end

    entry.local_factors isa Tuple && !isempty(entry.local_factors) ||
        throw(ArgumentError("fixture $(entry.id) local_factors must be a non-empty tuple"))
    entry.expected isa NamedTuple ||
        throw(ArgumentError("fixture $(entry.id) expected must be a NamedTuple"))
    entry.source_refs isa Tuple && !isempty(entry.source_refs) ||
        throw(ArgumentError("fixture $(entry.id) source_refs must be a non-empty tuple"))
    entry.consumer_issue_ids isa Tuple && !isempty(entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer_issue_ids must be a non-empty tuple"))
    all(id -> id isa AbstractString && startswith(id, "#"), entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must look like issue references"))
    "#63" in entry.consumer_issue_ids ||
        throw(ArgumentError("fixture $(entry.id) must be consumable by #63"))
    return true
end

function _quillen_assert_status(entry)
    status = _quillen_field(_quillen_field(entry, :expected), :current_status)
    R = _quillen_field(_quillen_field(entry, :ring), :object)
    if status == :passes
        global_correction = _quillen_field(entry.expected, :global_correction)
        _quillen_require_matrix_over(global_correction, R, entry.size, "expected global correction")
        det(global_correction) == one(R) ||
            throw(ArgumentError("fixture $(entry.id) expected global correction must have determinant one"))
    elseif status == :staged_fail
        missing = _quillen_field(entry.expected, :missing)
        missing isa Tuple && !isempty(missing) ||
            throw(ArgumentError("fixture $(entry.id) staged_fail must include missing witness metadata"))
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported current_status $(status)"))
    end
    return true
end

function _quillen_assert_denominator_cover(entry)
    R = _quillen_field(_quillen_field(entry, :ring), :object)
    total = sum(
        data.coverage_multiplier * data.denominator
        for data in entry.denominator_data;
        init = zero(R),
    )
    total == one(R) ||
        throw(ArgumentError("fixture $(entry.id) denominator coverage must sum to one"))
    return true
end

function _quillen_assert_local_factor(entry, local_factor)
    R = _quillen_field(_quillen_field(entry, :ring), :object)
    n = entry.size
    certificate = _quillen_field(local_factor, :certificate)
    indices = _quillen_field(certificate, :indices)
    denominators = _quillen_field(certificate, :denominators)
    indices isa AbstractVector && denominators isa AbstractVector && length(indices) == length(denominators) ||
        throw(ArgumentError("fixture $(entry.id) local certificate metadata is inconsistent"))

    denominator = _quillen_field(local_factor, :denominator)
    coverage_multiplier = _quillen_field(local_factor, :coverage_multiplier)
    parent(denominator) == R ||
        throw(ArgumentError("fixture $(entry.id) local denominator has wrong parent ring"))
    parent(coverage_multiplier) == R ||
        throw(ArgumentError("fixture $(entry.id) local coverage multiplier has wrong parent ring"))
    any(data -> data.denominator == denominator && data.coverage_multiplier == coverage_multiplier, entry.denominator_data) ||
        throw(ArgumentError("fixture $(entry.id) local factor denominator data is not cataloged"))

    correction = _quillen_field(local_factor, :correction)
    row = _quillen_field(correction, :row)
    col = _quillen_field(correction, :col)
    entry_value = _quillen_field(correction, :entry)
    row isa Int && col isa Int && 1 <= row <= n && 1 <= col <= n && row != col ||
        throw(ArgumentError("fixture $(entry.id) local correction indices are invalid"))
    row in indices && col in indices ||
        throw(ArgumentError("fixture $(entry.id) local certificate must include correction indices"))
    parent(entry_value) == R ||
        throw(ArgumentError("fixture $(entry.id) local correction entry has wrong parent ring"))
    for (index, certificate_denominator) in zip(indices, denominators)
        parent(certificate_denominator) == R ||
            throw(ArgumentError("fixture $(entry.id) certificate denominator has wrong parent ring"))
        if (index == row || index == col) && certificate_denominator != denominator
            throw(ArgumentError("fixture $(entry.id) certificate denominator does not match local denominator"))
        end
    end

    weighted_entry = coverage_multiplier * denominator * entry_value
    expected_correction = elementary_matrix(n, row, col, weighted_entry, R)
    _quillen_require_matrix_over(_quillen_field(local_factor, :expected_correction), R, n, "local expected correction")
    _quillen_require_matrix_over(_quillen_field(local_factor, :factor), R, n, "local factor")
    local_factor.expected_correction == expected_correction ||
        throw(ArgumentError("fixture $(entry.id) local expected correction is incorrect"))
    local_factor.factor == expected_correction ||
        throw(ArgumentError("fixture $(entry.id) local factor does not match expected correction"))
    return true
end

function _quillen_assert_local_factors(entry)
    for local_factor in entry.local_factors
        _quillen_assert_local_factor(entry, local_factor)
    end
    R = _quillen_field(_quillen_field(entry, :ring), :object)
    product = _quillen_factor_product([local_factor.factor for local_factor in entry.local_factors], R, entry.size)
    product == entry.expected.global_correction ||
        throw(ArgumentError("fixture $(entry.id) local factor product does not match expected global correction"))
    product == entry.target_matrix ||
        throw(ArgumentError("fixture $(entry.id) local factor product does not match target matrix"))
    return true
end

function _quillen_assert_patched_substitution_witness(entry)
    witness = entry.patched_substitution_witness
    witness === nothing && return true
    for field in (:matrix, :variable, :denominator, :exponent, :shift, :expected_matrix)
        _quillen_field(witness, field)
    end
    witness.exponent isa Integer ||
        throw(ArgumentError("fixture $(entry.id) patched substitution exponent must be an integer"))
    actual = Suslin.patched_substitution(
        witness.matrix,
        witness.variable,
        witness.denominator,
        witness.exponent,
        witness.shift,
    )
    actual == witness.expected_matrix ||
        throw(ArgumentError("fixture $(entry.id) patched substitution witness does not replay"))
    witness.variable == entry.substitution_variable ||
        throw(ArgumentError("fixture $(entry.id) patched substitution variable must match fixture variable"))
    return true
end

function validate_quillen_patch_fixture(entry)
    _quillen_assert_metadata(entry)
    _quillen_assert_status(entry)
    _quillen_assert_denominator_cover(entry)
    _quillen_assert_local_factors(entry)
    _quillen_assert_patched_substitution_witness(entry)
    return true
end

function validate_quillen_patch_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain valid cases"))
    case_ids = [entry.id for entry in catalog.cases]
    length(case_ids) == length(unique(case_ids)) ||
        throw(ArgumentError("catalog valid case ids must be unique"))
    control_ids = [entry.id for entry in catalog.negative_controls]
    length(control_ids) == length(unique(control_ids)) ||
        throw(ArgumentError("catalog negative control ids must be unique"))
    for entry in catalog.cases
        validate_quillen_patch_fixture(entry)
    end
    isempty(catalog.negative_controls) &&
        throw(ArgumentError("catalog must contain negative controls"))
    for entry in catalog.negative_controls
        try
            validate_quillen_patch_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "Quillen patch fixture catalog" begin
    include(QUILLEN_PATCH_CATALOG_PATH)
    catalog = QuillenPatchFixtureCatalog.catalog()
    validate_quillen_patch_fixture_catalog(catalog)
    entries = QuillenPatchFixtureCatalog.cases_by_id()
    negatives = Dict(entry.id => entry for entry in catalog.negative_controls)
    @test REQUIRED_QUILLEN_PATCH_IDS ⊆ Set(keys(entries))
    @test REQUIRED_QUILLEN_PATCH_NEGATIVE_IDS ⊆ Set(keys(negatives))
    @test length(entries) >= 5
    for entry in values(negatives)
        @test_throws ArgumentError validate_quillen_patch_fixture(entry)
    end

    nontrivial_entry = entries["quillen-nontrivial-multipliers-qq"]
    R_nontrivial = nontrivial_entry.ring.object
    @test any(data -> data.coverage_multiplier != one(R_nontrivial), nontrivial_entry.denominator_data)

    witness_entry = entries["quillen-patched-substitution-witness-qq"]
    witness = witness_entry.patched_substitution_witness
    @test witness.variable == witness_entry.substitution_variable
    @test witness.exponent == 2
    @test Suslin.patched_substitution(
        witness.matrix,
        witness.variable,
        witness.denominator,
        witness.exponent,
        witness.shift,
    ) == witness.expected_matrix

    cover_entry = entries["quillen-two-open-cover-qq"]
    R_cover = cover_entry.ring.object
    mutated_cover = merge(cover_entry, (;
        denominator_data = (
            cover_entry.denominator_data[1],
            merge(cover_entry.denominator_data[2], (;
                coverage_multiplier = zero(R_cover),
            )),
        ),
    ))
    @test_throws ArgumentError validate_quillen_patch_fixture(mutated_cover)

    factor_entry = entries["quillen-supplied-local-certificate-gf2"]
    R_factor = factor_entry.ring.object
    mutated_factor = merge(factor_entry, (;
        local_factors = (
            merge(factor_entry.local_factors[1], (;
                factor = factor_entry.local_factors[1].factor *
                    elementary_matrix(factor_entry.size, 1, 3, one(R_factor), R_factor),
            )),
            factor_entry.local_factors[2],
        ),
    ))
    @test_throws ArgumentError validate_quillen_patch_fixture(mutated_factor)
end
