using Test
using Oscar
using Suslin

const STEINBERG_OPTIMIZATION_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "steinberg_optimization_cases.jl")
const STEINBERG_SECTION_6_REF = "refs/arXiv-alg-geom9405003v1 Section 6"

const REQUIRED_STEINBERG_POSITIVE_IDS = Set([
    "steinberg-identity-removal-qq",
    "steinberg-same-position-merge-qq",
    "steinberg-inverse-cancellation-qq",
    "steinberg-commutator-forward-qq",
    "steinberg-commutator-reverse-qq",
    "steinberg-disjoint-commutator-identity-qq",
])

const REQUIRED_STEINBERG_NEGATIVE_IDS = Set([
    "steinberg-negative-mismatched-factor-rings",
    "steinberg-negative-stale-expected-product",
    "steinberg-negative-invalid-commutator-indices",
])

const REQUIRED_STEINBERG_RULE_NAMES = Set([
    :identity_removal,
    :same_position_merge,
    :inverse_cancellation,
    :commutator_forward,
    :commutator_reverse,
    :disjoint_commutator_identity,
])

const REQUIRED_STEINBERG_FIELDS = (
    :id,
    :rule_name,
    :description,
    :ring_constructor,
    :ring,
    :matrix_size,
    :factor_metadata,
    :factors,
    :expected_rewrite_factors,
    :original_product,
    :rewritten_product,
    :rewrite_span,
    :rule_metadata,
    :source_refs,
    :consumer_issue_ids,
)

if !isdefined(Main, :SteinbergOptimizationFixtureCatalog)
    include(STEINBERG_OPTIMIZATION_CATALOG_PATH)
end

function _steinberg_label(entry)
    return hasproperty(entry, :id) ? string(getproperty(entry, :id)) : "<missing id>"
end

function _steinberg_field(entry, field::Symbol)
    hasproperty(entry, field) ||
        throw(ArgumentError("fixture $(_steinberg_label(entry)) missing field $(field)"))
    return getproperty(entry, field)
end

function _steinberg_require_matrix_over(matrix_value, R, n::Int, label)
    matrix_value isa AbstractAlgebra.MatElem ||
        throw(ArgumentError("$(label) must be an Oscar matrix"))
    nrows(matrix_value) == n && ncols(matrix_value) == n ||
        throw(ArgumentError("$(label) must be a square matrix of fixture size"))
    base_ring(matrix_value) == R ||
        throw(ArgumentError("$(label) must be defined over the fixture ring"))
    return true
end

function _steinberg_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

const STEINBERG_EXPECTED_SPAN_LENGTHS = Dict(
    :identity_removal => 2,
    :same_position_merge => 2,
    :inverse_cancellation => 2,
    :commutator_forward => 4,
    :commutator_reverse => 4,
    :disjoint_commutator_identity => 4,
)

function _steinberg_assert_metadata_indices(entry, indices, row::Int, col::Int, rule_label::Symbol)
    _steinberg_field(indices, :i) == row ||
        throw(ArgumentError("fixture $(entry.id) $(rule_label) row index metadata mismatch"))
    _steinberg_field(indices, :j) == col ||
        throw(ArgumentError("fixture $(entry.id) $(rule_label) column index metadata mismatch"))
    return true
end

function _steinberg_assert_ring(entry)
    ring_constructor = _steinberg_field(entry, :ring_constructor)
    _steinberg_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use ordinary polynomial_ring metadata"))
    coefficient = _steinberg_field(ring_constructor, :coefficient)
    coefficient isa AbstractString && !isempty(coefficient) ||
        throw(ArgumentError("fixture $(entry.id) ring constructor coefficient must be a non-empty string"))
    variables = _steinberg_field(ring_constructor, :variables)
    variables isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be a tuple"))
    all(variable -> variable isa AbstractString && !isempty(variable), variables) ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be non-empty strings"))

    ring = _steinberg_field(entry, :ring)
    _steinberg_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    generator_names = _steinberg_field(ring, :generator_names)
    generators = _steinberg_field(ring, :generators)
    generator_names isa Tuple && generators isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring generator metadata must use tuples"))
    length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator metadata is inconsistent"))
    generator_names == variables ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must match ring metadata"))
    all(name -> name isa AbstractString && !isempty(name), generator_names) ||
        throw(ArgumentError("fixture $(entry.id) ring generator names must be non-empty strings"))

    R = _steinberg_field(ring, :object)
    R isa MPolyRing ||
        throw(ArgumentError("fixture $(entry.id) must use an exact field-backed MPolyRing"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("fixture $(entry.id) must use an exact coefficient type"))
    coefficient_ring(R) isa Field ||
        throw(ArgumentError("fixture $(entry.id) must be backed by a field coefficient ring"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator parent mismatch"))
    return R
end

function _steinberg_factor_metadata(entry, field::Symbol, metadata_field::Symbol)
    metadata_value = _steinberg_field(entry, metadata_field)
    if hasproperty(metadata_value, field)
        return getproperty(metadata_value, field)
    end
    return metadata_value
end

function _steinberg_assert_factor_sequence(entry, R, n::Int, field::Symbol, metadata_field::Symbol)
    factors = _steinberg_field(entry, field)
    factors isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) $(field) must be a tuple"))
    metadata = _steinberg_factor_metadata(entry, field, metadata_field)
    metadata isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) $(metadata_field) metadata for $(field) must be a tuple"))
    length(factors) == length(metadata) ||
        throw(ArgumentError("fixture $(entry.id) $(field) length must match $(metadata_field)"))

    for (index, (factor, factor_entry)) in enumerate(zip(factors, metadata))
        _steinberg_require_matrix_over(factor, R, n, "fixture $(entry.id) $(field)[$index]")
        row = _steinberg_field(factor_entry, :row)
        col = _steinberg_field(factor_entry, :col)
        coefficient = _steinberg_field(factor_entry, :coefficient)
        row isa Int && col isa Int ||
            throw(ArgumentError("fixture $(entry.id) $(field)[$index] indices must be integers"))
        1 <= row <= n && 1 <= col <= n && row != col ||
            throw(ArgumentError("fixture $(entry.id) $(field)[$index] indices are invalid"))
        try
            parent(coefficient) == R
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("fixture $(entry.id) $(field)[$index] coefficient has wrong parent ring"))
        end ||
            throw(ArgumentError("fixture $(entry.id) $(field)[$index] coefficient has wrong parent ring"))
        expected_factor = elementary_matrix(n, row, col, coefficient, R)
        factor == expected_factor ||
            throw(ArgumentError("fixture $(entry.id) $(field)[$index] does not match its elementary metadata"))
    end
    return metadata
end

function _steinberg_assert_rule_metadata(entry, R)
    rule_metadata = _steinberg_field(entry, :rule_metadata)
    indices = _steinberg_field(rule_metadata, :indices)
    factors = _steinberg_field(entry, :factors)
    factor_metadata = _steinberg_factor_metadata(entry, :factors, :factor_metadata)
    rewrite_factors = _steinberg_field(entry, :expected_rewrite_factors)
    rewrite_metadata = _steinberg_factor_metadata(entry, :expected_rewrite_factors, :factor_metadata)

    if entry.rule_name == :identity_removal
        length(factors) == 2 && length(rewrite_factors) == 1 ||
            throw(ArgumentError("fixture $(entry.id) identity_removal must use a two-factor rewrite window"))
        _steinberg_field(indices, :i) == factor_metadata[1].row ||
            throw(ArgumentError("fixture $(entry.id) identity_removal row index metadata mismatch"))
        _steinberg_field(indices, :j) == factor_metadata[1].col ||
            throw(ArgumentError("fixture $(entry.id) identity_removal column index metadata mismatch"))
        factor_metadata[1].coefficient == zero(R) ||
            throw(ArgumentError("fixture $(entry.id) identity_removal must remove a zero-coefficient factor"))
        rewrite_metadata[1] == factor_metadata[2] ||
            throw(ArgumentError("fixture $(entry.id) identity_removal rewrite must preserve the non-identity factor"))
    elseif entry.rule_name == :same_position_merge
        length(factors) == 2 && length(rewrite_factors) == 1 ||
            throw(ArgumentError("fixture $(entry.id) same_position_merge must use two factors and one rewrite factor"))
        left = factor_metadata[1]
        right = factor_metadata[2]
        merged = rewrite_metadata[1]
        _steinberg_assert_metadata_indices(entry, indices, left.row, left.col, :same_position_merge)
        left.row == right.row == merged.row ||
            throw(ArgumentError("fixture $(entry.id) same_position_merge row indices must agree"))
        left.col == right.col == merged.col ||
            throw(ArgumentError("fixture $(entry.id) same_position_merge column indices must agree"))
        merged.coefficient == left.coefficient + right.coefficient ||
            throw(ArgumentError("fixture $(entry.id) same_position_merge coefficient must equal the sum"))
    elseif entry.rule_name == :inverse_cancellation
        length(factors) == 2 && isempty(rewrite_factors) ||
            throw(ArgumentError("fixture $(entry.id) inverse_cancellation must rewrite to the empty factor sequence"))
        factor_metadata[1].row == factor_metadata[2].row &&
            factor_metadata[1].col == factor_metadata[2].col ||
            throw(ArgumentError("fixture $(entry.id) inverse_cancellation factors must share row and column indices"))
        _steinberg_assert_metadata_indices(
            entry,
            indices,
            factor_metadata[1].row,
            factor_metadata[1].col,
            :inverse_cancellation,
        )
        factor_metadata[1].coefficient + factor_metadata[2].coefficient == zero(R) ||
            throw(ArgumentError("fixture $(entry.id) inverse_cancellation coefficients must sum to zero"))
    elseif entry.rule_name == :commutator_forward
        length(factors) == 4 && length(rewrite_factors) == 1 ||
            throw(ArgumentError("fixture $(entry.id) commutator_forward must use four factors and one rewrite factor"))
        i = _steinberg_field(indices, :i)
        j = _steinberg_field(indices, :j)
        l = _steinberg_field(indices, :l)
        i != l || throw(ArgumentError("fixture $(entry.id) commutator_forward requires i != l"))
        a = factor_metadata[1].coefficient
        b = factor_metadata[2].coefficient
        factor_metadata[1].row == i && factor_metadata[1].col == j ||
            throw(ArgumentError("fixture $(entry.id) commutator_forward first factor must be E_ij(a)"))
        factor_metadata[2].row == j && factor_metadata[2].col == l ||
            throw(ArgumentError("fixture $(entry.id) commutator_forward second factor must be E_jl(b)"))
        factor_metadata[3].row == i && factor_metadata[3].col == j &&
            factor_metadata[3].coefficient == -a ||
            throw(ArgumentError("fixture $(entry.id) commutator_forward third factor must be E_ij(-a)"))
        factor_metadata[4].row == j && factor_metadata[4].col == l &&
            factor_metadata[4].coefficient == -b ||
            throw(ArgumentError("fixture $(entry.id) commutator_forward fourth factor must be E_jl(-b)"))
        rewrite_metadata[1].row == i && rewrite_metadata[1].col == l &&
            rewrite_metadata[1].coefficient == a * b ||
            throw(ArgumentError("fixture $(entry.id) commutator_forward rewrite must be E_il(ab)"))
    elseif entry.rule_name == :commutator_reverse
        length(factors) == 4 && length(rewrite_factors) == 1 ||
            throw(ArgumentError("fixture $(entry.id) commutator_reverse must use four factors and one rewrite factor"))
        l = _steinberg_field(indices, :l)
        i = _steinberg_field(indices, :i)
        j = _steinberg_field(indices, :j)
        j != l || throw(ArgumentError("fixture $(entry.id) commutator_reverse requires j != l"))
        a = factor_metadata[1].coefficient
        b = factor_metadata[2].coefficient
        factor_metadata[1].row == i && factor_metadata[1].col == j ||
            throw(ArgumentError("fixture $(entry.id) commutator_reverse first factor must be E_ij(a)"))
        factor_metadata[2].row == l && factor_metadata[2].col == i ||
            throw(ArgumentError("fixture $(entry.id) commutator_reverse second factor must be E_li(b)"))
        factor_metadata[3].row == i && factor_metadata[3].col == j &&
            factor_metadata[3].coefficient == -a ||
            throw(ArgumentError("fixture $(entry.id) commutator_reverse third factor must be E_ij(-a)"))
        factor_metadata[4].row == l && factor_metadata[4].col == i &&
            factor_metadata[4].coefficient == -b ||
            throw(ArgumentError("fixture $(entry.id) commutator_reverse fourth factor must be E_li(-b)"))
        rewrite_metadata[1].row == l && rewrite_metadata[1].col == j &&
            rewrite_metadata[1].coefficient == -(a * b) ||
            throw(ArgumentError("fixture $(entry.id) commutator_reverse rewrite must be E_lj(-ab)"))
    elseif entry.rule_name == :disjoint_commutator_identity
        length(factors) == 4 && isempty(rewrite_factors) ||
            throw(ArgumentError("fixture $(entry.id) disjoint_commutator_identity must rewrite to the empty factor sequence"))
        i = _steinberg_field(indices, :i)
        j = _steinberg_field(indices, :j)
        l = _steinberg_field(indices, :l)
        p = _steinberg_field(indices, :p)
        i != p || throw(ArgumentError("fixture $(entry.id) disjoint commutator requires i != p"))
        j != l || throw(ArgumentError("fixture $(entry.id) disjoint commutator requires j != l"))
        a = factor_metadata[1].coefficient
        b = factor_metadata[2].coefficient
        factor_metadata[1].row == i && factor_metadata[1].col == j ||
            throw(ArgumentError("fixture $(entry.id) disjoint commutator first factor must be E_ij(a)"))
        factor_metadata[2].row == l && factor_metadata[2].col == p ||
            throw(ArgumentError("fixture $(entry.id) disjoint commutator second factor must be E_lp(b)"))
        factor_metadata[3].row == i && factor_metadata[3].col == j &&
            factor_metadata[3].coefficient == -a ||
            throw(ArgumentError("fixture $(entry.id) disjoint commutator third factor must be E_ij(-a)"))
        factor_metadata[4].row == l && factor_metadata[4].col == p &&
            factor_metadata[4].coefficient == -b ||
            throw(ArgumentError("fixture $(entry.id) disjoint commutator fourth factor must be E_lp(-b)"))
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported rule_name $(entry.rule_name)"))
    end
    return true
end

function validate_steinberg_optimization_fixture(entry)
    for field in REQUIRED_STEINBERG_FIELDS
        _steinberg_field(entry, field)
    end

    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    entry.rule_name in REQUIRED_STEINBERG_RULE_NAMES ||
        throw(ArgumentError("fixture $(entry.id) must use a known Steinberg rule name"))
    entry.description isa AbstractString && !isempty(entry.description) ||
        throw(ArgumentError("fixture $(entry.id) description must be a non-empty string"))
    entry.matrix_size isa Int && entry.matrix_size >= 2 ||
        throw(ArgumentError("fixture $(entry.id) matrix_size must be an integer at least 2"))

    R = _steinberg_assert_ring(entry)
    n = entry.matrix_size
    _steinberg_assert_factor_sequence(entry, R, n, :factors, :factor_metadata)
    _steinberg_assert_factor_sequence(entry, R, n, :expected_rewrite_factors, :factor_metadata)

    rewrite_span = _steinberg_field(entry, :rewrite_span)
    start_idx = _steinberg_field(rewrite_span, :start)
    stop_idx = _steinberg_field(rewrite_span, :stop)
    start_idx isa Int && stop_idx isa Int ||
        throw(ArgumentError("fixture $(entry.id) rewrite_span indices must be integers"))
    1 <= start_idx <= stop_idx <= length(entry.factors) ||
        throw(ArgumentError("fixture $(entry.id) rewrite_span must be within the factor sequence"))
    stop_idx - start_idx + 1 == STEINBERG_EXPECTED_SPAN_LENGTHS[entry.rule_name] ||
        throw(ArgumentError("fixture $(entry.id) rewrite_span length must match its Steinberg rule"))

    source_refs = _steinberg_field(entry, :source_refs)
    source_refs isa Tuple && !isempty(source_refs) ||
        throw(ArgumentError("fixture $(entry.id) source_refs must be a non-empty tuple"))
    all(ref -> ref isa AbstractString && !isempty(ref), source_refs) ||
        throw(ArgumentError("fixture $(entry.id) source_refs must be non-empty strings"))
    STEINBERG_SECTION_6_REF in source_refs ||
        throw(ArgumentError("fixture $(entry.id) must cite Steinberg Section 6"))

    entry.consumer_issue_ids == ("#188",) ||
        throw(ArgumentError("fixture $(entry.id) consumer_issue_ids must equal (\"#188\",)"))

    _steinberg_require_matrix_over(entry.original_product, R, n, "fixture $(entry.id) original_product")
    _steinberg_require_matrix_over(entry.rewritten_product, R, n, "fixture $(entry.id) rewritten_product")
    entry.original_product == _steinberg_factor_product(entry.factors, R, n) ||
        throw(ArgumentError("fixture $(entry.id) original_product does not replay from its factors"))
    entry.rewritten_product == _steinberg_factor_product(entry.expected_rewrite_factors, R, n) ||
        throw(ArgumentError("fixture $(entry.id) rewritten_product does not replay from its rewrite factors"))
    entry.original_product == entry.rewritten_product ||
        throw(ArgumentError("fixture $(entry.id) original and rewritten products must match exactly"))

    _steinberg_assert_rule_metadata(entry, R)
    return true
end

function validate_steinberg_optimization_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain positive Steinberg cases"))

    case_ids = [entry.id for entry in catalog.cases]
    control_ids = [entry.id for entry in catalog.negative_controls]
    length(case_ids) == length(unique(case_ids)) ||
        throw(ArgumentError("catalog Steinberg positive ids must be unique"))
    length(control_ids) == length(unique(control_ids)) ||
        throw(ArgumentError("catalog Steinberg negative ids must be unique"))
    isempty(intersect(Set(case_ids), Set(control_ids))) ||
        throw(ArgumentError("catalog positive and negative Steinberg ids must be disjoint"))
    Set(case_ids) == REQUIRED_STEINBERG_POSITIVE_IDS ||
        throw(ArgumentError("catalog positive Steinberg ids do not match the required set"))
    Set(control_ids) == REQUIRED_STEINBERG_NEGATIVE_IDS ||
        throw(ArgumentError("catalog negative Steinberg ids do not match the required set"))

    seen_rules = Set{Symbol}()
    for entry in catalog.cases
        validate_steinberg_optimization_fixture(entry)
        push!(seen_rules, entry.rule_name)
    end
    seen_rules == REQUIRED_STEINBERG_RULE_NAMES ||
        throw(ArgumentError("catalog Steinberg positive cases must cover each required rule exactly once"))

    for entry in catalog.negative_controls
        hasproperty(entry, :base_case_id) && entry.base_case_id isa AbstractString &&
            !isempty(entry.base_case_id) ||
            throw(ArgumentError("negative control $(entry.id) must record base_case_id"))
        hasproperty(entry, :reason) && entry.reason isa AbstractString &&
            !isempty(entry.reason) ||
            throw(ArgumentError("negative control $(entry.id) must record reason"))
        failed_as_expected = false
        try
            validate_steinberg_optimization_fixture(entry)
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError || rethrow()
            failed_as_expected = true
        end
        failed_as_expected ||
            throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "Steinberg optimization fixture catalog" begin
    catalog = SteinbergOptimizationFixtureCatalog.catalog()
    entries = SteinbergOptimizationFixtureCatalog.cases_by_id()
    @test entries == Dict(entry.id => entry for entry in catalog.cases)
    @test validate_steinberg_optimization_fixture_catalog(catalog)
    for entry in catalog.negative_controls
        @test_throws ArgumentError validate_steinberg_optimization_fixture(entry)
    end

    unexpectedly_valid_negative = merge(first(catalog.cases), (;
        id = "steinberg-negative-mismatched-factor-rings",
        base_case_id = first(catalog.cases).id,
        reason = "negative controls must fail fixture validation",
    ))
    bad_catalog = (;
        cases = catalog.cases,
        negative_controls = (
            unexpectedly_valid_negative,
            catalog.negative_controls[2],
            catalog.negative_controls[3],
        ),
    )
    @test_throws ArgumentError validate_steinberg_optimization_fixture_catalog(bad_catalog)

    entries = SteinbergOptimizationFixtureCatalog.cases_by_id()
    stale_indices = merge(entries["steinberg-same-position-merge-qq"], (;
        rule_metadata = (; indices = (; i = 2, j = 1)),
    ))
    @test_throws ArgumentError validate_steinberg_optimization_fixture(stale_indices)

    stale_inverse_indices = merge(entries["steinberg-inverse-cancellation-qq"], (;
        rule_metadata = (; indices = (; i = 3, j = 2)),
    ))
    @test_throws ArgumentError validate_steinberg_optimization_fixture(stale_inverse_indices)

    stale_span = merge(entries["steinberg-commutator-forward-qq"], (;
        rewrite_span = (; start = 2, stop = 2),
    ))
    @test_throws ArgumentError validate_steinberg_optimization_fixture(stale_span)
end
