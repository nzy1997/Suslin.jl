using Test
using Suslin
using Oscar

const ECP_COLUMN_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "ecp_column_cases.jl")
const REQUIRED_ECP_COLUMN_FIELDS = (
    :id,
    :kind,
    :stage_coverage,
    :ring_constructor,
    :ring,
    :variable_order,
    :entries,
    :column_order,
    :monicity,
    :witnesses,
    :expected,
    :source_refs,
    :consumer_issue_ids,
)

const REQUIRED_ECP_COLUMN_IDS = Set([
    "ecp-unit-entry-gf2",
    "ecp-witness-unit-gf2",
    "ecp-variable-change-monic-gf2",
    "ecp-variable-change-permuted-gf2",
    "ecp-link-bezout-nonunit-witness-qq",
    "ecp-longer-embedded-block-gf2",
    "ecp-unsupported-unimodular-gf2",
    "ecp-non-unimodular-gf2",
    "ecp-monic-first-entry-qq",
])

const REQUIRED_ECP_COLUMN_NEGATIVE_IDS = Set([
    "ecp-corrupt-witness-control",
    "ecp-corrupt-monicity-control",
])

function _ecp_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _ecp_column(entry)
    entries = _ecp_field(entry, :entries)
    column_order = _ecp_field(entry, :column_order)
    return tuple((getproperty(entries, name) for name in column_order)...)
end

function _ecp_column_matrix(column, R)
    return matrix(R, length(column), 1, collect(column))
end

function _ecp_target_column(R, n::Int)
    target = zero_matrix(R, n, 1)
    target[n, 1] = one(R)
    return target
end

function _ecp_factor_product(factors, R, n::Int)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function _ecp_apply_factors(factors, column, R)
    return _ecp_factor_product(factors, R, length(column)) * _ecp_column_matrix(column, R)
end

function _ecp_variable_index(entry, variable_name::Symbol)
    variable_order = _ecp_field(entry, :variable_order)
    idx = findfirst(==(variable_name), variable_order)
    idx === nothing && throw(ArgumentError("fixture $(entry.id) missing variable $(variable_name)"))
    return idx
end

function _ecp_monic_in_variable(p, R, variable_name::Symbol)
    gens_R = collect(gens(R))
    idx = findfirst(g -> Symbol(string(g)) == variable_name, gens_R)
    idx === nothing && throw(ArgumentError("variable must be a generator of the ambient ring"))
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

function _ecp_substitution_values(entry, substitution)
    ring_generators = _ecp_field(_ecp_field(entry, :ring), :generators)
    values = Any[]
    for variable_name in _ecp_field(entry, :variable_order)
        if hasproperty(substitution, variable_name)
            push!(values, getproperty(substitution, variable_name))
        else
            idx = findfirst(==(variable_name), _ecp_field(entry, :variable_order))
            push!(values, ring_generators[idx])
        end
    end
    return values
end

function _ecp_assert_metadata(entry)
    for field in REQUIRED_ECP_COLUMN_FIELDS
        _ecp_field(entry, field)
    end
    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    entry.kind isa Symbol || throw(ArgumentError("fixture $(entry.id) kind must be a symbol"))
    entry.stage_coverage isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) stage_coverage must be a symbol"))

    ring_constructor = _ecp_field(entry, :ring_constructor)
    _ecp_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use polynomial_ring constructor metadata"))
    _ecp_field(ring_constructor, :coefficient) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring constructor coefficient must be a string"))
    _ecp_field(ring_constructor, :variables) isa Tuple ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must be a tuple"))

    ring = _ecp_field(entry, :ring)
    _ecp_field(ring, :description) isa AbstractString ||
        throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    R = _ecp_field(ring, :object)
    generator_names = _ecp_field(ring, :generator_names)
    generators = _ecp_field(ring, :generators)
    generator_names isa Tuple && generators isa Tuple && length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generators metadata is inconsistent"))
    ring_constructor.variables == generator_names ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables do not match ring metadata"))
    all(generator -> parent(generator) == R, generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator parent mismatch"))

    entry.variable_order isa Tuple && !isempty(entry.variable_order) ||
        throw(ArgumentError("fixture $(entry.id) variable_order must be a non-empty tuple"))
    Tuple(Symbol(name) for name in generator_names) == entry.variable_order ||
        throw(ArgumentError("fixture $(entry.id) variable_order does not match ring generator names"))
    entry.entries isa NamedTuple || throw(ArgumentError("fixture $(entry.id) entries must be a NamedTuple"))
    entry.column_order isa Tuple && !isempty(entry.column_order) ||
        throw(ArgumentError("fixture $(entry.id) column_order must be a non-empty tuple"))
    for name in entry.column_order
        hasproperty(entry.entries, name) ||
            throw(ArgumentError("fixture $(entry.id) column entry $(name) is missing"))
        parent(getproperty(entry.entries, name)) == R ||
            throw(ArgumentError("fixture $(entry.id) column entry $(name) has wrong parent ring"))
    end
    entry.witnesses isa Tuple || throw(ArgumentError("fixture $(entry.id) witnesses must be a tuple"))
    entry.expected isa NamedTuple || throw(ArgumentError("fixture $(entry.id) expected must be a NamedTuple"))
    entry.source_refs isa Tuple && !isempty(entry.source_refs) ||
        throw(ArgumentError("fixture $(entry.id) must include source refs"))
    entry.consumer_issue_ids isa Tuple && !isempty(entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) must include consumer issue ids"))
    all(id -> id isa AbstractString && startswith(id, "#"), entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must look like issue references"))
    if _ecp_field(entry.expected, :current_status) in (:passes, :staged_fail, :rejects_non_unimodular)
        "#62" in entry.consumer_issue_ids ||
            throw(ArgumentError("fixture $(entry.id) reducer-facing case must be consumable by #62"))
    end
    if entry.kind == :monic_first_entry
        "#87" in entry.consumer_issue_ids ||
            throw(ArgumentError("fixture $(entry.id) monic-first-entry case must be consumable by #87"))
    end
    if any(witness -> _ecp_field(witness, :kind) == :link_bezout, entry.witnesses)
        "#88" in entry.consumer_issue_ids ||
            throw(ArgumentError("fixture $(entry.id) link witness must be consumable by #88"))
    end
    if any(witness -> _ecp_field(witness, :kind) == :missing_link_witness, entry.witnesses)
        "#88" in entry.consumer_issue_ids ||
            throw(ArgumentError("fixture $(entry.id) missing link witness must be consumable by #88"))
    end
    return true
end

function _ecp_assert_unimodularity(entry)
    expected = _ecp_field(entry, :expected)
    current_status = _ecp_field(expected, :current_status)
    column = _ecp_column(entry)
    R = parent(column[1])
    is_unimodular = Suslin.is_unimodular_column(collect(column), R)
    if current_status == :rejects_non_unimodular
        is_unimodular &&
            throw(ArgumentError("fixture $(entry.id) column must be non-unimodular"))
    else
        is_unimodular ||
            throw(ArgumentError("fixture $(entry.id) column must be unimodular"))
    end
    return true
end

function _ecp_assert_current_status(entry)
    expected = _ecp_field(entry, :expected)
    current_status = _ecp_field(expected, :current_status)
    if current_status == :staged_fail
        hasproperty(expected, :message_substring) ||
            throw(ArgumentError("fixture $(entry.id) staged_fail must include message_substring"))
        err = try
            column = _ecp_column(entry)
            R = parent(column[1])
            Suslin.reduce_unimodular_column(collect(column), R)
            nothing
        catch caught
            caught
        end
        err isa ArgumentError ||
            throw(ArgumentError("fixture $(entry.id) expected staged_fail but reduction succeeded"))
        occursin(expected.message_substring, sprint(showerror, err)) ||
            throw(ArgumentError("fixture $(entry.id) staged failure message does not match"))
        return true
    elseif current_status == :passes
        column = _ecp_column(entry)
        R = parent(column[1])
        factors = Suslin.reduce_unimodular_column(collect(column), R)
        _ecp_apply_factors(factors, column, R) == _ecp_target_column(R, length(column)) ||
            throw(ArgumentError("fixture $(entry.id) reducer factors do not reduce to e_n"))
        return true
    elseif current_status == :rejects_non_unimodular
        hasproperty(expected, :message_substring) ||
            throw(ArgumentError("fixture $(entry.id) rejects_non_unimodular must include message_substring"))
        err = try
            column = _ecp_column(entry)
            R = parent(column[1])
            Suslin.reduce_unimodular_column(collect(column), R)
            nothing
        catch caught
            caught
        end
        err isa ArgumentError ||
            throw(ArgumentError("fixture $(entry.id) expected non-unimodular rejection but reduction succeeded"))
        occursin(expected.message_substring, sprint(showerror, err)) ||
            throw(ArgumentError("fixture $(entry.id) rejection message does not match"))
        return true
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported current_status $(current_status)"))
    end
end

function _ecp_assert_monicity(entry)
    monicity = _ecp_field(entry, :monicity)
    variable_name = _ecp_field(monicity, :variable_name)
    selected_entry_name = _ecp_field(monicity, :selected_entry)
    selected_entry = _ecp_field(entry, :entries)
    hasproperty(selected_entry, selected_entry_name) ||
        throw(ArgumentError("fixture $(entry.id) monicity selected entry $(selected_entry_name) is missing"))
    selected = getproperty(selected_entry, selected_entry_name)
    substitution_values = _ecp_substitution_values(entry, _ecp_field(monicity, :substitution))
    transformed = evaluate(selected, substitution_values)
    transformed_entry = _ecp_field(monicity, :transformed_entry)
    R = parent(transformed)
    _ecp_monic_in_variable(transformed_entry, R, variable_name) ||
        throw(ArgumentError("fixture $(entry.id) transformed entry is not monic in $(variable_name)"))
    transformed == transformed_entry ||
        throw(ArgumentError("fixture $(entry.id) transformed entry does not match monicity metadata"))
    return true
end

function _ecp_assert_witness(entry, witness)
    kind = _ecp_field(witness, :kind)
    column = _ecp_column(entry)

    if kind in (:ideal_membership, :link_bezout)
        coefficients = _ecp_field(witness, :coefficients)
        length(coefficients) == length(column) ||
            throw(ArgumentError("fixture $(entry.id) witness coefficients must match column length"))
        R = parent(column[1])
        total = sum(coefficients[idx] * column[idx] for idx in eachindex(column); init = zero(R))
        total == one(R) ||
            throw(ArgumentError("fixture $(entry.id) witness does not reconstruct the unit ideal"))
        get(witness, :require_nonunit_coefficients, false) == true && any(is_unit, coefficients) &&
            throw(ArgumentError("fixture $(entry.id) witness coefficients must be nonunits"))
        if kind == :link_bezout
            is_unit(_ecp_field(witness, :resultant)) ||
                throw(ArgumentError("fixture $(entry.id) link_bezout witness resultant must be a unit"))
            coverage = _ecp_field(witness, :coverage)
            _ecp_field(coverage, :covers_unit_ideal) == true ||
                throw(ArgumentError("fixture $(entry.id) link_bezout witness must cover the unit ideal"))
            path = _ecp_field(witness, :path)
            path isa Tuple && !isempty(path) ||
                throw(ArgumentError("fixture $(entry.id) link_bezout witness must include a path"))
        end
    elseif kind == :missing_link_witness
        _ecp_field(_ecp_field(entry, :expected), :current_status) == :staged_fail ||
            throw(ArgumentError("fixture $(entry.id) missing link witness requires staged_fail"))
        missing = _ecp_field(witness, :missing)
        missing isa Tuple && !isempty(missing) ||
            throw(ArgumentError("fixture $(entry.id) missing link witness must identify missing data"))
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported witness kind $(kind)"))
    end

    return true
end

function _ecp_assert_witnesses(entry)
    for witness in _ecp_field(entry, :witnesses)
        _ecp_assert_witness(entry, witness)
    end
    return true
end

function validate_ecp_column_fixture(entry)
    _ecp_assert_metadata(entry)
    _ecp_assert_unimodularity(entry)
    _ecp_assert_current_status(entry)
    _ecp_assert_monicity(entry)
    _ecp_assert_witnesses(entry)
    return true
end

function validate_ecp_column_fixture_catalog(catalog)
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
        validate_ecp_column_fixture(entry)
    end
    isempty(catalog.negative_controls) && throw(ArgumentError("catalog must contain negative controls"))
    for entry in catalog.negative_controls
        try
            validate_ecp_column_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "ECP column fixture catalog" begin
    include(ECP_COLUMN_CATALOG_PATH)
    catalog = ECPColumnFixtureCatalog.catalog()
    validate_ecp_column_fixture_catalog(catalog)
    entries = ECPColumnFixtureCatalog.cases_by_id()
    negatives = Dict(entry.id => entry for entry in catalog.negative_controls)
    @test REQUIRED_ECP_COLUMN_IDS ⊆ Set(keys(entries))
    @test REQUIRED_ECP_COLUMN_NEGATIVE_IDS ⊆ Set(keys(negatives))
    @test length(entries) >= 8
    for entry in values(negatives)
        @test_throws ArgumentError validate_ecp_column_fixture(entry)
    end

    witness_entry = entries["ecp-witness-unit-gf2"]
    witness = only(witness_entry.witnesses)
    mutated_witness = merge(
        witness,
        (;
            coefficients = (
                zero(parent(witness.coefficients[1])),
                zero(parent(witness.coefficients[2])),
                zero(parent(witness.coefficients[3])),
            ),
        ),
    )
    @test_throws ArgumentError _ecp_assert_witness(witness_entry, mutated_witness)

    monic_entry = entries["ecp-variable-change-monic-gf2"]
    mutated_monicity = merge(
        monic_entry.monicity,
        (;
            transformed_entry = zero(parent(_ecp_column(monic_entry)[1])),
        ),
    )
    mutated_entry = merge(monic_entry, (; monicity = mutated_monicity))
    @test_throws ArgumentError _ecp_assert_monicity(mutated_entry)
end
