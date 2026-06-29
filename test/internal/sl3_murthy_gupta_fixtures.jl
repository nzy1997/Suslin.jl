using Test
using Suslin
using Oscar

const SL3_MURTHY_GUPTA_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "sl3_murthy_gupta_cases.jl")
const REQUIRED_SL3_MURTHY_GUPTA_FIELDS = (
    :id,
    :branch,
    :ring_constructor,
    :ring,
    :variable,
    :entries,
    :target,
    :murthy_path,
    :expected_current_solver,
    :witnesses,
    :source_refs,
    :consumer_issue_ids,
)

const REQUIRED_SL3_MURTHY_GUPTA_BRANCHES = Set([
    :q_degree_normalization,
    :split_lemma,
    :q0_unit_recursion,
    :q0_nonunit_bezout_resultant,
    :open_slice_control,
])

const REQUIRED_SL3_MURTHY_GUPTA_LOCAL_CONTRACT_IDS = Set([
    "mg-local-q-degree-qq-u-x",
    "mg-local-q0-unit-at-u",
    "mg-local-q0-nonunit-bezout-at-u",
])

const REQUIRED_SL3_MURTHY_GUPTA_NEGATIVE_IDS = Set([
    "mg-negative-nonmonic-p",
    "mg-negative-determinant-not-one",
    "mg-negative-corrupted-split-witness",
    "mg-negative-corrupted-local-unit-witness",
    "mg-negative-corrupted-bezout-equality",
])

function _sl3_mg_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _sl3_mg_target(entry)
    return _sl3_mg_field(entry, :target)
end

function _sl3_mg_matrix(R, p, q, r, s)
    return matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
end

function _sl3_mg_product(factors, R)
    product = identity_matrix(R, 3)
    for factor in factors
        product *= factor
    end
    return product
end

function _sl3_mg_monic_in_variable(p, X)
    R = parent(p)
    parent(X) == R || throw(ArgumentError("p and variable must lie in same ring"))
    var_idx = findfirst(isequal(X), collect(gens(R)))
    var_idx === nothing && throw(ArgumentError("variable must be a generator of the ambient ring"))

    target_degree = degree(p, var_idx)
    target_degree < 0 && return false

    vars = collect(gens(R))
    total = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        exponents[var_idx] == target_degree || continue
        term = R(coeff)
        for idx in eachindex(vars)
            idx == var_idx && continue
            exponent = exponents[idx]
            exponent == 0 || (term *= vars[idx]^exponent)
        end
        total += term
    end

    return total == one(R)
end

function _sl3_mg_constant_coefficient(value)
    R = parent(value)
    const_coeff = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(value), AbstractAlgebra.exponent_vectors(value))
        all(exponent -> exponent == 0, exponents) || continue
        const_coeff += R(coeff)
    end
    return const_coeff
end

function _sl3_mg_constant_in_variable(value, X)
    R = parent(value)
    parent(X) == R || throw(ArgumentError("value and variable must lie in same ring"))
    var_idx = findfirst(isequal(X), collect(gens(R)))
    var_idx === nothing && throw(ArgumentError("variable must be a generator of the ambient ring"))

    vars = collect(gens(R))
    constant = zero(R)
    for (coeff, exponents) in zip(AbstractAlgebra.coefficients(value), AbstractAlgebra.exponent_vectors(value))
        exponents[var_idx] == 0 || continue
        term = R(coeff)
        for idx in eachindex(vars)
            idx == var_idx && continue
            exponent = exponents[idx]
            exponent == 0 || (term *= vars[idx]^exponent)
        end
        constant += term
    end
    return constant
end

function _sl3_mg_degree_in_variable(value, X)
    R = parent(value)
    parent(X) == R || throw(ArgumentError("value and variable must lie in same ring"))
    var_idx = findfirst(isequal(X), collect(gens(R)))
    var_idx === nothing && throw(ArgumentError("variable must be a generator of the ambient ring"))
    return degree(value, var_idx)
end

function _sl3_mg_assert_metadata(entry)
    for field in REQUIRED_SL3_MURTHY_GUPTA_FIELDS
        _sl3_mg_field(entry, field)
    end

    hasproperty(_sl3_mg_field(entry, :ring_constructor), :function_name) ||
        throw(ArgumentError("fixture $(entry.id) missing ring_constructor.function_name"))
    hasproperty(entry.ring_constructor, :coefficient) ||
        throw(ArgumentError("fixture $(entry.id) missing ring_constructor.coefficient"))
    hasproperty(entry.ring_constructor, :variables) ||
        throw(ArgumentError("fixture $(entry.id) missing ring_constructor.variables"))

    hasproperty(entry.ring, :description) ||
        throw(ArgumentError("fixture $(entry.id) missing ring.description"))
    hasproperty(entry.ring, :object) ||
        throw(ArgumentError("fixture $(entry.id) missing ring.object"))
    hasproperty(entry.ring, :generator_names) ||
        throw(ArgumentError("fixture $(entry.id) missing ring.generator_names"))
    hasproperty(entry.ring, :generators) ||
        throw(ArgumentError("fixture $(entry.id) missing ring.generators"))

    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be a non-empty string"))
    entry.branch in REQUIRED_SL3_MURTHY_GUPTA_BRANCHES ||
        throw(ArgumentError("fixture $(entry.id) has unsupported branch $(entry.branch)"))
    _sl3_mg_field(entry, :target) isa AbstractAlgebra.Generic.MatSpaceElem ||
        throw(ArgumentError("fixture $(entry.id) target must be a matrix"))
    _sl3_mg_field(entry, :variable) in entry.ring.generators ||
        throw(ArgumentError("fixture $(entry.id) variable is not a ring generator"))
    _sl3_mg_field(entry, :entries) isa NamedTuple ||
        throw(ArgumentError("fixture $(entry.id) entries must be a NamedTuple"))
    for name in (:p, :q, :r, :s)
        hasproperty(entry.entries, name) || throw(ArgumentError("fixture $(entry.id) missing entry $(name)"))
    end
    isfile(SL3_MURTHY_GUPTA_CATALOG_PATH) || throw(ArgumentError("fixture catalog path is invalid"))

    entry.source_refs isa Tuple && !isempty(entry.source_refs) ||
        throw(ArgumentError("fixture $(entry.id) must include at least one source reference"))
    entry.consumer_issue_ids isa Tuple &&
        !isempty(entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) must include at least one consumer issue id"))
    all(id -> isa(id, AbstractString), entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must be strings"))
    all(occursin.("#", entry.consumer_issue_ids)) ||
        throw(ArgumentError("fixture $(entry.id) consumer issue ids must look like issue references"))
    entry.murthy_path == true || throw(ArgumentError("fixture $(entry.id) must set murthy_path=true"))
    if hasproperty(entry, :local_contract) && entry.local_contract == true
        entry.expected_current_solver.status == :staged_fail ||
            throw(ArgumentError("fixture $(entry.id) local-contract entries must be staged_fail"))
        entry.consumer_issue_ids == ("#182", "#208", "#207", "#209", "#210") ||
            throw(ArgumentError("fixture $(entry.id) local-contract consumer issue ids are incorrect"))
    end
    return true
end

function _sl3_mg_assert_target(entry)
    R = entry.ring.object
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    target = _sl3_mg_matrix(R, p, q, r, s)
    target == entry.target ||
        throw(ArgumentError("fixture $(entry.id) target matrix does not match entries"))
    det(target) == one(R) || throw(ArgumentError("fixture $(entry.id) target matrix must have determinant one"))
    _sl3_mg_monic_in_variable(p, entry.variable) || throw(ArgumentError("fixture $(entry.id) p is not monic in variable"))
    return true
end

function _sl3_mg_assert_murthy_path(entry)
    entry.murthy_path || throw(ArgumentError("fixture $(entry.id) must have murthy_path=true"))
    return true
end

function _sl3_mg_assert_current_solver_status(entry)
    status = _sl3_mg_field(entry, :expected_current_solver)
    hasproperty(status, :status) || throw(ArgumentError("fixture $(entry.id) expected_current_solver missing status"))
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    X = entry.variable

    if status.status == :staged_fail
        hasproperty(status, :message_substring) ||
            throw(ArgumentError("fixture $(entry.id) staged_fail must set message_substring"))
        err = try
            Suslin.realize_sl3_local(p, q, r, s, X)
            nothing
        catch caught
            caught
        end
        err isa ArgumentError ||
            throw(ArgumentError("fixture $(entry.id) expected staged_fail but solver succeeded"))
        occursin(status.message_substring, sprint(showerror, err)) ||
            throw(ArgumentError("fixture $(entry.id) staged failure message does not contain expected text"))
        return true
    elseif status.status == :passes
        factors = Suslin.realize_sl3_local(p, q, r, s, X)
        Suslin.verify_factorization(entry.target, factors) || throw(
            ArgumentError("fixture $(entry.id) current solver factors do not verify"),
        )
        return true
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported expected_current_solver.status $(status.status)"))
    end
end

function _sl3_mg_assert_q_degree_witness(entry, witness)
    quotient = _sl3_mg_field(witness, :quotient)
    remainder = _sl3_mg_field(witness, :remainder)
    normalized_s = _sl3_mg_field(witness, :normalized_s)

    R = entry.ring.object
    X = entry.variable
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    q == p * quotient + remainder ||
        throw(ArgumentError("fixture $(entry.id) split relation q = p*quotient + remainder failed"))
    _sl3_mg_degree_in_variable(remainder, X) < _sl3_mg_degree_in_variable(p, X) ||
        throw(ArgumentError("fixture $(entry.id) remainder does not strictly reduce below p degree"))
    normalized_s == s - r * quotient ||
        throw(ArgumentError("fixture $(entry.id) q-degree normalization normalized_s incorrect"))
    lhs = _sl3_mg_matrix(R, p, remainder, r, normalized_s) * Suslin.elementary_matrix(3, 1, 2, quotient, R)
    lhs == entry.target ||
        throw(ArgumentError("fixture $(entry.id) q-degree normalization identity failed"))
    return true
end

function _sl3_mg_assert_split_lemma_witness(entry, witness; expected_target = entry.target)
    a = _sl3_mg_field(witness, :a)
    a_prime = _sl3_mg_field(witness, :a_prime)
    b = _sl3_mg_field(witness, :b)
    c = _sl3_mg_field(witness, :c)
    c1 = _sl3_mg_field(witness, :c1)
    c2 = _sl3_mg_field(witness, :c2)
    d1 = _sl3_mg_field(witness, :d1)
    d2 = _sl3_mg_field(witness, :d2)
    d = _sl3_mg_field(witness, :d)

    R = entry.ring.object
    M = _sl3_mg_matrix(R, a * a_prime, b, c, d)
    M == expected_target ||
        throw(ArgumentError("fixture $(entry.id) split-lemma witness does not reconstruct target"))
    first = _sl3_mg_matrix(R, a, b, c1, d1)
    second = _sl3_mg_matrix(R, a_prime, b, c2, d2)
    det(first) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) split-lemma first factor determinant is not one"))
    det(second) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) split-lemma second factor determinant is not one"))
    rhs =
        Suslin.elementary_matrix(3, 2, 1, c * d1 * d2 - d * (c2 + a_prime * c1 * d2), R) *
        Suslin.elementary_matrix(3, 2, 3, d2 - one(R), R) *
        Suslin.elementary_matrix(3, 3, 2, one(R), R) *
        Suslin.elementary_matrix(3, 2, 3, -one(R), R) *
        first *
        Suslin.elementary_matrix(3, 2, 3, one(R), R) *
        Suslin.elementary_matrix(3, 3, 2, -one(R), R) *
        Suslin.elementary_matrix(3, 2, 3, one(R), R) *
        second *
        Suslin.elementary_matrix(3, 2, 3, -one(R), R) *
        Suslin.elementary_matrix(3, 3, 2, one(R), R) *
        Suslin.elementary_matrix(3, 2, 3, a - one(R), R) *
        Suslin.elementary_matrix(3, 3, 1, -a_prime * c1, R) *
        Suslin.elementary_matrix(3, 3, 2, -d1, R)

    M == rhs || throw(ArgumentError("fixture $(entry.id) split-lemma elementary identity failed"))
    return true
end

function _sl3_mg_assert_q0_unit_witness(entry, witness)
    p0 = _sl3_mg_field(witness, :p0)
    q0 = _sl3_mg_field(witness, :q0)

    R = entry.ring.object
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    if hasproperty(witness, :local_unit_witness)
        p0 == _sl3_mg_constant_in_variable(p, entry.variable) ||
            throw(ArgumentError("fixture $(entry.id) q0-unit p(0) witness is incorrect"))
        q0 == _sl3_mg_constant_in_variable(q, entry.variable) ||
            throw(ArgumentError("fixture $(entry.id) q0-unit q(0) witness is incorrect"))
        _sl3_mg_assert_local_unit_witness(entry, witness.local_unit_witness, q0; label = "q0-unit")
        _sl3_mg_degree_in_variable(q, entry.variable) < _sl3_mg_degree_in_variable(p, entry.variable) ||
            throw(ArgumentError("fixture $(entry.id) q0-unit local witness must already satisfy q-degree guard"))
        hasproperty(witness, :formal_right_e21_coefficient) ||
            throw(ArgumentError("fixture $(entry.id) q0-unit local witness missing formal_right_e21_coefficient"))
        return true
    end

    q0_inverse = _sl3_mg_field(witness, :q0_inverse)
    right_e21_coefficient = _sl3_mg_field(witness, :right_e21_coefficient)
    normalized_p = _sl3_mg_field(witness, :normalized_p)
    normalized_r = _sl3_mg_field(witness, :normalized_r)
    normalized_s = _sl3_mg_field(witness, :normalized_s)
    split = _sl3_mg_field(witness, :split)

    p0 == _sl3_mg_constant_in_variable(p, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) q0-unit p(0) witness is incorrect"))
    q0 == _sl3_mg_constant_in_variable(q, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) q0-unit q(0) witness is incorrect"))
    q0 * q0_inverse == one(R) ||
        throw(ArgumentError("fixture $(entry.id) q0-unit q(0) inverse witness is incorrect"))
    right_e21_coefficient == -q0_inverse * p0 ||
        throw(ArgumentError("fixture $(entry.id) q0-unit E21 coefficient is inconsistent"))
    normalized_p == p + right_e21_coefficient * q ||
        throw(ArgumentError("fixture $(entry.id) q0-unit normalized_p inconsistent"))
    normalized_r == r + right_e21_coefficient * s ||
        throw(ArgumentError("fixture $(entry.id) q0-unit normalized_r inconsistent"))
    normalized_s == s ||
        throw(ArgumentError("fixture $(entry.id) q0-unit normalized_s inconsistent"))
    _sl3_mg_constant_in_variable(normalized_p, entry.variable) == zero(R) ||
        throw(ArgumentError("fixture $(entry.id) q0-unit normalization did not make p(0) zero"))
    target = _sl3_mg_target(entry)
    normalized_target = _sl3_mg_matrix(R, normalized_p, q, normalized_r, normalized_s)
    normalized_target.base_ring == R ||
        throw(ArgumentError("fixture $(entry.id) q0-unit target ring mismatch"))
    det(normalized_target) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) q0-unit normalized matrix has incorrect determinant"))
    target == normalized_target *
        Suslin.elementary_matrix(3, 2, 1, -right_e21_coefficient, R) ||
        throw(ArgumentError("fixture $(entry.id) q0-unit elementary identity failed"))
    _sl3_mg_assert_split_lemma_witness(entry, split; expected_target = normalized_target)
    return true
end

function _sl3_mg_assert_local_unit_witness(entry, witness, expected_unit; label)
    R = entry.ring.object
    context = _sl3_mg_field(witness, :context)
    unit = _sl3_mg_field(witness, :unit)
    residue_unit = _sl3_mg_field(witness, :residue_unit)
    residue_inverse = _sl3_mg_field(witness, :residue_inverse)
    generators = _sl3_mg_field(witness, :maximal_ideal_generators)
    coefficients = _sl3_mg_field(witness, :residue_difference_coefficients)
    global_unit = hasproperty(witness, :global_unit) ? witness.global_unit : is_unit(unit)

    context.kind == :localization_at_maximal_ideal ||
        throw(ArgumentError("fixture $(entry.id) $(label) local context kind is unsupported"))
    context.selected_variable == entry.variable ||
        throw(ArgumentError("fixture $(entry.id) $(label) local context variable mismatch"))
    context.maximal_ideal_generators == generators ||
        throw(ArgumentError("fixture $(entry.id) $(label) local context maximal_ideal_generators mismatch"))
    unit == expected_unit ||
        throw(ArgumentError("fixture $(entry.id) $(label) local unit does not match expected value"))
    length(generators) == length(coefficients) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local witness generator/coefficient length mismatch"))
    all(parent(generator) === R for generator in generators) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local witness generator ring mismatch"))
    all(parent(coefficient) === R for coefficient in coefficients) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local witness coefficient ring mismatch"))

    difference = zero(R)
    for (coefficient, generator) in zip(coefficients, generators)
        difference += coefficient * generator
    end
    unit - residue_unit == difference ||
        throw(ArgumentError("fixture $(entry.id) $(label) local residue equation failed"))
    residue_unit * residue_inverse == one(R) ||
        throw(ArgumentError("fixture $(entry.id) $(label) local residue inverse equation failed"))
    is_unit(unit) == global_unit ||
        throw(ArgumentError("fixture $(entry.id) $(label) global unit flag is incorrect"))
    return true
end

function _sl3_mg_assert_q0_nonunit_bezout_witness(entry, witness)
    p0 = _sl3_mg_field(witness, :p0)
    q0 = _sl3_mg_field(witness, :q0)
    p_prime = _sl3_mg_field(witness, :p_prime)
    q_prime = _sl3_mg_field(witness, :q_prime)
    resultant = _sl3_mg_field(witness, :resultant)
    p_prime_degree = _sl3_mg_field(witness, :p_prime_degree)
    q_prime_degree = _sl3_mg_field(witness, :q_prime_degree)
    branch_unit = _sl3_mg_field(witness, :branch_unit)
    case1_entries = _sl3_mg_field(witness, :case1_entries)

    R = entry.ring.object
    p, q, r, s = entry.entries.p, entry.entries.q, entry.entries.r, entry.entries.s
    p0 == _sl3_mg_constant_in_variable(p, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) q0-nonunit p(0) witness is incorrect"))
    q0 == _sl3_mg_constant_in_variable(q, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) q0-nonunit q(0) witness is incorrect"))
    !is_unit(q0) ||
        throw(ArgumentError("fixture $(entry.id) q0-nonunit witness requires q(0) to be nonunit"))
    p_prime * p - q_prime * q == resultant ||
        throw(ArgumentError("fixture $(entry.id) Bezout equality p_prime*p - q_prime*q failed"))
    if hasproperty(witness, :resultant_unit_witness)
        _sl3_mg_assert_local_unit_witness(entry, witness.resultant_unit_witness, resultant; label = "resultant")
    else
        is_unit(resultant) ||
            throw(ArgumentError("fixture $(entry.id) resultant witness must be a unit"))
    end
    _sl3_mg_degree_in_variable(p_prime, entry.variable) == p_prime_degree ||
        throw(ArgumentError("fixture $(entry.id) p_prime degree witness is incorrect"))
    _sl3_mg_degree_in_variable(q_prime, entry.variable) == q_prime_degree ||
        throw(ArgumentError("fixture $(entry.id) q_prime degree witness is incorrect"))
    p_prime_degree < _sl3_mg_degree_in_variable(q, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) p_prime degree does not satisfy resultant bound"))
    q_prime_degree < _sl3_mg_degree_in_variable(p, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) q_prime degree does not satisfy resultant bound"))
    branch_unit == q0 + _sl3_mg_constant_in_variable(p_prime, entry.variable) ||
        throw(ArgumentError("fixture $(entry.id) q0-nonunit branch unit witness is incorrect"))
    if hasproperty(witness, :branch_unit_witness)
        _sl3_mg_assert_local_unit_witness(entry, witness.branch_unit_witness, branch_unit; label = "branch-unit")
    else
        is_unit(branch_unit) ||
            throw(ArgumentError("fixture $(entry.id) q0-nonunit branch unit must be a unit"))
    end
    case1_entries.p == p + q_prime ||
        throw(ArgumentError("fixture $(entry.id) Case 2 reduction p entry is incorrect"))
    case1_entries.q == q + p_prime ||
        throw(ArgumentError("fixture $(entry.id) Case 2 reduction q entry is incorrect"))
    case1_entries.r == q_prime ||
        throw(ArgumentError("fixture $(entry.id) Case 2 reduction r entry is incorrect"))
    case1_entries.s == p_prime ||
        throw(ArgumentError("fixture $(entry.id) Case 2 reduction s entry is incorrect"))
    target = _sl3_mg_matrix(R, p, q, r, s)

    bezout_matrix = _sl3_mg_matrix(R, p, q, q_prime, p_prime)
    case1_matrix = _sl3_mg_matrix(R, case1_entries.p, case1_entries.q, case1_entries.r, case1_entries.s)
    det(case1_matrix) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) Case 2 q0-unit matrix determinant is not one"))
    _sl3_mg_constant_in_variable(case1_entries.q, entry.variable) == branch_unit ||
        throw(ArgumentError("fixture $(entry.id) Case 2 q0-unit branch constant is incorrect"))
    target ==
        Suslin.elementary_matrix(3, 2, 1, r * p_prime - s * q_prime, R) * bezout_matrix ||
        throw(ArgumentError("fixture $(entry.id) Bezout reduction first equality failed"))
    bezout_matrix ==
        Suslin.elementary_matrix(3, 1, 2, -one(R), R) * case1_matrix ||
        throw(ArgumentError("fixture $(entry.id) Bezout reduction q0-unit equality failed"))
    return true
end

function _sl3_mg_assert_witnesses(entry)
    witnesses = _sl3_mg_field(entry, :witnesses)
    witnesses_data = if witnesses isa NamedTuple
        (witnesses,)
    elseif witnesses isa Tuple
        witnesses
    elseif witnesses === nothing
        ()
    else
        throw(ArgumentError("fixture $(entry.id) witnesses must be a NamedTuple or Tuple"))
    end

    if entry.branch in (:open_slice_control,)
        isempty(witnesses_data) || throw(ArgumentError("fixture $(entry.id) should not define witnesses"))
        return true
    end

    !isempty(witnesses_data) ||
        throw(ArgumentError("fixture $(entry.id) must provide at least one witness"))

    for witness in witnesses_data
        if entry.branch == :q_degree_normalization
            _sl3_mg_assert_q_degree_witness(entry, witness)
        elseif entry.branch == :split_lemma
            _sl3_mg_assert_split_lemma_witness(entry, witness)
        elseif entry.branch == :q0_unit_recursion
            _sl3_mg_assert_q0_unit_witness(entry, witness)
        elseif entry.branch == :q0_nonunit_bezout_resultant
            _sl3_mg_assert_q0_nonunit_bezout_witness(entry, witness)
        else
            throw(ArgumentError("fixture $(entry.id) has unsupported branch $(entry.branch)"))
        end
    end

    return true
end

function validate_sl3_murthy_gupta_fixture(entry)
    _sl3_mg_assert_metadata(entry)
    _sl3_mg_assert_target(entry)
    _sl3_mg_assert_murthy_path(entry)
    _sl3_mg_assert_current_solver_status(entry)
    _sl3_mg_assert_witnesses(entry)
    return true
end

function validate_sl3_murthy_gupta_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("SL3 Murthy-Gupta fixture catalog missing cases"))
    isempty(catalog.cases) && throw(ArgumentError("SL3 Murthy-Gupta fixture catalog must not be empty"))

    ids = [entry.id for entry in catalog.cases]
    length(ids) == length(unique(ids)) || throw(ArgumentError("SL3 Murthy-Gupta fixture ids must be unique"))
    for entry in catalog.cases
        validate_sl3_murthy_gupta_fixture(entry)
    end
    if hasproperty(catalog, :negative_controls)
        negative_ids = [entry.id for entry in catalog.negative_controls]
        length(negative_ids) == length(unique(negative_ids)) ||
            throw(ArgumentError("SL3 Murthy-Gupta negative control ids must be unique"))
    end
    return true
end

@testset "Murthy-Gupta local SL3 fixture catalog" begin
    @test isfile(SL3_MURTHY_GUPTA_CATALOG_PATH)

    include(SL3_MURTHY_GUPTA_CATALOG_PATH)
    catalog = SL3MurthyGuptaFixtureCatalog.catalog()

    @test validate_sl3_murthy_gupta_fixture_catalog(catalog)

    by_id = Dict(entry.id => entry for entry in catalog.cases)
    @test haskey(by_id, "mg-q-degree-normalization")
    @test haskey(by_id, "mg-split-lemma-x-square")
    @test haskey(by_id, "mg-q0-unit-recursion")
    @test haskey(by_id, "mg-q0-nonunit-normalizes-to-q0-unit")
    @test haskey(by_id, "mg-q0-nonunit-normalized-bezout-resultant")
    @test haskey(by_id, "mg-q0-nonunit-extracted-bezout-resultant")
    @test haskey(by_id, "mg-open-slice-control")
    @test REQUIRED_SL3_MURTHY_GUPTA_LOCAL_CONTRACT_IDS ⊆ Set(keys(by_id))

    local_contract_entries = [by_id[id] for id in REQUIRED_SL3_MURTHY_GUPTA_LOCAL_CONTRACT_IDS]
    @test length(local_contract_entries) == 3
    for entry in local_contract_entries
        @test entry.local_contract == true
        @test entry.expected_current_solver.status == :staged_fail
        @test ("#182", "#208", "#207", "#209", "#210") == entry.consumer_issue_ids
    end

    @test hasproperty(catalog, :negative_controls)
    negative_by_id = Dict(entry.id => entry for entry in catalog.negative_controls)
    @test REQUIRED_SL3_MURTHY_GUPTA_NEGATIVE_IDS ⊆ Set(keys(negative_by_id))
    for entry in values(negative_by_id)
        @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(entry)
    end

    @test haskey(by_id, "mg-q-degree-normalization")
    @test by_id["mg-q-degree-normalization"].branch == :q_degree_normalization
    @test by_id["mg-split-lemma-x-square"].branch == :split_lemma
    @test by_id["mg-q0-unit-recursion"].branch == :q0_unit_recursion
    @test by_id["mg-q0-nonunit-normalizes-to-q0-unit"].branch == :q_degree_normalization
    @test by_id["mg-q0-nonunit-normalized-bezout-resultant"].branch == :q0_nonunit_bezout_resultant
    @test by_id["mg-q0-nonunit-extracted-bezout-resultant"].branch == :q0_nonunit_bezout_resultant
    @test by_id["mg-open-slice-control"].branch == :open_slice_control

    supported_nonunit_diagonal = [
        entry.id for entry in catalog.cases
        if entry.expected_current_solver.status == :passes && !is_unit(entry.target[1, 1]) &&
           !is_unit(entry.target[2, 2])
    ]
    @test "mg-q0-nonunit-normalized-bezout-resultant" in supported_nonunit_diagonal
    @test "mg-q0-nonunit-extracted-bezout-resultant" in supported_nonunit_diagonal
    @test length(supported_nonunit_diagonal) >= 2

    split_entry = by_id["mg-split-lemma-x-square"]
    split_witness = first(split_entry.witnesses)
    split_bad = merge(
        split_entry,
        (;
            witnesses = (merge(split_witness, (; a = split_witness.a + one(split_entry.ring.object))),),
        ),
    )
    @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(split_bad)

    q0_unit_entry = by_id["mg-q0-unit-recursion"]
    q0_unit_witness = first(q0_unit_entry.witnesses)
    q0_unit_bad = merge(
        q0_unit_entry,
        (;
            witnesses = (merge(q0_unit_witness, (; q0_inverse = q0_unit_witness.q0_inverse + one(q0_unit_entry.ring.object))),),
        ),
    )
    @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(q0_unit_bad)

    local_q0_unit_entry = by_id["mg-local-q0-unit-at-u"]
    local_q0_unit_witness = first(local_q0_unit_entry.witnesses)
    local_q0_unit_context_bad = merge(
        local_q0_unit_entry,
        (;
            witnesses = (merge(
                local_q0_unit_witness,
                (;
                    local_unit_witness = merge(
                        local_q0_unit_witness.local_unit_witness,
                        (;
                            context = merge(
                                local_q0_unit_witness.local_unit_witness.context,
                                (; maximal_ideal_generators = (local_q0_unit_entry.variable,)),
                            ),
                        ),
                    ),
                ),
            ),),
        ),
    )
    @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(local_q0_unit_context_bad)

    RU, (u, X) = Oscar.polynomial_ring(QQ, ["u", "X"])
    synthetic_q0_unit = (;
        id = "synthetic-q0-unit-selected-variable",
        branch = :q0_unit_recursion,
        ring_constructor = (;
            function_name = :polynomial_ring,
            coefficient = "QQ",
            variables = ("u", "X"),
        ),
        ring = (;
            description = "QQ[u, X]",
            object = RU,
            generator_names = ("u", "X"),
            generators = (u, X),
        ),
        variable = X,
        entries = (;
            p = X + u,
            q = RU(2),
            r = (X + u - one(RU)) * inv(RU(2)),
            s = one(RU),
        ),
        target = _sl3_mg_matrix(
            RU,
            X + u,
            RU(2),
            (X + u - one(RU)) * inv(RU(2)),
            one(RU),
        ),
        murthy_path = true,
        expected_current_solver = (; status = :passes),
        witnesses = ((
            p0 = u,
            q0 = RU(2),
            q0_inverse = inv(RU(2)),
            right_e21_coefficient = -inv(RU(2)) * u,
            normalized_p = X,
            normalized_r = (X - one(RU)) * inv(RU(2)),
            normalized_s = one(RU),
            split = (;
                a = X,
                a_prime = one(RU),
                b = RU(2),
                c = (X - one(RU)) * inv(RU(2)),
                c1 = (X - one(RU)) * inv(RU(2)),
                c2 = zero(RU),
                d1 = one(RU),
                d2 = one(RU),
                d = one(RU),
            ),
        ),),
        source_refs = ("synthetic regression",),
        consumer_issue_ids = ("#206",),
    )
    @test validate_sl3_murthy_gupta_fixture(synthetic_q0_unit)

    bezout_entry = by_id["mg-q0-nonunit-normalized-bezout-resultant"]
    bezout_witness = first(bezout_entry.witnesses)
    bezout_bad = merge(
        bezout_entry,
        (;
            witnesses = (merge(bezout_witness, (; p_prime = bezout_witness.p_prime + one(bezout_entry.ring.object))),),
        ),
    )
    @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(bezout_bad)

    extracted_entry = by_id["mg-q0-nonunit-extracted-bezout-resultant"]
    extracted_witness = first(extracted_entry.witnesses)
    extracted_bad = merge(
        extracted_entry,
        (;
            witnesses = (merge(extracted_witness, (; p_prime = extracted_witness.p_prime + one(extracted_entry.ring.object))),),
        ),
    )
    @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(extracted_bad)

    branch_unit_bad = merge(
        bezout_entry,
        (;
            witnesses = (merge(bezout_witness, (; branch_unit = bezout_witness.branch_unit + bezout_entry.variable)),),
        ),
    )
    @test_throws ArgumentError validate_sl3_murthy_gupta_fixture(branch_unit_bad)
end
