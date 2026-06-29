using Test
using Oscar
using Suslin

const PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "polynomial_normality_cases.jl")

const REQUIRED_PW_POLYNOMIAL_NORMALITY_IDS = Set([
    "pw-section2-cohn-type-qq",
    "pw-section2-orthogonal-rank-one-qq",
    "pw-section2-conjugated-elementary-qq",
])

const REQUIRED_PW_POLYNOMIAL_NORMALITY_NEGATIVE_IDS = Set([
    "pw-section2-cohn-type-tampered-target-control",
    "pw-section2-rank-one-bad-orthogonality-control",
    "pw-section2-conjugated-elementary-tampered-target-control",
])

const REQUIRED_NORMALITY_CASE_FIELDS = (
    :id,
    :section2_layer,
    :ring_constructor,
    :ring,
    :inputs,
    :target_matrix,
    :expected_convention,
    :consumer_issue_ids,
)

const EXPECTED_SECTION2_LAYER_BY_ID = Dict(
    "pw-section2-cohn-type-qq" => :cohn_type,
    "pw-section2-orthogonal-rank-one-qq" => :orthogonal_rank_one,
    "pw-section2-conjugated-elementary-qq" => :conjugated_elementary,
)

const SECTION2_LAYERS = Set([
    :cohn_type,
    :orthogonal_rank_one,
    :conjugated_elementary,
])

const EXPECTED_RING_COEFFICIENT = "QQ"
const EXPECTED_RING_VARIABLES = ("x", "y")

function _pn_name_tuple(values)
    values isa Tuple || throw(ArgumentError("expected tuple for name metadata"))
    return map(
        name -> (name isa AbstractString ? name : name isa Symbol ? String(name) : throw(ArgumentError("ring names must be strings or symbols"))),
        values,
    )
end

function _pn_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture entry missing field $(field)"))
    return getproperty(entry, field)
end

function _pn_field_optional(entry, fields::Vector{Symbol})
    for field in fields
        if hasproperty(entry, field)
            return getproperty(entry, field)
        end
    end
    return nothing
end

function _pn_as_vector(values, n::Int, R, label::AbstractString)
    if values isa AbstractVector
        length(values) == n ||
            throw(ArgumentError("$(label) must be a length-$(n) vector"))
        all(value -> parent(value) == R, values) ||
            throw(ArgumentError("$(label) must use the fixture ring"))
        return values
    elseif values isa AbstractAlgebra.MatElem
        nrows(values) == n && ncols(values) == 1 ||
            throw(ArgumentError("$(label) must be an n×1 vector matrix"))
        all(row -> parent(values[row, 1]) == R, 1:n) ||
            throw(ArgumentError("$(label) entries must use the fixture ring"))
        return [values[row, 1] for row in 1:n]
    else
        throw(ArgumentError("$(label) must be a vector or column matrix"))
    end
end

function _pn_target_matrix(entry, R)
    layer = _pn_field(entry, :section2_layer)
    inputs = _pn_field(entry, :inputs)
    if layer == :cohn_type
        i = _pn_field(inputs, :i)
        j = _pn_field(inputs, :j)
        a = R(_pn_field(inputs, :a))
        v = _pn_as_vector(_pn_field(inputs, :v), nrows(_pn_field(entry, :target_matrix)), R, "$(entry.id) input v")
        n = length(v)
        i == j && throw(ArgumentError("$(entry.id) must use distinct indices for cohn-type"))
        1 <= i <= n || throw(ArgumentError("$(entry.id) cohn index i out of range"))
        1 <= j <= n || throw(ArgumentError("$(entry.id) cohn index j out of range"))

        target = identity_matrix(R, n)
        vi, vj = v[i], v[j]
        for row in 1:n
            target[row, i] += a * v[row] * vj
            target[row, j] -= a * v[row] * vi
        end
        return target
    elseif layer == :orthogonal_rank_one
        v = _pn_as_vector(
            _pn_field(inputs, :v),
            nrows(_pn_field(entry, :target_matrix)),
            R,
            "$(entry.id) input v",
        )
        w = _pn_as_vector(
            _pn_field(inputs, :w),
            nrows(_pn_field(entry, :target_matrix)),
            R,
            "$(entry.id) input w",
        )
        n = length(v)
        target = identity_matrix(R, n)
        for row in 1:n, col in 1:n
            target[row, col] += v[row] * w[col]
        end
        return target
    elseif layer == :conjugated_elementary
        B = _pn_field(inputs, :B)
        i = _pn_field(inputs, :i)
        j = _pn_field(inputs, :j)
        a = R(_pn_field(inputs, :a))
        n = nrows(B)
        E = elementary_matrix(n, i, j, a, R)
        return B * E * inv(B)
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported section2 layer $(layer)"))
    end
end

function _pn_dot(values_a::AbstractVector, values_b::AbstractVector, R)
    length(values_a) == length(values_b) || throw(ArgumentError("dot product inputs have different lengths"))
    total = zero(R)
    for idx in eachindex(values_a, values_b)
        total += values_a[idx] * values_b[idx]
    end
    return total
end

function _pn_require_ordinary_field_polynomial_ring(entry)
    ring_meta = _pn_field(entry, :ring)
    ring_name = _pn_field(ring_meta, :description)
    ring_name isa AbstractString || throw(ArgumentError("fixture $(entry.id) ring description must be a string"))
    R = _pn_field(ring_meta, :object)
    R isa MPolyRing || throw(ArgumentError("fixture $(entry.id) must use an ordinary polynomial ring"))
    ring_ctor = _pn_field(entry, :ring_constructor)
    _pn_field(ring_ctor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use polynomial_ring metadata"))
    ring_coefficient = _pn_field(ring_ctor, :coefficient)
    ring_coefficient == EXPECTED_RING_COEFFICIENT ||
        throw(ArgumentError("fixture $(entry.id) must use $(EXPECTED_RING_COEFFICIENT) coefficient in metadata"))
    ring_vars = _pn_field(ring_ctor, :variables)
    ring_gen_names = _pn_field(ring_meta, :generator_names)
    ring_vars = _pn_name_tuple(ring_vars)
    ring_gen_names = _pn_name_tuple(ring_gen_names)
    ring_generators = _pn_field(ring_meta, :generators)
    ring_vars == _pn_name_tuple(EXPECTED_RING_VARIABLES) ||
        throw(ArgumentError("fixture $(entry.id) must use ring variables $(EXPECTED_RING_VARIABLES)"))
    length(ring_vars) == length(ring_gen_names) ||
        throw(ArgumentError("fixture $(entry.id) ring constructor and metadata generator names must have same length"))
    ring_vars == ring_gen_names ||
        throw(ArgumentError("fixture $(entry.id) ring metadata generator names do not match constructor metadata"))
    ring_generators isa Tuple || throw(ArgumentError("fixture $(entry.id) ring generators must be a tuple"))
    length(ring_generators) == length(ring_gen_names) ||
        throw(ArgumentError("fixture $(entry.id) ring metadata must match constructor arity"))
    all(parent(generator) == R for generator in ring_generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generators must use the ring object"))
    Oscar.is_exact_type(typeof(zero(coefficient_ring(R)))) ||
        throw(ArgumentError("fixture $(entry.id) coefficient ring must be exact"))
    coefficient_ring(R) isa Field ||
        throw(ArgumentError("fixture $(entry.id) must use field-backed coefficients"))
    Suslin._is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("fixture $(entry.id) must use ordinary (non-Laurent) polynomial ring"))
    return (R, ring_name)
end

function validate_polynomial_normality_fixture(entry)
    for field in REQUIRED_NORMALITY_CASE_FIELDS
        _pn_field(entry, field)
    end

    entry.id isa AbstractString && !isempty(entry.id) ||
        throw(ArgumentError("fixture id must be non-empty"))
    entry.section2_layer isa Symbol || throw(ArgumentError("fixture $(entry.id) section2_layer must be a symbol"))
    entry.section2_layer in SECTION2_LAYERS ||
        throw(ArgumentError("fixture $(entry.id) has unsupported section2_layer $(entry.section2_layer)"))
    entry.expected_convention isa AbstractString || entry.expected_convention isa Symbol ||
        throw(ArgumentError("fixture $(entry.id) expected_convention must be string or symbol"))
    if haskey(EXPECTED_SECTION2_LAYER_BY_ID, entry.id)
        EXPECTED_SECTION2_LAYER_BY_ID[entry.id] == entry.section2_layer ||
            throw(ArgumentError("fixture $(entry.id) has wrong section2 layer"))
    end

    source = _pn_field_optional(entry, [:source, :source_refs, :provenance])
    source === nothing && throw(ArgumentError("fixture $(entry.id) missing source metadata"))
    if source isa Tuple
        !isempty(source) || throw(ArgumentError("fixture $(entry.id) source metadata tuple must be non-empty"))
    elseif source isa AbstractString
        !isempty(source) || throw(ArgumentError("fixture $(entry.id) source metadata must be non-empty"))
    elseif source isa NamedTuple
        isempty(keys(source)) && throw(ArgumentError("fixture $(entry.id) source metadata must be non-empty"))
    else
        throw(ArgumentError("fixture $(entry.id) source metadata has unsupported format $(typeof(source))"))
    end

    entry.consumer_issue_ids isa Tuple || throw(ArgumentError("fixture $(entry.id) consumer_issue_ids must be a tuple"))
    all(id -> id isa AbstractString && startswith(id, "#"), entry.consumer_issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) consumer_issue_ids must look like issue references"))

    R, _ = _pn_require_ordinary_field_polynomial_ring(entry)

    target = _pn_field(entry, :target_matrix)
    target isa AbstractAlgebra.MatElem || throw(ArgumentError("fixture $(entry.id) target_matrix must be a matrix"))
    base_ring(target) == R || throw(ArgumentError("fixture $(entry.id) target_matrix must use ring object"))
    nrows(target) == ncols(target) || throw(ArgumentError("fixture $(entry.id) target must be square"))

    computed_target = try
        _pn_target_matrix(entry, R)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("fixture $(entry.id) cannot reconstruct target from stored inputs"))
    end
    computed_target == target || throw(ArgumentError("fixture $(entry.id) target_matrix does not match computed section 2 convention"))
    det(target) == one(R) || throw(ArgumentError("fixture $(entry.id) target_matrix determinant must be one"))

    inputs = _pn_field(entry, :inputs)
    if entry.section2_layer == :cohn_type
        v = _pn_as_vector(_pn_field(inputs, :v), nrows(target), R, "$(entry.id) v")
        length(v) == nrows(target) || throw(ArgumentError("fixture $(entry.id) Cohn layer must provide n-length vector"))
    elseif entry.section2_layer == :orthogonal_rank_one
        v = _pn_as_vector(_pn_field(inputs, :v), nrows(target), R, "$(entry.id) v")
        w = _pn_as_vector(_pn_field(inputs, :w), nrows(target), R, "$(entry.id) w")
        g = _pn_as_vector(_pn_field(inputs, :g), nrows(target), R, "$(entry.id) g")
        _pn_dot(w, v, R) == zero(R) ||
            throw(ArgumentError("fixture $(entry.id) must satisfy w*v = 0"))
        _pn_dot(g, v, R) == one(R) ||
            throw(ArgumentError("fixture $(entry.id) must satisfy g*v = 1"))
    elseif entry.section2_layer == :conjugated_elementary
        B = _pn_field(inputs, :B)
        nrows(B) == nrows(target) &&
            ncols(B) == nrows(target) ||
            throw(ArgumentError("fixture $(entry.id) B must be square and same size as target"))
        base_ring(B) == R || throw(ArgumentError("fixture $(entry.id) B must use fixture ring"))

        i = _pn_field(inputs, :i)
        j = _pn_field(inputs, :j)
        1 <= i <= nrows(target) || throw(ArgumentError("fixture $(entry.id) conjugated i out of range"))
        1 <= j <= nrows(target) || throw(ArgumentError("fixture $(entry.id) conjugated j out of range"))
        i != j || throw(ArgumentError("fixture $(entry.id) requires i != j"))

        Binv = inv(B)
        v = [B[row, i] for row in 1:nrows(B)]
        w = [R(_pn_field(inputs, :a)) * Binv[j, col] for col in 1:ncols(Binv)]
        g = [Binv[i, col] for col in 1:ncols(Binv)]
        _pn_dot(v, w, R) == zero(R) ||
            throw(ArgumentError("fixture $(entry.id) extracted conjugated cohn decomposition must be orthogonal"))
        _pn_dot(g, v, R) == one(R) ||
            throw(ArgumentError("fixture $(entry.id) extracted conjugated cohn decomposition must include unimodular g"))

        if hasproperty(inputs, :v) || hasproperty(inputs, :w) || hasproperty(inputs, :g)
            stored_v = _pn_as_vector(_pn_field(inputs, :v), nrows(target), R, "$(entry.id) stored v")
            stored_w = _pn_as_vector(_pn_field(inputs, :w), nrows(target), R, "$(entry.id) stored w")
            stored_g = _pn_as_vector(_pn_field(inputs, :g), nrows(target), R, "$(entry.id) stored g")
            stored_v == v || throw(ArgumentError("fixture $(entry.id) stored v does not match extracted conjugated v"))
            stored_w == w || throw(ArgumentError("fixture $(entry.id) stored w does not match extracted conjugated w"))
            stored_g == g || throw(ArgumentError("fixture $(entry.id) stored g does not match extracted conjugated g"))
        end
    else
        throw(ArgumentError("fixture $(entry.id) has unsupported section2 layer $(entry.section2_layer)"))
    end

    return true
end

function validate_polynomial_normality_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("fixture catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("fixture catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("fixture catalog must contain positive cases"))
    case_ids = [entry.id for entry in catalog.cases]
    length(case_ids) == length(unique(case_ids)) ||
        throw(ArgumentError("fixture catalog ids must be unique"))

    for entry in catalog.cases
        validate_polynomial_normality_fixture(entry)
    end

    isempty(catalog.negative_controls) && throw(ArgumentError("fixture catalog must contain negative controls"))
    control_ids = [entry.id for entry in catalog.negative_controls]
    length(control_ids) == length(unique(control_ids)) ||
        throw(ArgumentError("negative control ids must be unique"))
    for entry in catalog.negative_controls
        try
            validate_polynomial_normality_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end

    section2_counts = Dict{Symbol,Int}(layer => 0 for layer in SECTION2_LAYERS)
    for entry in catalog.cases
        section2_counts[entry.section2_layer] = section2_counts[entry.section2_layer] + 1
    end
    for layer in SECTION2_LAYERS
        layer_count = section2_counts[layer]
        layer_count == 1 ||
            throw(ArgumentError("fixture catalog must contain exactly one positive case with section2_layer $(layer), found $(layer_count)"))
    end

    return true
end

@testset "Park-Woodburn polynomial normality fixture catalog" begin
    if !isdefined(Main, :ParkWoodburnPolynomialNormalityFixtureCatalog)
        include(PARK_WOODBURN_POLYNOMIAL_NORMALITY_CATALOG_PATH)
    end
    catalog = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.catalog()

    @test validate_polynomial_normality_fixture_catalog(catalog)

    cases = Main.ParkWoodburnPolynomialNormalityFixtureCatalog.cases_by_id()
    @test REQUIRED_PW_POLYNOMIAL_NORMALITY_IDS ⊆ Set(keys(cases))

    negatives = if isdefined(Main.ParkWoodburnPolynomialNormalityFixtureCatalog, :negative_controls_by_id)
        Main.ParkWoodburnPolynomialNormalityFixtureCatalog.negative_controls_by_id()
    else
        Dict(entry.id => entry for entry in catalog.negative_controls)
    end
    @test REQUIRED_PW_POLYNOMIAL_NORMALITY_NEGATIVE_IDS ⊆ Set(keys(negatives))

    for id in REQUIRED_PW_POLYNOMIAL_NORMALITY_NEGATIVE_IDS
        @test haskey(negatives, id)
    end

    for entry in values(negatives)
        @test_throws ArgumentError validate_polynomial_normality_fixture(entry)
    end
end
