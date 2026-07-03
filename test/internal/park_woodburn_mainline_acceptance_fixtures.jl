using Test
using Oscar
using Suslin

const PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "park_woodburn_mainline_acceptance_cases.jl")

const PW_MAINLINE_SECTION_REFS = Dict(
    :section2 => "refs/arXiv-alg-geom9405003v1 Section 2",
    :section3 => "refs/arXiv-alg-geom9405003v1 Section 3",
    :section4 => "refs/arXiv-alg-geom9405003v1 Section 4",
    :section5 => "refs/arXiv-alg-geom9405003v1 Section 5",
)

if !isdefined(Main, :ParkWoodburnMainlineAcceptanceFixtureCatalog)
    include(PARK_WOODBURN_MAINLINE_ACCEPTANCE_CATALOG_PATH)
end

const REQUIRED_PW_MAINLINE_CASE_IDS = Set([
    "pw-mainline-sl3-multivariate-issue184-qq",
    "pw-mainline-sln-recursive-issue185-186-gf2",
    "pw-mainline-readme-ordinary-polynomial-qq",
    "pw-mainline-staged-missing-evidence-qq",
])

const REQUIRED_PW_MAINLINE_NEGATIVE_IDS = Set([
    "pw-mainline-negative-det-not-one",
    "pw-mainline-negative-unsupported-coefficient-ring",
    "pw-mainline-negative-missing-evidence",
])

const REQUIRED_PW_MAINLINE_FIELDS = (
    :id,
    :entry_class,
    :expected_status,
    :public_route,
    :ring_constructor,
    :ring,
    :matrix,
    :determinant_metadata,
    :source_refs,
    :upstream_issue_ids,
    :upstream_evidence,
    :acceptance_metadata,
)

const ALLOWED_PW_MAINLINE_ENTRY_CLASSES = Set([
    :issue184_sl3_multivariate,
    :issue185_186_sln_recursive,
    :readme_public_example,
    :staged_missing_evidence_boundary,
])

const REQUIRED_PW_MAINLINE_CLASS_SECTIONS = Dict(
    :issue184_sl3_multivariate => (:section5,),
    :issue185_186_sln_recursive => (:section3, :section4, :section5),
    :readme_public_example => (:section2, :section5),
    :staged_missing_evidence_boundary => (:section2, :section4),
)

const REQUIRED_PW_MAINLINE_CLASS_ISSUES = Dict(
    :issue184_sl3_multivariate => Set(["#184"]),
    :issue185_186_sln_recursive => Set(["#184", "#185", "#186"]),
    :readme_public_example => Set(["#184"]),
    :staged_missing_evidence_boundary => Set(["#186"]),
)

const REQUIRED_PW_MAINLINE_CLASS_EVIDENCE = Dict(
    :issue184_sl3_multivariate => (:driver_case_id,),
    :issue185_186_sln_recursive => (:sln_case_id, :ecp_case_id, :final_sl3_case_id),
    :readme_public_example => (:driver_case_id,),
    :staged_missing_evidence_boundary => (:sln_case_id,),
)

function _pwma_label(entry)
    return hasproperty(entry, :id) ? string(getproperty(entry, :id)) : "<missing id>"
end

function _pwma_field(entry, field::Symbol)
    hasproperty(entry, field) || throw(ArgumentError("fixture $(_pwma_label(entry)) missing field $(field)"))
    return getproperty(entry, field)
end

function _pwma_required_string(value, label)
    value isa AbstractString && !isempty(value) ||
        throw(ArgumentError("$(label) must be a non-empty string"))
    return value
end

function _pwma_required_symbol(value, label)
    value isa Symbol || throw(ArgumentError("$(label) must be a symbol"))
    return value
end

function _pwma_required_tuple(value, label)
    value isa Tuple || throw(ArgumentError("$(label) must be a tuple"))
    return value
end

function _pwma_catalog_module()
    return getfield(Main, :ParkWoodburnMainlineAcceptanceFixtureCatalog)
end

function _pwma_nested_module(name::Symbol)
    catalog_module = _pwma_catalog_module()
    if isdefined(catalog_module, name)
        return getfield(catalog_module, name)
    end
    isdefined(Main, name) || throw(ArgumentError("validator requires upstream catalog module $(name)"))
    return getfield(Main, name)
end

function _pwma_sl3_cases_by_id()
    catalog_module = _pwma_nested_module(:ParkWoodburnSL3DriverFixtureCatalog)
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function _pwma_ecp_cases_by_id()
    catalog_module = _pwma_nested_module(:ECPMainlineFixtureCatalog)
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function _pwma_sln_cases_by_id()
    catalog_module = _pwma_nested_module(:ParkWoodburnSLnDriverFixtureCatalog)
    return Base.invokelatest(getfield(catalog_module, :cases_by_id))
end

function _pwma_assert_public_route(entry)
    route = _pwma_field(entry, :public_route)
    _pwma_field(route, :entrypoint) == :elementary_factorization ||
        throw(ArgumentError("fixture $(entry.id) public route must use elementary_factorization"))
    _pwma_field(route, :issue_id) == "#187" ||
        throw(ArgumentError("fixture $(entry.id) public route must record parent issue #187"))
    _pwma_required_symbol(_pwma_field(route, :route_marker), "fixture $(entry.id) public route marker")
    _pwma_field(route, :status) == _pwma_field(entry, :expected_status) ||
        throw(ArgumentError("fixture $(entry.id) public route status must match expected_status"))
    return true
end

function _pwma_assert_ring_metadata(entry)
    ring_constructor = _pwma_field(entry, :ring_constructor)
    _pwma_field(ring_constructor, :function_name) == :polynomial_ring ||
        throw(ArgumentError("fixture $(entry.id) must use ordinary polynomial_ring metadata"))
    _pwma_required_string(_pwma_field(ring_constructor, :coefficient), "fixture $(entry.id) ring coefficient")
    variables = _pwma_required_tuple(_pwma_field(ring_constructor, :variables), "fixture $(entry.id) ring variables")
    all(name -> name isa AbstractString, variables) ||
        throw(ArgumentError("fixture $(entry.id) ring variables must be strings"))

    ring = _pwma_field(entry, :ring)
    _pwma_required_string(_pwma_field(ring, :description), "fixture $(entry.id) ring description")
    generator_names = _pwma_required_tuple(_pwma_field(ring, :generator_names), "fixture $(entry.id) ring generator_names")
    generators = _pwma_required_tuple(_pwma_field(ring, :generators), "fixture $(entry.id) ring generators")
    generator_names == variables ||
        throw(ArgumentError("fixture $(entry.id) ring constructor variables must match generator metadata"))
    length(generator_names) == length(generators) ||
        throw(ArgumentError("fixture $(entry.id) ring generator metadata is inconsistent"))
    all(name -> name isa AbstractString, generator_names) ||
        throw(ArgumentError("fixture $(entry.id) ring generator names must be strings"))

    R = _pwma_field(ring, :object)
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

function _pwma_assert_matrix_and_determinant(entry, R)
    matrix_value = _pwma_field(entry, :matrix)
    matrix_value isa AbstractAlgebra.MatElem ||
        throw(ArgumentError("fixture $(entry.id) matrix must be an Oscar matrix"))
    nrows(matrix_value) == ncols(matrix_value) ||
        throw(ArgumentError("fixture $(entry.id) matrix must be square"))
    nrows(matrix_value) >= 3 ||
        throw(ArgumentError("fixture $(entry.id) matrix must have rank at least 3"))
    base_ring(matrix_value) == R ||
        throw(ArgumentError("fixture $(entry.id) matrix must be defined over the fixture ring"))

    determinant_metadata = _pwma_field(entry, :determinant_metadata)
    _pwma_field(determinant_metadata, :expected) == :one ||
        throw(ArgumentError("fixture $(entry.id) determinant metadata must expect one"))
    _pwma_field(determinant_metadata, :certified) === true ||
        throw(ArgumentError("fixture $(entry.id) determinant metadata must be certified"))
    _pwma_field(determinant_metadata, :unit) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) determinant metadata unit must match the fixture ring"))
    _pwma_field(determinant_metadata, :value) == det(matrix_value) ||
        throw(ArgumentError("fixture $(entry.id) determinant metadata value must replay"))
    det(matrix_value) == one(R) ||
        throw(ArgumentError("fixture $(entry.id) matrix determinant must be one"))
    return true
end

function _pwma_assert_source_refs(entry)
    source_refs = _pwma_required_tuple(_pwma_field(entry, :source_refs), "fixture $(entry.id) source_refs")
    !isempty(source_refs) || throw(ArgumentError("fixture $(entry.id) source_refs must not be empty"))
    all(ref -> ref isa AbstractString && !isempty(ref), source_refs) ||
        throw(ArgumentError("fixture $(entry.id) source_refs must be non-empty strings"))
    entry_class = _pwma_field(entry, :entry_class)
    for section in REQUIRED_PW_MAINLINE_CLASS_SECTIONS[entry_class]
        PW_MAINLINE_SECTION_REFS[section] in source_refs ||
            throw(ArgumentError("fixture $(entry.id) must include Park-Woodburn $(section) source ref"))
    end
    return true
end

function _pwma_assert_upstream_issue_ids(entry)
    issue_ids = _pwma_required_tuple(_pwma_field(entry, :upstream_issue_ids), "fixture $(entry.id) upstream_issue_ids")
    all(id -> id isa AbstractString && occursin(r"^#[0-9]+$", id), issue_ids) ||
        throw(ArgumentError("fixture $(entry.id) upstream issue ids must be #number references"))
    required = REQUIRED_PW_MAINLINE_CLASS_ISSUES[_pwma_field(entry, :entry_class)]
    issubset(required, Set(issue_ids)) ||
        throw(ArgumentError("fixture $(entry.id) missing required upstream issue ids"))
    return true
end

function _pwma_assert_evidence_field(entry, field::Symbol)
    evidence = _pwma_field(entry, :upstream_evidence)
    hasproperty(evidence, field) ||
        throw(ArgumentError("fixture $(entry.id) upstream_evidence missing $(field)"))
    value = getproperty(evidence, field)
    value isa AbstractString && !isempty(value) ||
        throw(ArgumentError("fixture $(entry.id) upstream evidence $(field) must be a non-empty string"))
    return value
end

function _pwma_assert_upstream_evidence(entry)
    entry_class = _pwma_field(entry, :entry_class)
    for field in REQUIRED_PW_MAINLINE_CLASS_EVIDENCE[entry_class]
        _pwma_assert_evidence_field(entry, field)
    end

    if entry_class in (:issue184_sl3_multivariate, :readme_public_example)
        driver_case_id = _pwma_assert_evidence_field(entry, :driver_case_id)
        haskey(_pwma_sl3_cases_by_id(), driver_case_id) ||
            throw(ArgumentError("fixture $(entry.id) upstream SL_3 driver case id $(driver_case_id) does not exist"))
    elseif entry_class == :issue185_186_sln_recursive
        sln_case_id = _pwma_assert_evidence_field(entry, :sln_case_id)
        ecp_case_id = _pwma_assert_evidence_field(entry, :ecp_case_id)
        final_sl3_case_id = _pwma_assert_evidence_field(entry, :final_sl3_case_id)
        haskey(_pwma_sln_cases_by_id(), sln_case_id) ||
            throw(ArgumentError("fixture $(entry.id) upstream SL_n case id $(sln_case_id) does not exist"))
        haskey(_pwma_ecp_cases_by_id(), ecp_case_id) ||
            throw(ArgumentError("fixture $(entry.id) upstream ECP case id $(ecp_case_id) does not exist"))
        haskey(_pwma_sl3_cases_by_id(), final_sl3_case_id) ||
            throw(ArgumentError("fixture $(entry.id) upstream final SL_3 case id $(final_sl3_case_id) does not exist"))
    elseif entry_class == :staged_missing_evidence_boundary
        sln_case_id = _pwma_assert_evidence_field(entry, :sln_case_id)
        haskey(_pwma_sln_cases_by_id(), sln_case_id) ||
            throw(ArgumentError("fixture $(entry.id) upstream staged SL_n case id $(sln_case_id) does not exist"))
    end
    return true
end

function _pwma_assert_acceptance_metadata(entry)
    metadata = _pwma_field(entry, :acceptance_metadata)
    _pwma_required_symbol(_pwma_field(metadata, :boundary), "fixture $(entry.id) acceptance boundary")
    _pwma_required_string(_pwma_field(metadata, :note), "fixture $(entry.id) acceptance note")
    coverage = _pwma_field(metadata, :parent_issue_coverage)
    for field in (:issue181, :issue182, :issue183, :issue184, :issue185, :issue186)
        _pwma_required_symbol(_pwma_field(coverage, field), "fixture $(entry.id) parent coverage $(field)")
    end

    if _pwma_field(entry, :expected_status) == :mainline_accepted
        metadata.boundary != :staged_missing_evidence ||
            throw(ArgumentError("fixture $(entry.id) accepted entries must not use the staged boundary"))
        for issue_id in _pwma_field(entry, :upstream_issue_ids)
            field = Symbol("issue", issue_id[2:end])
            hasproperty(coverage, field) && getproperty(coverage, field) == :covered ||
                throw(ArgumentError("fixture $(entry.id) accepted entries must mark upstream $(issue_id) as covered"))
        end
    end
    return true
end

function _pwma_assert_status_boundary(entry)
    expected_status = _pwma_required_symbol(_pwma_field(entry, :expected_status), "fixture $(entry.id) expected_status")
    expected_status in Set([:mainline_accepted, :staged]) ||
        throw(ArgumentError("fixture $(entry.id) expected_status must be :mainline_accepted or :staged"))

    if expected_status == :mainline_accepted
        hasproperty(entry, :missing_evidence) &&
            throw(ArgumentError("fixture $(entry.id) accepted entries must not carry missing_evidence"))
        hasproperty(entry, :staged_reason) &&
            throw(ArgumentError("fixture $(entry.id) accepted entries must not carry staged_reason"))
    elseif expected_status == :staged
        missing_evidence = _pwma_required_tuple(_pwma_field(entry, :missing_evidence), "fixture $(entry.id) missing_evidence")
        !isempty(missing_evidence) ||
            throw(ArgumentError("fixture $(entry.id) staged missing_evidence must not be empty"))
        all(item -> item isa Symbol, missing_evidence) ||
            throw(ArgumentError("fixture $(entry.id) staged missing_evidence entries must be symbols"))
        _pwma_required_string(_pwma_field(entry, :staged_reason), "fixture $(entry.id) staged_reason")
        :final_sl3_case_id in missing_evidence ||
            throw(ArgumentError("fixture $(entry.id) staged boundary must name missing final SL_3 evidence"))
    end
    return true
end

function validate_park_woodburn_mainline_acceptance_fixture(entry)
    for field in REQUIRED_PW_MAINLINE_FIELDS
        _pwma_field(entry, field)
    end
    _pwma_required_string(_pwma_field(entry, :id), "fixture id")
    entry_class = _pwma_required_symbol(_pwma_field(entry, :entry_class), "fixture $(entry.id) entry_class")
    entry_class in ALLOWED_PW_MAINLINE_ENTRY_CLASSES ||
        throw(ArgumentError("fixture $(entry.id) entry_class must be recognized"))

    _pwma_assert_status_boundary(entry)
    _pwma_assert_public_route(entry)
    R = _pwma_assert_ring_metadata(entry)
    _pwma_assert_matrix_and_determinant(entry, R)
    _pwma_assert_source_refs(entry)
    _pwma_assert_upstream_issue_ids(entry)
    _pwma_assert_upstream_evidence(entry)
    _pwma_assert_acceptance_metadata(entry)
    return true
end

function validate_park_woodburn_mainline_acceptance_fixture_catalog(catalog)
    hasproperty(catalog, :cases) || throw(ArgumentError("catalog missing cases"))
    hasproperty(catalog, :negative_controls) || throw(ArgumentError("catalog missing negative_controls"))
    isempty(catalog.cases) && throw(ArgumentError("catalog must contain valid cases"))
    isempty(catalog.negative_controls) && throw(ArgumentError("catalog must contain negative controls"))

    case_ids = [entry.id for entry in catalog.cases]
    control_ids = [entry.id for entry in catalog.negative_controls]
    all_ids = vcat(case_ids, control_ids)
    length(all_ids) == length(unique(all_ids)) ||
        throw(ArgumentError("catalog case and negative control ids must be unique"))
    issubset(REQUIRED_PW_MAINLINE_CASE_IDS, Set(case_ids)) ||
        throw(ArgumentError("catalog missing required mainline acceptance case ids"))
    issubset(REQUIRED_PW_MAINLINE_NEGATIVE_IDS, Set(control_ids)) ||
        throw(ArgumentError("catalog missing required mainline acceptance negative-control ids"))

    catalog_refs = Set{String}()
    for entry in catalog.cases
        validate_park_woodburn_mainline_acceptance_fixture(entry)
        union!(catalog_refs, Set(entry.source_refs))
    end
    for ref in values(PW_MAINLINE_SECTION_REFS)
        ref in catalog_refs ||
            throw(ArgumentError("catalog must cover Park-Woodburn source ref $(ref)"))
    end

    valid_case_ids = Set(case_ids)
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
            validate_park_woodburn_mainline_acceptance_fixture(entry)
        catch err
            err isa ArgumentError || rethrow()
            continue
        end
        throw(ArgumentError("negative control $(entry.id) unexpectedly validated"))
    end
    return true
end

@testset "Park-Woodburn mainline acceptance fixture catalog" begin
    catalog_module = _pwma_catalog_module()
    catalog = Base.invokelatest(getfield(catalog_module, :catalog))
    @test validate_park_woodburn_mainline_acceptance_fixture_catalog(catalog)

    entries = Base.invokelatest(getfield(catalog_module, :cases_by_id))
    @test Set(keys(entries)) == Set(entry.id for entry in catalog.cases)
    @test issubset(REQUIRED_PW_MAINLINE_CASE_IDS, Set(keys(entries)))
    @test issubset(REQUIRED_PW_MAINLINE_NEGATIVE_IDS, Set(entry.id for entry in catalog.negative_controls))

    for entry in values(entries)
        @test validate_park_woodburn_mainline_acceptance_fixture(entry)
    end
    for entry in catalog.negative_controls
        @test_throws ArgumentError validate_park_woodburn_mainline_acceptance_fixture(entry)
    end

    @testset "validator mutation gates" begin
        sl3_entry = entries["pw-mainline-sl3-multivariate-issue184-qq"]
        missing_section_ref = merge(sl3_entry, (;
            source_refs = Tuple(ref for ref in sl3_entry.source_refs if ref != PW_MAINLINE_SECTION_REFS[:section5]),
        ))
        @test_throws ArgumentError validate_park_woodburn_mainline_acceptance_fixture(missing_section_ref)

        sln_entry = entries["pw-mainline-sln-recursive-issue185-186-gf2"]
        missing_final_evidence = merge(sln_entry, (;
            upstream_evidence = (;
                sln_case_id = sln_entry.upstream_evidence.sln_case_id,
                ecp_case_id = sln_entry.upstream_evidence.ecp_case_id,
            ),
        ))
        @test_throws ArgumentError validate_park_woodburn_mainline_acceptance_fixture(missing_final_evidence)

        staged_entry = entries["pw-mainline-staged-missing-evidence-qq"]
        staged_without_reason = NamedTuple{Tuple(k for k in keys(staged_entry) if k != :staged_reason)}(
            Tuple(v for (k, v) in pairs(staged_entry) if k != :staged_reason),
        )
        @test_throws ArgumentError validate_park_woodburn_mainline_acceptance_fixture(staged_without_reason)

        duplicate_catalog = merge(catalog, (; cases = [catalog.cases; catalog.cases[1]]))
        @test_throws ArgumentError validate_park_woodburn_mainline_acceptance_fixture_catalog(duplicate_catalog)
    end
end
