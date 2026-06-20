module ToricBuilderLaurentProblemCatalog

using Oscar
using Suslin

include("toricbuilder_issue38_cases.jl")
include("toricbuilder_factor_toric_block_3.jl")
include("laurent_large_acceptance_cases.jl")

_matrix_dimensions(A) = (nrows(A), ncols(A))

function _ring_metadata_from_matrix(A; description = nothing, generators = nothing)
    R = base_ring(A)
    return (;
        description = description === nothing ? string(R) : description,
        object = R,
        generators,
    )
end

function _determinant_profile(A)
    profile = Suslin.classify_laurent_determinant(A)
    return (;
        expected_class = profile.classification,
        expected_determinant = profile.determinant,
        monomial_exponents = hasproperty(profile, :monomial_exponents) ? profile.monomial_exponents : nothing,
        monomial_coefficient = hasproperty(profile, :monomial_coefficient) ? profile.monomial_coefficient : nothing,
    )
end

function _verifier(path::String, scenario::Symbol)
    return (;
        path,
        scenario,
    )
end

function _consumers(; milestone::Int, issues, tests)
    return (;
        milestone,
        issues = Tuple(issues),
        tests = Tuple(tests),
    )
end

function _issue38_problem()
    issue38 = only(ToricBuilderIssue38Cases.catalog().cases)
    return (;
        id = issue38.id,
        kind = :issue38_q_block,
        source_fixture = :toricbuilder_issue38_cases,
        ring = issue38.ring,
        dimensions = issue38.dimensions,
        matrix = issue38.inputs.matrix,
        determinant_profile = (;
            expected_class = issue38.determinant_profile.expected_class,
            expected_determinant = issue38.determinant_profile.expected_determinant,
            monomial_exponents = issue38.determinant_profile.monomial_exponents,
            monomial_coefficient = issue38.determinant_profile.monomial_coefficient,
        ),
        expected_current_status = :unsupported_now,
        expected_suslin_path = nothing,
        verifier = _verifier("internal/toricbuilder_issue38_fixture.jl", :issue38_fixture_validator),
        provenance = merge(issue38.provenance, (;
            fixture_id = issue38.id,
            source_fixture = "test/fixtures/toricbuilder_issue38_cases.jl",
        )),
        consumers = _consumers(
            milestone = 4,
            issues = ("#38", "#39", "#40"),
            tests = ("test/internal/toricbuilder_issue38_fixture.jl", "test/internal/toricbuilder_problem_catalog.jl"),
        ),
    )
end

function _contract_problem(entry)
    id_suffix = lowercase(entry.toricbuilder_role)
    return (;
        id = "toricbuilder-factor-toric-block-3-$(id_suffix)",
        kind = :toricbuilder_contract,
        source_fixture = :toricbuilder_factor_toric_block_3,
        ring = _ring_metadata_from_matrix(entry.matrix; description = entry.ring),
        dimensions = (; matrix = entry.size),
        matrix = entry.matrix,
        determinant_profile = (;
            expected_class = Symbol(entry.determinant_classification),
            expected_determinant = det(entry.matrix),
            monomial_exponents = nothing,
            monomial_coefficient = nothing,
        ),
        expected_current_status = :supported_column_peel,
        expected_suslin_path = entry.expected_suslin_path,
        verifier = _verifier("public/toricbuilder_factor_toric_block_acceptance.jl", :toricbuilder_factor_toric_block_column_peel_acceptance),
        provenance = merge(entry.provenance, (;
            fixture_id = entry.name,
            toricbuilder_role = entry.toricbuilder_role,
            source_fixture = "test/fixtures/toricbuilder_factor_toric_block_3.jl",
        )),
        consumers = _consumers(
            milestone = 4,
            issues = ("#19", "#40", "#58"),
            tests = ("test/public/toricbuilder_factor_toric_block_acceptance.jl", "test/internal/toricbuilder_problem_catalog.jl"),
        ),
    )
end

function _large_acceptance_problem(entry, status::Symbol)
    return (;
        id = entry.id,
        kind = :synthetic_block_local_acceptance,
        source_fixture = :laurent_large_acceptance_cases,
        ring = entry.ring,
        dimensions = (; matrix = entry.size),
        matrix = entry.matrix,
        determinant_profile = _determinant_profile(entry.matrix),
        expected_current_status = status,
        expected_suslin_path = nothing,
        verifier = _verifier("public/laurent_large_acceptance.jl", :large_laurent_acceptance),
        provenance = merge(entry.provenance, (;
            fixture_id = entry.id,
            source_fixture = "test/fixtures/laurent_large_acceptance_cases.jl",
        )),
        consumers = _consumers(
            milestone = 4,
            issues = ("#17", "#40"),
            tests = ("test/public/laurent_large_acceptance.jl", "test/internal/toricbuilder_problem_catalog.jl"),
        ),
    )
end

function catalog()
    contract_fixture = ToricBuilderFactorToricBlock3Fixture.fixture()
    contract_cases = Dict(entry.name => entry for entry in contract_fixture.cases)
    large_cases = Dict(entry.id => entry for entry in LaurentLargeAcceptanceCases.acceptance_catalog().cases)

    return (;
        cases = [
            _issue38_problem(),
            _contract_problem(contract_cases["factor_toric_block_3_qinv"]),
            _contract_problem(contract_cases["factor_toric_block_3_pinv"]),
            _large_acceptance_problem(large_cases["laurent-block-local-40x40"], :supported_block_local),
            _large_acceptance_problem(large_cases["laurent-block-local-48x48"], :target_acceptance),
        ],
    )
end

end
