# Issue 217 Park-Woodburn Substitution Chain Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an internal replayable Park-Woodburn substitution-chain certificate for Quillen induction denominator covers.

**Architecture:** Extend `src/algorithm/quillen_induction.jl` beside the #216 denominator-cover solver with a scaled-variable matrix substitution helper, chain step records, chain replay verification, and a constructor that consumes a verified solver result. Focused expert tests build a two-open chain and tamper one field at a time to prove the verifier rejects corrupted records.

**Tech Stack:** Julia, Oscar ordinary polynomial rings and exact matrix arithmetic, existing Suslin Quillen solver records, `Test`.

## Global Constraints

- Input is an ordinary-polynomial matrix `A`, selected variable `X`, and a verified `QuillenDenominatorCoverSolverResult`.
- Output is a `QuillenPatchSubstitutionChain` whose replay proves final substitution `X -> 0`, base term `A(0)`, and bracket telescope `A(X) * bracket_1 * ... * bracket_n == A(0)`.
- The chain supports substitutions `X -> c*X` and `X -> (c + r_i^l*g_i)*X` through exact scaled-variable matrix evaluation, but the Park-Woodburn chain constructor uses the explicit sign convention `:park_woodburn_minus`.
- Record each cumulative coefficient, intermediate substituted matrix, bracket target, selected variable, exponent, multiplier, sign convention, and replay metadata.
- Do not force the new chain through the older `patched_substitution(A, X, r, l, g)` `X + r^l*g` interface.
- Do not solve local realizability, normalize local factors, factor `A(0)`, assemble global factors, or claim public `elementary_factorization` support.
- Keep new names expert/internal and do not export them from `src/Suslin.jl`.
- Focused chain command is `julia --project=. -e 'include("test/expert/quillen_patch_substitution_chain.jl")'`.
- Full package command is `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## File Structure

- Modify `src/algorithm/quillen_induction.jl`: add the scaled-variable helper, step and chain records, replay verification, and constructor.
- Create `test/expert/quillen_patch_substitution_chain.jl`: positive two-open chain replay plus required negative controls.
- Modify `test/runtests.jl`: register the new expert test after `expert/quillen_denominator_cover_solver.jl`.

### Task 1: Add Red Chain Tests

**Files:**
- Create: `test/expert/quillen_patch_substitution_chain.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes future `Suslin.QuillenPatchSubstitutionChain`, `Suslin.quillen_patch_substitution_chain`, `Suslin.replay_quillen_patch_substitution_chain`, `Suslin.verify_quillen_patch_substitution_chain`, and `Suslin._quillen_substitute_matrix_scaled_variable`.
- Produces failing behavioral coverage for the Park-Woodburn chain implementation.

- [ ] **Step 1: Write the failing expert test**

Create `test/expert/quillen_patch_substitution_chain.jl`:

```julia
using Test
using Suslin
using Oscar

const QUILLEN_CHAIN_CATALOG_PATH =
    joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_CHAIN_CATALOG_PATH)
end

function chain_product(factors, R, n)
    product = identity_matrix(R, n)
    for factor in factors
        product *= factor
    end
    return product
end

function chain_rebuild_step(step; kwargs...)
    fields = merge((
        step_index = step.step_index,
        selected_variable = step.selected_variable,
        raw_denominator = step.raw_denominator,
        exponent = step.exponent,
        powered_denominator = step.powered_denominator,
        coverage_multiplier = step.coverage_multiplier,
        sign_convention = step.sign_convention,
        previous_coefficient = step.previous_coefficient,
        next_coefficient = step.next_coefficient,
        previous_matrix = step.previous_matrix,
        next_matrix = step.next_matrix,
        bracket_target = step.bracket_target,
        replay_metadata = step.replay_metadata,
    ), kwargs)
    return Suslin.QuillenPatchSubstitutionStep(
        fields.step_index,
        fields.selected_variable,
        fields.raw_denominator,
        fields.exponent,
        fields.powered_denominator,
        fields.coverage_multiplier,
        fields.sign_convention,
        fields.previous_coefficient,
        fields.next_coefficient,
        fields.previous_matrix,
        fields.next_matrix,
        fields.bracket_target,
        fields.replay_metadata,
    )
end

function chain_rebuild(chain; kwargs...)
    fields = merge((
        original_matrix = chain.original_matrix,
        ring = chain.ring,
        size = chain.size,
        selected_variable = chain.selected_variable,
        sign_convention = chain.sign_convention,
        solver_result = chain.solver_result,
        cumulative_coefficients = chain.cumulative_coefficients,
        intermediate_matrices = chain.intermediate_matrices,
        steps = chain.steps,
        bracket_matrices = chain.bracket_matrices,
        base_term = chain.base_term,
        metadata = chain.metadata,
        replay_metadata = chain.replay_metadata,
        verification = chain.verification,
    ), kwargs)
    return Suslin.QuillenPatchSubstitutionChain(
        fields.original_matrix,
        fields.ring,
        fields.size,
        fields.selected_variable,
        fields.sign_convention,
        fields.solver_result,
        fields.cumulative_coefficients,
        fields.intermediate_matrices,
        fields.steps,
        fields.bracket_matrices,
        fields.base_term,
        fields.metadata,
        fields.replay_metadata,
        fields.verification,
    )
end

@testset "Park-Woodburn substitution chain" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    entry = entries["quillen-two-open-cover-qq"]
    R = entry.ring.object
    X = entry.substitution_variable
    raw = [data.denominator for data in entry.denominator_data]
    solver = Suslin.solve_quillen_denominator_cover(R, raw; max_exponent = 2)

    chain = Suslin.quillen_patch_substitution_chain(
        entry.base_matrix,
        X,
        solver;
        metadata = (; fixture_id = entry.id, consumer_issue_id = "#217"),
    )

    @test chain isa Suslin.QuillenPatchSubstitutionChain
    @test Suslin.verify_quillen_denominator_cover_solver_result(solver)
    @test Suslin.verify_quillen_patch_substitution_chain(chain)
    replay = Suslin.replay_quillen_patch_substitution_chain(chain)
    @test replay.overall_ok
    @test replay.final_coefficient_ok
    @test replay.base_term_ok
    @test replay.telescope_ok

    r = raw[1]
    expected_coefficients = [one(R), one(R) - r, zero(R)]
    @test chain.cumulative_coefficients == expected_coefficients
    @test chain.sign_convention == :park_woodburn_minus
    @test chain.selected_variable == X
    @test chain.solver_result === solver
    @test chain.metadata == (; fixture_id = entry.id, consumer_issue_id = "#217")
    @test chain.replay_metadata.sign_convention == :park_woodburn_minus
    @test chain.replay_metadata.denominator_count == length(raw)
    @test chain.replay_metadata.coverage_sum == one(R)

    expected_matrices = [
        Suslin._quillen_substitute_matrix_scaled_variable(entry.base_matrix, X, coefficient)
        for coefficient in expected_coefficients
    ]
    @test chain.intermediate_matrices == expected_matrices
    @test chain.base_term ==
          Suslin._quillen_substitute_matrix_scaled_variable(entry.base_matrix, X, zero(R))
    @test chain.base_term == last(chain.intermediate_matrices)

    @test length(chain.steps) == length(raw)
    @test length(chain.bracket_matrices) == length(raw)
    for (idx, step) in enumerate(chain.steps)
        @test step.step_index == idx
        @test step.selected_variable == X
        @test step.raw_denominator == solver.raw_denominators[idx]
        @test step.exponent == solver.exponent
        @test step.powered_denominator == solver.powered_denominators[idx]
        @test step.coverage_multiplier == solver.coverage_multipliers[idx]
        @test step.sign_convention == :park_woodburn_minus
        @test step.previous_coefficient == chain.cumulative_coefficients[idx]
        @test step.next_coefficient == chain.cumulative_coefficients[idx + 1]
        @test step.previous_matrix == chain.intermediate_matrices[idx]
        @test step.next_matrix == chain.intermediate_matrices[idx + 1]
        @test step.bracket_target == chain.bracket_matrices[idx]
        @test step.previous_matrix * step.bracket_target == step.next_matrix
        @test step.replay_metadata.step_index == idx
        @test step.replay_metadata.exponent == solver.exponent
        @test step.replay_metadata.sign_convention == :park_woodburn_minus
    end

    telescope = entry.base_matrix
    for bracket in chain.bracket_matrices
        telescope *= bracket
    end
    @test telescope == chain.base_term
    @test chain_product(chain.bracket_matrices, R, chain.size) ==
          inv(entry.base_matrix) * chain.base_term

    corrupted_exponent_steps = copy(chain.steps)
    corrupted_exponent_steps[1] = chain_rebuild_step(
        chain.steps[1];
        exponent = chain.steps[1].exponent + 1,
    )
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; steps = corrupted_exponent_steps),
    )

    corrupted_multiplier_steps = copy(chain.steps)
    corrupted_multiplier_steps[1] = chain_rebuild_step(
        chain.steps[1];
        coverage_multiplier = chain.steps[1].coverage_multiplier + one(R),
    )
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; steps = corrupted_multiplier_steps),
    )

    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; sign_convention = :park_woodburn_plus),
    )

    corrupted_intermediates = copy(chain.intermediate_matrices)
    corrupted_intermediates[2] =
        corrupted_intermediates[2] * elementary_matrix(chain.size, 1, 2, one(R), R)
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; intermediate_matrices = corrupted_intermediates),
    )

    wrong_variable = collect(gens(R))[2]
    @test wrong_variable != X
    @test !Suslin.verify_quillen_patch_substitution_chain(
        chain_rebuild(chain; selected_variable = wrong_variable),
    )
end
```

- [ ] **Step 2: Register the red expert test**

In `test/runtests.jl`, add:

```julia
        "expert/quillen_patch_substitution_chain.jl",
```

immediately after:

```julia
        "expert/quillen_denominator_cover_solver.jl",
```

- [ ] **Step 3: Run the focused red test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_substitution_chain.jl")'
```

Expected: FAIL with `UndefVarError` for `quillen_patch_substitution_chain` or `QuillenPatchSubstitutionChain`.

- [ ] **Step 4: Commit the red tests**

Run:

```bash
git add test/expert/quillen_patch_substitution_chain.jl test/runtests.jl
git commit -m "test: cover park woodburn substitution chain"
```

### Task 2: Implement Replayable Chain Certificate

**Files:**
- Modify: `src/algorithm/quillen_induction.jl`

**Interfaces:**
- Consumes: `_require_quillen_denominator_cover_ring`, `_require_substitution_generator`, `_require_square_matrix`, `_coerce_into_ring`, and `verify_quillen_denominator_cover_solver_result`.
- Produces: `_quillen_substitute_matrix_scaled_variable`, `QuillenPatchSubstitutionStep`, `QuillenPatchSubstitutionChainVerification`, `QuillenPatchSubstitutionChain`, `quillen_patch_substitution_chain`, `replay_quillen_patch_substitution_chain`, and `verify_quillen_patch_substitution_chain`.

- [ ] **Step 1: Add substitution helper and records**

In `src/algorithm/quillen_induction.jl`, after `patched_substitution`, add:

```julia
function _quillen_substitute_matrix_scaled_variable(A, X, coefficient)
    R = base_ring(A)
    selected = _require_substitution_generator(R, X)
    scaled_coefficient = _coerce_into_ring(R, coefficient, "substitution coefficient")
    ring_gens = collect(gens(R))
    variable_idx = findfirst(gen -> gen == selected, ring_gens)
    variable_idx === nothing &&
        throw(ArgumentError("X must be a generator of the matrix base ring"))

    values = copy(ring_gens)
    values[variable_idx] = scaled_coefficient * ring_gens[variable_idx]
    entries = [
        _coerce_into_ring(R, evaluate(A[row, col], values), "scaled substitution entry")
        for col in 1:ncols(A), row in 1:nrows(A)
    ]
    return matrix(R, nrows(A), ncols(A), vec(entries))
end

struct QuillenPatchSubstitutionStep
    step_index::Int
    selected_variable
    raw_denominator
    exponent::Int
    powered_denominator
    coverage_multiplier
    sign_convention::Symbol
    previous_coefficient
    next_coefficient
    previous_matrix
    next_matrix
    bracket_target
    replay_metadata
end

struct QuillenPatchSubstitutionChainVerification
    solver_result_ok::Bool
    ring_ok::Bool
    matrix_ok::Bool
    selected_variable_ok::Bool
    sign_convention_ok::Bool
    coefficient_count::Int
    step_count::Int
    bracket_count::Int
    cumulative_coefficients::Vector
    cumulative_coefficients_ok::Bool
    intermediate_matrices::Vector
    intermediate_matrices_ok::Bool
    expected_steps::Vector{QuillenPatchSubstitutionStep}
    steps_ok::Bool
    bracket_matrices::Vector
    bracket_matrices_ok::Bool
    final_coefficient
    final_coefficient_ok::Bool
    base_term
    base_term_ok::Bool
    telescope_product
    telescope_ok::Bool
    replay_metadata
    replay_metadata_ok::Bool
    overall_ok::Bool
end

struct QuillenPatchSubstitutionChain
    original_matrix
    ring
    size::Int
    selected_variable
    sign_convention::Symbol
    solver_result::QuillenDenominatorCoverSolverResult
    cumulative_coefficients::Vector
    intermediate_matrices::Vector
    steps::Vector{QuillenPatchSubstitutionStep}
    bracket_matrices::Vector
    base_term
    metadata
    replay_metadata
    verification::QuillenPatchSubstitutionChainVerification
end
```

- [ ] **Step 2: Add expected-chain helpers**

After the #216 solver functions and before `extract_quillen_denominator_cover_candidate`, add:

```julia
function _quillen_patch_substitution_step_metadata(
    solver_result::QuillenDenominatorCoverSolverResult,
    selected_variable,
    sign_convention::Symbol,
    step_index::Int,
)
    return (;
        source = :park_woodburn_substitution_chain,
        step_index = step_index,
        selected_variable = selected_variable,
        raw_denominator = solver_result.raw_denominators[step_index],
        exponent = solver_result.exponent,
        powered_denominator = solver_result.powered_denominators[step_index],
        coverage_multiplier = solver_result.coverage_multipliers[step_index],
        coverage_term = solver_result.coverage_terms[step_index],
        sign_convention = sign_convention,
    )
end

function _quillen_patch_substitution_chain_metadata(
    solver_result::QuillenDenominatorCoverSolverResult,
    selected_variable,
    sign_convention::Symbol,
    metadata,
)
    return (;
        source = :park_woodburn_substitution_chain,
        selected_variable = selected_variable,
        exponent = solver_result.exponent,
        denominator_count = length(solver_result.raw_denominators),
        coverage_sum = solver_result.coverage_sum,
        sign_convention = sign_convention,
        metadata = metadata,
    )
end

function _quillen_patch_next_coefficient(
    previous,
    powered_denominator,
    coverage_multiplier,
    sign_convention::Symbol,
)
    sign_convention == :park_woodburn_minus ||
        throw(ArgumentError("unsupported Park-Woodburn substitution sign convention"))
    return previous - coverage_multiplier * powered_denominator
end

function _same_quillen_patch_substitution_step(
    left::QuillenPatchSubstitutionStep,
    right::QuillenPatchSubstitutionStep,
)::Bool
    return left.step_index == right.step_index &&
           left.selected_variable == right.selected_variable &&
           left.raw_denominator == right.raw_denominator &&
           left.exponent == right.exponent &&
           left.powered_denominator == right.powered_denominator &&
           left.coverage_multiplier == right.coverage_multiplier &&
           left.sign_convention == right.sign_convention &&
           left.previous_coefficient == right.previous_coefficient &&
           left.next_coefficient == right.next_coefficient &&
           left.previous_matrix == right.previous_matrix &&
           left.next_matrix == right.next_matrix &&
           left.bracket_target == right.bracket_target &&
           left.replay_metadata == right.replay_metadata
end

function _same_quillen_patch_substitution_steps(left, right)::Bool
    length(left) == length(right) || return false
    return all(
        _same_quillen_patch_substitution_step(left[idx], right[idx])
        for idx in eachindex(left)
    )
end

function _quillen_patch_expected_substitution_chain(
    A,
    selected_variable,
    solver_result::QuillenDenominatorCoverSolverResult,
    sign_convention::Symbol,
    metadata,
)
    R = base_ring(A)
    cumulative_coefficients = Any[one(R)]
    intermediate_matrices = Any[
        _quillen_substitute_matrix_scaled_variable(A, selected_variable, one(R)),
    ]
    steps = QuillenPatchSubstitutionStep[]
    bracket_matrices = Any[]

    for idx in eachindex(solver_result.raw_denominators)
        previous_coefficient = last(cumulative_coefficients)
        next_coefficient = _quillen_patch_next_coefficient(
            previous_coefficient,
            solver_result.powered_denominators[idx],
            solver_result.coverage_multipliers[idx],
            sign_convention,
        )
        previous_matrix = last(intermediate_matrices)
        next_matrix = _quillen_substitute_matrix_scaled_variable(
            A,
            selected_variable,
            next_coefficient,
        )
        bracket_target = inv(previous_matrix) * next_matrix
        push!(cumulative_coefficients, next_coefficient)
        push!(intermediate_matrices, next_matrix)
        push!(bracket_matrices, bracket_target)
        push!(
            steps,
            QuillenPatchSubstitutionStep(
                idx,
                selected_variable,
                solver_result.raw_denominators[idx],
                solver_result.exponent,
                solver_result.powered_denominators[idx],
                solver_result.coverage_multipliers[idx],
                sign_convention,
                previous_coefficient,
                next_coefficient,
                previous_matrix,
                next_matrix,
                bracket_target,
                _quillen_patch_substitution_step_metadata(
                    solver_result,
                    selected_variable,
                    sign_convention,
                    idx,
                ),
            ),
        )
    end

    base_term = _quillen_substitute_matrix_scaled_variable(A, selected_variable, zero(R))
    telescope_product = A
    for bracket in bracket_matrices
        telescope_product *= bracket
    end

    return (;
        cumulative_coefficients = cumulative_coefficients,
        intermediate_matrices = intermediate_matrices,
        steps = steps,
        bracket_matrices = bracket_matrices,
        final_coefficient = last(cumulative_coefficients),
        base_term = base_term,
        telescope_product = telescope_product,
        replay_metadata = _quillen_patch_substitution_chain_metadata(
            solver_result,
            selected_variable,
            sign_convention,
            metadata,
        ),
    )
end
```

- [ ] **Step 3: Add replay, verify, and constructor**

Continue in `src/algorithm/quillen_induction.jl` with:

```julia
function _same_quillen_patch_substitution_chain_verification(
    left::QuillenPatchSubstitutionChainVerification,
    right::QuillenPatchSubstitutionChainVerification,
)::Bool
    return left.solver_result_ok == right.solver_result_ok &&
           left.ring_ok == right.ring_ok &&
           left.matrix_ok == right.matrix_ok &&
           left.selected_variable_ok == right.selected_variable_ok &&
           left.sign_convention_ok == right.sign_convention_ok &&
           left.coefficient_count == right.coefficient_count &&
           left.step_count == right.step_count &&
           left.bracket_count == right.bracket_count &&
           left.cumulative_coefficients == right.cumulative_coefficients &&
           left.cumulative_coefficients_ok == right.cumulative_coefficients_ok &&
           left.intermediate_matrices == right.intermediate_matrices &&
           left.intermediate_matrices_ok == right.intermediate_matrices_ok &&
           _same_quillen_patch_substitution_steps(left.expected_steps, right.expected_steps) &&
           left.steps_ok == right.steps_ok &&
           left.bracket_matrices == right.bracket_matrices &&
           left.bracket_matrices_ok == right.bracket_matrices_ok &&
           left.final_coefficient == right.final_coefficient &&
           left.final_coefficient_ok == right.final_coefficient_ok &&
           left.base_term == right.base_term &&
           left.base_term_ok == right.base_term_ok &&
           left.telescope_product == right.telescope_product &&
           left.telescope_ok == right.telescope_ok &&
           left.replay_metadata == right.replay_metadata &&
           left.replay_metadata_ok == right.replay_metadata_ok &&
           left.overall_ok == right.overall_ok
end

function _quillen_patch_substitution_chain_verification(
    A,
    R,
    size::Int,
    selected_variable,
    sign_convention::Symbol,
    solver_result::QuillenDenominatorCoverSolverResult,
    cumulative_coefficients,
    intermediate_matrices,
    steps,
    bracket_matrices,
    base_term,
    metadata,
    replay_metadata,
)
    solver_result_ok = verify_quillen_denominator_cover_solver_result(solver_result)
    ring_ok = solver_result_ok && solver_result.ring == R
    matrix_ok = nrows(A) == size && ncols(A) == size && base_ring(A) == R
    selected = _require_substitution_generator(R, selected_variable)
    selected_variable_ok = selected == selected_variable
    sign_convention_ok = sign_convention == :park_woodburn_minus
    expected = _quillen_patch_expected_substitution_chain(
        A,
        selected,
        solver_result,
        sign_convention,
        metadata,
    )

    coefficient_count = length(cumulative_coefficients)
    step_count = length(steps)
    bracket_count = length(bracket_matrices)
    cumulative_coefficients_ok =
        collect(cumulative_coefficients) == expected.cumulative_coefficients
    intermediate_matrices_ok =
        collect(intermediate_matrices) == expected.intermediate_matrices
    steps_ok = _same_quillen_patch_substitution_steps(steps, expected.steps)
    bracket_matrices_ok = collect(bracket_matrices) == expected.bracket_matrices
    final_coefficient_ok = expected.final_coefficient == zero(R)
    base_term_ok = base_term == expected.base_term &&
                   base_term == last(expected.intermediate_matrices)
    telescope_ok = expected.telescope_product == expected.base_term
    replay_metadata_ok = replay_metadata == expected.replay_metadata
    overall_ok =
        solver_result_ok &&
        ring_ok &&
        matrix_ok &&
        selected_variable_ok &&
        sign_convention_ok &&
        coefficient_count == length(expected.cumulative_coefficients) &&
        step_count == length(expected.steps) &&
        bracket_count == length(expected.bracket_matrices) &&
        cumulative_coefficients_ok &&
        intermediate_matrices_ok &&
        steps_ok &&
        bracket_matrices_ok &&
        final_coefficient_ok &&
        base_term_ok &&
        telescope_ok &&
        replay_metadata_ok

    return QuillenPatchSubstitutionChainVerification(
        solver_result_ok,
        ring_ok,
        matrix_ok,
        selected_variable_ok,
        sign_convention_ok,
        coefficient_count,
        step_count,
        bracket_count,
        expected.cumulative_coefficients,
        cumulative_coefficients_ok,
        expected.intermediate_matrices,
        intermediate_matrices_ok,
        expected.steps,
        steps_ok,
        expected.bracket_matrices,
        bracket_matrices_ok,
        expected.final_coefficient,
        final_coefficient_ok,
        expected.base_term,
        base_term_ok,
        expected.telescope_product,
        telescope_ok,
        expected.replay_metadata,
        replay_metadata_ok,
        overall_ok,
    )
end

function replay_quillen_patch_substitution_chain(
    chain::QuillenPatchSubstitutionChain,
)
    R = _require_quillen_denominator_cover_ring(chain.ring)
    _require_square_matrix(chain.original_matrix, "substitution-chain original matrix") ==
        chain.size || throw(DimensionMismatch("substitution-chain size must match original matrix"))
    return _quillen_patch_substitution_chain_verification(
        chain.original_matrix,
        R,
        chain.size,
        chain.selected_variable,
        chain.sign_convention,
        chain.solver_result,
        chain.cumulative_coefficients,
        chain.intermediate_matrices,
        chain.steps,
        chain.bracket_matrices,
        chain.base_term,
        chain.metadata,
        chain.replay_metadata,
    )
end

function verify_quillen_patch_substitution_chain(chain)::Bool
    try
        replay = replay_quillen_patch_substitution_chain(chain)
        return replay.overall_ok &&
               _same_quillen_patch_substitution_chain_verification(
                   chain.verification,
                   replay,
               )
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function quillen_patch_substitution_chain(
    A,
    selected_variable,
    solver_result::QuillenDenominatorCoverSolverResult;
    sign_convention::Symbol = :park_woodburn_minus,
    metadata = (;),
)
    R = _require_quillen_denominator_cover_ring(base_ring(A))
    solver_result.ring == R ||
        throw(ArgumentError("substitution-chain solver result ring must match matrix ring"))
    verify_quillen_denominator_cover_solver_result(solver_result) ||
        throw(ArgumentError("substitution-chain solver result must replay"))
    n = _require_square_matrix(A, "substitution-chain original matrix")
    selected = _require_substitution_generator(R, selected_variable)
    sign_convention == :park_woodburn_minus ||
        throw(ArgumentError("unsupported Park-Woodburn substitution sign convention"))

    expected = _quillen_patch_expected_substitution_chain(
        A,
        selected,
        solver_result,
        sign_convention,
        metadata,
    )
    verification = _quillen_patch_substitution_chain_verification(
        A,
        R,
        n,
        selected,
        sign_convention,
        solver_result,
        expected.cumulative_coefficients,
        expected.intermediate_matrices,
        expected.steps,
        expected.bracket_matrices,
        expected.base_term,
        metadata,
        expected.replay_metadata,
    )
    verification.overall_ok ||
        throw(ArgumentError("Park-Woodburn substitution chain does not replay"))
    return QuillenPatchSubstitutionChain(
        A,
        R,
        n,
        selected,
        sign_convention,
        solver_result,
        expected.cumulative_coefficients,
        expected.intermediate_matrices,
        expected.steps,
        expected.bracket_matrices,
        expected.base_term,
        metadata,
        expected.replay_metadata,
        verification,
    )
end
```

- [ ] **Step 4: Run the focused green test**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_substitution_chain.jl")'
```

Expected: PASS with all tests in `Park-Woodburn substitution chain`.

- [ ] **Step 5: Commit the implementation**

Run:

```bash
git add src/algorithm/quillen_induction.jl
git commit -m "feat: add park woodburn substitution chain"
```

### Task 3: Verify Package Integration

**Files:**
- No new implementation files. This task verifies the full branch state.

**Interfaces:**
- Consumes the registered expert test and package test entrypoint.
- Produces final verification evidence before PR creation.

- [ ] **Step 1: Run diff hygiene**

Run:

```bash
git diff --check
```

Expected: exit 0.

- [ ] **Step 2: Run focused chain verification**

Run:

```bash
julia --project=. -e 'include("test/expert/quillen_patch_substitution_chain.jl")'
```

Expected: exit 0.

- [ ] **Step 3: Run full package verification**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: exit 0.

- [ ] **Step 4: Inspect final git state**

Run:

```bash
git status --short
git log --oneline --decorate -5
```

Expected: clean worktree except no untracked build artifacts; branch includes design, plan, test, and implementation commits.

## Self-Review

- The plan covers every file from the design spec.
- The red test fails before production implementation because all new chain names are absent.
- Negative controls cover exponent, multiplier, sign convention, intermediate matrix, and selected variable corruption.
- The implementation records coefficients, substituted matrices, bracket targets, base term, sign convention, and replay metadata.
- No public exports or out-of-scope Quillen factor assembly are included.
