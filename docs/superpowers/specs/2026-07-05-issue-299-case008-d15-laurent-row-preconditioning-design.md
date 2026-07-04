# Issue 299 Case008 D15 Laurent Row Preconditioning Design

## Goal

Add production reducer support for the `case_008 d=15` Laurent boundary column so `Suslin.reduce_unimodular_column` reduces it to `e_15` and `Suslin.ecp_column_reduction_certificate` verifies through replayable elementary factors.

## Context

Issue #295 added the compact `case_008 d=15` column fixture and internal validator. Issue #297 recorded the pre-repair unsupported diagnostic profile. Issue #298 then found a replay-verified candidate: a bounded row-side synthesis on the `15 x 15` matrix whose effect on the target column is a sequence of left elementary row additions targeting row `1`. The transformed target column has a direct unit in row `1`, so the existing Laurent base reducer supports it by `:unit_entry`.

The production reducer already has a certified `:laurent_elementary_row_preconditioning` stage for a narrow d16 case. That stage runs only after Laurent base reduction fails, constructs row preconditioning from a finite spec, reduces the transformed column through the base Laurent reducer, and verifies all metadata in `_ecp_replay_stage`.

No `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md` file is present in this worktree. GitHub issue body fetches were blocked by the sandbox proxy, so this design uses the issue body supplied by Agent Desk plus local merged PR metadata for #295, #297, and #298.

## Approaches Considered

Recommended: generalize the existing Laurent row-preconditioning spec shape so it can describe either fixed coefficients or an algebraic row-synthesis rule. Keep d16 as a fixed single-source spec. Add a d15 spec with fixed target row `1`, source rows `2:15`, and coefficient strategy `:target_unit_laurent_linear_synthesis`. The strategy recomputes coefficients from the input column by solving `sum(coeff_i * column[source_i]) == one(R) - column[target]`, drops zero coefficients, builds one elementary row-addition factor per nonzero coefficient, and accepts only if the transformed column is reducible by the existing Laurent base reducer.

Alternative: hard-code the fourteen d15 coefficients found by #298. This would be fast but violates the issue requirement that production support not replay a hard-coded factor list or recognize the fixture.

Alternative: add a broad fallback that tries row-synthesis for many target rows, source sets, and lengths. This is more general than the issue asks for and risks claiming arbitrary Laurent support.

Chosen approach: finite column-only row-synthesis spec. It is algebraic over the column and ring, not fixture id; it translates #298's matrix-side candidate into a certified row-preconditioning stage; and it keeps the stage bounded to `length(column) == 15` and two-generator Laurent rings.

## Production Rule

`_laurent_row_preconditioning_specs(column, R)` will return two finite spec families:

- d16 fixed spec: target row `1`, source row `10`, coefficient `one(R)`.
- d15 synthesis spec: target row `1`, source rows `2:15`, coefficient strategy `:target_unit_laurent_linear_synthesis`, maximum nonzero coefficients `14`.

`_laurent_row_preconditioning_candidate(column, R)` will resolve each spec into concrete row factors:

1. For fixed specs, coerce the configured coefficients into `R`.
2. For synthesis specs, call `solve_laurent_linear` on the selected source entries and target residual `one(R) - column[target]`.
3. Reject failed solves, zero-only coefficient vectors, target/source aliasing, out-of-range indices, and coefficient vectors over the nonzero limit.
4. Build preconditioning factors as `elementary_matrix(n, target, source, coeff, R)` for each nonzero coefficient.
5. Compute the transformed column from the product of those factors.
6. Run `_reduce_laurent_unimodular_column_base_certificate` on the transformed column, not the full row-preconditioning reducer, so the stage cannot recurse through itself.
7. Accept only if the transformed certificate verifies and the composed factors reduce the original column to `e_n`.

## Certificate Replay

The `:laurent_elementary_row_preconditioning` stage will store enough data for replay:

- `input_column`;
- `target_index`;
- legacy summary `source_index` and `coefficient`;
- full `source_indices` and `coefficients`;
- `coefficient_strategy`;
- `precondition_factors` and their product `precondition_factor`;
- `transformed_column`;
- `transformed_certificate`;
- composed `factors`;
- `output_column`.

Replay will recompute the accepted spec candidates from the input column and ring, require that the stored indices, coefficients, strategy, factors, transformed column, nested certificate, composed factors, and output column match, and return `false` on tampered metadata.

## Diagnostics

`diagnose_unimodular_column_reduction` will continue to attempt Laurent row-preconditioning only after the base Laurent path fails. For d15, it should report:

- `status == :supported`;
- `failure_code === nothing`;
- `:laurent_elementary_row_preconditioning` in `attempted_stages`;
- stage detail with `outcome == :supported`, target/source metadata, coefficient strategy, coefficient count, and transformed stage `:unit_entry`.

The d15 internal boundary validator and historical witness-profile expert test must stop asserting the permanent unsupported profile. They should now assert the supported profile while still preserving negative controls for non-unimodular inputs.

## Tests

Add `test/expert/case008_d15_laurent_column_reduction.jl` modeled on the d16 expert test. It must assert:

- `reduce_unimodular_column` sends the fixture column exactly to `e_15`;
- replacing one returned factor with identity no longer reduces to `e_15`;
- `ecp_column_reduction_certificate` verifies;
- tampering a certificate factor returns `false`;
- tampering stage coefficient metadata returns `false`;
- tampering stage source metadata returns `false`;
- the certificate contains `:laurent_elementary_row_preconditioning`;
- no stage named `:case008_special_case` appears;
- diagnostics are supported and record the preconditioning detail;
- the non-unimodular negative control fails before reducer stages.

Update:

- `test/internal/toricbuilder_case008_d15_column_boundary.jl`;
- `test/expert/case008_d15_laurent_witness_profile.jl`;
- `test/expert/laurent_column_reduction_diagnostics.jl`;
- existing d16 row-preconditioning assertions as needed for the generalized stage metadata.

Verification commands:

```bash
julia --project=. -e 'include("test/expert/case008_d15_laurent_column_reduction.jl")'
julia --project=. -e 'include("test/internal/toricbuilder_case008_d15_column_boundary.jl")'
julia --project=. -e 'include("test/expert/laurent_column_reduction_diagnostics.jl")'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Out Of Scope

Do not add arbitrary Laurent unimodular-column support. Do not recognize `case_008` fixture ids. Do not hard-code the #298 coefficient list. Do not require full `case_008` certificate success beyond this column-reduction boundary. Do not change the public `elementary_factorization` contract for original Laurent `GL_n` inputs.

## Automatic Decisions

- Visual companion skipped because the task is algebraic reducer work with no visual design question.
- Clarifying questions skipped because Agent Desk is non-interactive and the issue gives the acceptance contract.
- Design approval auto-approved under the Standing Answer Policy.
- Recommended approach selected: finite column-only row-synthesis spec, because it is the narrowest production translation of #298's replay-verified candidate.
