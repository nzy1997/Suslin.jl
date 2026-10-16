### Task 2: Q-Degree Helper And Replay Branch

**Files:**
- Modify: `src/algorithm/sl3_local.jl`
- Modify: `test/expert/sl3_local_q_degree_normalization.jl`

**Interfaces:**
- Consumes: `_sl3_local_target_entries(target)`, `_sl3_local_factor_sequences_equal(left, right)`, and `elementary_matrix`.
- Produces: `_recognize_sl3_local_special_form_matrix(A, X; check_monic=true)`.
- Produces: `_recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic=true)`.
- Produces: `normalize_sl3_local_q_degree(A, X; check_monic::Bool=true)::SL3LocalRealizationCertificate`.
- Produces: `normalize_sl3_local_q_degree(p, q, r, s, X; check_monic::Bool=true)::SL3LocalRealizationCertificate`.
- Extends: `verify_sl3_local_realization(cert)::Bool` for `cert.branch == :q_degree_normalization`.

- [ ] **Step 1: Factor out general special-form recognition**

In `src/algorithm/sl3_local.jl`, replace the bodies of
`_recognize_sl3_local_matrix` and `_recognize_sl3_local_parameters` with a
general special-form validator plus the existing supported-family selection:

```julia
function _recognize_sl3_local_matrix(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_matrix(A, X; check_monic)
    return _recognize_sl3_local_supported_family(form)
end

function _recognize_sl3_local_parameters(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic)
    return _recognize_sl3_local_supported_family(form)
end

function _recognize_sl3_local_special_form_matrix(A, X; check_monic::Bool=true)
    nrows(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))
    ncols(A) == 3 || throw(ArgumentError("local SL_3 special-form recognition failed: A must be 3x3"))

    R = base_ring(A)
    parent(X) === R || throw(ArgumentError("A and X must lie in the same polynomial or Laurent polynomial ring"))
    if A[1, 3] != zero(R) || A[2, 3] != zero(R) ||
            A[3, 1] != zero(R) || A[3, 2] != zero(R) ||
            A[3, 3] != one(R)
        throw(ArgumentError("local SL_3 special-form recognition failed: A must be an embedded 2x2 block with trailing identity"))
    end

    return _recognize_sl3_local_special_form_parameters(A[1, 1], A[1, 2], A[2, 1], A[2, 2], X; check_monic)
end

function _recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic::Bool=true)
    R = parent(p)
    parent(q) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(r) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(s) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))
    parent(X) === R || throw(ArgumentError("p, q, r, s, and X must lie in the same polynomial or Laurent polynomial ring"))

    ring_gens = collect(gens(R))
    var_idx = findfirst(isequal(X), ring_gens)
    var_idx === nothing && throw(ArgumentError("X must be one of the polynomial or Laurent polynomial ring generators"))

    if check_monic
        if _is_laurent_polynomial_ring(R)
            throw(ArgumentError("p monicity check is only supported for ordinary polynomial local SL_3 inputs; pass check_monic=false only when the caller has discharged the local monicity assumption"))
        end
        _is_monic_in_variable(p, var_idx, R) || throw(ArgumentError("p must be monic in X"))
    end

    target = matrix(R, [
        p q zero(R);
        r s zero(R);
        zero(R) zero(R) one(R)
    ])
    det(target) == one(R) || throw(ArgumentError("constructed matrix must have determinant 1"))

    return (; R, p, q, r, s, X, target)
end

function _recognize_sl3_local_supported_family(form)
    if form.s == one(form.R) && form.p == one(form.R) + form.q * form.r
        return merge(form, (; family = :open_s_one))
    end

    if form.p == one(form.R) && form.s == one(form.R) + form.q * form.r
        return merge(form, (; family = :open_p_one))
    end

    s_inverse = _unit_inverse_or_nothing(form.s)
    if s_inverse !== nothing
        return merge(form, (; family = :s_unit, pivot_inverse = s_inverse))
    end

    p_inverse = _unit_inverse_or_nothing(form.p)
    if p_inverse !== nothing
        return merge(form, (; family = :p_unit, pivot_inverse = p_inverse))
    end

    _throw_staged_sl3_local_failure("supported families require one open unipotent slice or a unit diagonal pivot p or s")
end
```

- [ ] **Step 2: Add q-degree helper methods**

In `src/algorithm/sl3_local.jl`, insert these methods after the existing
`realize_sl3_local_certificate` methods:

```julia
function normalize_sl3_local_q_degree(A, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_matrix(A, X; check_monic)
    return _normalize_sl3_local_q_degree_form(form)
end

function normalize_sl3_local_q_degree(p, q, r, s, X; check_monic::Bool=true)
    form = _recognize_sl3_local_special_form_parameters(p, q, r, s, X; check_monic)
    return _normalize_sl3_local_q_degree_form(form)
end
```

- [ ] **Step 3: Add normalization construction helpers**

In `src/algorithm/sl3_local.jl`, insert these helpers after
`_realize_sl3_local_form(form)`:

```julia
function _normalize_sl3_local_q_degree_form(form)
    _sl3_local_require_ordinary_polynomial_ring(form.R, "q-degree normalization")
    quotient, remainder = divrem(form.q, form.p)
    p_degree = degree(form.p, _sl3_local_variable_index(form.R, form.X))
    remainder_degree = _sl3_local_degree_in_variable(remainder, form.R, form.X)
    remainder_degree < p_degree ||
        throw(ArgumentError("q-degree normalization failed to reduce q below p in X"))

    normalized_s = form.s - quotient * form.r
    normalized_target = matrix(form.R, [
        form.p remainder zero(form.R);
        form.r normalized_s zero(form.R);
        zero(form.R) zero(form.R) one(form.R)
    ])
    correction = elementary_matrix(3, 1, 2, quotient, form.R)
    normalized_target * correction == form.target ||
        error("internal q-degree normalization reconstruction failed")

    certificate = SL3LocalRealizationCertificate(
        form.target,
        :q_degree_normalization,
        [correction],
        form.X,
        (;
            p = form.p,
            q = form.q,
            r = form.r,
            s = form.s,
            quotient,
            remainder,
            normalized_target,
            correction,
            p_degree,
            remainder_degree,
        ),
    )
    verify_sl3_local_realization(certificate) ||
        error("internal q-degree normalization certificate verification failed")
    return certificate
end

function _sl3_local_require_ordinary_polynomial_ring(R, label::AbstractString)
    _is_laurent_polynomial_ring(R) &&
        throw(ArgumentError("local SL_3 $(label) is only supported for ordinary polynomial rings"))
    return nothing
end

function _sl3_local_variable_index(R, X)
    var_idx = findfirst(isequal(X), collect(gens(R)))
    var_idx === nothing &&
        throw(ArgumentError("X must be one of the polynomial ring generators"))
    return var_idx
end

function _sl3_local_degree_in_variable(value, R, X)
    return degree(value, _sl3_local_variable_index(R, X))
end
```

- [ ] **Step 4: Extend verifier for partial q-degree certificates**

In `_sl3_local_realization_verification(certificate)`, replace:

```julia
factors_ok = size_ok && verify_factorization(target, certificate.factors)
```

with:

```julia
factors_ok = size_ok && (
    certificate.branch == :q_degree_normalization ?
    _sl3_local_q_degree_reconstruction_ok(certificate) :
    verify_factorization(target, certificate.factors)
)
```

In `_sl3_local_certificate_expected_factors(certificate)`, add this branch before
the final `end`:

```julia
elseif certificate.branch == :q_degree_normalization
    return [witness.correction]
```

In `_sl3_local_branch_witness_ok(certificate)`, add this branch before the final
`end`:

```julia
elseif certificate.branch == :q_degree_normalization
    return _sl3_local_q_degree_witness_ok(certificate)
```

In `_sl3_local_witness_keys_ok(branch, witness)`, add this branch before the
final `return false`:

```julia
elseif branch == :q_degree_normalization
    return keys(witness) == (
        :p,
        :q,
        :r,
        :s,
        :quotient,
        :remainder,
        :normalized_target,
        :correction,
        :p_degree,
        :remainder_degree,
    )
```

Then add these helpers after `_sl3_local_witness_keys_ok`:

```julia
function _sl3_local_q_degree_witness_ok(certificate)
    entries = _sl3_local_target_entries(certificate.target)
    entries === nothing && return false

    R = base_ring(certificate.target)
    witness = certificate.witness
    witness.p == entries.p || return false
    witness.q == entries.q || return false
    witness.r == entries.r || return false
    witness.s == entries.s || return false
    parent(witness.quotient) === R || return false
    parent(witness.remainder) === R || return false
    witness.q == witness.quotient * witness.p + witness.remainder || return false
    witness.p_degree == _sl3_local_degree_in_variable(witness.p, R, certificate.selected_variable) ||
        return false
    witness.remainder_degree ==
        _sl3_local_degree_in_variable(witness.remainder, R, certificate.selected_variable) ||
        return false
    witness.remainder_degree < witness.p_degree || return false

    expected_normalized = matrix(R, [
        witness.p witness.remainder zero(R);
        witness.r witness.s - witness.quotient * witness.r zero(R);
        zero(R) zero(R) one(R)
    ])
    witness.normalized_target == expected_normalized || return false
    witness.correction == elementary_matrix(3, 1, 2, witness.quotient, R) || return false
    return _sl3_local_q_degree_reconstruction_ok(certificate)
end

function _sl3_local_q_degree_reconstruction_ok(certificate)
    witness = certificate.witness
    return witness.normalized_target * witness.correction == certificate.target &&
        _sl3_local_factor_sequences_equal(certificate.factors, [witness.correction])
end
```

- [ ] **Step 5: Run focused test to verify GREEN**

Run:

```bash
julia --project=. -e 'include("test/expert/sl3_local_q_degree_normalization.jl")'
```

Expected: PASS with the `local SL3 q-degree normalization certificate` testset.

---

