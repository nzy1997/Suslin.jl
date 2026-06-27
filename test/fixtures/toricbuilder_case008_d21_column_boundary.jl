module ToricBuilderCase008D21ColumnBoundary

using Oscar
using Suslin

include(joinpath(@__DIR__, "toricbuilder_cache_q_blocks.jl"))

const CASE_ID = "case_008"
const EXPECTED_RING_DESCRIPTION = "GF(2)[u^+/-1, v^+/-1]"
const EXPECTED_PASSED_PEEL_DIMENSIONS = (30, 29, 28, 27, 26, 25, 24, 23, 22)
const FIRST_FAILING_PEEL_DIMENSION = 21
const REQUIRED_BOUNDARY_FIELDS = (
    :case_id,
    :source_entry,
    :original_matrix,
    :source_matrix_dimensions,
    :source_sparse_entry_count,
    :normalization_provenance,
    :passed_peel_dimensions,
    :first_failing_peel_dimension,
    :failing_input_matrix,
    :failing_column,
    :ring,
    :ring_description,
    :expected_diagnostic,
)

const FAILING_INPUT_SPARSE_ENTRIES = (
    (1, 1, "u^6*v^-1 + u^5 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3 + u^3*v^-3 + u^2*v^2 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v^-1"),
    (1, 2, "u^3*v^-1 + u^2 + u^2*v^-3 + u*v^-1 + u*v^-2 + u*v^-3"),
    (1, 3, "u^2 + u^2*v^-1 + u^2*v^-2 + u*v + u + u*v^-2"),
    (1, 4, "u^2*v^-2 + u*v^-1 + u*v^-2"),
    (1, 5, "u^3*v^-1 + u^2 + u^2*v^-3 + u*v^-1 + u*v^-2 + u*v^-3"),
    (1, 6, "u^3*v^-1 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^-3"),
    (1, 7, "u^3*v^-1 + u^2*v^-2 + u*v + u*v^-2"),
    (1, 8, "u^2*v^-2 + u*v^-1 + u*v^-2"),
    (1, 9, "u^2*v^-1 + u^2*v^-3 + u + u*v^-1 + u*v^-2 + u*v^-3"),
    (1, 10, "u^3*v^-1 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^-3"),
    (1, 11, "u^2 + u^2*v^-1 + u^2*v^-2 + u*v + u + u*v^-2"),
    (1, 12, "u^6*v^-1 + u^5 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-2 + u^3*v^2 + u^3 + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v^-1 + u^2*v^-3 + u*v^2 + u*v^-1 + v + 1 + v^-2 + v^-3"),
    (1, 14, "u^3*v^-1 + u^2 + u^2*v^-1 + u"),
    (1, 15, "u^2*v^-2 + u^2*v^-3 + u*v^-2 + u*v^-3 + 1 + v^-2"),
    (1, 16, "u^6*v^-1 + u^5 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-2 + u^3*v^2 + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v^-1 + u*v^-2 + u*v^-3 + v + v^-3"),
    (1, 17, "u^6*v^-1 + u^5 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3*v^-3 + u^2*v + u^2*v^-3 + u*v^2 + u*v^-1 + u*v^-3 + v + v^-3"),
    (1, 18, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-1 + u^3*v^-2 + u^2*v^2 + u^2*v^-1 + u^2*v^-3 + u*v + u"),
    (1, 19, "u^6*v^-1 + u^5 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-2 + u^3*v^2 + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-1 + u*v^-2 + u*v^-3 + v + v^-3"),
    (1, 20, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-1 + u^3*v^-2 + u^2*v^2 + u^2*v^-2 + u^2*v^-3 + u*v + u*v^-2"),
    (1, 21, "u^6*v^-1 + u^5 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3*v^-3 + u^2*v + u^2*v^-1 + u^2*v^-3 + u*v^2 + u + u*v^-1 + u*v^-2 + u*v^-3 + v + 1 + v^-2 + v^-3"),
    (2, 1, "u^3 + u^2*v^-1 + u + v^-1"),
    (2, 2, "u^-1*v^-1"),
    (2, 5, "u^-1*v^-1"),
    (2, 6, "u^-1*v^-1"),
    (2, 9, "u^-1*v^-1"),
    (2, 10, "u^-1*v^-1"),
    (2, 12, "u^3 + u^2 + u^2*v^-1 + u*v^-1 + u^-1 + u^-1*v^-1 + u^-2 + u^-2*v^-1"),
    (2, 15, "1 + u^-1 + u^-1*v^-1"),
    (2, 16, "u^3 + u^2 + u^2*v^-1 + u*v^-1 + 1 + u^-2 + u^-2*v^-1"),
    (2, 17, "u^3 + u^2*v^-1 + u + v^-1 + u^-1 + u^-2 + u^-2*v^-1"),
    (2, 18, "u^2 + u + u*v^-1 + v^-1"),
    (2, 19, "u^3 + u^2 + u^2*v^-1 + u*v^-1 + 1 + u^-2 + u^-2*v^-1"),
    (2, 20, "u^2 + u + u*v^-1 + 1 + v^-1 + u^-1"),
    (2, 21, "u^3 + u^2*v^-1 + u + 1 + v^-1 + u^-2 + u^-2*v^-1"),
    (3, 1, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-1 + u^3*v^-3 + u^2*v^2 + u^2 + u^2*v^-3 + u*v^2 + u*v^-2 + u*v^-3 + v^2 + 1 + v^-1"),
    (3, 2, "u^2*v^-1 + u + u*v^-3 + v^-1 + v^-2 + v^-3"),
    (3, 3, "u + u*v^-1 + u*v^-2 + v + v^-2"),
    (3, 4, "u*v^-2 + v^-1 + v^-2"),
    (3, 5, "u^2*v^-1 + u + u*v^-3 + v^-1 + v^-2 + v^-3"),
    (3, 6, "u^2*v^-1 + u + u*v^-2 + u*v^-3 + v^-3"),
    (3, 7, "u^2*v^-1 + u*v^-2 + v + v^-2"),
    (3, 8, "u*v^-2 + v^-1 + v^-2"),
    (3, 9, "u*v^-1 + u*v^-3 + 1 + v^-1 + v^-2 + v^-3"),
    (3, 10, "u^2*v^-1 + u + u*v^-2 + u*v^-3 + v^-3"),
    (3, 11, "u + u*v^-1 + u*v^-2 + v + 1 + v^-2"),
    (3, 12, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-2 + u^2*v^2 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v^-1 + u*v^-3 + v^2 + v^-1 + u^-1*v + u^-1 + u^-1*v^-2 + u^-1*v^-3"),
    (3, 14, "u^2*v^-1 + u + u*v^-1"),
    (3, 15, "u*v^-2 + u*v^-3 + v^-2 + v^-3 + u^-1 + u^-1*v^-2"),
    (3, 16, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-2 + u^2*v^2 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v + u*v^-1 + u*v^-2 + v^2 + v^-1 + v^-2 + v^-3 + u^-1*v + u^-1*v^-3"),
    (3, 17, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-1 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2*v^-3 + u*v + u*v^-3 + v^2 + v^-1 + v^-3 + u^-1*v + u^-1*v^-3"),
    (3, 18, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v^-1 + u*v^-3 + v"),
    (3, 19, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-2 + u^2*v^2 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v + u + u*v^-1 + u*v^-2 + v^2 + v + v^-1 + v^-2 + v^-3 + u^-1*v + u^-1*v^-3"),
    (3, 20, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v^-2 + u*v^-3 + v + v^-2"),
    (3, 21, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-1 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2*v^-3 + u*v + u*v^-1 + u*v^-3 + v^2 + 1 + v^-1 + v^-2 + v^-3 + u^-1*v + u^-1 + u^-1*v^-2 + u^-1*v^-3"),
    (4, 1, "u^3*v^-1 + u^3*v^-2 + u^2 + u^2*v^-2 + u*v^-1 + u*v^-2 + v^-2"),
    (4, 2, "v^-2 + u^-1*v^-2"),
    (4, 4, "v^-1"),
    (4, 5, "v^-1 + v^-2 + u^-1*v^-2"),
    (4, 6, "v^-1 + v^-2 + u^-1*v^-2"),
    (4, 8, "v^-1"),
    (4, 9, "v^-2 + u^-1*v^-2"),
    (4, 10, "v^-1 + v^-2 + u^-1*v^-2"),
    (4, 11, "v^-1"),
    (4, 12, "u^3*v^-1 + u^3*v^-2 + u^2 + u^2*v^-1 + u + u*v^-1 + u*v^-2 + 1 + v^-2 + u^-1 + u^-1*v^-1 + u^-2*v^-1 + u^-2*v^-2"),
    (4, 15, "v^-1 + v^-2 + u^-1*v^-1 + u^-1*v^-2"),
    (4, 16, "u^3*v^-1 + u^3*v^-2 + u^2 + u^2*v^-1 + u + u*v^-1 + u*v^-2 + 1 + v^-1 + u^-1 + u^-1*v^-2 + u^-2*v^-1 + u^-2*v^-2"),
    (4, 17, "u^3*v^-1 + u^3*v^-2 + u^2 + u^2*v^-2 + u*v^-2 + v^-2 + u^-1 + u^-1*v^-1 + u^-1*v^-2 + u^-2*v^-1 + u^-2*v^-2"),
    (4, 18, "u^2*v^-1 + u^2*v^-2 + u + u*v^-1 + 1 + v^-2"),
    (4, 19, "u^3*v^-1 + u^3*v^-2 + u^2 + u^2*v^-1 + u + u*v^-1 + u*v^-2 + 1 + v^-1 + u^-1 + u^-1*v^-2 + u^-2*v^-1 + u^-2*v^-2"),
    (4, 20, "u^2*v^-1 + u^2*v^-2 + u + u*v^-1 + 1 + v^-1 + v^-2 + u^-1*v^-1"),
    (4, 21, "u^3*v^-1 + u^3*v^-2 + u^2 + u^2*v^-2 + u*v^-2 + v^-2 + u^-1 + u^-1*v^-2 + u^-2*v^-1 + u^-2*v^-2"),
    (5, 1, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^3*v + u^3 + u^2*v^2 + u^2*v^-2 + u*v^2 + u + u*v^-1 + v^2"),
    (5, 2, "u^2*v^-1 + u + u*v^-2 + 1"),
    (5, 3, "u + u*v^-1 + v + 1"),
    (5, 4, "u*v^-1 + 1"),
    (5, 5, "u^2*v^-1 + u + u*v^-2 + 1"),
    (5, 6, "u^2*v^-1 + u + u*v^-1 + u*v^-2"),
    (5, 7, "u^2*v^-1 + v + 1"),
    (5, 8, "u*v^-1 + 1"),
    (5, 9, "u*v^-1 + u*v^-2 + 1"),
    (5, 10, "u^2*v^-1 + u + u*v^-1 + u*v^-2"),
    (5, 11, "u + u*v^-1 + v + 1"),
    (5, 12, "u^5*v^-1 + u^4 + u^4*v^-2 + u^3*v + u^3*v^-2 + u^2*v^2 + u^2 + u*v^2 + u*v^-2 + v^2 + v + v^-1 + v^-2"),
    (5, 14, "u^2*v^-1 + u + u*v^-1"),
    (5, 15, "u*v^-1 + u*v^-2 + 1 + u^-1"),
    (5, 16, "u^5*v^-1 + u^4 + u^4*v^-2 + u^3*v + u^3*v^-2 + u^2*v^2 + u*v^2 + u*v + u*v^-1 + v^2 + v + 1 + v^-1 + v^-2 + u^-1"),
    (5, 17, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^3*v + u^3 + u^2*v^2 + u^2*v + u^2 + u^2*v^-2 + u*v + u*v^-1 + v^2 + v + 1 + v^-1 + v^-2 + u^-1"),
    (5, 18, "u^4*v^-1 + u^3 + u^3*v^-2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + v + 1"),
    (5, 19, "u^5*v^-1 + u^4 + u^4*v^-2 + u^3*v + u^3*v^-2 + u^2*v^2 + u*v^2 + u*v + u + u*v^-1 + v^2 + v^-1 + v^-2 + u^-1"),
    (5, 20, "u^4*v^-1 + u^3 + u^3*v^-2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v^-1 + v + 1"),
    (5, 21, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^3*v + u^3 + u^2*v^2 + u^2*v + u^2 + u^2*v^-2 + u*v + v^2 + v + v^-1 + v^-2"),
    (6, 1, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-3 + u^2*v + u^2*v^-1 + u^2*v^-3 + u*v^2 + u + u*v^-1 + u*v^-2 + u*v^-3 + v^2 + v + u^-1*v^2"),
    (6, 2, "u^2*v^-1 + u*v^-3 + v^-2 + v^-3 + u^-1 + u^-1*v^-1"),
    (6, 3, "u + u*v^-1 + u*v^-2 + v^-2 + u^-1*v + u^-1"),
    (6, 4, "u*v^-2 + v^-1 + v^-2 + u^-1"),
    (6, 5, "u^2*v^-1 + u*v^-3 + 1 + v^-2 + v^-3 + u^-1 + u^-1*v^-1"),
    (6, 6, "u^2*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-1 + v^-3 + u^-1*v^-1"),
    (6, 7, "u^2*v^-1 + u + u*v^-2 + 1 + v^-2 + u^-1*v + u^-1"),
    (6, 8, "u*v^-2 + v^-1 + v^-2 + u^-1"),
    (6, 9, "u*v^-1 + u*v^-3 + v^-2 + v^-3 + u^-1 + u^-1*v^-1"),
    (6, 10, "u^2*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-1 + v^-3 + u^-1*v^-1"),
    (6, 11, "u + u*v^-1 + u*v^-2 + v^-2 + u^-1*v + u^-1"),
    (6, 12, "u^5*v^-1 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^2*v + u^2*v^-2 + u^2*v^-3 + u*v^2 + u + u*v^-3 + v^2 + v + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2 + u^-2*v^-1"),
    (6, 13, "1"),
    (6, 14, "u^2*v^-1 + u*v^-1"),
    (6, 15, "u*v^-2 + u*v^-3 + 1 + v^-1 + v^-2 + v^-3 + u^-1*v^-1 + u^-1*v^-2 + u^-2"),
    (6, 16, "u^5*v^-1 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^2*v + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u + u*v^-2 + v^2 + 1 + v^-1 + v^-2 + v^-3 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-1 + u^-1*v^-3 + u^-2*v^-1"),
    (6, 17, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-3 + u^2 + u^2*v^-1 + u^2*v^-3 + u*v^2 + u*v^-1 + u*v^-3 + 1 + v^-3 + u^-1*v^2 + u^-1*v + u^-1*v^-1 + u^-1*v^-3 + u^-2*v^-1"),
    (6, 18, "u^4*v^-1 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-2 + u*v^-1 + u*v^-3 + v^2 + 1 + v^-1 + u^-1*v + u^-1"),
    (6, 19, "u^5*v^-1 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^2*v + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v^-2 + v^2 + v^-1 + v^-2 + v^-3 + u^-1*v^2 + u^-1*v^-1 + u^-1*v^-3 + u^-2*v^-1"),
    (6, 20, "u^4*v^-1 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-2 + u*v^-2 + u*v^-3 + v^2 + 1 + v^-1 + v^-2 + u^-1*v"),
    (6, 21, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-3 + u^2 + u^2*v^-1 + u^2*v^-3 + u*v^2 + u*v^-3 + 1 + v^-2 + v^-3 + u^-1*v^2 + u^-1*v + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2 + u^-2*v^-1"),
    (7, 1, "u*v^2 + v + 1 + v^-1 + u^-1*v^3 + u^-2*v^3 + u^-2*v + u^-2*v^-1 + u^-3*v^3 + u^-3*v"),
    (7, 2, "u^-2*v^2 + u^-3*v^2 + u^-3*v + u^-3*v^-1"),
    (7, 3, "u^-3*v + u^-3"),
    (7, 4, "u^-3*v^2 + u^-3"),
    (7, 5, "u^-3*v^2 + u^-3*v + u^-3*v^-1"),
    (7, 6, "u^-2*v^2 + u^-3*v + u^-3 + u^-3*v^-1"),
    (7, 7, "u^-2*v + u^-3*v + u^-3"),
    (7, 8, "u^-3*v^2 + u^-3"),
    (7, 9, "u^-2*v + u^-3*v^2 + u^-3*v + u^-3*v^-1"),
    (7, 10, "u^-3*v + u^-3 + u^-3*v^-1"),
    (7, 11, "u^-3*v + u^-3"),
    (7, 12, "u*v^2 + v^2 + v + 1 + v^-1 + u^-1*v^3 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v^2 + u^-2*v + u^-2 + u^-3*v + u^-3 + u^-3*v^-1 + u^-4*v^2 + u^-4*v^-1"),
    (7, 14, "u^-2*v^2 + u^-2*v"),
    (7, 15, "u^-2*v^2 + u^-2*v + u^-3 + u^-3*v^-1 + u^-4*v + u^-4"),
    (7, 16, "u*v^2 + v^2 + v + 1 + v^-1 + u^-1*v^3 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v^2 + u^-2 + u^-3*v + u^-4*v^2 + u^-4*v + u^-4 + u^-4*v^-1"),
    (7, 17, "u*v^2 + v + 1 + v^-1 + u^-1*v^3 + u^-2*v^3 + u^-2*v + u^-2*v^-1 + u^-3 + u^-4*v^2 + u^-4*v + u^-4 + u^-4*v^-1"),
    (7, 18, "v^2 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v^3 + u^-2 + u^-2*v^-1"),
    (7, 19, "u*v^2 + v^2 + v + 1 + v^-1 + u^-1*v^3 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v^2 + u^-2 + u^-3*v + u^-4*v^2 + u^-4*v + u^-4 + u^-4*v^-1"),
    (7, 20, "v^2 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v^3 + u^-2*v^2 + u^-2 + u^-2*v^-1 + u^-3*v + u^-3"),
    (7, 21, "u*v^2 + v + 1 + v^-1 + u^-1*v^3 + u^-2*v^3 + u^-2*v + u^-2*v^-1 + u^-3*v + u^-3 + u^-4*v^2 + u^-4*v^-1"),
    (8, 1, "u^3*v^-1 + u^2*v + u*v^-1"),
    (8, 2, "1 + v^-1"),
    (8, 5, "1 + v^-1"),
    (8, 6, "1 + v^-1"),
    (8, 8, "1"),
    (8, 9, "v^-1"),
    (8, 10, "1 + v^-1"),
    (8, 11, "1"),
    (8, 12, "u^3*v^-1 + u^2*v + u^2*v^-1 + u*v + u + v + v^-1 + u^-1*v + u^-1*v^-1"),
    (8, 15, "v^-1"),
    (8, 16, "u^3*v^-1 + u^2*v + u^2*v^-1 + u*v + u + v + u^-1*v + u^-1*v^-1"),
    (8, 17, "u^3*v^-1 + u^2*v + u + u*v^-1 + u^-1*v + u^-1*v^-1"),
    (8, 18, "u^2*v^-1 + u*v + u*v^-1 + v"),
    (8, 19, "u^3*v^-1 + u^2*v + u^2*v^-1 + u*v + u + v + u^-1*v + u^-1*v^-1"),
    (8, 20, "u^2*v^-1 + u*v + u*v^-1 + v"),
    (8, 21, "u^3*v^-1 + u^2*v + u + u*v^-1 + 1 + u^-1*v + u^-1*v^-1"),
    (9, 1, "u^3 + u^3*v^-1 + u + u*v^-1"),
    (9, 2, "1 + v^-1"),
    (9, 5, "1 + v^-1"),
    (9, 6, "1 + v^-1"),
    (9, 9, "1 + v^-1"),
    (9, 10, "1 + v^-1"),
    (9, 12, "u^3 + u^3*v^-1 + u^2 + u^2*v^-1 + u + v + v^-1 + u^-1*v + u^-1*v^-1"),
    (9, 15, "1 + v^-1"),
    (9, 16, "u^3 + u^3*v^-1 + u^2 + u^2*v^-1 + u + v + 1 + u^-1*v + u^-1*v^-1"),
    (9, 17, "u^3 + u^3*v^-1 + u*v^-1 + v + 1 + u^-1*v + u^-1*v^-1"),
    (9, 18, "u^2 + u^2*v^-1 + u + u*v^-1"),
    (9, 19, "u^3 + u^3*v^-1 + u^2 + u^2*v^-1 + u + v + 1 + u^-1*v + u^-1*v^-1"),
    (9, 20, "u^2 + u^2*v^-1 + u + u*v^-1"),
    (9, 21, "u^3 + u^3*v^-1 + u*v^-1 + v + 1 + u^-1*v + u^-1*v^-1"),
    (10, 1, "u^6*v^-1 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-2 + u^3 + u^2 + u^2*v^-1 + u^2*v^-2 + u*v + u*v^-2 + u*v^-3 + 1 + v^-1"),
    (10, 2, "u^3*v^-1 + u^2*v^-1 + u^2*v^-3 + u*v^-2 + v^-2 + v^-3"),
    (10, 3, "u^2 + u^2*v^-1 + u^2*v^-2 + u + u*v^-1 + v^-2"),
    (10, 4, "u^2*v^-2 + u*v^-1 + v^-1 + v^-2"),
    (10, 5, "u^3*v^-1 + u^2*v^-1 + u^2*v^-3 + u*v^-2 + v^-2 + v^-3"),
    (10, 6, "u^3*v^-1 + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^-1 + u*v^-2 + v^-1 + v^-3"),
    (10, 7, "u^3*v^-1 + u^2 + u^2*v^-1 + u^2*v^-2 + u + v^-2"),
    (10, 8, "u^2*v^-2 + u*v^-1 + v^-1 + v^-2"),
    (10, 9, "u^2*v^-1 + u^2*v^-3 + u*v^-1 + u*v^-2 + v^-2 + v^-3"),
    (10, 10, "u^3*v^-1 + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^-1 + u*v^-2 + v^-1 + v^-3"),
    (10, 11, "u^2 + u^2*v^-1 + u^2*v^-2 + u + u*v^-1 + v^-2"),
    (10, 12, "u^6*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-3 + u^3 + u^3*v^-3 + u^2 + u^2*v^-2 + u*v + u + u*v^-1 + u*v^-3 + v^-1 + v^-2 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3"),
    (10, 13, "u"),
    (10, 14, "u^3*v^-1 + u*v^-1"),
    (10, 15, "u^2*v^-2 + u^2*v^-3 + u*v^-1 + v^-1 + v^-3 + u^-1*v^-2"),
    (10, 16, "u^6*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-3 + u^3*v^-3 + u^2*v^-3 + u*v + u + u*v^-3 + v^-2 + u^-1 + u^-1*v^-1 + u^-1*v^-3"),
    (10, 17, "u^6*v^-1 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-2 + u^3*v + u^2*v + u^2 + u^2*v^-1 + u*v + u + u^-1 + u^-1*v^-1 + u^-1*v^-3"),
    (10, 18, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-1 + u^3*v^-3 + u^2*v + u^2*v^-2 + u^2*v^-3 + u + u*v^-3 + 1"),
    (10, 19, "u^6*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-3 + u^3*v^-3 + u^2 + u^2*v^-3 + u*v + u + u*v^-3 + 1 + v^-2 + u^-1 + u^-1*v^-1 + u^-1*v^-3"),
    (10, 20, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-1 + u^3*v^-3 + u^2*v + u^2*v^-1 + u^2*v^-3 + u + u*v^-1 + u*v^-3 + 1 + v^-2"),
    (10, 21, "u^6*v^-1 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-2 + u^3*v + u^2*v + u^2 + u*v + u*v^-1 + u*v^-2 + u^-1 + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3"),
    (11, 1, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v^-3 + v^2 + 1 + v^-2 + u^-1*v + u^-1"),
    (11, 2, "u^2*v^-1 + u + u*v^-3 + v^-1 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (11, 3, "u + u*v^-1 + u*v^-2 + v + 1 + v^-1 + v^-2 + u^-1*v^-1"),
    (11, 4, "u*v^-2 + v^-2 + u^-1 + u^-1*v^-1"),
    (11, 5, "u^2*v^-1 + u + u*v^-3 + v^-1 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (11, 6, "u^2*v^-1 + u + u*v^-2 + u*v^-3 + v^-1 + v^-2 + v^-3 + u^-1*v^-2"),
    (11, 7, "u^2*v^-1 + u*v^-2 + v + v^-1 + v^-2 + u^-1*v^-1"),
    (11, 8, "u*v^-2 + v^-2 + u^-1 + u^-1*v^-1"),
    (11, 9, "u*v^-1 + u*v^-3 + 1 + v^-1 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (11, 10, "u^2*v^-1 + u + u*v^-2 + u*v^-3 + v^-1 + v^-2 + v^-3 + u^-1*v^-2"),
    (11, 11, "u + u*v^-1 + u*v^-2 + v + v^-1 + v^-2 + u^-1*v^-1"),
    (11, 12, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^2*v^2 + u^2*v + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v + u*v^-2 + u*v^-3 + v^2 + v^-1 + v^-2 + u^-1*v + u^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v^-1 + u^-2*v^-2"),
    (11, 14, "u^2*v^-1 + u + u*v^-1 + 1"),
    (11, 15, "u*v^-2 + u*v^-3 + v^-1 + v^-3 + u^-1*v^-1 + u^-2*v^-1"),
    (11, 16, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^2*v^2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^2 + v^2 + v^-2 + v^-3 + u^-1*v + u^-1 + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v^-2"),
    (11, 17, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v + u + u*v^-2 + u*v^-3 + v^2 + v + 1 + v^-1 + v^-2 + v^-3 + u^-1*v + u^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v^-2"),
    (11, 18, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u*v^2 + u*v + u + u*v^-3 + 1 + v^-2 + u^-1*v"),
    (11, 19, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^2*v^2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u + v^2 + v + v^-2 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v^-2"),
    (11, 20, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u*v^2 + u*v + u + u*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-1 + u^-1*v + u^-1*v^-1"),
    (11, 21, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v + u + u*v^-1 + u*v^-2 + u*v^-3 + v^2 + v + v^-1 + v^-3 + u^-1*v + u^-1 + u^-1*v^-1 + u^-1*v^-3 + u^-2*v^-1 + u^-2*v^-2"),
    (12, 1, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-2 + u*v^2 + u + u*v^-1 + u*v^-3 + v^2 + v + v^-1 + v^-2 + u^-1*v^2 + u^-1*v + u^-1"),
    (12, 2, "u*v^-1 + 1 + v^-1 + v^-3 + u^-1 + u^-1*v^-2"),
    (12, 3, "1 + v^-1 + v^-2 + u^-1*v + u^-1 + u^-1*v^-1"),
    (12, 4, "v^-2 + u^-1*v^-1"),
    (12, 5, "u*v^-1 + 1 + v^-1 + v^-3 + u^-1 + u^-1*v^-2"),
    (12, 6, "u*v^-1 + 1 + v^-1 + v^-2 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (12, 7, "u*v^-1 + v^-1 + v^-2 + u^-1*v + u^-1 + u^-1*v^-1"),
    (12, 8, "v^-2 + u^-1*v^-1"),
    (12, 9, "v^-1 + v^-3 + u^-1 + u^-1*v^-2"),
    (12, 10, "u*v^-1 + 1 + v^-1 + v^-2 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (12, 11, "1 + v^-1 + v^-2 + u^-1*v + u^-1 + u^-1*v^-1"),
    (12, 12, "u^4*v^-1 + u^3 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2 + u^2*v^-3 + u*v^2 + u + u*v^-1 + v^2 + v + 1 + v^-1 + v^-2 + v^-3 + u^-1*v^2 + u^-1*v + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v + u^-2*v^-2"),
    (12, 14, "u*v^-1 + 1"),
    (12, 15, "v^-2 + v^-3 + u^-2 + u^-2*v^-1"),
    (12, 16, "u^4*v^-1 + u^3 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2 + u^2*v^-3 + u*v^2 + u*v^-1 + v^2 + v^-1 + u^-1*v^2 + u^-1*v + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v + u^-2 + u^-2*v^-1 + u^-2*v^-2"),
    (12, 17, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-2 + u*v^2 + u*v + u*v^-1 + u*v^-3 + 1 + v^-1 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v + u^-2 + u^-2*v^-1 + u^-2*v^-2"),
    (12, 18, "u^3*v^-1 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v + u*v^-1 + u*v^-3 + v^2 + 1 + v^-1 + v^-2 + u^-1*v"),
    (12, 19, "u^4*v^-1 + u^3 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2 + u^2*v^-3 + u*v^2 + u*v^-1 + v^2 + 1 + v^-1 + u^-1*v^2 + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v + u^-2 + u^-2*v^-1 + u^-2*v^-2"),
    (12, 20, "u^3*v^-1 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v + u*v^-1 + u*v^-3 + v^2 + 1 + u^-1*v + u^-1 + u^-1*v^-1"),
    (12, 21, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-2 + u*v^2 + u*v + u*v^-1 + u*v^-3 + 1 + u^-1*v^2 + u^-1*v + u^-1*v^-1 + u^-1*v^-3 + u^-2*v + u^-2*v^-2"),
    (13, 1, "u^5 + u^5*v^-1 + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3 + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u + u*v^-3 + v + 1 + v^-1"),
    (13, 2, "u^2 + u^2*v^-1 + u*v^-2 + u*v^-3 + v^-1 + v^-3"),
    (13, 3, "u*v + u*v^-2 + v^-1 + v^-2"),
    (13, 4, "u*v^-1 + u*v^-2 + 1 + v^-2"),
    (13, 5, "u^2 + u^2*v^-1 + u*v^-2 + u*v^-3 + v^-1 + v^-3"),
    (13, 6, "u^2 + u^2*v^-1 + u*v^-1 + u*v^-3 + 1 + v^-1 + v^-2 + v^-3"),
    (13, 7, "u^2 + u^2*v^-1 + u*v + u + u*v^-1 + u*v^-2 + v^-1 + v^-2"),
    (13, 8, "u*v^-1 + u*v^-2 + 1 + v^-2"),
    (13, 9, "u + u*v^-1 + u*v^-2 + u*v^-3 + v^-1 + v^-3"),
    (13, 10, "u^2 + u^2*v^-1 + u*v^-1 + u*v^-3 + 1 + v^-1 + v^-2 + v^-3"),
    (13, 11, "u*v + u*v^-2 + v^-1 + v^-2"),
    (13, 12, "u^5 + u^5*v^-1 + u^4 + u^4*v^-3 + u^3*v^2 + u^3*v^-2 + u^2*v^2 + u^2*v + u^2*v^-3 + u*v^2 + u*v^-1 + u*v^-2 + u*v^-3 + u^-1*v + u^-1 + u^-1*v^-3"),
    (13, 13, "v"),
    (13, 14, "u^2 + u^2*v^-1 + u + u*v^-1"),
    (13, 15, "u*v^-1 + u*v^-3 + v^-3 + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (13, 16, "u^5 + u^5*v^-1 + u^4 + u^4*v^-3 + u^3*v^2 + u^3*v^-2 + u^2*v^2 + u^2 + u^2*v^-3 + u*v^2 + u*v^-2 + v^-3 + u^-1*v + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3"),
    (13, 17, "u^5 + u^5*v^-1 + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v + u*v^-1 + u*v^-2 + u*v^-3 + v + 1 + v^-1 + v^-2 + v^-3 + u^-1*v + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3"),
    (13, 18, "u^4 + u^4*v^-1 + u^3 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2*v^-1 + u^2*v^-2 + u + u*v^-2 + u*v^-3 + v"),
    (13, 19, "u^5 + u^5*v^-1 + u^4 + u^4*v^-3 + u^3*v^2 + u^3*v^-2 + u^2*v^2 + u^2 + u^2*v^-3 + u*v^2 + u*v + u + u*v^-2 + v + v^-3 + u^-1*v + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3"),
    (13, 20, "u^4 + u^4*v^-1 + u^3 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2*v^-1 + u^2*v^-2 + u*v^-3 + v + v^-1 + v^-2"),
    (13, 21, "u^5 + u^5*v^-1 + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v + u + u*v^-2 + u*v^-3 + 1 + v^-3 + u^-1*v + u^-1 + u^-1*v^-3"),
    (14, 1, "u^3 + u^3*v^-1 + u^2*v + u^2 + u + u*v^-1"),
    (14, 2, "1 + v^-1"),
    (14, 5, "1 + v^-1"),
    (14, 6, "1 + v^-1"),
    (14, 9, "1 + v^-1"),
    (14, 10, "1 + v^-1"),
    (14, 12, "u^3 + u^3*v^-1 + u^2*v + u^2*v^-1 + u*v + u + v + v^-1 + u^-1*v + u^-1*v^-1"),
    (14, 14, "v"),
    (14, 15, "1 + v^-1"),
    (14, 16, "u^3 + u^3*v^-1 + u^2*v + u^2*v^-1 + u*v + u + v + 1 + u^-1*v + u^-1*v^-1"),
    (14, 17, "u^3 + u^3*v^-1 + u^2*v + u^2 + u + u*v^-1 + u^-1*v + u^-1*v^-1"),
    (14, 18, "u^2 + u^2*v^-1 + u*v + u*v^-1"),
    (14, 19, "u^3 + u^3*v^-1 + u^2*v + u^2*v^-1 + u*v + u + v + 1 + u^-1*v + u^-1*v^-1"),
    (14, 20, "u^2 + u^2*v^-1 + u*v + u*v^-1 + v + 1"),
    (14, 21, "u^3 + u^3*v^-1 + u^2*v + u^2 + u + u*v^-1 + u^-1*v + u^-1*v^-1"),
    (15, 1, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^2*v^2 + u^2*v + u^2 + u^2*v^-3 + u*v^2 + u*v + u*v^-2 + v^2 + v + 1"),
    (15, 2, "u^2*v^-1 + u + u*v^-3 + 1 + v^-1 + v^-2"),
    (15, 3, "u + u*v^-1 + u*v^-2 + v + 1 + v^-1"),
    (15, 4, "u*v^-2 + v^-1"),
    (15, 5, "u^2*v^-1 + u + u*v^-3 + 1 + v^-1 + v^-2"),
    (15, 6, "u^2*v^-1 + u + u*v^-2 + u*v^-3 + 1 + v^-2"),
    (15, 7, "u^2*v^-1 + u*v^-1 + u*v^-2 + v + 1 + v^-1"),
    (15, 8, "u*v^-2 + v^-1"),
    (15, 9, "u*v^-3 + 1 + v^-1 + v^-2"),
    (15, 10, "u^2*v^-1 + u + u*v^-2 + u*v^-3 + 1 + v^-2"),
    (15, 11, "u + u*v^-1 + u*v^-2 + v + 1 + v^-1"),
    (15, 12, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-3 + u^2*v^2 + u^2*v + u*v^2 + u + u*v^-2 + u*v^-3 + v^2 + 1 + v^-1 + v^-2 + v^-3 + u^-1*v + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (15, 14, "u^2*v^-1 + u"),
    (15, 15, "u*v^-1 + u*v^-2 + u*v^-3 + v^-1 + u^-1 + u^-1*v^-1"),
    (15, 16, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-3 + u^2*v^2 + u^2*v + u^2 + u*v^2 + u*v + u*v^-1 + v^2 + 1 + v^-2 + v^-3 + u^-1*v + u^-1*v^-2"),
    (15, 17, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^2*v^2 + u^2*v^-3 + v^2 + 1 + v^-1 + v^-2 + v^-3 + u^-1*v + u^-1*v^-2"),
    (15, 18, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-3 + u*v^2 + u*v + u + u*v^-2 + v"),
    (15, 19, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3*v^-3 + u^2*v^2 + u^2*v + u^2 + u*v^2 + u*v + u + u*v^-1 + v^2 + v + 1 + v^-2 + v^-3 + u^-1*v + u^-1*v^-2"),
    (15, 20, "u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-3 + u*v^2 + u*v + u*v^-1 + v + v^-1"),
    (15, 21, "u^5*v^-1 + u^4 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^2*v^2 + u^2*v^-3 + u + u*v^-1 + v^2 + v^-3 + u^-1*v + u^-1 + u^-1*v^-1 + u^-1*v^-2"),
    (16, 1, "u^3*v^-1 + u^3*v^-3 + u^2*v^-3 + u + u*v^-2 + u*v^-3 + v^-2 + v^-3 + u^-1 + u^-1*v^-1"),
    (16, 2, "v^-2 + v^-3 + u^-1*v^-2 + u^-1*v^-3"),
    (16, 3, "v^-2 + u^-1*v^-2"),
    (16, 4, "v^-1 + v^-2 + u^-1*v^-1 + u^-1*v^-2"),
    (16, 5, "v^-2 + v^-3 + u^-1*v^-2 + u^-1*v^-3"),
    (16, 6, "v^-1 + v^-3 + u^-1*v^-1 + u^-1*v^-3"),
    (16, 7, "v^-2 + u^-1*v^-2"),
    (16, 8, "v^-1 + v^-2 + u^-1*v^-1 + u^-1*v^-2"),
    (16, 9, "v^-2 + v^-3 + u^-1*v^-2 + u^-1*v^-3"),
    (16, 10, "v^-1 + v^-3 + u^-1*v^-1 + u^-1*v^-3"),
    (16, 11, "v^-2 + u^-1*v^-2"),
    (16, 12, "u^3*v^-1 + u^3*v^-3 + u^2*v^-1 + u + u*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-1 + v^-2 + v^-3 + u^-1*v^-1 + u^-1*v^-2 + u^-2 + u^-2*v^-1 + u^-2*v^-2 + u^-2*v^-3"),
    (16, 15, "v^-1 + v^-3 + u^-1*v^-1 + u^-1*v^-2 + u^-1*v^-3 + u^-2*v^-2"),
    (16, 16, "u^3*v^-1 + u^3*v^-3 + u^2*v^-1 + u + u*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-2 + u^-1*v^-3 + u^-2 + u^-2*v^-1 + u^-2*v^-3"),
    (16, 17, "u^3*v^-1 + u^3*v^-3 + u^2*v^-3 + u + u*v^-2 + u*v^-3 + v^-3 + u^-1*v^-2 + u^-1*v^-3 + u^-2 + u^-2*v^-1 + u^-2*v^-3"),
    (16, 18, "u^2*v^-1 + u^2*v^-3 + u*v^-1 + 1 + v^-3 + u^-1"),
    (16, 19, "u^3*v^-1 + u^3*v^-3 + u^2*v^-1 + u + u*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-2 + u^-1 + u^-1*v^-3 + u^-2 + u^-2*v^-1 + u^-2*v^-3"),
    (16, 20, "u^2*v^-1 + u^2*v^-3 + u*v^-1 + 1 + v^-2 + v^-3 + u^-1 + u^-1*v^-2"),
    (16, 21, "u^3*v^-1 + u^3*v^-3 + u^2*v^-3 + u + u*v^-2 + u*v^-3 + v^-3 + u^-1*v^-3 + u^-2 + u^-2*v^-1 + u^-2*v^-2 + u^-2*v^-3"),
    (17, 1, "u^5*v^-1 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2 + u^2*v^-3 + u*v^-1 + u*v^-3 + v + 1"),
    (17, 2, "u^2*v^-1 + u*v^-1 + u*v^-3 + v^-3"),
    (17, 3, "u + u*v^-1 + u*v^-2 + 1 + v^-1 + v^-2"),
    (17, 4, "u*v^-2 + v^-2"),
    (17, 5, "u^2*v^-1 + u*v^-1 + u*v^-3 + v^-3"),
    (17, 6, "u^2*v^-1 + u*v^-1 + u*v^-2 + u*v^-3 + v^-2 + v^-3"),
    (17, 7, "u^2*v^-1 + u + u*v^-1 + u*v^-2 + 1 + v^-2"),
    (17, 8, "u*v^-2 + v^-2"),
    (17, 9, "u*v^-1 + u*v^-3 + v^-1 + v^-3"),
    (17, 10, "u^2*v^-1 + u*v^-1 + u*v^-2 + u*v^-3 + v^-2 + v^-3"),
    (17, 11, "u + u*v^-1 + u*v^-2 + 1 + v^-1 + v^-2"),
    (17, 12, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^2 + u^2*v^-1 + u^2*v^-3 + u*v^-1 + u*v^-3 + v + v^-2 + u^-1*v^-3"),
    (17, 13, "1"),
    (17, 14, "u^2*v^-1 + v^-1"),
    (17, 15, "u*v^-2 + u*v^-3 + v^-3 + u^-1*v^-2"),
    (17, 16, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^2*v^-1 + u^2*v^-3 + u + u*v^-1 + u*v^-2 + v + v^-2 + v^-3 + u^-1*v^-2 + u^-1*v^-3"),
    (17, 17, "u^5*v^-1 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-3 + u*v + u*v^-1 + u*v^-2 + u*v^-3 + v + 1 + v^-3 + u^-1*v^-2 + u^-1*v^-3"),
    (17, 18, "u^4*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u*v + u*v^-1 + u*v^-2 + u*v^-3"),
    (17, 19, "u^5*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^2*v^-1 + u^2*v^-3 + u*v^-1 + u*v^-2 + v + v^-2 + v^-3 + u^-1*v^-2 + u^-1*v^-3"),
    (17, 20, "u^4*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u*v + u*v^-3 + v^-1 + v^-2"),
    (17, 21, "u^5*v^-1 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v + u^2*v^-3 + u*v + u*v^-2 + u*v^-3 + v + v^-1 + v^-2 + v^-3 + u^-1*v^-3"),
    (18, 14, "1"),
    (18, 18, "1"),
    (19, 1, "u^6*v^-1 + u^5 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-3 + u^3*v^2 + u^3 + u^3*v^-3 + u^2*v^2 + u^2 + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u + u*v^-1"),
    (19, 2, "u^3*v^-1 + u^2 + u^2*v^-3 + u + u*v^-2 + u*v^-3"),
    (19, 3, "u^2 + u^2*v^-1 + u^2*v^-2 + u*v + u + u*v^-2"),
    (19, 4, "u^2*v^-2 + u + u*v^-1 + u*v^-2"),
    (19, 5, "u^3*v^-1 + u^2 + u^2*v^-3 + u + u*v^-2 + u*v^-3"),
    (19, 6, "u^3*v^-1 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^-1 + u*v^-3"),
    (19, 7, "u^3*v^-1 + u^2*v^-2 + u*v + u + u*v^-2"),
    (19, 8, "u^2*v^-2 + u + u*v^-1 + u*v^-2"),
    (19, 9, "u^2*v^-1 + u^2*v^-3 + u + u*v^-2 + u*v^-3"),
    (19, 10, "u^3*v^-1 + u^2 + u^2*v^-2 + u^2*v^-3 + u*v^-1 + u*v^-3"),
    (19, 11, "u^2 + u^2*v^-1 + u^2*v^-2 + u*v + u + u*v^-2"),
    (19, 12, "u^6*v^-1 + u^5 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4*v^-1 + u^4*v^-2 + u^3*v^2 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v^-1 + u^2*v^-3 + u*v^2 + u*v + 1 + v^-1 + v^-2 + v^-3"),
    (19, 14, "u^3*v^-1 + u^2 + u^2*v^-1 + u"),
    (19, 15, "u^2*v^-2 + u^2*v^-3 + u + u*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-2"),
    (19, 16, "u^6*v^-1 + u^5 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4*v^-1 + u^4*v^-2 + u^3*v^2 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u + u*v^-1 + u*v^-2 + u*v^-3 + v^-1 + v^-3"),
    (19, 17, "u^6*v^-1 + u^5 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3*v^-3 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-3 + u*v^2 + u*v + u*v^-1 + u*v^-3 + v^-1 + v^-3"),
    (19, 18, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-2 + u^2*v^2 + u^2 + u^2*v^-3 + u*v"),
    (19, 19, "u^6*v^-1 + u^5 + u^5*v^-1 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4*v^-1 + u^4*v^-2 + u^3*v^2 + u^3 + u^3*v^-1 + u^3*v^-2 + u^3*v^-3 + u^2*v^2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u + u*v^-1 + u*v^-2 + u*v^-3 + v^-1 + v^-3"),
    (19, 20, "u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-2 + u^4*v^-3 + u^3*v + u^3 + u^3*v^-2 + u^2*v^2 + u^2 + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v + u + u*v^-2"),
    (19, 21, "u^6*v^-1 + u^5 + u^5*v^-2 + u^5*v^-3 + u^4*v + u^4 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3*v^-3 + u^2*v + u^2 + u^2*v^-3 + u*v^2 + u*v + u + u*v^-1 + u*v^-2 + u*v^-3 + 1 + v^-1 + v^-2 + v^-3"),
    (20, 1, "u^5 + u^5*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3 + u^3*v^-1 + u^2*v^2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u + u*v^-1"),
    (20, 2, "u^2 + u^2*v^-1 + u + u*v^-1 + u*v^-2 + u*v^-3"),
    (20, 3, "u*v + u*v^-2"),
    (20, 4, "u*v^-1 + u*v^-2"),
    (20, 5, "u^2 + u^2*v^-1 + u*v^-1 + u*v^-2 + u*v^-3"),
    (20, 6, "u^2 + u^2*v^-1 + u + u*v^-3"),
    (20, 7, "u^2 + u^2*v^-1 + u*v + u*v^-2"),
    (20, 8, "u*v^-1 + u*v^-2"),
    (20, 9, "u + u*v^-1 + u*v^-2 + u*v^-3"),
    (20, 10, "u^2 + u^2*v^-1 + u*v^-3"),
    (20, 11, "u*v + u + u*v^-2"),
    (20, 12, "u^5 + u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3*v^-3 + u^2*v^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-1 + u*v^-3 + v + 1 + v^-2 + v^-3"),
    (20, 13, "u"),
    (20, 14, "u^2 + u^2*v^-1"),
    (20, 15, "u + u*v^-1 + u*v^-3 + v^-2"),
    (20, 16, "u^5 + u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3*v^-3 + u^2*v^2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + v + 1 + v^-3"),
    (20, 17, "u^5 + u^5*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3 + u^3*v^-1 + u^2*v + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v + u*v^-2 + v + 1 + v^-3"),
    (20, 18, "u^4 + u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-3 + u^2*v^2 + u^2*v^-1 + u^2*v^-3 + u*v + v + 1"),
    (20, 19, "u^5 + u^5*v^-1 + u^4 + u^4*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3*v^-3 + u^2*v^2 + u^2*v + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + v^-3"),
    (20, 20, "u^4 + u^4*v^-1 + u^3 + u^3*v^-1 + u^3*v^-3 + u^2*v^2 + u^2*v^-1 + u^2*v^-3 + u*v + u + u*v^-2 + v + 1"),
    (20, 21, "u^5 + u^5*v^-1 + u^4*v^-3 + u^3*v^2 + u^3*v + u^3 + u^3*v^-1 + u^2*v + u^2*v^-1 + u^2*v^-2 + u^2*v^-3 + u*v^2 + u*v + u*v^-1 + u*v^-2 + v + 1 + v^-2 + v^-3"),
    (21, 1, "u^3 + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-1 + u^-1 + u^-1*v^-1 + u^-2*v^2 + u^-2*v + u^-2*v^-1 + u^-2*v^-2 + u^-3*v + u^-3"),
    (21, 2, "1 + u^-1 + u^-1*v^-2 + u^-2*v^-1 + u^-3*v^-1 + u^-3*v^-2"),
    (21, 3, "u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v + u^-2 + u^-3*v^-1"),
    (21, 4, "u^-1*v^-1 + u^-2 + u^-3 + u^-3*v^-1"),
    (21, 5, "1 + u^-1 + u^-1*v^-2 + u^-2*v^-1 + u^-3*v^-1 + u^-3*v^-2"),
    (21, 6, "1 + u^-1 + u^-1*v^-1 + u^-1*v^-2 + u^-2 + u^-2*v^-1 + u^-3 + u^-3*v^-2"),
    (21, 7, "1 + u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v + u^-3*v^-1"),
    (21, 8, "u^-1*v^-1 + u^-2 + u^-3 + u^-3*v^-1"),
    (21, 9, "u^-1 + u^-1*v^-2 + u^-2 + u^-2*v^-1 + u^-3*v^-1 + u^-3*v^-2"),
    (21, 10, "1 + u^-1 + u^-1*v^-1 + u^-1*v^-2 + u^-2 + u^-2*v^-1 + u^-3 + u^-3*v^-2"),
    (21, 11, "u^-1*v + u^-1 + u^-1*v^-1 + u^-2*v + u^-2 + u^-3*v^-1"),
    (21, 12, "u^3 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-2 + v^-2 + u^-1*v + u^-1*v^-1 + u^-2*v^2 + u^-2*v + u^-2 + u^-2*v^-2 + u^-3 + u^-3*v^-1 + u^-3*v^-2 + u^-4*v + u^-4 + u^-4*v^-1 + u^-4*v^-2"),
    (21, 13, "u^-2*v"),
    (21, 14, "1 + u^-2"),
    (21, 15, "u^-1*v^-1 + u^-1*v^-2 + u^-2 + u^-3*v + u^-3 + u^-3*v^-2 + u^-4*v^-1"),
    (21, 16, "u^3 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-2 + v + v^-2 + u^-1*v^-2 + u^-2*v^2 + u^-2*v + u^-2*v^-2 + u^-3*v + u^-3*v^-1 + u^-4*v + u^-4 + u^-4*v^-2"),
    (21, 17, "u^3 + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-1 + v^2 + v + u^-1*v^2 + u^-1 + u^-2*v^2 + u^-2*v + u^-4*v + u^-4 + u^-4*v^-2"),
    (21, 18, "u^2 + u*v^-1 + u*v^-2 + v^2 + 1 + v^-2 + u^-1*v^2 + u^-1*v + u^-1*v^-1 + u^-1*v^-2 + u^-2*v + u^-2*v^-2 + u^-3*v"),
    (21, 19, "u^3 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-2 + v + v^-2 + u^-1*v + u^-1*v^-2 + u^-2*v^2 + u^-2*v + u^-2*v^-2 + u^-3*v^-1 + u^-4*v + u^-4 + u^-4*v^-2"),
    (21, 20, "u^2 + u*v^-1 + u*v^-2 + v^2 + 1 + v^-2 + u^-1*v^2 + u^-1*v + u^-1 + u^-1*v^-2 + u^-2*v + u^-2 + u^-2*v^-2 + u^-3*v^-1"),
    (21, 21, "u^3 + u^2 + u^2*v^-1 + u^2*v^-2 + u*v^2 + u*v + u*v^-1 + v^2 + v + u^-1*v^2 + u^-2*v^2 + u^-2 + u^-2*v^-1 + u^-3*v + u^-4*v + u^-4 + u^-4*v^-1 + u^-4*v^-2"),
)
function _case008_entry()
    matches = filter(entry -> entry.id == CASE_ID, ToricBuilderCacheQBlocks.catalog().cases)
    length(matches) == 1 ||
        throw(ArgumentError("expected exactly one ToricBuilder cache Q-block entry for " * CASE_ID))
    return only(matches)
end

function _ring_generator_names(R)
    return Tuple(string.(gens(R)))
end

function _is_expected_uv_laurent_ring(R)::Bool
    try
        return Suslin._is_laurent_polynomial_ring(R) &&
            base_ring(R) == GF(2) &&
            _ring_generator_names(R) == ("u", "v")
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function _snapshot_matrix()
    R, (u, v) = Suslin.suslin_laurent_polynomial_ring(GF(2), ["u", "v"])
    return ToricBuilderCacheQBlocks._sparse_laurent_matrix(
        R,
        u,
        v,
        FIRST_FAILING_PEEL_DIMENSION,
        FIRST_FAILING_PEEL_DIMENSION,
        FAILING_INPUT_SPARSE_ENTRIES,
    )
end

function _matrix_size(M)
    return (nrows(M), ncols(M))
end

function _last_column(M)
    d = ncols(M)
    return [M[row, d] for row in 1:nrows(M)]
end

function _column_entries_are_over_ring(column, R)::Bool
    try
        return all(entry -> R(entry) == entry, column)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

function boundary_fixture()
    entry = _case008_entry()
    original_matrix = ToricBuilderCacheQBlocks.materialize_matrix(entry)
    failing_input_matrix = _snapshot_matrix()
    R = base_ring(failing_input_matrix)

    return (;
        case_id = entry.id,
        source_entry = entry,
        original_matrix,
        source_matrix_dimensions = entry.dimensions.matrix,
        source_sparse_entry_count = entry.sparse_entry_count,
        normalization_provenance = (;
            source_issue = 131,
            method = :normalize_laurent_gl_matrix,
            determinant_classification = :one,
            normalized_matrix_dimensions = entry.dimensions.matrix,
        ),
        passed_peel_dimensions = EXPECTED_PASSED_PEEL_DIMENSIONS,
        first_failing_peel_dimension = FIRST_FAILING_PEEL_DIMENSION,
        failing_input_matrix,
        failing_column = _last_column(failing_input_matrix),
        ring = R,
        ring_description = entry.ring.description,
        expected_diagnostic = (;
            status = :unsupported,
            failure_code = :unsupported_laurent_column_family,
        ),
    )
end

function _has_required_boundary_fields(fixture)::Bool
    return all(field -> hasproperty(fixture, field), REQUIRED_BOUNDARY_FIELDS)
end

function _diagnostic_matches_expected(column, R, expected)::Bool
    diagnostic = Suslin.diagnose_unimodular_column_reduction(column, R)
    return diagnostic.status == expected.status &&
        diagnostic.failure_code == expected.failure_code &&
        diagnostic.column_length == FIRST_FAILING_PEEL_DIMENSION
end

function validate_boundary_fixture(fixture)::Symbol
    _has_required_boundary_fields(fixture) || return :missing_metadata
    fixture.case_id == CASE_ID || return :wrong_case
    fixture.source_matrix_dimensions == (30, 30) || return :wrong_source_matrix_dimensions
    fixture.first_failing_peel_dimension == FIRST_FAILING_PEEL_DIMENSION || return :wrong_peel_dimension
    fixture.passed_peel_dimensions == EXPECTED_PASSED_PEEL_DIMENSIONS || return :wrong_passed_peel_dimensions
    _matrix_size(fixture.failing_input_matrix) == (FIRST_FAILING_PEEL_DIMENSION, FIRST_FAILING_PEEL_DIMENSION) ||
        return :wrong_failing_input_dimension
    length(fixture.failing_column) == FIRST_FAILING_PEEL_DIMENSION || return :wrong_column_length
    fixture.ring_description == EXPECTED_RING_DESCRIPTION || return :wrong_ring

    R = fixture.ring
    _is_expected_uv_laurent_ring(R) || return :wrong_ring
    base_ring(fixture.failing_input_matrix) == R || return :wrong_ring
    fixture.failing_column == _last_column(fixture.failing_input_matrix) || return :wrong_column
    _column_entries_are_over_ring(fixture.failing_column, R) || return :wrong_ring
    Suslin.is_unimodular_column(fixture.failing_column, R) || return :not_unimodular
    _diagnostic_matches_expected(fixture.failing_column, R, fixture.expected_diagnostic) ||
        return :wrong_diagnostic

    return :ok
end

function non_unimodular_negative_control(fixture = boundary_fixture())
    _, v = gens(fixture.ring)
    nonunit = v + one(fixture.ring)
    corrupted_matrix = copy(fixture.failing_input_matrix)
    d = FIRST_FAILING_PEEL_DIMENSION
    for row in 1:d
        corrupted_matrix[row, d] = nonunit * fixture.failing_column[row]
    end
    return merge(
        fixture,
        (;
            failing_input_matrix = corrupted_matrix,
            failing_column = _last_column(corrupted_matrix),
        ),
    )
end

end
