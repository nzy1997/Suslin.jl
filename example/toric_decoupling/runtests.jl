using Test

include("try_column_block_decoupling.jl")

@testset "optional ToricBuilder cache smoke gate" begin
    @testset "missing local cache reports stable CACHE_ERROR" begin
        io = IOBuffer()
        rows = run_smoke(
            ["case_013"];
            toricbuilder_dir = "/definitely/missing",
            cache_dir = "/definitely/missing",
            io,
        )

        @test only(rows).case == "case_013"
        @test only(rows).status == "FAIL"
        @test only(rows).failure == "CACHE_ERROR"
        @test occursin("TORIC_SMOKE case=case_013", String(take!(io)))
    end

    @testset "live local cache exercises selected blocks" begin
        if local_toricbuilder_available()
            io = IOBuffer()
            rows = run_smoke(["case_001", "case_004"]; io)

            @test any(row -> row.block == "column_Q", rows)
            @test any(row -> row.block == "pair_mix_2_1", rows)
            pair_mix_rows = filter(row -> row.block == "pair_mix_2_1", rows)
            @test !isempty(pair_mix_rows)
            @test all(row -> row.sl_core == "SL_CORE_PASS", pair_mix_rows)
            @test all(row -> row.verified == "true", pair_mix_rows)
            @test occursin("TORIC_SMOKE case=case_001", String(take!(io)))
        end
    end
end
