using Test
using Suslin
using Oscar

const QUILLEN_PATCH_CATALOG_PATH = joinpath(@__DIR__, "..", "fixtures", "quillen_patch_cases.jl")
if !isdefined(Main, :QuillenPatchFixtureCatalog)
    include(QUILLEN_PATCH_CATALOG_PATH)
end

function _cover_inputs(entry)
    return (
        [data.denominator for data in entry.denominator_data],
        [data.coverage_multiplier for data in entry.denominator_data],
    )
end

@testset "Quillen denominator cover certificates" begin
    entries = Main.QuillenPatchFixtureCatalog.cases_by_id()
    for id in ("quillen-two-open-cover-qq", "quillen-nontrivial-multipliers-qq")
        entry = entries[id]
        R = entry.ring.object
        denominators, multipliers = _cover_inputs(entry)
        certificate = Suslin.quillen_denominator_cover_certificate(R, denominators, multipliers)
        replay = Suslin.replay_quillen_denominator_cover(certificate)

        @test Suslin.verify_quillen_denominator_cover(certificate)
        @test certificate.coverage_sum == one(R)
        @test replay.coverage_sum == certificate.coverage_sum
        @test replay.coverage_ok
        @test certificate.denominators == [R(denominator) for denominator in denominators]
        @test certificate.coverage_multipliers == [R(multiplier) for multiplier in multipliers]
        @test all(denominator -> parent(denominator) == R, certificate.denominators)
        @test all(multiplier -> parent(multiplier) == R, certificate.coverage_multipliers)
    end

    two_open = entries["quillen-two-open-cover-qq"]
    R_two = two_open.ring.object
    r = two_open.denominator_data[1].denominator
    @test [data.denominator for data in two_open.denominator_data] == [r, one(R_two) - r]

    nontrivial = entries["quillen-nontrivial-multipliers-qq"]
    R_nontrivial = nontrivial.ring.object
    @test any(data -> data.coverage_multiplier != one(R_nontrivial), nontrivial.denominator_data)

    denominators, multipliers = _cover_inputs(two_open)
    @test_throws ArgumentError Suslin.quillen_denominator_cover_certificate(
        R_two,
        denominators[1:1],
        multipliers[1:1],
    )
    @test_throws ArgumentError Suslin.quillen_denominator_cover_certificate(
        R_two,
        denominators,
        [one(R_two), zero(R_two)],
    )

    RR, (Y, s) = Oscar.polynomial_ring(RealField(), ["Y", "s"])
    @test_throws ArgumentError Suslin.quillen_denominator_cover_certificate(
        RR,
        [s, one(RR) - s],
        [one(RR), one(RR)],
    )

    certificate = Suslin.quillen_denominator_cover_certificate(R_two, denominators, multipliers)
    tampered = Suslin.QuillenDenominatorCoverCertificate(
        certificate.ring,
        certificate.denominators,
        [one(R_two), zero(R_two)],
        certificate.coverage_sum,
        certificate.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover(tampered)

    malformed_replay = Suslin.QuillenDenominatorCoverCertificate(
        QQ,
        certificate.denominators,
        certificate.coverage_multipliers,
        certificate.coverage_sum,
        certificate.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover(malformed_replay)

    R_exact, (_, r_exact) = Oscar.polynomial_ring(QQ, ["X", "r"])
    certificate_exact = Suslin.quillen_denominator_cover_certificate(
        R_exact,
        [r_exact, one(R_exact) - r_exact],
        [one(R_exact), one(R_exact)],
    )
    tampered_exact = Suslin.QuillenDenominatorCoverCertificate(
        certificate_exact.ring,
        [2 * r_exact, one(R_exact) - r_exact],
        [R_exact(1//2), one(R_exact)],
        certificate_exact.coverage_sum,
        certificate_exact.verification,
    )
    @test !Suslin.verify_quillen_denominator_cover(tampered_exact)
end
