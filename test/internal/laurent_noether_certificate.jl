using Test
using Suslin
using Oscar

function _laurent_noether_test_certificate()
    R, (x, y) = suslin_laurent_polynomial_ring(QQ, ["x", "y"])
    column = [x^2 + x * y^-1 + y^2, one(R) + x^-1 * y]
    return R, x, y, column, Suslin._laurent_noether_certificate(column, 1, x)
end

function _laurent_noether_tamper(certificate; kwargs...)
    fields = propertynames(certificate)
    values = map(field -> get(kwargs, field, getproperty(certificate, field)), fields)
    return Suslin.LaurentNoetherCertificate(values...)
end

function _assert_laurent_noether_replay(R, column, certificate)
    @test certificate.ring === R
    @test certificate.original_column == tuple(column...)
    @test certificate.validation_status == :ok
    @test Suslin._validate_laurent_noether_certificate(certificate) == :ok

    forward_values = Suslin._laurent_noether_substitution_values(certificate.forward_substitution)
    inverse_values = Suslin._laurent_noether_substitution_values(certificate.inverse_substitution)
    for entry in column
        @test evaluate(evaluate(entry, collect(forward_values)), collect(inverse_values)) == entry
        @test evaluate(evaluate(entry, collect(inverse_values)), collect(forward_values)) == entry
    end
    for entry in certificate.transformed_column
        @test evaluate(evaluate(entry, collect(inverse_values)), collect(forward_values)) == entry
    end
    @test certificate.transformed_column == tuple((evaluate(entry, collect(forward_values)) for entry in column)...)
    @test certificate.replayed_selected_entry == certificate.transformed_column[certificate.selected_entry_index]
    @test certificate.leading_coefficient ==
        Suslin._laurent_noether_leading_coefficient(certificate.replayed_selected_entry, R, certificate.selected_generator_index)
    @test certificate.trailing_coefficient ==
        Suslin._laurent_noether_trailing_coefficient(certificate.replayed_selected_entry, R, certificate.selected_generator_index)
    @test certificate.leading_coefficient_is_unit
    @test certificate.trailing_coefficient_is_unit
end

@testset "Laurent Noether variable-change certificate" begin
    R, x, y, column, certificate = _laurent_noether_test_certificate()

    @test certificate.selected_entry_index == 1
    @test certificate.selected_generator_index == 1
    @test certificate.selected_generator == x
    @test certificate.other_generator_index == 2
    _assert_laurent_noether_replay(R, column, certificate)

    varied_column = [x^-3 * y^4 + x * y^-2 + y^3, one(R) + x^2 * y^-3]
    varied_certificate = Suslin._laurent_noether_certificate(varied_column, 1, x)
    @test varied_certificate.noether_power == 5
    _assert_laurent_noether_replay(R, varied_column, varied_certificate)

    S5, (u5, v5) = suslin_laurent_polynomial_ring(GF(5), ["u", "v"])
    second_generator_column = [v5^2 + u5^-1 * v5 + u5^3 * v5^-1, one(S5) + u5 * v5^-2]
    second_generator_certificate =
        Suslin._laurent_noether_certificate(second_generator_column, 1, v5)
    @test second_generator_certificate.selected_generator_index == 2
    @test second_generator_certificate.selected_generator == v5
    @test second_generator_certificate.other_generator_index == 1
    _assert_laurent_noether_replay(S5, second_generator_column, second_generator_certificate)

    @test Suslin._validate_laurent_noether_certificate(_laurent_noether_tamper(
        certificate;
        noether_power = certificate.noether_power + 1,
    )) != :ok
    @test Suslin._validate_laurent_noether_certificate(_laurent_noether_tamper(
        certificate;
        inverse_substitution = certificate.forward_substitution,
    )) != :ok
    noninvertible_forward = (
        certificate.forward_substitution[1],
        merge(certificate.forward_substitution[2], (; value = certificate.forward_substitution[1].value)),
    )
    @test Suslin._validate_laurent_noether_certificate(_laurent_noether_tamper(
        certificate;
        forward_substitution = noninvertible_forward,
    )) != :ok
    @test Suslin._validate_laurent_noether_certificate(_laurent_noether_tamper(
        certificate;
        transformed_column = tuple(certificate.original_column...),
    )) != :ok
    @test Suslin._validate_laurent_noether_certificate(_laurent_noether_tamper(
        certificate;
        selected_generator_index = 2,
        selected_generator = y,
        other_generator_index = 1,
    )) != :ok

    S, (u, v) = suslin_laurent_polynomial_ring(QQ, ["u", "v"])
    @test Suslin._validate_laurent_noether_certificate(_laurent_noether_tamper(
        certificate;
        ring = S,
    )) != :ok
    @test Suslin._validate_laurent_noether_certificate(_laurent_noether_tamper(
        certificate;
        selected_entry_index = length(column) + 1,
    )) != :ok

    @test_throws ArgumentError Suslin._laurent_noether_certificate(column, 0, x)
    @test_throws ArgumentError Suslin._laurent_noether_certificate([zero(R), one(R)], 1, x)
    @test_throws ArgumentError Suslin._laurent_noether_certificate(column, 1, y + one(R))
    P, (t,) = suslin_polynomial_ring(QQ, ["t"])
    @test_throws ArgumentError Suslin._laurent_noether_certificate([t], 1, t)
    T, (a, b, c) = suslin_laurent_polynomial_ring(QQ, ["a", "b", "c"])
    @test_throws ArgumentError Suslin._laurent_noether_certificate([a + b + c], 1, a)
end
