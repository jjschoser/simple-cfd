#include "catch2/catch.hpp"
#include "helpers.H"
#include "FluxSolver.H"

#include <cmath>

// Helper: test that two flux arrays agree component-by-component.
static void requireFluxEq(const std::array<REAL, Euler::NVARS>& got,
                           const std::array<REAL, Euler::NVARS>& expected,
                           REAL tol = 1e-12)
{
    for (int v = 0; v < Euler::NVARS; ++v)
    {
        INFO("variable index " << v);
        REQUIRE(got[v] == Approx(expected[v]).margin(tol));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared fixture: γ=1.4 EOS, Euler physics, HLLC solver.
// dx and dt are needed by the interface but HLLC ignores them.
// ─────────────────────────────────────────────────────────────────────────────
static const IdealGas g_eos(1.4);
static const Euler    g_euler(&g_eos);
static const HLLCSolver g_hllc(g_euler);
static const std::array<REAL, GRIDDIM> g_dx = {GRIDDIM_DECL(0.1, 0.1, 0.1)};
static const REAL g_dt = 0.01;


// ─────────────────────────────────────────────────────────────────────────────
// Property 1: consistency – F_HLLC(U, U) == F_Euler(U)
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("HLLC: consistency F(U,U) = F_Euler(U)", "[flux][hllc]")
{
    // Test several physically distinct states and all three flux dimensions.
    const auto states = {
        makeState(g_eos, 1.0,  0.3,  0.0, 0.0, 1.0),   // subsonic x-flow
        makeState(g_eos, 2.0, -0.5,  0.2, 0.0, 3.0),   // subsonic mixed flow
        makeState(g_eos, 0.125, 0.0, 0.0, 0.0, 0.1),   // Sod right state
        makeState(g_eos, 1.0,  0.0,  0.5, 0.0, 2.0),   // y-directed flow
        makeState(g_eos, 1.0,  0.0,  0.0, 0.7, 1.5),   // z-directed flow
    };

    for (int dim = 0; dim < GRIDDIM; ++dim)
    {
        for (const auto& U : states)
        {
            std::array<REAL, Euler::NVARS> F_hllc, F_euler;
            g_hllc(F_hllc, U, U, g_dx, g_dt, dim);
            g_euler.evalFluxFunc(F_euler, U, dim);
            INFO("dim=" << dim);
            requireFluxEq(F_hllc, F_euler);
        }
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// Property 2: contact-discontinuity preservation
//
// For two states with equal pressure and equal velocity the HLLC solver is
// exact.  Specifically, the contact speed equals the flow velocity:
//   • if u > 0 the interface is "left of contact" → F = F(U_Lo)
//   • if u < 0 the interface is "right of contact" → F = F(U_Hi)
//   • if u = 0 the contact is stationary → F = [0, p, 0…, 0]
//
// Derivation: equal p and u give s* = u from the contact-speed formula,
// and the resulting star state collapses to the corresponding side state.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("HLLC: contact preservation in x-direction", "[flux][hllc]")
{
    // States differ only in density; pressure and velocity are equal.
    const REAL u = 0.5, p = 1.0;
    auto ULo = makeState(g_eos, 1.0, u, 0.0, 0.0, p);
    auto UHi = makeState(g_eos, 2.0, u, 0.0, 0.0, p);

    SECTION("right-moving contact: F = F(U_Lo)")
    {
        std::array<REAL, Euler::NVARS> F_hllc, F_expected;
        g_hllc(F_hllc, ULo, UHi, g_dx, g_dt, 0);
        g_euler.evalFluxFunc(F_expected, ULo, 0);
        requireFluxEq(F_hllc, F_expected);
    }

    SECTION("left-moving contact: F = F(U_Hi)")
    {
        auto ULo2 = makeState(g_eos, 1.0, -u, 0.0, 0.0, p);
        auto UHi2 = makeState(g_eos, 2.0, -u, 0.0, 0.0, p);

        std::array<REAL, Euler::NVARS> F_hllc, F_expected;
        g_hllc(F_hllc, ULo2, UHi2, g_dx, g_dt, 0);
        g_euler.evalFluxFunc(F_expected, UHi2, 0);
        requireFluxEq(F_hllc, F_expected);
    }

    SECTION("stationary contact: mass flux = 0, normal pressure flux = p, energy flux = 0")
    {
        auto ULo3 = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 2.0);
        auto UHi3 = makeState(g_eos, 3.0, 0.0, 0.0, 0.0, 2.0);

        std::array<REAL, Euler::NVARS> F;
        g_hllc(F, ULo3, UHi3, g_dx, g_dt, 0);

        REQUIRE(F[Euler::RHO]    == Approx(0.0).margin(1e-12)); // no mass flux
        REQUIRE(F[Euler::MOM[0]] == Approx(2.0).margin(1e-12)); // pressure = 2
        REQUIRE(F[Euler::ENE]    == Approx(0.0).margin(1e-12)); // no energy flux
    }
}

TEST_CASE("HLLC: contact preservation in y-direction", "[flux][hllc]")
{
    const REAL v = 0.4, p = 1.5;
    auto ULo = makeState(g_eos, 1.0, 0.0, v, 0.0, p);
    auto UHi = makeState(g_eos, 2.5, 0.0, v, 0.0, p);

    std::array<REAL, Euler::NVARS> F_hllc, F_expected;
    g_hllc(F_hllc, ULo, UHi, g_dx, g_dt, 1);   // dim=1 (y-direction)
    g_euler.evalFluxFunc(F_expected, ULo, 1);
    requireFluxEq(F_hllc, F_expected);
}

TEST_CASE("HLLC: contact preservation in z-direction", "[flux][hllc]")
{
    const REAL w = 0.6, p = 2.0;
    auto ULo = makeState(g_eos, 1.0, 0.0, 0.0, w, p);
    auto UHi = makeState(g_eos, 3.0, 0.0, 0.0, w, p);

    std::array<REAL, Euler::NVARS> F_hllc, F_expected;
    g_hllc(F_hllc, ULo, UHi, g_dx, g_dt, 2);   // dim=2 (z-direction)
    g_euler.evalFluxFunc(F_expected, ULo, 2);
    requireFluxEq(F_hllc, F_expected);
}


// ─────────────────────────────────────────────────────────────────────────────
// Property 3: Sod initial states – flux is in the correct direction
//
// For the Sod left/right states at t=0 the flow is at rest (u=0) and the
// contact speed s* > 0, so the HLLC flux is evaluated from the left-star state.
// The mass flux must be positive (left→right), and the pressure flux must
// exceed both side pressures (driven by the shock).
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("HLLC: Sod initial states – qualitative checks", "[flux][hllc]")
{
    auto ULo = makeState(g_eos, 1.0,   0.0, 0.0, 0.0, 1.0);
    auto UHi = makeState(g_eos, 0.125, 0.0, 0.0, 0.0, 0.1);

    std::array<REAL, Euler::NVARS> F;
    g_hllc(F, ULo, UHi, g_dx, g_dt, 0);

    // Mass flux is positive (net flow from high- to low-pressure side)
    REQUIRE(F[Euler::RHO] > 0.0);

    // Normal-momentum flux lies between the two side pressures.
    // The contact pressure p* ≈ 0.30, so the numerical flux is above the
    // right-state pressure (0.1) and below the left-state pressure (1.0).
    REQUIRE(F[Euler::MOM[0]] > 0.1);
    REQUIRE(F[Euler::MOM[0]] < 1.0);

    // Energy flux is positive
    REQUIRE(F[Euler::ENE] > 0.0);
}


// ─────────────────────────────────────────────────────────────────────────────
// Lax-Friedrichs: consistency (same as HLLC above but for the other solver)
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Lax-Friedrichs: consistency F(U,U) = F_Euler(U)", "[flux][lax]")
{
    const LaxFriedrichsSolver lax(g_euler);

    auto U = makeState(g_eos, 1.0, 0.3, 0.0, 0.0, 1.0);
    std::array<REAL, Euler::NVARS> F_lax, F_euler;
    lax(F_lax, U, U, g_dx, g_dt, 0);
    g_euler.evalFluxFunc(F_euler, U, 0);
    requireFluxEq(F_lax, F_euler);
}
