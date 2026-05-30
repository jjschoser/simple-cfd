#include "catch2/catch.hpp"
#include "helpers.H"
#include "Reconstruction.H"

#include <cmath>

static const IdealGas g_eos(1.4);
static const Euler    g_euler(&g_eos);

// Standard slab geometry: unit cells in all directions.
static const std::array<REAL, GRIDDIM> g_dx = {GRIDDIM_DECL(1.0, 1.0, 1.0)};
static const REAL g_dt = 0.1;
static const int  g_dim = 0;   // reconstruct in x-direction

// ─────────────────────────────────────────────────────────────────────────────
// NoReconstruction: always returns U on both sides regardless of neighbours
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("NoReconstruction: always outputs U on both faces", "[recon][norecon]")
{
    NoReconstruction recon;

    auto U   = makeState(g_eos, 1.0, 0.5, 0.0, 0.0, 2.0);
    auto Nbr = makeState(g_eos, 2.0, 0.0, 0.0, 0.0, 0.5);  // very different

    std::array<REAL, Euler::NVARS> Lo, Hi;
    recon(Lo, Hi, Nbr, Nbr, U, g_dx, g_dt, g_dim);

    for (int v = 0; v < Euler::NVARS; ++v)
    {
        REQUIRE(Lo[v] == Approx(U[v]).margin(1e-15));
        REQUIRE(Hi[v] == Approx(U[v]).margin(1e-15));
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// MUSCLHancock: spatial reconstruction properties
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("MUSCLHancock: constant field → no slope, faces equal cell centre",
          "[recon][muscl]")
{
    MUSCLHancock recon(g_euler);
    auto U = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 1.0);

    std::array<REAL, Euler::NVARS> Lo, Hi;
    recon(Lo, Hi, U, U, U, g_dx, g_dt, g_dim);

    for (int v = 0; v < Euler::NVARS; ++v)
    {
        REQUIRE(Lo[v] == Approx(U[v]).margin(1e-14));
        REQUIRE(Hi[v] == Approx(U[v]).margin(1e-14));
    }
}

TEST_CASE("MUSCLHancock: linear density profile → exact face values", "[recon][muscl]")
{
    // rho varies linearly (0.9, 1.0, 1.1) with step dx=1.
    // Pressure is uniform so the energy variable is also constant (E = p/(γ-1)).
    // With r = 1 (uniform gradient) the limiter gives φ=1, so:
    //   delta_rho = 0.5 * 1 * (1.1 - 0.9) = 0.1
    //   face_lo   = 1.0 - 0.05 = 0.95   (exact for a linear field)
    //   face_hi   = 1.0 + 0.05 = 1.05
    // The half-step predictor leaves the faces unchanged here because with u=0
    // the flux difference (ΔF) across a uniform-pressure, zero-velocity field
    // evaluates to zero.
    MUSCLHancock recon(g_euler);

    const REAL p = 1.0;
    auto UNbrLo = makeState(g_eos, 0.9, 0.0, 0.0, 0.0, p);
    auto U      = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, p);
    auto UNbrHi = makeState(g_eos, 1.1, 0.0, 0.0, 0.0, p);

    std::array<REAL, Euler::NVARS> Lo, Hi;
    recon(Lo, Hi, UNbrLo, UNbrHi, U, g_dx, g_dt, g_dim);

    // Density face values are exact for a linear profile.
    REQUIRE(Lo[Euler::RHO] == Approx(0.95).margin(1e-14));
    REQUIRE(Hi[Euler::RHO] == Approx(1.05).margin(1e-14));

    // Pressure is recoverable from conserved variables: p = (γ-1)·ρ·e
    const REAL p_lo = g_eos.getPressure(Lo[Euler::RHO],
                                        g_euler.getSpecificInternalEnergy(Lo));
    const REAL p_hi = g_eos.getPressure(Hi[Euler::RHO],
                                        g_euler.getSpecificInternalEnergy(Hi));
    REQUIRE(p_lo == Approx(p).epsilon(1e-12));
    REQUIRE(p_hi == Approx(p).epsilon(1e-12));
}

TEST_CASE("MUSCLHancock: local density maximum → limiter suppresses slope",
          "[recon][muscl]")
{
    // U is a local maximum in density: both neighbours are lower.
    // The ratio r = (U - UNbrLo) / (UNbrHi - U) has the wrong sign → φ = 0.
    // The half-step then also produces zero delta (F(U) - F(U) = 0).
    MUSCLHancock recon(g_euler);

    auto UNbrLo = makeState(g_eos, 0.5, 0.0, 0.0, 0.0, 1.0);
    auto U      = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 1.0);
    auto UNbrHi = makeState(g_eos, 0.5, 0.0, 0.0, 0.0, 1.0);

    std::array<REAL, Euler::NVARS> Lo, Hi;
    recon(Lo, Hi, UNbrLo, UNbrHi, U, g_dx, g_dt, g_dim);

    for (int v = 0; v < Euler::NVARS; ++v)
    {
        REQUIRE(Lo[v] == Approx(U[v]).margin(1e-14));
        REQUIRE(Hi[v] == Approx(U[v]).margin(1e-14));
    }
}

TEST_CASE("MUSCLHancock: step in density → limiter suppresses slope",
          "[recon][muscl]")
{
    // The left neighbour equals the centre cell (r=0 numerator → φ=0).
    // This models the first interior cell next to a step discontinuity.
    MUSCLHancock recon(g_euler);

    auto UNbrLo = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 1.0);
    auto U      = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 1.0);
    auto UNbrHi = makeState(g_eos, 0.125, 0.0, 0.0, 0.0, 0.1);

    std::array<REAL, Euler::NVARS> Lo, Hi;
    recon(Lo, Hi, UNbrLo, UNbrHi, U, g_dx, g_dt, g_dim);

    for (int v = 0; v < Euler::NVARS; ++v)
    {
        REQUIRE(Lo[v] == Approx(U[v]).margin(1e-14));
        REQUIRE(Hi[v] == Approx(U[v]).margin(1e-14));
    }
}

TEST_CASE("MUSCLHancock: slope limiter is symmetric (φ(r) = φ(1/r))", "[recon][muscl]")
{
    // For a symmetric limiter φ(r) = φ(1/r) the reconstruction of a mirrored
    // profile gives the same magnitude of delta.  Here the density gradient has
    // ratio r=2 in one case and r=0.5 in the mirrored case; the face values
    // must satisfy |Lo-U| = |Hi-U| in each case.
    MUSCLHancock recon(g_euler);

    // r = 2  (gradient steepens toward Hi neighbour)
    auto UNbrLo_a = makeState(g_eos, 0.0, 0.0, 0.0, 0.0, 1.0);
    auto U_a      = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 1.0);
    auto UNbrHi_a = makeState(g_eos, 1.5, 0.0, 0.0, 0.0, 1.0);

    // r = 0.5  (mirror: gradient steepens toward Lo neighbour)
    auto UNbrLo_b = makeState(g_eos, 0.5, 0.0, 0.0, 0.0, 1.0);
    auto U_b      = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 1.0);
    auto UNbrHi_b = makeState(g_eos, 2.0, 0.0, 0.0, 0.0, 1.0);

    std::array<REAL, Euler::NVARS> Lo_a, Hi_a, Lo_b, Hi_b;
    recon(Lo_a, Hi_a, UNbrLo_a, UNbrHi_a, U_a, g_dx, g_dt, g_dim);
    recon(Lo_b, Hi_b, UNbrLo_b, UNbrHi_b, U_b, g_dx, g_dt, g_dim);

    // |face - centre| in density should be equal in both configurations
    REAL delta_a = std::fabs(Hi_a[Euler::RHO] - U_a[Euler::RHO]);
    REAL delta_b = std::fabs(Hi_b[Euler::RHO] - U_b[Euler::RHO]);
    REQUIRE(delta_a == Approx(delta_b).epsilon(1e-12));
}

TEST_CASE("MUSCLHancock: reconstructed faces bracket the cell average", "[recon][muscl]")
{
    // For any monotone stencil the TVD property guarantees that the
    // face values don't introduce new extrema.  Face values must lie between
    // the minimum and maximum of the three-cell stencil.
    MUSCLHancock recon(g_euler);

    auto UNbrLo = makeState(g_eos, 0.8, 0.0, 0.0, 0.0, 1.0);
    auto U      = makeState(g_eos, 1.0, 0.0, 0.0, 0.0, 1.0);
    auto UNbrHi = makeState(g_eos, 1.3, 0.0, 0.0, 0.0, 1.0);

    std::array<REAL, Euler::NVARS> Lo, Hi;
    recon(Lo, Hi, UNbrLo, UNbrHi, U, g_dx, g_dt, g_dim);

    const REAL rho_min = 0.8, rho_max = 1.3;
    REQUIRE(Lo[Euler::RHO] >= rho_min - 1e-12);
    REQUIRE(Lo[Euler::RHO] <= rho_max + 1e-12);
    REQUIRE(Hi[Euler::RHO] >= rho_min - 1e-12);
    REQUIRE(Hi[Euler::RHO] <= rho_max + 1e-12);
}

TEST_CASE("MUSCLHancock: half-step predictor advances states in time", "[recon][muscl]")
{
    // For a uniform state with non-zero velocity the half-step should produce
    // zero update (ΔF = 0 since F(Lo)=F(Hi) for equal reconstructed states).
    // For a non-trivial gradient the predictor must differ from the spatial-only
    // result.  Here we simply verify that the half-step does not corrupt
    // positivity: density and internal energy must remain positive.
    MUSCLHancock recon(g_euler);

    auto UNbrLo = makeState(g_eos, 0.9, 0.5, 0.0, 0.0, 1.0);
    auto U      = makeState(g_eos, 1.0, 0.5, 0.0, 0.0, 1.0);
    auto UNbrHi = makeState(g_eos, 1.1, 0.5, 0.0, 0.0, 1.0);

    std::array<REAL, Euler::NVARS> Lo, Hi;
    recon(Lo, Hi, UNbrLo, UNbrHi, U, g_dx, g_dt, g_dim);

    REQUIRE(Lo[Euler::RHO] > 0.0);
    REQUIRE(Hi[Euler::RHO] > 0.0);
    REQUIRE(g_euler.getSpecificInternalEnergy(Lo) > 0.0);
    REQUIRE(g_euler.getSpecificInternalEnergy(Hi) > 0.0);
}
