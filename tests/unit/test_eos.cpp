#include "catch2/catch.hpp"
#include "EquationOfState.H"

#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// IdealGas: p = (γ-1) ρ e,  e = p / ((γ-1) ρ),  c = √(γ p / ρ)
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("IdealGas: getPressure formula", "[eos]")
{
    const REAL gamma = 1.4;
    IdealGas eos(gamma);

    // p = (γ-1)·ρ·e
    REQUIRE(eos.getPressure(1.0, 1.0) == Approx((gamma - 1.0) * 1.0 * 1.0));
    REQUIRE(eos.getPressure(2.0, 0.5) == Approx((gamma - 1.0) * 2.0 * 0.5));
    REQUIRE(eos.getPressure(0.125, 0.8) == Approx((gamma - 1.0) * 0.125 * 0.8));
    REQUIRE(eos.getPressure(3.0, 2.0) == Approx((gamma - 1.0) * 3.0 * 2.0));
}

TEST_CASE("IdealGas: getSpecificInternalEnergy formula", "[eos]")
{
    const REAL gamma = 1.4;
    IdealGas eos(gamma);

    // e = p / ((γ-1) ρ)
    REQUIRE(eos.getSpecificInternalEnergy(1.0, 1.0) == Approx(1.0 / (gamma - 1.0)));
    REQUIRE(eos.getSpecificInternalEnergy(2.0, 0.4) == Approx(0.4 / ((gamma - 1.0) * 2.0)));
    REQUIRE(eos.getSpecificInternalEnergy(0.5, 2.0) == Approx(2.0 / ((gamma - 1.0) * 0.5)));
}

TEST_CASE("IdealGas: getSoundSpeed formula", "[eos]")
{
    const REAL gamma = 1.4;
    IdealGas eos(gamma);

    // c = √(γ p / ρ)
    REQUIRE(eos.getSoundSpeed(1.0, 1.0) == Approx(std::sqrt(gamma)));
    REQUIRE(eos.getSoundSpeed(1.0, 0.1) == Approx(std::sqrt(gamma * 0.1)));
    REQUIRE(eos.getSoundSpeed(2.0, 3.0) == Approx(std::sqrt(gamma * 3.0 / 2.0)));
    REQUIRE(eos.getSoundSpeed(0.125, 0.1) == Approx(std::sqrt(gamma * 0.1 / 0.125)));
}

TEST_CASE("IdealGas: pressure–energy round-trip", "[eos]")
{
    // getPressure(ρ, getSpecificInternalEnergy(ρ, p)) == p
    for (REAL gamma : {1.2, 1.4, 5.0 / 3.0})
    {
        IdealGas eos(gamma);
        for (REAL rho : {0.125, 0.5, 1.0, 2.0})
        {
            for (REAL p : {0.1, 1.0, 5.0, 101325.0})
            {
                REAL e = eos.getSpecificInternalEnergy(rho, p);
                REQUIRE(eos.getPressure(rho, e) == Approx(p).epsilon(1e-12));
            }
        }
    }
}

TEST_CASE("IdealGas: getGamma returns constructor value", "[eos]")
{
    REQUIRE(IdealGas(1.4).getGamma()        == Approx(1.4));
    REQUIRE(IdealGas(1.2).getGamma()        == Approx(1.2));
    REQUIRE(IdealGas(5.0 / 3.0).getGamma() == Approx(5.0 / 3.0));
}

TEST_CASE("IdealGas: sound speed consistent with pressure", "[eos]")
{
    // c² = γ p / ρ = γ (γ-1) e, which means c = √(γ(γ-1)e)
    IdealGas eos(1.4);
    const REAL rho = 1.5, p = 2.0;
    const REAL e   = eos.getSpecificInternalEnergy(rho, p);
    const REAL c   = eos.getSoundSpeed(rho, p);
    REQUIRE(c * c == Approx(1.4 * p / rho).epsilon(1e-12));
    REQUIRE(c * c == Approx(1.4 * (1.4 - 1.0) * e).epsilon(1e-12));
}
