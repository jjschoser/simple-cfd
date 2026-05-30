#pragma once
#include "Euler.H"
#include "EquationOfState.H"

// Build a conserved-variable state from primitive variables.
// Works for any GRIDDIM / SPACEDIM because it always fills all SPACEDIM
// momentum slots and uses Euler::MOM / Euler::ENE indices.
inline std::array<REAL, Euler::NVARS>
makeState(const IdealGas& eos, REAL rho, REAL ux, REAL uy, REAL uz, REAL p)
{
    std::array<REAL, Euler::NVARS> U{};
    U[Euler::RHO]    = rho;
    U[Euler::MOM[0]] = rho * ux;
    #if SPACEDIM >= 2
        U[Euler::MOM[1]] = rho * uy;
    #endif
    #if SPACEDIM >= 3
        U[Euler::MOM[2]] = rho * uz;
    #endif
    U[Euler::ENE] = eos.getSpecificInternalEnergy(rho, p) * rho
                  + 0.5 * rho * (ux*ux + uy*uy + uz*uz);
    return U;
}
