#include "Euler.H"

#include <cmath>

Euler::Euler(const EquationOfState* const eos) : m_eos(eos)
{
    
}

const EquationOfState* Euler::getEquationOfState() const
{
    return m_eos;
}

REAL Euler::getSpecificInternalEnergy(const std::array<REAL, NVARS>& U) const
{
    const REAL rho = U[RHO];
    REAL momMag2 = 0.0;
    for(int d = 0; d < SPACEDIM; ++d)
    {
        momMag2 += U[MOM[d]] * U[MOM[d]];
    }
    return (U[ENE] - 0.5 * momMag2 / rho) / rho;
}

REAL Euler::getTotalEnergy(const REAL rho, const std::array<REAL, SPACEDIM>& vel, const REAL p) const
{
    REAL velMag2 = 0.0;
    for(int d = 0; d < SPACEDIM; ++d)
    {
        velMag2 += vel[d] * vel[d];
    }
    const REAL e = m_eos->getSpecificInternalEnergy(rho, p);
    return rho * (e + 0.5 * velMag2);
}

void Euler::evalFluxFunc(std::array<REAL, NVARS>& F, const std::array<REAL, NVARS>& U, const int dim) const
{
    const REAL rho = U[RHO];
    const REAL e = getSpecificInternalEnergy(U);
    const REAL p = m_eos->getPressure(rho, e);
    const REAL velDim = U[MOM[dim]] / rho;
    F[RHO] = U[MOM[dim]];
    for(int d = 0; d < SPACEDIM; ++d)
    {
        F[MOM[d]] = U[MOM[d]] * velDim;
    }
    F[MOM[dim]] += p;
    F[ENE] = (U[ENE] + p) * velDim;
}

REAL Euler::getMaxWaveSpeed(const std::array<REAL, NVARS>& U) const
{
    const REAL rho = U[RHO];
    const REAL e = getSpecificInternalEnergy(U);
    const REAL p = m_eos->getPressure(rho, e);
    REAL velMag2 = 0.0;
    for(int d = 0; d < SPACEDIM; ++d)
    {
        velMag2 += U[MOM[d]] * U[MOM[d]] / (rho * rho);
    }
    return m_eos->getSoundSpeed(rho, p) + std::sqrt(velMag2);
}

#ifdef USE_GRAVITY
void Euler::getGravitySource(std::array<REAL, NVARS>& S, const std::array<REAL, NVARS>& U, const std::array<REAL, SPACEDIM>& g)
{
    S[ENE] = 0.0;
    for(int d = 0; d < SPACEDIM; ++d)
    {
        S[MOM[d]] = U[RHO] * g[d];
        S[ENE] += U[MOM[d]] * g[d];
    }
}

void Euler::addGravityToGhostState(std::array<REAL, NVARS>& U, const std::array<REAL, GRIDDIM>& probePos, 
                                   const std::array<REAL, GRIDDIM>& ghostPos, const std::array<REAL, SPACEDIM>& g) const
{
    const int dim = std::min(GRIDDIM, SPACEDIM);
    REAL gDotD = 0.0;
    for(int d = 0; d < dim; ++d)
    {
        gDotD += g[d] * (ghostPos[d] - probePos[d]);
    }
    const REAL rho = U[RHO];
    const REAL e = getSpecificInternalEnergy(U);
    const REAL p = m_eos->getPressure(rho, e);
    const REAL deltaP = rho * gDotD;
    U[ENE] += rho * (m_eos->getSpecificInternalEnergy(rho, p + deltaP) - e);
}
#endif

#if GRIDDIM == 3 && defined(USE_VDB)
void Euler::getOutState(std::array<std::vector<REAL>, NOUTVARS>& outVars, const std::array<REAL, NVARS>& U) const
{
    const REAL rho = U[RHO];
    const REAL e = getSpecificInternalEnergy(U);
    outVars[0] = {rho};
    outVars[1].resize(SPACEDIM);
    for(int d = 0; d < SPACEDIM; ++d)
    {
        outVars[1][d] = U[MOM[d]] / rho;
    }
    outVars[2] = {m_eos->getPressure(rho, e)};
}
#endif
