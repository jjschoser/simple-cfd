#include "FluxSolver.H"

LaxFriedrichsSolver::LaxFriedrichsSolver(const Euler& euler) : m_euler(euler)
{

}

void LaxFriedrichsSolver::operator()(std::array<REAL, Euler::NVARS>& F, const std::array<REAL, Euler::NVARS>& ULo, 
                                     const std::array<REAL, Euler::NVARS>& UHi, const std::array<REAL, GRIDDIM>& dx, 
                                     const REAL dt, const int dim) const
{
    std::array<REAL, Euler::NVARS> FLo, FHi;
    m_euler.evalFluxFunc(FLo, ULo, dim);
    m_euler.evalFluxFunc(FHi, UHi, dim);
    for(int v = 0; v < Euler::NVARS; ++v)
    {
        F[v] = 0.5 * (FLo[v] + FHi[v]) - 0.5 * dx[dim] / dt * (UHi[v] - ULo[v]);
    }
}

HLLCSolver::HLLCSolver(const Euler& euler) : m_euler(euler)
{
    
}

void HLLCSolver::operator()(std::array<REAL, Euler::NVARS>& F, const std::array<REAL, Euler::NVARS>& ULo, 
                      const std::array<REAL, Euler::NVARS>& UHi, const std::array<REAL, GRIDDIM>&, 
                      const REAL, const int dim) const
{
    const REAL rhoLo = ULo[Euler::RHO];
    const REAL rhoHi = UHi[Euler::RHO];
    const REAL velLo = ULo[Euler::MOM[dim]] / rhoLo;
    const REAL velHi = UHi[Euler::MOM[dim]] / rhoHi;
    const REAL eLo = m_euler.getSpecificInternalEnergy(ULo);
    const REAL eHi = m_euler.getSpecificInternalEnergy(UHi);
    const REAL pLo = m_euler.getEquationOfState()->getPressure(rhoLo, eLo);
    const REAL pHi = m_euler.getEquationOfState()->getPressure(rhoHi, eHi);
    const REAL cLo = m_euler.getEquationOfState()->getSoundSpeed(rhoLo, pLo);
    const REAL cHi = m_euler.getEquationOfState()->getSoundSpeed(rhoHi, pHi);

    REAL sLo, sHi;
    getFastWaveSpeeds(sLo, sHi, velLo, velHi, cLo, cHi);

    if(sLo >= 0.0)
    {
        m_euler.evalFluxFunc(F, ULo, dim);
    }
    else if(sHi <= 0.0)
    {
        m_euler.evalFluxFunc(F, UHi, dim);
    }
    else
    {
        const REAL sStar = getContactSpeed(rhoLo, rhoHi, velLo, velHi, pLo, pHi, sLo, sHi);

        if(sStar >= 0.0)
        {
            std::array<REAL, Euler::NVARS> UStarLo;
            getStarState(UStarLo, ULo, rhoLo, velLo, pLo, sLo, sStar, dim);
            m_euler.evalFluxFunc(F, ULo, dim);
            for(int v = 0; v < Euler::NVARS; ++v)
            {
                F[v] += sLo * (UStarLo[v] - ULo[v]);
            }
        }
        else
        {
            std::array<REAL, Euler::NVARS> UStarHi;
            getStarState(UStarHi, UHi, rhoHi, velHi, pHi, sHi, sStar, dim);
            m_euler.evalFluxFunc(F, UHi, dim);
            for(int v = 0; v < Euler::NVARS; ++v)
            {
                F[v] += sHi * (UStarHi[v] - UHi[v]);
            }
        }
    }
}

void HLLCSolver::getFastWaveSpeeds(REAL& sLo, REAL& sHi, const REAL velLo, const REAL velHi, const REAL cLo, const REAL cHi) const
{
    const REAL sMax = std::max(std::fabs(velLo) + cLo, std::fabs(velHi) + cHi);
    sLo = -sMax;
    sHi = sMax;
}

REAL HLLCSolver::getContactSpeed(const REAL rhoLo, const REAL rhoHi, const REAL velLo, const REAL velHi, 
                           const REAL pLo, const REAL pHi, const REAL sLo, const REAL sHi) const
{
    return (pHi - pLo + rhoLo * velLo * (sLo - velLo) - rhoHi * velHi * (sHi - velHi)) / (rhoLo * (sLo - velLo) - rhoHi * (sHi - velHi));
}

void HLLCSolver::getStarState(std::array<REAL, Euler::NVARS>& UStar, const std::array<REAL, Euler::NVARS>& U, const REAL rho, 
                        const REAL vel, const REAL p, const REAL s, const REAL sStar, const int dim) const
{
    UStar[Euler::RHO] = 1.0;
    for(int d = 0; d < SPACEDIM; ++d)
    {
        UStar[Euler::MOM[d]] = U[Euler::MOM[d]] / rho;
    }
    UStar[Euler::MOM[dim]] = sStar;
    UStar[Euler::ENE] = U[Euler::ENE] / rho + (sStar - vel) * (sStar + p / (rho * (s - vel)));

    const REAL factor = rho * (s - vel) / (s - sStar);
    for(int v = 0; v < Euler::NVARS; ++v)
    {
        UStar[v] *= factor;
    }
}
