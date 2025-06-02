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
    for (int v = 0; v < Euler::NVARS; ++v)
    {
        F[v] = 0.5 * (FLo[v] + FHi[v]) - 0.5 * dx[dim] / dt * (UHi[v] - ULo[v]);
    }
}