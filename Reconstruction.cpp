#include "Reconstruction.H"

int Reconstruction::getStencilSize() const
{
    return 1;
}

void NoReconstruction::operator()(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi,
                                  const std::array<REAL, Euler::NVARS>&, const std::array<REAL, Euler::NVARS>&,
                                  const std::array<REAL, Euler::NVARS>& U, const std::array<REAL, GRIDDIM>,
                                  const REAL, const int) const
{
    UReconLo = U;
    UReconHi = U;
}

MUSCLHancock::MUSCLHancock(const Euler& euler) : m_euler(euler)
{
    
}

void MUSCLHancock::operator()(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi,
                              const std::array<REAL, Euler::NVARS>& UNbrLo, const std::array<REAL, Euler::NVARS>& UNbrHi,
                              const std::array<REAL, Euler::NVARS>& U, const std::array<REAL, GRIDDIM> dx,
                              const REAL dt, const int dim) const
{
    reconLinear(UReconLo, UReconHi, UNbrLo, UNbrHi, U);
    doHalfStep(UReconLo, UReconHi, dx, dt, dim);
}

REAL MUSCLHancock::slopeLimiter(const REAL r) const
{
    return std::max(0.0, std::min(2.0 * r / (1.0 + r), 2.0 / (1.0 + r)));
}

void MUSCLHancock::reconLinear(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi,
                               const std::array<REAL, Euler::NVARS>& UNbrLo, const std::array<REAL, Euler::NVARS>& UNbrHi,
                               const std::array<REAL, Euler::NVARS>& U) const
{
    REAL denom, r, phi, delta;
    for (int v = 0; v < Euler::NVARS; ++v)
    {
        denom = UNbrHi[v] - U[v];
        r = (std::fabs(denom) > 1e-16) ? (U[v] - UNbrLo[v]) / denom : 0.0;
        phi = slopeLimiter(r);
        delta = 0.5 * phi * (UNbrHi[v] - UNbrLo[v]);
        UReconLo[v] = U[v] - 0.5 * delta;
        UReconHi[v] = U[v] + 0.5 * delta;
    }
}

void MUSCLHancock::doHalfStep(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi,
                              const std::array<REAL, GRIDDIM> dx, const REAL dt, const int dim) const
{
    std::array<REAL, Euler::NVARS> FRecLo, FRecHi;
    m_euler.evalFluxFunc(FRecLo, UReconLo, dim);
    m_euler.evalFluxFunc(FRecHi, UReconHi, dim);

    REAL delta;
    for (int v = 0; v < Euler::NVARS; ++v)
    {
        delta = 0.5 * dt / dx[dim] * (FRecHi[v] - FRecLo[v]);
        UReconLo[v] -= delta;
        UReconHi[v] -= delta;
    }
}