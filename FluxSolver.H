#ifndef FLUX_SOLVER_H
#define FLUX_SOLVER_H

#include "Euler.H"
#include "Macros.H"

class FluxSolver
{
    public:
        virtual void operator()(std::array<REAL, Euler::NVARS>& F, const std::array<REAL, Euler::NVARS>& ULo, 
                                const std::array<REAL, Euler::NVARS>& UHi, const std::array<REAL, GRIDDIM>& dx, 
                                const REAL dt, const int dim) const = 0;
};

class LaxFriedrichsSolver : public FluxSolver
{
    public:
        LaxFriedrichsSolver(const Euler& euler);
        void operator()(std::array<REAL, Euler::NVARS>& F, const std::array<REAL, Euler::NVARS>& ULo, 
                        const std::array<REAL, Euler::NVARS>& UHi, const std::array<REAL, GRIDDIM>& dx, 
                        const REAL dt, const int dim) const override;

    private:
        const Euler m_euler;
};

class HLLCSolver : public FluxSolver
{
    public:
        HLLCSolver(const Euler& euler);
        void operator()(std::array<REAL, Euler::NVARS>& F, const std::array<REAL, Euler::NVARS>& ULo, 
                        const std::array<REAL, Euler::NVARS>& UHi, const std::array<REAL, GRIDDIM>& dx, 
                        const REAL dt, const int dim) const override;

    private:
        const Euler m_euler;
        void getFastWaveSpeeds(REAL& sLo, REAL& sHi, const REAL velLo, const REAL velHi, const REAL cLo, const REAL cHi) const;
        REAL getContactSpeed(const REAL rhoLo, const REAL rhoHi, const REAL velLo, const REAL velHi, 
                             const REAL pLo, const REAL pHi, const REAL sLo, const REAL sHi) const;
        void getStarState(std::array<REAL, Euler::NVARS>& UStar, const std::array<REAL, Euler::NVARS>& U, const REAL rho, 
                          const REAL vel, const REAL p, const REAL s, const REAL sStar, const int dim) const;
};

#endif