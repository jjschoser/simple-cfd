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

#endif