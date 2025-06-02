#ifndef EULER_H
#define EULER_H

#include <array>

#include "EquationOfState.H"
#include "Macros.H"

class Euler
{
    public:
        static constexpr int NVARS = 1 + SPACEDIM + 1;  // Density, momentum, total energy
        static constexpr int RHO = 0;  // Index for density
        static constexpr std::array<int, SPACEDIM> MOM = {SPACEDIM_DECL(RHO + 1, RHO + 2, RHO + 3)};  // Indices for momentum components
        static constexpr int ENE = MOM[SPACEDIM - 1] + 1;  // Index for total energy

        Euler(const EquationOfState* const eos);
        const EquationOfState* getEquationOfState() const;
        REAL getSpecificInternalEnergy(const std::array<REAL, NVARS>& U) const;
        REAL getTotalEnergy(const REAL rho, const std::array<REAL, SPACEDIM>& vel, const REAL p) const;
        void evalFluxFunc(std::array<REAL, NVARS>& F, const std::array<REAL, NVARS>& U, const int dim) const;
        REAL getMaxWaveSpeed(const std::array<REAL, NVARS>& U) const;

    private:
        const EquationOfState* const m_eos;
};

#endif