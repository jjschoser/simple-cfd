#ifndef EQUATION_OF_STATE_H
#define EQUATION_OF_STATE_H

#include <cmath>

#include "Macros.H"

class EquationOfState
{
    public:
        virtual REAL getPressure(const REAL rho, const REAL e) const = 0;
        virtual REAL getSpecificInternalEnergy(const REAL rho, const REAL p) const = 0;
        virtual REAL getSoundSpeed(const REAL rho, const REAL p) const = 0;
};

class IdealGas : public EquationOfState
{
    public:
        IdealGas(const REAL gamma=1.4);
        REAL getPressure(const REAL rho, const REAL e) const override;
        REAL getSpecificInternalEnergy(const REAL rho, const REAL p) const override;
        REAL getSoundSpeed(const REAL rho, const REAL p) const override;

    private:
        const REAL m_gamma;
};

#endif