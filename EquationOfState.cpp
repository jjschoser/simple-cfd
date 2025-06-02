#include "EquationOfState.H"

IdealGas::IdealGas(const REAL gamma) : m_gamma(gamma)
{

}

REAL IdealGas::getPressure(const REAL rho, const REAL e) const
{
    return (m_gamma - 1.0) * rho * e;
}

REAL IdealGas::getSpecificInternalEnergy(const REAL rho, const REAL p) const
{
    return p / ((m_gamma - 1.0) * rho);
}

REAL IdealGas::getSoundSpeed(const REAL rho, const REAL p) const
{
    return std::sqrt(m_gamma * p / rho);
}