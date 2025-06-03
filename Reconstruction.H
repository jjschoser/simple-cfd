#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include "Euler.H"
#include "Macros.H"

class Reconstruction
{
    public:
        virtual void operator()(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi, 
                                const std::array<REAL, Euler::NVARS>& UNbrLo, const std::array<REAL, Euler::NVARS>& UNbrHi, 
                                const std::array<REAL, Euler::NVARS>& U, const std::array<REAL, GRIDDIM> dx, 
                                const REAL dt, const int dim) const = 0;
        int getStencilSize() const;
};

class NoReconstruction : public Reconstruction
{
    public:
        void operator()(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi, 
                        const std::array<REAL, Euler::NVARS>& UNbrLo, const std::array<REAL, Euler::NVARS>& UNbrHi, 
                        const std::array<REAL, Euler::NVARS>& U, const std::array<REAL, GRIDDIM> dx, 
                        const REAL dt, const int dim) const override;
};

class MUSCLHancock : public Reconstruction
{
    public:
        MUSCLHancock(const Euler& euler);
        void operator()(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi, 
                        const std::array<REAL, Euler::NVARS>& UNbrLo, const std::array<REAL, Euler::NVARS>& UNbrHi, 
                        const std::array<REAL, Euler::NVARS>& U, const std::array<REAL, GRIDDIM> dx, 
                        const REAL dt, const int dim) const override;
    
    private:
        const Euler m_euler;
        REAL slopeLimiter(const REAL r) const;
        void reconLinear(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi, 
                         const std::array<REAL, Euler::NVARS>& UNbrLo, const std::array<REAL, Euler::NVARS>& UNbrHi, 
                         const std::array<REAL, Euler::NVARS>& U) const;
        void doHalfStep(std::array<REAL, Euler::NVARS>& UReconLo, std::array<REAL, Euler::NVARS>& UReconHi, 
                        const std::array<REAL, GRIDDIM> dx, const REAL dt, const int dim) const;
};

#endif