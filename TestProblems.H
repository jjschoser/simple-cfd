#ifndef TEST_PROBLEMS_H
#define TEST_PROBLEMS_H

#include <string>

#include "Euler.H"
#include "FluxSolver.H"
#include "Reconstruction.H"

void runSimpleTest(const Euler& euler, 
                   const FluxSolver* const fluxSolver, 
                   const Reconstruction* const recon, 
                   const std::array<int, GRIDDIM>& res);

#if GRIDDIM == 2
void runKelvinHelmholtzTest(const Euler& euler, 
                            const FluxSolver* const fluxSolver, 
                            const Reconstruction* const recon, 
                            const std::array<int, GRIDDIM>& res);
#endif

#ifdef USE_RIGID
std::string getSDFFilename(const std::string& filename);

#if GRIDDIM == 2
void runShockReflectionTest(const Euler& euler, 
                            const FluxSolver* const fluxSolver, 
                            const Reconstruction* const recon, 
                            const std::array<int, GRIDDIM>& res);
#endif

#if GRIDDIM == 3
void runHypersonicSphereTest(const Euler& euler, 
                             const FluxSolver* const fluxSolver, 
                             const Reconstruction* const recon, 
                             const std::array<int, GRIDDIM>& res,
                             const bool useSTL);

void runSpaceShuttleTest(const Euler& euler, 
                         const FluxSolver* const fluxSolver, 
                         const Reconstruction* const recon, 
                         const std::array<int, GRIDDIM>& res);
#endif

#endif

#endif