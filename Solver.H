#ifndef SOLVER_H
#define SOLVER_H

#include "Euler.H"
#include "FluxSolver.H"
#include "Macros.H"
#include "Mesh.H"
#include "Reconstruction.H"

int solve(const Euler& euler, const REAL finalTime, Mesh<Euler::NVARS>& mesh, 
          const std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, 
          const FluxSolver* const fluxSolver, const Reconstruction* recon, 
          const REAL cfl=0.9, const int startStep=0, const REAL startTime=0.0);

REAL calcDt(const Euler& euler, const Mesh<Euler::NVARS>& mesh, const REAL cfl);

void doReconstruction(const Reconstruction* const recon, 
                      const Mesh<Euler::NVARS>& mesh, 
                      DataArray<Euler::NVARS>& reconDataLo, 
                      DataArray<Euler::NVARS>& reconDataHi, 
                      const REAL dt, const int dim);

void calcFlux(const FluxSolver* const fluxSolver, const Geometry& geom,
              const DataArray<Euler::NVARS>& reconDataLo, 
              const DataArray<Euler::NVARS>& reconDataHi, 
              DataArray<Euler::NVARS>& fluxData,  REAL dt, const int dim);

void updateMesh(Mesh<Euler::NVARS>& mesh, 
                const DataArray<Euler::NVARS>& fluxData, 
                const REAL dt, const int dim);

#endif