#ifndef SOLVER_H
#define SOLVER_H

#include "Euler.H"
#include "FluxSolver.H"
#include "Macros.H"
#include "Mesh.H"

int solve(const Euler& euler, const REAL finalTime, Mesh<Euler::NVARS>& mesh, 
          const std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, 
          const FluxSolver* const fluxSolver, const REAL cfl=0.9, 
          const int startStep=0, const REAL startTime=0.0);

REAL calcDt(const Euler& euler, const Mesh<Euler::NVARS>& mesh, const REAL cfl);

void calcFlux(const FluxSolver* const fluxSolver, const Mesh<Euler::NVARS>& mesh, 
              Mesh<Euler::NVARS>& fluxMesh,  REAL dt, const int dim);

void updateMesh(Mesh<Euler::NVARS>& mesh, const Mesh<Euler::NVARS>& fluxMesh, 
                const REAL dt, const int dim);

#endif