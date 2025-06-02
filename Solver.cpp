#include <iostream>

#include "Solver.H"

int solve(const Euler& euler, const REAL finalTime, Mesh<Euler::NVARS>& mesh, 
          const std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, 
          const FluxSolver* const fluxSolver, const REAL cfl, 
          const int startStep, const REAL startTime)
{
    int step = startStep;
    REAL t = startTime;
    std::array<Mesh<Euler::NVARS>, GRIDDIM> fluxMeshes = {GRIDDIM_DECL(mesh.createFaceMesh(0), 
                                                                       mesh.createFaceMesh(1), 
                                                                       mesh.createFaceMesh(2))};
    assert(SPACEDIM >= GRIDDIM);
    const std::vector<std::array<int, GRIDDIM>> vecIdx = {{GRIDDIM_DECL(Euler::MOM[0], 
                                                                        Euler::MOM[1], 
                                                                        Euler::MOM[2])}};
    REAL dt;
    while(t < finalTime)
    {
        dt = std::min(calcDt(euler, mesh, cfl), finalTime - t);
        std::cout << "step = " << step << ", time = " << t << " (" << (t - startTime) / (finalTime - startTime) * 100 << "%), dt = " << dt << std::endl;
        for(int d = 0; d < GRIDDIM; ++d)
        {
            mesh.fillGhost(bc, vecIdx);
            calcFlux(fluxSolver, mesh, fluxMeshes[d], dt, d);
            updateMesh(mesh, fluxMeshes[d], dt, d);
        }
        t += dt;
        ++step;
    }
    return step;
}

REAL calcDt(const Euler& euler, const Mesh<Euler::NVARS>& mesh, const REAL cfl)
{
    const std::array<int, GRIDDIM>& res = mesh.getGeometry().getRes();
    const std::array<REAL, GRIDDIM>& dx = mesh.getGeometry().getDx();

    REAL maxWaveSpeed = 0.0;
    #if GRIDDIM == 3
        for(int k = 0; k < res[2]; ++k)
    #endif
        {
    #if GRIDDIM >= 2
            for(int j = 0; j < res[1]; ++j)
    #endif
            {
                for(int i = 0; i < res[0]; ++i)
                {
                    const std::array<REAL, Euler::NVARS>& U = mesh(GRIDDIM_DECL(i, j, k));
                    maxWaveSpeed = std::max(maxWaveSpeed, euler.getMaxWaveSpeed(U));
                }
            }
        }
    
    REAL minDx = dx[0];
    #if GRIDDIM >= 2
        for(int d = 1; d < GRIDDIM; ++d)
        {
            minDx = std::min(minDx, dx[d]);
        }
    #endif

    return cfl * minDx / maxWaveSpeed;
}

void calcFlux(const FluxSolver* const fluxSolver, const Mesh<Euler::NVARS>& mesh, 
              Mesh<Euler::NVARS>& fluxMesh,  REAL dt, const int dim)
{
    const std::array<int, GRIDDIM>& res = mesh.getGeometry().getRes();
    const std::array<REAL, GRIDDIM>& dx = mesh.getGeometry().getDx();

    std::array<int, GRIDDIM> offset = {GRIDDIM_DECL(0, 0, 0)};
    offset[dim] = 1;

    #if GRIDDIM == 3
        for(int k = 0; k < res[2] + offset[2]; ++k)
    #endif
        {
    #if GRIDDIM >= 2
            for(int j = 0; j < res[1] + offset[1]; ++j)
    #endif
            {
                for(int i = 0; i < res[0] + offset[0]; ++i)
                {
                    std::array<REAL, Euler::NVARS>& F = fluxMesh(GRIDDIM_DECL(i, j, k));
                    const std::array<REAL, Euler::NVARS>& ULo = mesh(GRIDDIM_DECL(i - offset[0], j - offset[1], k - offset[2]));
                    const std::array<REAL, Euler::NVARS>& UHi = mesh(GRIDDIM_DECL(i, j, k));
                    (*fluxSolver)(F, ULo, UHi, dx, dt, dim);
                }
            }
        }
}

void updateMesh(Mesh<Euler::NVARS>& mesh, const Mesh<Euler::NVARS>& fluxMesh, 
                const REAL dt, const int dim)
{
    const std::array<int, GRIDDIM>& res = mesh.getGeometry().getRes();
    const std::array<REAL, GRIDDIM>& dx = mesh.getGeometry().getDx();

    std::array<int, GRIDDIM> offset = {GRIDDIM_DECL(0, 0, 0)};
    offset[dim] = 1;

    #if GRIDDIM == 3
        for(int k = 0; k < res[2]; ++k)
    #endif
        {
    #if GRIDDIM >= 2
            for(int j = 0; j < res[1]; ++j)
    #endif
            {
                for(int i = 0; i < res[0]; ++i)
                {
                    std::array<REAL, Euler::NVARS>& U = mesh(GRIDDIM_DECL(i, j, k));
                    const std::array<REAL, Euler::NVARS>& FLo = fluxMesh(GRIDDIM_DECL(i, j, k));
                    const std::array<REAL, Euler::NVARS>& FHi = fluxMesh(GRIDDIM_DECL(i + offset[0], j + offset[1], k + offset[2]));
                    for(int v = 0; v < Euler::NVARS; ++v)
                    {
                        U[v] -= dt / dx[dim] * (FHi[v] - FLo[v]);
                    }
                }
            }
        }
}
