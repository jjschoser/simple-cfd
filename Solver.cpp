#include <iostream>

#include "Solver.H"

int solve(const Euler& euler, const REAL finalTime, Mesh<Euler::NVARS>& mesh, 
          const std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, 
          const FluxSolver* const fluxSolver, const Reconstruction* const recon, 
          const REAL cfl, const int startStep, const REAL startTime)
{
    assert(startStep >= 0);
    assert(startTime >= 0.0);
    assert(finalTime > startTime);
    assert(cfl > 0.0 && cfl < 1.0);
    assert(mesh.getNGhost() >= 1 + recon->getStencilSize());

    int step = startStep;
    REAL t = startTime;

    std::array<int, GRIDDIM> reconRes = mesh.getRes();
    std::array<int, GRIDDIM> fluxRes = mesh.getRes();
    for(int d = 0; d < GRIDDIM; ++d)
    {
        fluxRes[d] += 1;
    }
    DataArray<Euler::NVARS> reconDataLo(reconRes, mesh.getNGhost() - recon->getStencilSize());
    DataArray<Euler::NVARS> reconDataHi(reconRes, mesh.getNGhost() - recon->getStencilSize());
    DataArray<Euler::NVARS> fluxData(fluxRes);

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
            doReconstruction(recon, mesh, reconDataLo, reconDataHi, dt, d);
            calcFlux(fluxSolver, mesh.getGeometry(), reconDataLo, reconDataHi, fluxData, dt, d);
            updateMesh(mesh, fluxData, dt, d);
        }
        t += dt;
        ++step;
    }
    return step;
}

REAL calcDt(const Euler& euler, const Mesh<Euler::NVARS>& mesh, const REAL cfl)
{
    const std::array<int, GRIDDIM>& res = mesh.getRes();
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

void doReconstruction(const Reconstruction* const recon, 
                      const Mesh<Euler::NVARS>& mesh, 
                      DataArray<Euler::NVARS>& reconDataLo, 
                      DataArray<Euler::NVARS>& reconDataHi, 
                      const REAL dt, const int dim)
{
    const int nGhost = mesh.getNGhost() - recon->getStencilSize();
    const std::array<int, GRIDDIM>& res = mesh.getRes();
    const std::array<REAL, GRIDDIM>& dx = mesh.getGeometry().getDx();

    std::array<int, GRIDDIM> offset = {GRIDDIM_DECL(0, 0, 0)};
    offset[dim] = 1;

    #if GRIDDIM == 3
        for(int k = -nGhost * offset[2]; k < res[2] + nGhost * offset[2]; ++k)
    #endif
        {
    #if GRIDDIM >= 2
            for(int j = -nGhost * offset[1]; j < res[1] + nGhost * offset[1]; ++j)
    #endif
            {
                for(int i = -nGhost * offset[0]; i < res[0] + nGhost * offset[0]; ++i)
                {
                    std::array<REAL, Euler::NVARS>& UReconLo = reconDataLo(GRIDDIM_DECL(i, j, k));
                    std::array<REAL, Euler::NVARS>& UReconHi = reconDataHi(GRIDDIM_DECL(i, j, k));
                    const std::array<REAL, Euler::NVARS>& UNbrLo = mesh(GRIDDIM_DECL(i - offset[0], j - offset[1], k - offset[2]));
                    const std::array<REAL, Euler::NVARS>& UNbrHi = mesh(GRIDDIM_DECL(i + offset[0], j + offset[1], k + offset[2]));
                    const std::array<REAL, Euler::NVARS>& U = mesh(GRIDDIM_DECL(i, j, k));
                    (*recon)(UReconLo, UReconHi, UNbrLo, UNbrHi, U, dx, dt, dim);
                }
            }
        }
}

void calcFlux(const FluxSolver* const fluxSolver, const Geometry& geom,
              const DataArray<Euler::NVARS>& reconDataLo,
              const DataArray<Euler::NVARS>& reconDataHi,
              DataArray<Euler::NVARS>& fluxData,  REAL dt, const int dim)
{
    const std::array<int, GRIDDIM>& res = geom.getRes();
    const std::array<REAL, GRIDDIM>& dx = geom.getDx();

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
                    std::array<REAL, Euler::NVARS>& F = fluxData(GRIDDIM_DECL(i, j, k));
                    const std::array<REAL, Euler::NVARS>& ULo = reconDataHi(GRIDDIM_DECL(i - offset[0], j - offset[1], k - offset[2]));
                    const std::array<REAL, Euler::NVARS>& UHi = reconDataLo(GRIDDIM_DECL(i, j, k));
                    (*fluxSolver)(F, ULo, UHi, dx, dt, dim);
                }
            }
        }
}

void updateMesh(Mesh<Euler::NVARS>& mesh, const DataArray<Euler::NVARS>& fluxData, 
                const REAL dt, const int dim)
{
    const std::array<int, GRIDDIM>& res = mesh.getRes();
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
                    const std::array<REAL, Euler::NVARS>& FLo = fluxData(GRIDDIM_DECL(i, j, k));
                    const std::array<REAL, Euler::NVARS>& FHi = fluxData(GRIDDIM_DECL(i + offset[0], j + offset[1], k + offset[2]));
                    for(int v = 0; v < Euler::NVARS; ++v)
                    {
                        U[v] -= dt / dx[dim] * (FHi[v] - FLo[v]);
                    }
                }
            }
        }
}
