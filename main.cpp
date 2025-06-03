#include "EquationOfState.H"
#include "Euler.H"
#include "FluxSolver.H"
#include "Mesh.H"
#include "Reconstruction.H"
#include "Solver.H"

#ifdef USE_OMP
    #include <omp.h>
#endif

void runSimpleTest(const Euler& euler, 
                   const FluxSolver* const fluxSolver, 
                   const Reconstruction* const recon, 
                   const std::array<int, GRIDDIM>& res)
{
    const std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
    const std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(1.0, 1.0, 1.0)};
    const REAL finalTime = 0.25;

    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for(int s = 0; s < 2; ++s)
    {
        for(int d = 0; d < GRIDDIM; ++d)
        {
            if(s == 0)
            {
                bc[s][d] = BoundaryCondition::REFLECTIVE;
            }
            else
            {
                bc[s][d] = BoundaryCondition::TRANSMISSIVE;
            }
        }
    }

    const Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 2);

    #if GRIDDIM == 1
        const REAL rInter = 0.5;
    #else
        const REAL rInter = 0.4;
    #endif

    #ifdef USE_OMP
    #pragma omp parallel for default(none) shared(res, geom, mesh, euler) schedule(static)
    #endif
    for(int i = 0; i < res[0]; ++i)
    {
        #if GRIDDIM >= 2
        for(int j = 0; j < res[1]; ++j)
        #endif
        {
            #if GRIDDIM == 3
            for(int k = 0; k < res[2]; ++k)
            #endif
            {
                const std::array<int, GRIDDIM> idx = {GRIDDIM_DECL(i, j, k)};
                std::array<REAL, GRIDDIM> pos;
                geom.getPos(pos, idx);
                const REAL r = std::sqrt(GRIDDIM_TERM(pos[0] * pos[0], + pos[1] * pos[1], + pos[2] * pos[2]));
                REAL rho, p;
                std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(0.0, 0.0, 0.0)};
                if(r < rInter)
                {
                    rho = 1.0;
                    p = 1.0;
                }
                else
                {
                    rho = 0.125;
                    p = 0.1;
                }
                mesh(idx)[euler.RHO] = rho;
                for(int d = 0; d < SPACEDIM; ++d)
                {
                    mesh(idx)[euler.MOM[d]] = rho * vel[d];
                }
                mesh(idx)[euler.ENE] = euler.getTotalEnergy(rho, vel, p);
            }
        }
    }
    
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon);

    #if GRIDDIM == 1
        mesh.writeToFile("SodTest.txt", finalStep, finalTime);
    #elif GRIDDIM == 2
        mesh.writeToFile("CylindricalExplosion.txt", finalStep, finalTime);
    #else  // GRIDDIM == 3
        mesh.writeToFile("SphericalExplosion.txt", finalStep, finalTime);
    #endif
}

void runKelvinHelmholtzTest(const Euler& euler, 
                            const FluxSolver* const fluxSolver, 
                            const Reconstruction* const recon, 
                            const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM == 2);

    const std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(-0.5, -0.5, -0.5)};
    const std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(0.5, 0.5, 0.5)};
    const REAL finalTime = 1.5;

    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for(int s = 0; s < 2; ++s)
    {
        for(int d = 0; d < GRIDDIM; ++d)
        {
            bc[s][d] = BoundaryCondition::PERIODIC;
        }
    }

    const Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 2);

    #ifdef USE_OMP
    #pragma omp parallel for default(none) shared(res, geom, mesh, euler) schedule(static)
    #endif
    for(int i = 0; i < res[0]; ++i)
    {
        #if GRIDDIM >= 2
        for(int j = 0; j < res[1]; ++j)
        #endif
        {
            #if GRIDDIM == 3
            for(int k = 0; k < res[2]; ++k)
            #endif
            {
                const std::array<int, GRIDDIM> idx = {GRIDDIM_DECL(i, j, k)};
                std::array<REAL, GRIDDIM> pos;
                geom.getPos(pos, idx);
                REAL rho;
                std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(0.0, 0.0, 0.0)};
                const REAL p = 2.5;
                if(std::fabs(pos[1]) < 0.25)
                {
                    rho = 2.0;
                    vel[0] = -0.5;
                }
                else
                {
                    rho = 1.0;
                    vel[0] = 0.5;
                }
                vel[1] = 0.01 * std::sin(2.0 * M_PI * pos[0]); 
                mesh(idx)[euler.RHO] = rho;
                for(int d = 0; d < SPACEDIM; ++d)
                {
                    mesh(idx)[euler.MOM[d]] = rho * vel[d];
                }
                mesh(idx)[euler.ENE] = euler.getTotalEnergy(rho, vel, p);
            }
        }
    }
    
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon);
    mesh.writeToFile("KelvinHelmholtz.txt", finalStep, finalTime);
}

int main()
{
    const IdealGas eos(1.4);
    const Euler euler(&eos);
    const HLLCSolver fluxSolver(euler);
    const MUSCLHancock recon(euler);

    #if GRIDDIM == 1
        const std::array<int, GRIDDIM> res = {2048};
    #elif GRIDDIM == 2
        const std::array<int, GRIDDIM> res = {512, 512};
    #else  // GRIDDIM == 3
        const std::array<int, GRIDDIM> res = {128, 128, 128};
    #endif

    runSimpleTest(euler, &fluxSolver, &recon, res);
    #if GRIDDIM == 2
        runKelvinHelmholtzTest(euler, &fluxSolver, &recon, res);
    #endif

    return 0;
}