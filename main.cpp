#include "EquationOfState.H"
#include "Euler.H"
#include "FluxSolver.H"
#include "Mesh.H"
#include "Solver.H"

int main()
{
    const std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
    const std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(1.0, 1.0, 1.0)};
    const REAL finalTime = 0.25;

    const IdealGas eos(1.4);
    const Euler euler(&eos);
    const HLLCSolver fluxSolver(euler);

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

    #if GRIDDIM == 1
        const std::array<int, GRIDDIM> res = {2048};
    #elif GRIDDIM == 2
        const std::array<int, GRIDDIM> res = {512, 512};
    #else  // GRIDDIM == 3
        const std::array<int, GRIDDIM> res = {128, 128, 128};
    #endif

    const Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 1);

    REAL rho, p;
    std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(0.0, 0.0, 0.0)};
    REAL r;
    std::array<int, GRIDDIM> idx;
    std::array<REAL, GRIDDIM> pos;

    #if GRIDDIM == 1
        const REAL rInter = 0.5;
    #else
        const REAL rInter = 0.4;
    #endif

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
                    idx = {GRIDDIM_DECL(i, j, k)};
                    geom.getPos(pos, idx);
                    r = std::sqrt(GRIDDIM_TERM(pos[0] * pos[0], + pos[1] * pos[1], + pos[2] * pos[2]));
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
    
    const int finalStep = solve(euler, finalTime, mesh, bc, &fluxSolver);
    mesh.writeToFile("out.txt", finalStep, finalTime);

    return 0;
}