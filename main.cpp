#include <iostream>

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
    #ifdef DEBUG
        assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));
    #endif
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
        mesh.writeToFile("SodTest.txt", "SodTest.dat", finalStep, finalTime);
    #elif GRIDDIM == 2
        mesh.writeToFile("CylindricalExplosion.txt", "CylindricalExplosion.dat", finalStep, finalTime);
    #else  // GRIDDIM == 3
        mesh.writeToFile("SphericalExplosion.txt", "SphericalExplosion.dat", finalStep, finalTime);
    #endif
}

void runKelvinHelmholtzTest(const Euler& euler, 
                            const FluxSolver* const fluxSolver, 
                            const Reconstruction* const recon, 
                            const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM == 2);
    #ifdef DEBUG
        assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));
    #endif

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
    mesh.writeToFile("KelvinHelmholtz.txt", "KelvinHelmholtz.dat", finalStep, finalTime);
}

#ifdef USE_RIGID
void runShockReflectionTest(const Euler& euler, 
                            const FluxSolver* const fluxSolver, 
                            const Reconstruction* const recon, 
                            const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM == 2);
    #ifdef DEBUG
        assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));
    #endif

    const REAL xShock = 4e-3;
    const REAL xWedge = 4.96e-3;
    const REAL alphaWedge = 25.0 * M_PI / 180.0;

    const REAL MShock = 1.7;
    const REAL rhoInf = 1.225;
    const REAL velInf = 0.0;
    const REAL pInf = 101325.0;

    const IdealGas* const eos = dynamic_cast<const IdealGas* const>(euler.getEquationOfState());
    const REAL gamma = eos->getGamma();
    const REAL cInf = eos->getSoundSpeed(rhoInf, pInf);
    const REAL MInf = velInf / cInf;
    const REAL velShock = MShock * cInf;
    const REAL rhoStar = rhoInf * ((gamma + 1) * (MInf - MShock) * (MInf - MShock)) / ((gamma - 1) * (MInf - MShock) * (MInf - MShock) + 2);
    const REAL velStar = (1 - rhoInf / rhoStar) * (velInf + velShock);
    const REAL pStar = pInf * ((2 * gamma * (MInf - MShock) * (MInf - MShock) - (gamma - 1)) / (gamma + 1));

    const std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(-4e-3, 0.0, 0.0)};
    const std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(29e-3, 16.5e-3, 0.0)};
    const REAL finalTime = 35e-6;

    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for(int s = 0; s < 2; ++s)
    {
        for(int d = 0; d < GRIDDIM; ++d)
        {
            bc[s][d] = BoundaryCondition::TRANSMISSIVE;
        }
    }
    bc[0][1] = BoundaryCondition::REFLECTIVE;

    const Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 2);

    #ifdef USE_OMP
    #pragma omp parallel for default(none) shared(res, geom, mesh, euler, rhoStar, velStar, pStar) schedule(static)
    #endif
    for(int i = -mesh.SDFNGHOST; i < res[0] + mesh.SDFNGHOST; ++i)
    {
        #if GRIDDIM >= 2
        for(int j = -mesh.SDFNGHOST; j < res[1] + mesh.SDFNGHOST; ++j)
        #endif
        {
            #if GRIDDIM == 3
            for(int k = -mesh.SDFNGHOST; k < res[2] + mesh.SDFNGHOST; ++k)
            #endif
            {
                const std::array<int, GRIDDIM> idx = {GRIDDIM_DECL(i, j, k)};
                std::array<REAL, GRIDDIM> pos;
                geom.getPos(pos, idx);
                REAL rho, p;
                std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(0.0, 0.0, 0.0)};
                if(pos[0] < xShock)
                {
                    rho = rhoStar;
                    vel[0] = velStar;
                    p = pStar;
                }
                else
                {
                    rho = rhoInf;
                    vel[0] = velInf;
                    p = pInf;
                }
                mesh(idx)[euler.RHO] = rho;
                for(int d = 0; d < SPACEDIM; ++d)
                {
                    mesh(idx)[euler.MOM[d]] = rho * vel[d];
                }
                mesh(idx)[euler.ENE] = euler.getTotalEnergy(rho, vel, p);
                mesh.setSDF(std::cos(alphaWedge) * pos[1] - std::sin(alphaWedge) * (pos[0] - xWedge), idx);
            }
        }
    }
    
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon);
    mesh.writeToFile("ShockReflection.txt", "ShockReflection.dat", finalStep, finalTime);
    mesh.writeSDFToFile("ShockReflectionSDF.txt", "ShockReflectionSDF.dat");
}

void runHypersonicSphereTest(const Euler& euler, 
                             const FluxSolver* const fluxSolver, 
                             const Reconstruction* const recon, 
                             const std::array<int, GRIDDIM>& res,
                             const bool useSTL)
{
    assert(GRIDDIM == 3);
    #ifdef DEBUG
        assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));
    #endif

    std::string name = "HypersonicSphere";

    #if GRIDDIM == 3
        if(useSTL)
        {
            name += "FromSTL";
        }
    #endif

    const REAL rSphere = 10e-3;
    const REAL rhoInf = 0.0798;
    const REAL velInf = 1002.25 / 2.0;  // Reduce true velocity from test problem because solver crashes otherwise
    const REAL pInf = 2290.85;

    const std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(-25e-3, -25e-3, -25e-3)};
    const std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(0.0, 25e-3, 25e-3)};
    const REAL finalTime = 100e-6;

    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for(int s = 0; s < 2; ++s)
    {
        for(int d = 0; d < GRIDDIM; ++d)
        {
            bc[s][d] = BoundaryCondition::TRANSMISSIVE;
        }
    }

    const Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 2);

    #ifdef USE_OMP
    #pragma omp parallel for default(none) shared(res, geom, mesh, euler, useSTL, rhoInf, velInf, pInf) schedule(static)
    #endif
    for(int i = -mesh.SDFNGHOST; i < res[0] + mesh.SDFNGHOST; ++i)
    {
        #if GRIDDIM >= 2
        for(int j = -mesh.SDFNGHOST; j < res[1] + mesh.SDFNGHOST; ++j)
        #endif
        {
            #if GRIDDIM == 3
            for(int k = -mesh.SDFNGHOST; k < res[2] + mesh.SDFNGHOST; ++k)
            #endif
            {
                const std::array<int, GRIDDIM> idx = {GRIDDIM_DECL(i, j, k)};
                std::array<REAL, GRIDDIM> pos;
                geom.getPos(pos, idx);
                std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(velInf, 0.0, 0.0)};
                mesh(idx)[euler.RHO] = rhoInf;
                for(int d = 0; d < SPACEDIM; ++d)
                {
                    mesh(idx)[euler.MOM[d]] = rhoInf * vel[d];
                }
                mesh(idx)[euler.ENE] = euler.getTotalEnergy(rhoInf, vel, pInf);
                #if GRIDDIM == 3
                if(!useSTL)
                #endif
                {
                    mesh.setSDF(std::sqrt(GRIDDIM_TERM(pos[0]*pos[0], + pos[1]*pos[1], + pos[2]*pos[2])) - rSphere, idx);
                }
            }
        }
    }

    #if GRIDDIM == 3
    if(useSTL)
    {
        mesh.readSDFFromSTL("sphere.stl");
    }
    #endif
    
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon);
    mesh.writeToFile(name + ".txt", name + ".dat", finalStep, finalTime);
    mesh.writeSDFToFile(name + "SDF.txt", name + "SDF.dat");
}

void runWingTest(const Euler& euler, 
                 const FluxSolver* const fluxSolver, 
                 const Reconstruction* const recon, 
                 const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM == 3);
    #ifdef DEBUG
        assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));
    #endif

    std::string name = "Wing";

    const REAL rhoInf = 1.225;
    const REAL velInf = 315.81;
    const REAL pInf = 101325.0;

    const std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(-200e-3, -200e-3, 0.0)};
    const std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(600e-3, 200e-3, 400e-3)};
    const REAL finalTime = 5e-3;

    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for(int s = 0; s < 2; ++s)
    {
        for(int d = 0; d < GRIDDIM; ++d)
        {
            bc[s][d] = BoundaryCondition::TRANSMISSIVE;
        }
    }
    bc[0][2] = BoundaryCondition::REFLECTIVE;

    const Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 2);

    #ifdef USE_OMP
    #pragma omp parallel for default(none) shared(res, geom, mesh, euler, rhoInf, velInf, pInf) schedule(static)
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
                std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(velInf, 0.0, 0.0)};
                mesh(idx)[euler.RHO] = rhoInf;
                for(int d = 0; d < SPACEDIM; ++d)
                {
                    mesh(idx)[euler.MOM[d]] = rhoInf * vel[d];
                }
                mesh(idx)[euler.ENE] = euler.getTotalEnergy(rhoInf, vel, pInf);
            }
        }
    }

    #if GRIDDIM == 3
    mesh.readSDFFromSTL("wing.stl");
    #endif
    
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon);
    mesh.writeToFile(name + ".txt", name + ".dat", finalStep, finalTime);
    mesh.writeSDFToFile(name + "SDF.txt", name + "SDF.dat");
}
#endif

int main(int argc, char *argv[])
{
    std::string initFileName, finalFileName;
    REAL finalTime;
    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    REAL gamma = 1.4;
    #ifdef USE_RIGID
        std::string sdfFileName;
    #endif

    if(argc >= 2)
    {
        const std::string settingsFileName = argv[1];
        std::ifstream file(settingsFileName + ".txt");
        assert(file.is_open());
        std::string finalTimeLine, loBCLine, hiBCLine, gammaLine;
        std::getline(file, initFileName);
        std::getline(file, finalFileName);
        std::getline(file, finalTimeLine);
        std::getline(file, loBCLine);
        std::getline(file, hiBCLine);
        std::istringstream finalTimeISS(finalTimeLine);
        std::istringstream loBCISS(loBCLine);
        std::istringstream hiBCISS(hiBCLine);
        finalTimeISS >> finalTime;
        int loBC, hiBC;
        for(int d = 0; d < GRIDDIM; ++d)
        {
            loBCISS >> loBC;
            hiBCISS >> hiBC;
            bc[0][d] = static_cast<BoundaryCondition>(loBC);
            bc[1][d] = static_cast<BoundaryCondition>(hiBC);
        }
        if(std::getline(file, gammaLine))
        {
            std::istringstream gammaISS(gammaLine);
            gammaISS >> gamma;
        }
        #ifdef USE_RIGID
            std::getline(file, sdfFileName);
        #endif
        file.close();
    }

    const IdealGas eos(gamma);
    const Euler euler(&eos);
    const HLLCSolver fluxSolver(euler);
    const MUSCLHancock recon(euler);

    if(!initFileName.empty())
    {
        int startStep;
        REAL startTime;
        Mesh<Euler::NVARS> mesh = Mesh<Euler::NVARS>::createFromFile(initFileName + ".txt", startStep, startTime, 2);
        #ifdef USE_RIGID
            if(!sdfFileName.empty())
            {
                #if GRIDDIM == 3
                    if(!mesh.readSDFFromSTL(sdfFileName + ".stl"))
                #endif
                {
                    mesh.readSDFFromFile(sdfFileName + ".txt");
                }
            }
        #endif
        const int finalStep = solve(euler, finalTime, mesh, bc, &fluxSolver, &recon, 0.9, startStep, startTime);
        mesh.writeToFile(finalFileName + ".txt", finalFileName + ".dat", finalStep, finalTime);
    }
    else
    {
        std::cout << "Running test problems..." << std::endl;
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
            #ifdef USE_RIGID
                const std::array<int, GRIDDIM> shockReflectionRes = {res[0], res[1] / 2};
                runShockReflectionTest(euler, &fluxSolver, &recon, shockReflectionRes);
            #endif
        #endif

        #if GRIDDIM == 3
            #ifdef USE_RIGID
                const std::array<int, GRIDDIM> hypersonicSphereRes = {GRIDDIM_DECL(res[0] / 2, res[1], res[2])};
                runHypersonicSphereTest(euler, &fluxSolver, &recon, hypersonicSphereRes, false);
                runHypersonicSphereTest(euler, &fluxSolver, &recon, hypersonicSphereRes, true);

                const std::array<int, GRIDDIM> wingRes = {GRIDDIM_DECL(2 * res[0], res[1], res[2])};
                runWingTest(euler, &fluxSolver, &recon, wingRes);
            #endif
        #endif
    }

    return 0;
}