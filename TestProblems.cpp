#include "TestProblems.H"

#include <chrono>

#ifdef USE_OMP
    #include <omp.h>
#endif

#ifdef USE_RIGID
    #include <filesystem>
#endif

#include "Mesh.H"
#include "Solver.H"

const std::string testOutDir = "test-output/";

void runSimpleTest(const Euler& euler, 
                   const FluxSolver* const fluxSolver, 
                   const Reconstruction* const recon, 
                   const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));

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
                const std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(0.0, 0.0, 0.0)};
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

    #if GRIDDIM == 1
        const std::string name = "SodTest";
    #elif GRIDDIM == 2
        const std::string name = "CylindricalExplosion";
    #else  // GRIDDIM == 3
        const std::string name = "SphericalExplosion";
    #endif

    const std::string filenameWPath = testOutDir + name + ".txt";

    const auto start = std::chrono::high_resolution_clock::now();
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon, filenameWPath);
    const auto stop = std::chrono::high_resolution_clock::now();

    mesh.save(addStepCounter(filenameWPath, finalStep), finalStep, finalTime);

    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Ran " << name << " in " << duration << " ms" << std::endl;
}

#if GRIDDIM == 2
void runKelvinHelmholtzTest(const Euler& euler, 
                            const FluxSolver* const fluxSolver, 
                            const Reconstruction* const recon, 
                            const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));

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

    const std::string name = "KelvinHelmholtz";
    const std::string filenameWPath = testOutDir + name + ".txt";
    
    const auto start = std::chrono::high_resolution_clock::now();
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon, filenameWPath);
    const auto stop = std::chrono::high_resolution_clock::now();

    mesh.save(addStepCounter(filenameWPath, finalStep), finalStep, finalTime);

    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Ran " << name << " in " << duration << " ms" << std::endl;
}
#endif

#ifdef USE_RIGID
std::string getSDFFilename(const std::string& filename)
{
    std::filesystem::path path(filename);
    std::string dir = path.parent_path().string();
    std::string stem = path.stem().string();
    std::string extension = path.extension().string();
    return dir + "/" + stem + "SDF" + extension;
}

#if GRIDDIM == 2
void runShockReflectionTest(const Euler& euler, 
                            const FluxSolver* const fluxSolver, 
                            const Reconstruction* const recon, 
                            const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));

    const std::string name = "ShockReflection";
    const std::string filename = testOutDir + name + ".txt";

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
    
    mesh.saveSDF(getSDFFilename(filename));

    const auto start = std::chrono::high_resolution_clock::now();
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon, filename);
    const auto stop = std::chrono::high_resolution_clock::now();

    mesh.save(addStepCounter(filename, finalStep), finalStep, finalTime);

    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Ran " << name << " in " << duration << " ms" << std::endl;
}
#endif

#if GRIDDIM == 3
void runHypersonicSphereTest(const Euler& euler, 
                             const FluxSolver* const fluxSolver, 
                             const Reconstruction* const recon, 
                             const std::array<int, GRIDDIM>& res,
                             const bool useSTL)
{
    assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));

    std::string name = "HypersonicSphere";
    if(useSTL)
    {
        name += "FromSTL";
    }
    const std::string filename = testOutDir + name + ".txt";

    const REAL rSphere = 10e-3;
    const REAL rhoInf = 0.0798;
    const REAL velInf = 1002.25;
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
                const std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(velInf, 0.0, 0.0)};
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

    if(useSTL)
    {
        if(!mesh.loadSDF(getSDFFilename(filename)))
        {
            assert(mesh.loadSDF("sphere.stl"));
            mesh.saveSDF(getSDFFilename(filename));
        }
    }
    else
    {
        mesh.saveSDF(getSDFFilename(filename));
    }
    
    const auto start = std::chrono::high_resolution_clock::now();
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon, filename);
    const auto stop = std::chrono::high_resolution_clock::now();

    mesh.save(addStepCounter(filename, finalStep), finalStep, finalTime);

    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Ran " << name << " in " << duration << " ms" << std::endl;
}

void runWingTest(const Euler& euler, 
                 const FluxSolver* const fluxSolver, 
                 const Reconstruction* const recon, 
                 const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));

    const std::string name = "Wing";
    const std::string filename = testOutDir + name + ".txt";

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
                const std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(velInf, 0.0, 0.0)};
                mesh(idx)[euler.RHO] = rhoInf;
                for(int d = 0; d < SPACEDIM; ++d)
                {
                    mesh(idx)[euler.MOM[d]] = rhoInf * vel[d];
                }
                mesh(idx)[euler.ENE] = euler.getTotalEnergy(rhoInf, vel, pInf);
            }
        }
    }

    if(!mesh.loadSDF(getSDFFilename(filename)))
    {
        assert(mesh.loadSDF("wing.stl"));
        mesh.saveSDF(getSDFFilename(filename));
    }
    
    const auto start = std::chrono::high_resolution_clock::now();
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon, filename);
    const auto stop = std::chrono::high_resolution_clock::now();

    mesh.save(addStepCounter(filename, finalStep), finalStep, finalTime);

    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Ran " << name << " in " << duration << " ms" << std::endl;
}

void runSpaceShuttleTest(const Euler& euler, 
                         const FluxSolver* const fluxSolver, 
                         const Reconstruction* const recon, 
                         const std::array<int, GRIDDIM>& res)
{
    assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));

    const std::string name = "SpaceShuttle";
    const std::string filename = testOutDir + name + ".txt";

    const REAL rhoStar = 0.01841;
    const REAL velStar = -543.027;
    const REAL pStar = 1197.0;
    const REAL xShock = 30.0;

    const IdealGas* const eos = dynamic_cast<const IdealGas* const>(euler.getEquationOfState());
    const REAL gamma = eos->getGamma();
    const REAL a = 1.0/(gamma-1.0)*rhoStar*velStar;
    const REAL b = -0.5*(3.0-gamma)/(gamma-1.0)*rhoStar*velStar*velStar;
    const REAL c = -velStar*(0.5*rhoStar*velStar*velStar + gamma/(gamma-1.0)*pStar);
    const REAL S = (-b + std::sqrt(b*b - 4*a*c)) / (2*a);
    const REAL rhoInf = rhoStar*(S - velStar) / S;
    const REAL pInf = -S*rhoStar*velStar + rhoStar*velStar*velStar + pStar;

    const std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(-65.0, -43.0, -50.0)};
    const std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(35.0, 57.0, 50.0)};
    const REAL finalTime = (hi[0] - lo[0]) / std::fabs(velStar);

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
    #pragma omp parallel for default(none) shared(res, geom, mesh, euler, rhoInf, pInf, rhoStar, velStar, pStar) schedule(static)
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
                REAL rho, p;
                std::array<REAL, SPACEDIM> vel = {SPACEDIM_DECL(0.0, 0.0, 0.0)};
                if(pos[0] < xShock)
                {
                    rho = rhoInf;
                    p = pInf;
                }
                else
                {
                    rho = rhoStar;
                    vel[0] = velStar;
                    p = pStar;
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

    if(!mesh.loadSDF(getSDFFilename(filename)))
    {
        assert(mesh.loadSDF("space-shuttle.stl"));
        mesh.saveSDF(getSDFFilename(filename));
    }

    const REAL cfl = 0.9;
    const REAL outInterval = 0.25 * finalTime;
    const int startStep = 0;
    const REAL startTime = 0.0;
    #ifdef USE_VDB
        const std::string vdbBaseName = testOutDir + name + ".vdb";
        const REAL vdbInterval = 0.01 * finalTime;
    #endif
    
    const auto start = std::chrono::high_resolution_clock::now();
    const int finalStep = solve(euler, finalTime, mesh, bc, fluxSolver, recon, filename, 
                                cfl, outInterval, startStep, startTime
                                #ifdef USE_VDB
                                    , vdbBaseName, vdbInterval
                                #endif
                                );
    const auto stop = std::chrono::high_resolution_clock::now();

    mesh.save(addStepCounter(filename, finalStep), finalStep, finalTime);

    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Ran " << name << " in " << duration << " ms" << std::endl;
}
#endif

#endif