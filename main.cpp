#include <chrono>
#include <filesystem>
#include <iostream>
#include <sstream>

#include "EquationOfState.H"
#include "Euler.H"
#include "FluxSolver.H"
#include "FileHandler.H"
#include "Mesh.H"
#include "Reconstruction.H"
#include "Solver.H"
#include "TestProblems.H"

void parseSettings(std::string& settingsName, std::string& initHeaderName, REAL& finalTime, 
                   std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, REAL& gamma, REAL& cfl,
                   #ifdef USE_RIGID
                       std::string& sdfName,
                   #endif
                   std::string& outHeaderBaseName, REAL& outInterval
                   #ifdef USE_GRAVITY
                       , std::array<REAL, SPACEDIM>& g
                   #endif
                   #if GRIDDIM == 3 && defined(USE_VDB)
                       , std::string& vdbBaseName, REAL& vdbInterval, int& vdbStartIdx
                   #endif
                )
{
    std::ifstream file(settingsName);
    assert(file.is_open());

    std::string line;
    
    while (std::getline(file, line))
    {
        std::stringstream lineStream(line);
        std::string segment;
        std::vector<std::string> segList;
        while(std::getline(lineStream, segment, '='))
        {
            segment.erase(std::remove_if(segment.begin(), segment.end(),
                          [](char c) {return c == '\n' || c == '\r';}), 
                          segment.end());
            auto start = std::find_if_not(segment.begin(), segment.end(),
                                          [](unsigned char ch) {return std::isspace(ch);});
            auto end = std::find_if_not(segment.rbegin(), segment.rend(),
                                        [](unsigned char ch) {return std::isspace(ch);}).base();
            segList.push_back((start < end) ? std::string(start, end) : "");
        }
        if(segList.size() == 2)
        {
            if(segList[0] == "init_header_fname")
            {
                initHeaderName = segList[1];
            }
            else if(segList[0] == "final_time")
            {
                std::istringstream iss(segList[1]);
                iss >> finalTime;
            }
            else if(segList[0] == "bc_lo")
            {
                std::istringstream iss(segList[1]);
                int bcLo;
                for(int d = 0; d < GRIDDIM; ++d)
                {
                    iss >> bcLo;
                    assert(bcLo == BoundaryCondition::TRANSMISSIVE 
                           || bcLo == BoundaryCondition::REFLECTIVE 
                           || bcLo == BoundaryCondition::PERIODIC);
                    bc[0][d] = static_cast<BoundaryCondition>(bcLo);
                }
            }
            else if(segList[0] == "bc_hi")
            {
                std::istringstream iss(segList[1]);
                int bcHi;
                for(int d = 0; d < GRIDDIM; ++d)
                {
                    iss >> bcHi;
                    assert(bcHi == BoundaryCondition::TRANSMISSIVE 
                           || bcHi == BoundaryCondition::REFLECTIVE 
                           || bcHi == BoundaryCondition::PERIODIC);
                    bc[1][d] = static_cast<BoundaryCondition>(bcHi);
                }
            }
            else if(segList[0] == "gamma")
            {
                std::istringstream iss(segList[1]);
                iss >> gamma;
            }
            else if(segList[0] == "cfl")
            {
                std::istringstream iss(segList[1]);
                iss >> cfl;
            }
            #ifdef USE_RIGID
                else if(segList[0] == "sdf_header_fname")
                {
                    sdfName = segList[1];
                }
            #endif
            else if(segList[0] == "out_header_base_fname")
            {
                outHeaderBaseName = segList[1];
            }
            else if(segList[0] == "out_interval")
            {
                std::istringstream iss(segList[1]);
                iss >> outInterval;
            }
            #ifdef USE_GRAVITY
                else if(segList[0] == "g")
                {
                    std::istringstream iss(segList[1]);
                    for(int d = 0; d < GRIDDIM; ++d)
                    {
                        iss >> g[d];
                    }
                }
            #endif
            #if GRIDDIM == 3 && defined(USE_VDB)
                else if(segList[0] == "vdb_base_fname")
                {
                    vdbBaseName = segList[1];
                }
                else if(segList[0] == "vdb_interval")
                {
                    std::istringstream iss(segList[1]);
                    iss >> vdbInterval;
                }
                else if(segList[0] == "vdb_start_idx")
                {
                    std::istringstream iss(segList[1]);
                    iss >> vdbStartIdx;
                }
            #endif
        }
    }

    assert(!initHeaderName.empty());
    assert(finalTime > 1e-16);
    if(outHeaderBaseName.empty())
    {
        outHeaderBaseName = removeStepCounter(initHeaderName);
    }
    #if GRIDDIM == 3 && defined(USE_VDB)
        if(vdbBaseName.empty())
        {
            vdbBaseName = setExtension(removeStepCounter(initHeaderName), ".vdb");
        }

    #endif
    
    file.close();
}

int main(int argc, char *argv[])
{
    std::string settingsName, initHeaderName, outHeaderBaseName;
    REAL finalTime = 0.0;
    REAL gamma = 1.4;
    REAL cfl = 0.9;
    REAL outInterval = 0.0;
    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for(int d = 0; d < GRIDDIM; ++d)
    {
        bc[0][d] = BoundaryCondition::TRANSMISSIVE;
        bc[1][d] = BoundaryCondition::TRANSMISSIVE;
    }
    #ifdef USE_RIGID
        std::string sdfName;
    #endif
    #ifdef USE_GRAVITY
        std::array<REAL, SPACEDIM> g = {SPACEDIM_DECL(0.0, 0.0, 0.0)};
    #endif
    #if GRIDDIM == 3 && defined(USE_VDB)
        std::string vdbBaseName;
        int vdbStartIdx = 0;
        REAL vdbInterval = 0.0;
    #endif

    if(argc >= 2)
    {
        settingsName = argv[1];
        parseSettings(settingsName, initHeaderName, finalTime, bc, gamma, cfl,
                      #ifdef USE_RIGID
                          sdfName,
                      #endif
                      outHeaderBaseName, outInterval
                      #ifdef USE_GRAVITY
                         , g
                      #endif
                      #if GRIDDIM == 3 && defined(USE_VDB)
                          , vdbBaseName, vdbInterval, vdbStartIdx
                      #endif
                      );
    }

    const IdealGas eos(gamma);
    const Euler euler(&eos);
    const HLLCSolver fluxSolver(euler);
    const MUSCLHancock recon(euler);

    if(!settingsName.empty())
    {
        const int nGhost = 1 + recon.getStencilSize();
        int startStep;
        REAL startTime;
        Mesh<Euler::NVARS> mesh = Mesh<Euler::NVARS>::createFromFile(addPath(settingsName, initHeaderName), startStep, startTime, nGhost);
        #ifdef USE_RIGID
            if(!sdfName.empty())
            {
                std::cout << addPath(settingsName, sdfName) << std::endl;
                assert(mesh.loadSDF(addPath(settingsName, sdfName)));
            }
        #endif
        const std::string outHeaderBaseNameWPath = addPath(settingsName, outHeaderBaseName);
        const auto start = std::chrono::high_resolution_clock::now();
        const int finalStep = solve(euler, finalTime, mesh, bc, &fluxSolver, &recon, outHeaderBaseNameWPath, cfl, outInterval, startStep, startTime
                                    #ifdef USE_GRAVITY
                                        , &g
                                    #endif
                                    #if GRIDDIM == 3 && defined(USE_VDB)
                                        , vdbBaseName, vdbInterval, vdbStartIdx
                                    #endif
                                    );
        const auto stop = std::chrono::high_resolution_clock::now();
        mesh.save(addStepCounter(outHeaderBaseNameWPath, finalStep), finalStep, finalTime);
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        const std::filesystem::path path(outHeaderBaseNameWPath);
        const std::string name = path.stem().string();
        std::cout << "Ran " << name << " in " << duration << " ms" << std::endl;
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
                const std::array<int, GRIDDIM> shockReflectionRes = {2 * res[0], res[1]};
                runShockReflectionTest(euler, &fluxSolver, &recon, shockReflectionRes);
            #endif
        #endif

        #if GRIDDIM == 3 && defined(USE_RIGID)
            const std::array<int, GRIDDIM> hypersonicSphereRes = {GRIDDIM_DECL(res[0], 2 * res[1], 2 * res[2])};
            runHypersonicSphereTest(euler, &fluxSolver, &recon, hypersonicSphereRes, false);
            runHypersonicSphereTest(euler, &fluxSolver, &recon, hypersonicSphereRes, true);

            const std::array<int, GRIDDIM> wingRes = {GRIDDIM_DECL(2 * res[0], res[1], res[2])};
            runWingTest(euler, &fluxSolver, &recon, wingRes);

            const int spaceShuttleRes1D = std::max(512, res[0]);
            const std::array<int, GRIDDIM> spaceShuttleRes = {GRIDDIM_DECL(spaceShuttleRes1D, spaceShuttleRes1D, spaceShuttleRes1D)};
            runSpaceShuttleTest(euler, &fluxSolver, &recon, spaceShuttleRes);
        #endif
    }

    return 0;
}