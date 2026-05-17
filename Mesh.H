#ifndef MESH_H
#define MESH_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Macros.H"
#include "FileHandler.H"

#if GRIDDIM == 3 && defined(USE_RIGID)
    #include "STLReader.H"
#endif

#ifdef USE_OMP
    #include <omp.h>
#endif

#if GRIDDIM == 3 && defined(USE_VDB)
    #include <openvdb/openvdb.h>
    #include "Euler.H"
#endif

enum BoundaryCondition
{
    TRANSMISSIVE = 0,
    REFLECTIVE = 1,
    PERIODIC = 2
};

class Geometry
{
    public:
        Geometry(const std::array<REAL, GRIDDIM>& lo, 
                 const std::array<REAL, GRIDDIM>& hi, 
                 const std::array<int, GRIDDIM>& res);
        const std::array<REAL, GRIDDIM>& getLo() const;
        const std::array<REAL, GRIDDIM>& getHi() const;
        const std::array<int, GRIDDIM>& getRes() const;
        const std::array<REAL, GRIDDIM>& getDx() const;
        void getPos(std::array<REAL, GRIDDIM>& pos, const std::array<int, GRIDDIM>& idx) const;
        void getIdx(std::array<REAL, GRIDDIM>& idx, const std::array<REAL, GRIDDIM>& pos) const;
        
    private:
        const std::array<REAL, GRIDDIM> m_lo;
        const std::array<REAL, GRIDDIM> m_hi;
        const std::array<int, GRIDDIM> m_res;
        const std::array<REAL, GRIDDIM> m_dx;
};

template<int NVARS>
class DataArray
{
    public:
        DataArray(const std::array<int, GRIDDIM>& res, const int nGhost=0);

        const std::array<int, GRIDDIM>& getRes() const;
        int getNGhost() const;

        void clampIdx(std::array<REAL, GRIDDIM>& idx, const REAL padding=1e-14) const;

        std::array<REAL, NVARS>& operator()(GRIDDIM_DECL(const int i, const int j, const int k));
        const std::array<REAL, NVARS>& operator()(GRIDDIM_DECL(const int i, const int j, const int k)) const;
        std::array<REAL, NVARS> operator()(GRIDDIM_DECL(const REAL i, const REAL j, const REAL k)) const;

        std::array<REAL, NVARS>& operator()(const std::array<int, GRIDDIM>& idx);
        const std::array<REAL, NVARS>& operator()(const std::array<int, GRIDDIM>& idx) const;
        std::array<REAL, NVARS> operator()(const std::array<REAL, GRIDDIM>& idx) const;

        void fillGhost(const std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, 
                       const std::vector<std::array<int, GRIDDIM>>& vecIdx);
        
        void saveData(const std::string& filename, const int nGhost=0);
        bool loadData(const std::string& filename, const int nGhost=0);

    private:
        const std::array<int, GRIDDIM> m_res;
        const int m_nGhost;

        #if GRIDDIM == 1
            std::vector<std::array<REAL, NVARS>> m_data;
        #elif GRIDDIM == 2
            std::vector<std::vector<std::array<REAL, NVARS>>> m_data;
        #else  // GRIDDIM == 3
            std::vector<std::vector<std::vector<std::array<REAL, NVARS>>>> m_data;
        #endif

        std::vector<REAL> m_buffer;
};

template <int NVARS>
class Mesh : public DataArray<NVARS>
{
    public:
        #ifdef USE_RIGID
            static constexpr int SDFNGHOST = 1;
        #endif

        Mesh(const Geometry& geom, const int nGhost=0);

        const Geometry& getGeometry() const;

        void save(const std::string& headerFilename, const int step, const REAL time);
        void writeHeader(const std::string& headerFilename, const int step, const REAL time) const;

        #if GRIDDIM == 3 && defined(USE_VDB)
            void saveToVDB(const std::string& filename, const Euler& euler) const;
        #endif

        static Mesh<NVARS> createFromFile(const std::string& headerFilename, int& step, REAL& time, const int nGhost);

        #ifdef USE_RIGID
            void setSDF(const REAL sdf, const std::array<int, GRIDDIM>& idx);

            bool isRigid(GRIDDIM_DECL(const int i, const int j, const int k)) const;
            bool isRigid(const std::array<int, GRIDDIM>& idx) const;

            std::array<REAL, GRIDDIM> getRigidNormal(GRIDDIM_DECL(const int i, const int j, const int k)) const;

            void fillInternalGhost(const std::vector<std::array<int, GRIDDIM>>& vecIdx, const REAL multiplier=1.1);

            void saveSDF(const std::string& headerFilename);
            void writeSDFHeader(const std::string& headerFilename) const;
            bool loadSDF(const std::string& headerFilename);
            #if GRIDDIM == 3
                bool loadSDFFromSTL(const std::string& stlFilename);
            #endif
        #endif

    private:
        const Geometry m_geom;

        #ifdef USE_RIGID
            DataArray<1> m_sdf;
        #endif
};

template<int NVARS>
DataArray<NVARS>::DataArray(const std::array<int, GRIDDIM>& res, const int nGhost) : m_res(res), m_nGhost(nGhost)
{
    assert(GRIDDIM_TERM(res[0] > 0, && res[1] > 0, && res[2] > 0));
    assert(nGhost >= 0);

    #if GRIDDIM == 1
        m_data.resize(res[0] + 2 * nGhost);
    #elif GRIDDIM == 2
        m_data.resize(res[0] + 2 * nGhost);
        for(int i = 0; i < res[0] + 2 * nGhost; ++i)
        {
            m_data[i].resize(res[1] + 2 * nGhost);
        }
    #else  // GRIDDIM == 3
        m_data.resize(res[0] + 2 * nGhost);
        for(int i = 0; i < res[0] + 2 * nGhost; ++i)
        {
            m_data[i].resize(res[1] + 2 * nGhost);
            for(int j = 0; j < res[1] + 2 * nGhost; ++j)
            {
                m_data[i][j].resize(res[2] + 2 * nGhost);
            }
        }
    #endif
}

template<int NVARS>
const std::array<int, GRIDDIM>& DataArray<NVARS>::getRes() const 
{
    return m_res;
}

template<int NVARS>
int DataArray<NVARS>::getNGhost() const 
{
    return m_nGhost;
}

template<int NVARS>
void DataArray<NVARS>::clampIdx(std::array<REAL, GRIDDIM>& idx, const REAL padding) const
{
    for(int d = 0; d < GRIDDIM; ++d)
    {
        idx[d] = std::clamp(idx[d], (REAL) -m_nGhost + padding, (REAL) m_res[d] + m_nGhost - 1 - padding);
    }
}

template<int NVARS>
std::array<REAL, NVARS>& DataArray<NVARS>::operator()(GRIDDIM_DECL(const int i, const int j, const int k))
{
    #if GRIDDIM == 1
        return m_data[i + m_nGhost];
    #elif GRIDDIM == 2
        return m_data[i + m_nGhost][j + m_nGhost];
    #else  // GRIDDIM == 3
        return m_data[i + m_nGhost][j + m_nGhost][k + m_nGhost];
    #endif
}

template<int NVARS>
const std::array<REAL, NVARS>& DataArray<NVARS>::operator()(GRIDDIM_DECL(const int i, const int j, const int k)) const
{
    #if GRIDDIM == 1
        return m_data[i + m_nGhost];
    #elif GRIDDIM == 2
        return m_data[i + m_nGhost][j + m_nGhost];
    #else  // GRIDDIM == 3
        return m_data[i + m_nGhost][j + m_nGhost][k + m_nGhost];
    #endif
}

template<int NVARS>
std::array<REAL, NVARS> DataArray<NVARS>::operator()(GRIDDIM_DECL(const REAL i, const REAL j, const REAL k)) const
{
    return (*this)({GRIDDIM_DECL(i, j, k)});
}

template<int NVARS>
std::array<REAL, NVARS>& DataArray<NVARS>::operator()(const std::array<int, GRIDDIM>& idx)
{
    return (*this)(GRIDDIM_DECL(idx[0], idx[1], idx[2]));
}

template<int NVARS>
const std::array<REAL, NVARS>& DataArray<NVARS>::operator()(const std::array<int, GRIDDIM>& idx) const
{
    return (*this)(GRIDDIM_DECL(idx[0], idx[1], idx[2]));
}

template<int NVARS>
std::array<REAL, NVARS> DataArray<NVARS>::operator()(const std::array<REAL, GRIDDIM>& idx) const
{
    std::array<int, GRIDDIM> idxLo;
    std::array<std::array<REAL, 2>, GRIDDIM> fac;
    for(int d = 0; d < GRIDDIM; ++d)
    {
        idxLo[d] = static_cast<int>(std::floor(idx[d]));
        fac[d][1] = idx[d] - idxLo[d];
        fac[d][0] = 1.0 - fac[d][1];
    }

    std::array<REAL, NVARS> interpData;
    for(int v = 0; v < NVARS; ++v)
    {
        interpData[v] = 0.0;
    }
    
    for(int iOffset = 0; iOffset < 2; ++iOffset)
    {
        #if GRIDDIM >= 2
        for(int jOffset = 0; jOffset < 2; ++jOffset)
        #endif
        {
            #if GRIDDIM == 3
            for(int kOffset = 0; kOffset < 2; ++kOffset)
            #endif
            {
                const std::array<REAL, NVARS> data = (*this)(GRIDDIM_DECL(idxLo[0] + iOffset, 
                                                                            idxLo[1] + jOffset, 
                                                                            idxLo[2] + kOffset));
                for(int v = 0; v < NVARS; ++v)
                {
                    interpData[v] += data[v] * GRIDDIM_TERM(fac[0][iOffset], 
                                                            * fac[1][jOffset], 
                                                            * fac[2][kOffset]);
                }
            }
        }
    }

    return interpData;
}

template<int NVARS>
void DataArray<NVARS>::fillGhost(const std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, 
                                 const std::vector<std::array<int, GRIDDIM>>& vecIdx)
{
    std::array<int, GRIDDIM> idxSource, idxTarget;
    std::array<int, GRIDDIM> nIter;
    for(int d = 0; d < GRIDDIM; ++d)
    {
        nIter = m_res;
        nIter[d] = m_nGhost;
        for(int s = 0; s < 2; ++s)
        {
            for(int i = 0; i < nIter[0]; ++i)
            {
                #if GRIDDIM >= 2
                for(int j = 0; j < nIter[1]; ++j)
                #endif
                {
                    #if GRIDDIM == 3
                    for(int k = 0; k < nIter[2]; ++k)
                    #endif
                    {
                        idxTarget = {GRIDDIM_DECL(i, j, k)};
                        const int g = idxTarget[d];
                        idxTarget[d] = (s == 0 ? -1 - g : m_res[d] + g);
                        idxSource = idxTarget;
                        if(bc[s][d] == TRANSMISSIVE)
                        {
                            idxSource[d] = (s == 0 ? 0 : m_res[d] - 1);
                        }
                        else if(bc[s][d] == REFLECTIVE)
                        {
                            idxSource[d] = (s == 0 ? g : m_res[d] - 1 - g);
                        }
                        else if(bc[s][d] == PERIODIC)
                        {
                            idxSource[d] =  (s == 0 ? m_res[d] - 1 - g : g);
                        }
                        (*this)(idxTarget) = (*this)(idxSource);
                        if(bc[s][d] == REFLECTIVE)
                        {
                            for(size_t v = 0; v < vecIdx.size(); ++v)
                            {
                                (*this)(idxTarget)[vecIdx[v][d]] *= -1.0;
                            }
                        }
                    }
                }
            }
        }
    }
}

template<int NVARS>
void DataArray<NVARS>::saveData(const std::string& filename, const int nGhost)
{
    assert(0 <= nGhost && nGhost <= m_nGhost);
    const int nCellsInclGhost = GRIDDIM_TERM((m_res[0] + 2 * nGhost), * (m_res[1] + 2 * nGhost), * (m_res[2] + 2 * nGhost));
    m_buffer.resize(nCellsInclGhost * NVARS);
    size_t idx = 0;
    for(int i = -nGhost; i < m_res[0] + nGhost; ++i)
    {
        #if GRIDDIM >= 2
        for(int j = -nGhost; j < m_res[1] + nGhost; ++j)
        #endif
        {
            #if GRIDDIM == 3
            for(int k = -nGhost; k < m_res[2] + nGhost; ++k)
            #endif
            {
                const std::array<REAL, NVARS>& data = (*this)(GRIDDIM_DECL(i, j, k));
                std::copy_n(data.begin(), NVARS, &m_buffer[idx]);
                idx += NVARS;
            }
        }
    }
    createDirForFile(filename);
    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<const char*>(m_buffer.data()), m_buffer.size() * sizeof(REAL));
    file.close();
}

template<int NVARS>
bool DataArray<NVARS>::loadData(const std::string& filename, const int nGhost)
{
    assert(0 <= nGhost && nGhost <= m_nGhost);
    const int nCellsInclGhost = GRIDDIM_TERM((m_res[0] + 2 * nGhost), * (m_res[1] + 2 * nGhost), * (m_res[2] + 2 * nGhost));
    m_buffer.resize(nCellsInclGhost * NVARS);
    std::ifstream in(filename, std::ios::binary);
    if(!in.is_open())
    {
        return false;
    }
    in.read(reinterpret_cast<char*>(m_buffer.data()), m_buffer.size() * sizeof(REAL));
    in.close();
    size_t idx = 0;
    for(int i = -nGhost; i < m_res[0] + nGhost; ++i)
    {
        #if GRIDDIM >= 2
        for(int j = -nGhost; j < m_res[1] + nGhost; ++j)
        #endif
        {
            #if GRIDDIM == 3
            for(int k = -nGhost; k < m_res[2] + nGhost; ++k)
            #endif
            {
                std::array<REAL, NVARS>& data = (*this)(GRIDDIM_DECL(i, j, k));
                std::copy_n(&m_buffer[idx], NVARS, data.begin());
                idx += NVARS;
            }
        }
    }
    return true;
}

template<int NVARS>
Mesh<NVARS>::Mesh(const Geometry& geom, const int nGhost) : DataArray<NVARS>(geom.getRes(), nGhost), m_geom(geom) 
#ifdef USE_RIGID
    , m_sdf(geom.getRes(), SDFNGHOST)
#endif
{
    assert(nGhost >= 0);
}

template<int NVARS>
const Geometry& Mesh<NVARS>::getGeometry() const 
{
    return m_geom;
}

template<int NVARS>
void Mesh<NVARS>::save(const std::string& headerFilename, const int step, const REAL time)
{
    writeHeader(headerFilename, step, time);
    this->saveData(getDataFilename(headerFilename));
}

template<int NVARS>
void Mesh<NVARS>::writeHeader(const std::string& headerFilename, const int step, const REAL time) const
{
    createDirForFile(headerFilename);
    std::ofstream file(headerFilename);
    assert(file.is_open());
    const std::array<REAL, GRIDDIM>& lo = m_geom.getLo();
    const std::array<REAL, GRIDDIM>& hi = m_geom.getHi();
    const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
    file << step << std::endl;
    file << time << std::endl;
    file << GRIDDIM_TERM(lo[0], << " " << lo[1], << " " << lo[2]) << std::endl;
    file << GRIDDIM_TERM(hi[0], << " " << hi[1], << " " << hi[2]) << std::endl;
    file << GRIDDIM_TERM(res[0], << " " << res[1], << " " << res[2]) << std::endl;
    file << NVARS << std::endl;
    file << removePath(getDataFilename(headerFilename)) << std::endl;
    file.close();
}

#if GRIDDIM == 3 && defined(USE_VDB)
template<int NVARS>
void Mesh<NVARS>::saveToVDB(const std::string& filename, const Euler& euler) const
{
    openvdb::initialize();

    const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
    const std::array<REAL, GRIDDIM>& dx = m_geom.getDx();
    
    openvdb::GridPtrVec grids(Euler::NOUTVARS);
    std::vector<openvdb::FloatGrid::Accessor> floatAccessors;
    std::vector<openvdb::Vec3SGrid::Accessor> vecAccessors;
    floatAccessors.reserve(Euler::NOUTVARS);
    vecAccessors.reserve(Euler::NOUTVARS);
    std::array<int, Euler::NOUTVARS> accessorIdx;
    for(int n = 0; n < Euler::NOUTVARS; ++n)
    {
        if(Euler::ISVEC[n])
        {
            auto grid = openvdb::Vec3SGrid::create();
            grids[n] = grid;
            vecAccessors.push_back(grid->getAccessor());
            accessorIdx[n] = vecAccessors.size() - 1;
        }
        else
        {
            auto grid = openvdb::FloatGrid::create();
            grids[n] = grid;
            auto accessor = grid->getAccessor();
            floatAccessors.push_back(grid->getAccessor());
            accessorIdx[n] = floatAccessors.size() - 1;
        }
        grids[n]->setTransform(openvdb::math::Transform::createLinearTransform(dx[0]));  // Assumes that voxels have equal side lengths
        grids[n]->setName(std::string(Euler::OUTVARNAMES[n]));
    }

    openvdb::Coord ijk;
    int& i = ijk[0];
    int& j = ijk[1];
    int& k = ijk[2];

    for(i = 0; i < res[0]; ++i)
    {
        for(j = 0; j < res[1]; ++j)
        {
            for(k = 0; k < res[2]; ++k)
            {
                #ifdef USE_RIGID
                if(!isRigid(GRIDDIM_DECL(i, j, k)))
                #endif
                {
                    const std::array<REAL, NVARS>& data = (*this)(GRIDDIM_DECL(i, j, k));
                    std::array<std::vector<REAL>, Euler::NOUTVARS> outData;
                    euler.getOutState(outData, data);
                    for(int n = 0; n < Euler::NOUTVARS; ++n)
                    {
                        if(Euler::ISVEC[n])
                        {
                            #ifdef DEBUG
                                assert(outData[n].size() == SPACEDIM);
                            #endif
                            std::array<REAL, GRIDDIM> vec = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
                            for(int d = 0; d < SPACEDIM; ++d)
                            {
                                vec[d] = outData[n][d];
                            }
                            vecAccessors[accessorIdx[n]].setValue(ijk, openvdb::Vec3s(GRIDDIM_DECL(vec[0], vec[1], vec[2])));
                        }
                        else
                        {
                            #ifdef DEBUG
                                assert(outData[n].size() == 1);
                            #endif
                            floatAccessors[accessorIdx[n]].setValue(ijk, outData[n][0]);
                        }
                    }
                }
            }
        }
    }

    createDirForFile(filename);
    openvdb::io::File file(filename);
 
    file.write(grids);
    file.close();
}
#endif

template<int NVARS>
Mesh<NVARS> Mesh<NVARS>::createFromFile(const std::string& headerFilename, int& step, REAL& time, const int nGhost)
{
    assert(nGhost >= 0);

    std::ifstream file(headerFilename);
    assert(file.is_open());
    std::string stepLine, timeLine, loLine, hiLine, resLine, NVARSLine, dataFilename;
    std::getline(file, stepLine);
    std::getline(file, timeLine);
    std::getline(file, loLine);
    std::getline(file, hiLine);
    std::getline(file, resLine);
    std::getline(file, NVARSLine);
    std::getline(file, dataFilename);
    file.close();
    std::istringstream stepISS(stepLine);
    std::istringstream timeISS(timeLine);
    std::istringstream loISS(loLine);
    std::istringstream hiISS(hiLine);
    std::istringstream resISS(resLine);
    std::istringstream NVARSISS(NVARSLine);
    stepISS >> step;
    timeISS >> time;
    std::array<REAL, GRIDDIM> lo, hi;
    std::array<int, GRIDDIM> res;
    for(int d = 0; d < GRIDDIM; ++d)
    {
        loISS >> lo[d];
        hiISS >> hi[d];
        resISS >> res[d];
    }
    int NVARSFromFile;
    NVARSISS >> NVARSFromFile;
    assert(NVARSFromFile == NVARS);
    const Geometry geom(lo, hi, res);
    Mesh<NVARS> mesh(geom, nGhost);
    assert(mesh.loadData(addPath(headerFilename, dataFilename)));
    return mesh;
}

#ifdef USE_RIGID
template<int NVARS>
void Mesh<NVARS>::setSDF(const REAL sdf, const std::array<int, GRIDDIM>& idx)
{
    m_sdf(idx)[0] = sdf;
}

template<int NVARS>
bool Mesh<NVARS>::isRigid(GRIDDIM_DECL(const int i, const int j, const int k)) const
{
    return m_sdf(GRIDDIM_DECL(i, j, k))[0] < 0;
}

template<int NVARS>
bool Mesh<NVARS>::isRigid(const std::array<int, GRIDDIM>& idx) const
{
    return m_sdf(idx)[0] < 0;
}

template<int NVARS>
std::array<REAL, GRIDDIM> Mesh<NVARS>::getRigidNormal(GRIDDIM_DECL(const int i, const int j, const int k)) const
{
    const std::array<REAL, GRIDDIM>& dx = m_geom.getDx();
    std::array<REAL, GRIDDIM> normal;
    std::array<int, GRIDDIM> offset;
    for(int d = 0; d < GRIDDIM; ++d)
    {
        offset = {GRIDDIM_DECL(0, 0, 0)};
        offset[d] = 1;
        const REAL sdfLo = m_sdf(GRIDDIM_DECL(i - offset[0],
                                                j - offset[1],
                                                k - offset[2]))[0];
        const REAL sdfHi = m_sdf(GRIDDIM_DECL(i + offset[0],
                                                j + offset[1],
                                                k + offset[2]))[0];
        normal[d] = 0.5 * (sdfHi - sdfLo) / dx[d];
    }
    const REAL normalMag = std::sqrt(GRIDDIM_TERM(normal[0]*normal[0], + normal[1]*normal[1], + normal[2]*normal[2]));
    #ifdef DEBUG
        if(normalMag <= 1e-16)
        {
            std::cout << "Invalid normal with length = " << normalMag << std::endl;
        }
    #endif
    for(int d = 0; d < GRIDDIM; ++d)
    {
        normal[d] /= normalMag;
    }
    return normal;
}

template<int NVARS>
void Mesh<NVARS>::fillInternalGhost(const std::vector<std::array<int, GRIDDIM>>& vecIdx, const REAL multiplier)
{
    const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
    const int nGhost = DataArray<NVARS>::getNGhost();
    const std::array<REAL, GRIDDIM>& dx = m_geom.getDx();
    const REAL cellSize = std::sqrt(GRIDDIM_TERM(dx[0] * dx[0],
                                                 + dx[1] * dx[1],
                                                 + dx[2] * dx[2]));
    #ifdef USE_OMP
    #pragma omp parallel for default(none) shared(vecIdx, multiplier, res, nGhost, dx, cellSize) schedule(dynamic)
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
                const REAL sdf = m_sdf(GRIDDIM_DECL(i, j, k))[0];
                if(isRigid(idx) && sdf > -multiplier * nGhost * cellSize)
                {
                    std::array<REAL, GRIDDIM> pos, probePos;
                    m_geom.getPos(pos, idx);
                    const std::array<REAL, GRIDDIM> normal = getRigidNormal(GRIDDIM_DECL(i, j, k));
                    for(int d = 0; d < GRIDDIM; ++d)
                    {
                        probePos[d] = pos[d] - sdf * normal[d] * multiplier;
                    }
                    std::array<REAL, GRIDDIM> probeIdx;
                    m_geom.getIdx(probeIdx, probePos);
                    this->clampIdx(probeIdx);
                    std::array<REAL, NVARS> probeData = (*this)(probeIdx);
                    std::array<REAL*, GRIDDIM> vec;
                    for(size_t v = 0; v < vecIdx.size(); ++v)
                    {
                        vec = {GRIDDIM_DECL(&probeData[vecIdx[v][0]], 
                                            &probeData[vecIdx[v][1]], 
                                            &probeData[vecIdx[v][2]])};
                        const REAL dot = GRIDDIM_TERM(normal[0] * (*vec[0]), 
                                                        + normal[1] * (*vec[1]), 
                                                        + normal[2] * (*vec[2]));
                        for(int d = 0; d < GRIDDIM; ++d)
                        {
                            (*vec[d]) -= 2 * normal[d] * dot;
                        }
                    }
                    for(int v = 0; v < NVARS; ++v)
                    {
                        (*this)(GRIDDIM_DECL(i, j, k))[v] = probeData[v];
                    }
                }
            }
        }
    }
}

template<int NVARS>
void Mesh<NVARS>::saveSDF(const std::string& headerFilename)
{
    writeSDFHeader(headerFilename);
    m_sdf.saveData(getDataFilename(headerFilename), SDFNGHOST);
}

template<int NVARS>
void Mesh<NVARS>::writeSDFHeader(const std::string& headerFilename) const
{
    createDirForFile(headerFilename);
    std::ofstream file(headerFilename);
    assert(file.is_open());
    const std::array<REAL, GRIDDIM>& lo = m_geom.getLo();
    const std::array<REAL, GRIDDIM>& hi = m_geom.getHi();
    const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
    file << GRIDDIM_TERM(lo[0], << " " << lo[1], << " " << lo[2]) << std::endl;
    file << GRIDDIM_TERM(hi[0], << " " << hi[1], << " " << hi[2]) << std::endl;
    file << GRIDDIM_TERM(res[0], << " " << res[1], << " " << res[2]) << std::endl;
    file << removePath(getDataFilename(headerFilename)) << std::endl;
    file.close();
}

template<int NVARS>
bool Mesh<NVARS>::loadSDF(const std::string& headerFilename)
{
    const std::array<REAL, GRIDDIM>& lo = m_geom.getLo();
    const std::array<REAL, GRIDDIM>& hi = m_geom.getHi();
    const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
    std::ifstream file(headerFilename);
    if(!file.is_open())
    {
        return false;
    }
    std::string loLine, hiLine, resLine, dataFilename;
    std::getline(file, loLine);
    std::getline(file, hiLine);
    std::getline(file, resLine);
    std::getline(file, dataFilename);
    file.close();
    std::istringstream loISS(loLine);
    std::istringstream hiISS(hiLine);
    std::istringstream resISS(resLine);
    REAL loFromFile, hiFromFile;
    int resFromFile;
    bool validFile = true;
    for(int d = 0; d < GRIDDIM; ++d)
    {
        loISS >> loFromFile;
        hiISS >> hiFromFile;
        resISS >> resFromFile;
        validFile = validFile && std::fabs(loFromFile - lo[d]) < 1e-16;
        validFile = validFile && std::fabs(hiFromFile - hi[d]) < 1e-16;
        validFile = validFile && resFromFile == res[d];
    }
    if(validFile && m_sdf.loadData(addPath(headerFilename, dataFilename), SDFNGHOST))
    {
        return true;
    }
    else
    {
        #if GRIDDIM == 3
            return loadSDFFromSTL(headerFilename);
        #else
            return false;
        #endif
    }
}

#if GRIDDIM == 3
template<int NVARS>
bool Mesh<NVARS>::loadSDFFromSTL(const std::string& stlFilename)
{
    const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
    std::vector<std::array<std::array<REAL, GRIDDIM>, 3>> verts;
    std::vector<std::array<REAL, GRIDDIM>> norms;
    if(!loadSTL(stlFilename, verts, norms))
    {
        return false;
    }
    std::cout << "Computing signed distance function..." << std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for default(none) shared(verts, norms, res) schedule(static)
    #endif
    for(int i = -SDFNGHOST; i < res[0] + SDFNGHOST; ++i)
    {
        #if GRIDDIM >= 2
        for(int j = -SDFNGHOST; j < res[1] + SDFNGHOST; ++j)
        #endif
        {
            #if GRIDDIM == 3
            for(int k = -SDFNGHOST; k < res[2] + SDFNGHOST; ++k)
            #endif
            {
                const std::array<int, GRIDDIM> idx = {GRIDDIM_DECL(i, j, k)};
                std::array<REAL, GRIDDIM> pos;
                m_geom.getPos(pos, idx);
                setSDF(computeMeshSignedDistance(pos, verts, norms), idx);
            }
        }
    }
    return true;
}
#endif

#endif

#endif