#include "Mesh.H"

Geometry::Geometry(const std::array<REAL, GRIDDIM>& lo, 
                   const std::array<REAL, GRIDDIM>& hi, 
                   const std::array<int, GRIDDIM>& res)
    : m_lo(lo), m_hi(hi), m_res(res), m_dx({GRIDDIM_DECL((hi[0] - lo[0]) / (REAL) res[0], 
                                                         (hi[1] - lo[1]) / (REAL) res[1], 
                                                         (hi[2] - lo[2]) / (REAL) res[2])})
{

}

const std::array<REAL, GRIDDIM>& Geometry::getLo() const 
{
    return m_lo;
}

const std::array<REAL, GRIDDIM>& Geometry::getHi() const 
{
    return m_hi;
}

const std::array<int, GRIDDIM>& Geometry::getRes() const 
{
    return m_res;
}

const std::array<REAL, GRIDDIM>& Geometry::getDx() const 
{
    return m_dx;
}

void Geometry::getPos(std::array<REAL, GRIDDIM>& pos, const std::array<int, GRIDDIM>& idx) const 
{
    for (int d = 0; d < GRIDDIM; ++d)
    {
        pos[d] = m_lo[d] + ((REAL) idx[d] + 0.5) * m_dx[d];
    }
}

void Geometry::getIdx(std::array<REAL, GRIDDIM>& idx, const std::array<REAL, GRIDDIM>& pos) const 
{
    for (int d = 0; d < GRIDDIM; ++d)
    {
        idx[d] = (pos[d] - m_lo[d]) / m_dx[d] - 0.5;
    }
}
