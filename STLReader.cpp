#include "STLReader.H"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>

// Contents of this file were written with the help of ChatGPT

#if GRIDDIM == 3 && defined(USE_RIGID)

bool loadSTL(const std::string& filename, 
             std::vector<std::array<std::array<REAL, GRIDDIM>, 3>>& verts, 
             std::vector<std::array<REAL, GRIDDIM>>& norms)
{
    std::ifstream file(filename, std::ios::binary);
    if(!file.is_open())
    {
        return false;
    }

    char header[80];
    file.read(header, 80);

    uint32_t numTriangles;
    file.read(reinterpret_cast<char*>(&numTriangles), sizeof(uint32_t));

    verts.resize(numTriangles);
    norms.resize(numTriangles);

    std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max())};
    std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(std::numeric_limits<REAL>::min(), std::numeric_limits<REAL>::min(), std::numeric_limits<REAL>::min())};

    std::array<float, GRIDDIM> buf;
    for (uint32_t i = 0; i < numTriangles; ++i)
    {
        file.read(reinterpret_cast<char*>(buf.data()), GRIDDIM * sizeof(float));
        for(int d = 0; d < GRIDDIM; ++d)
        {
            norms[i][d] = static_cast<REAL>(buf[d]);
        }
        for (int v = 0; v < 3; ++v)
        {
            file.read(reinterpret_cast<char*>(buf.data()), GRIDDIM * sizeof(float));
            for(int d = 0; d < GRIDDIM; ++d)
            {
                verts[i][v][d] = static_cast<REAL>(buf[d]);
                lo[d] = std::min(verts[i][v][d], lo[d]);
                hi[d] = std::max(verts[i][v][d], hi[d]);
            }
        }
        file.ignore(2);
    }

    std::cout << "Loaded STL file " << filename << " with bounding box (" << GRIDDIM_TERM(lo[0], << ", " << lo[1], << ", " << lo[2]) << ") (" << GRIDDIM_TERM(hi[0], << ", " << hi[1], << ", " << hi[2]) << ")" << std::endl;

    file.close();
    return true;
}

REAL computeMeshSignedDistance(const std::array<REAL, GRIDDIM>& p,
                               const std::vector<std::array<std::array<REAL, GRIDDIM>, 3>>& verts, 
                               const std::vector<std::array<REAL, GRIDDIM>>& norms)
{
    assert(verts.size() == norms.size());
    REAL signedDistance = std::numeric_limits<REAL>::max();
    for(size_t i = 0; i < verts.size(); ++i)
    {
        const REAL triSignedDistance = computeTriangleSignedDistance(p, verts[i], norms[i]);
        if(std::fabs(triSignedDistance) < std::fabs(signedDistance))
        {
            signedDistance = triSignedDistance;
        }
    }
    return signedDistance;
}

REAL computeTriangleSignedDistance(const std::array<REAL, GRIDDIM>& p, 
                                   const std::array<std::array<REAL, GRIDDIM>, 3>& verts, 
                                   const std::array<REAL, GRIDDIM>& norm)
{
    const std::array<REAL, GRIDDIM>& a = verts[0];
    const std::array<REAL, GRIDDIM>& b = verts[1];
    const std::array<REAL, GRIDDIM>& c = verts[2];

    // Project point onto triangle plane
    std::array<REAL, GRIDDIM> ap = sub(p, a);
    REAL distToPlane = dot(ap, norm);
    std::array<REAL, GRIDDIM> projected = sub(p, scale(norm, distToPlane));

    // Check if projected point lies inside triangle using barycentric coordinates
    std::array<REAL, GRIDDIM> v0 = sub(b, a);
    std::array<REAL, GRIDDIM> v1 = sub(c, a);
    std::array<REAL, GRIDDIM> v2 = sub(projected, a);

    REAL d00 = dot(v0, v0);
    REAL d01 = dot(v0, v1);
    REAL d11 = dot(v1, v1);
    REAL d20 = dot(v2, v0);
    REAL d21 = dot(v2, v1);
    REAL denom = d00 * d11 - d01 * d01;

    #ifdef DEBUG
        assert(std::fabs(denom) > 1e-16);
    #endif

    REAL v = (d11 * d20 - d01 * d21) / denom;
    REAL w = (d00 * d21 - d01 * d20) / denom;
    REAL u = 1.0 - v - w;

    std::array<REAL, GRIDDIM> closest;
    if (u >= 0 && v >= 0 && w >= 0)  // Inside triangle
    {
        closest = projected;
    } 
    else  // Clamp to nearest edge or vertex 
    {
        std::array<REAL, GRIDDIM> cp1 = findClosestPointOnEdge(p, a, b);
        std::array<REAL, GRIDDIM> cp2 = findClosestPointOnEdge(p, b, c);
        std::array<REAL, GRIDDIM> cp3 = findClosestPointOnEdge(p, c, a);

        REAL d1 = len(sub(cp1, p));
        REAL d2 = len(sub(cp2, p));
        REAL d3 = len(sub(cp3, p));

        if (d1 < d2 && d1 < d3)
        {
            closest = cp1;
        }
        else if (d2 < d3)
        {
            closest = cp2;
        }
        else
        {
            closest = cp3;
        }
    }

    // Signed distance
    std::array<REAL, GRIDDIM> dir = sub(p, closest);
    REAL signedDist = len(dir);
    if (dot(dir, norm) < 0)
    {
        signedDist *= -1;
    }

    return signedDist;
}

std::array<REAL, GRIDDIM> findClosestPointOnEdge(const std::array<REAL, GRIDDIM>& p, 
                                                 const std::array<REAL, GRIDDIM>& vert1,
                                                 const std::array<REAL, GRIDDIM>& vert2)
{
    const std::array<REAL, GRIDDIM> ab = sub(vert2, vert1);
    const std::array<REAL, GRIDDIM> ap = sub(p, vert1);
    #ifdef DEBUG
        assert(std::fabs(dot(ab, ab)) > 1e-16);
    #endif
    REAL t = dot(ap, ab) / dot(ab, ab);
    t = std::clamp(t, 0.0, 1.0);
    return add(vert1, scale(ab, t));
}

REAL dot(const std::array<REAL, GRIDDIM>& v1, const std::array<REAL, GRIDDIM>& v2)
{
    return GRIDDIM_TERM(v1[0]*v2[0], + v1[1]*v2[1], + v1[2]*v2[2]);
}

REAL len(const std::array<REAL, GRIDDIM>& v)
{
    return std::sqrt(dot(v, v));
}

std::array<REAL, GRIDDIM> scale(const std::array<REAL, GRIDDIM>& v, const REAL s)
{
    return {GRIDDIM_DECL(v[0]*s, v[1]*s, v[2]*s)};
}

std::array<REAL, GRIDDIM> add(const std::array<REAL, GRIDDIM>& v1, const std::array<REAL, GRIDDIM>& v2)
{
    return {GRIDDIM_DECL(v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])};
}

std::array<REAL, GRIDDIM> sub(const std::array<REAL, GRIDDIM>& v1, const std::array<REAL, GRIDDIM>& v2)
{
    return {GRIDDIM_DECL(v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])};
}

#endif
