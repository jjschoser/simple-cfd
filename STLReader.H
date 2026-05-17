#ifndef STLReader_H
#define STLReader_H

#include <array>
#include <string>
#include <vector>

#include "Macros.H"

#if GRIDDIM == 3 && defined(USE_RIGID)

bool loadSTL(const std::string& filename, 
             std::vector<std::array<std::array<REAL, GRIDDIM>, 3>>& verts, 
             std::vector<std::array<REAL, GRIDDIM>>& norms);

REAL computeMeshSignedDistance(const std::array<REAL, GRIDDIM>& p,
                               const std::vector<std::array<std::array<REAL, GRIDDIM>, 3>>& verts, 
                               const std::vector<std::array<REAL, GRIDDIM>>& norms);

REAL computeTriangleSignedDistance(const std::array<REAL, GRIDDIM>& p, 
                                   const std::array<std::array<REAL, GRIDDIM>, 3>& verts, 
                                   const std::array<REAL, GRIDDIM>& norm);

std::array<REAL, GRIDDIM> findClosestPointOnEdge(const std::array<REAL, GRIDDIM>& p, 
                                                 const std::array<REAL, GRIDDIM>& vert1,
                                                 const std::array<REAL, GRIDDIM>& vert2);

REAL dot(const std::array<REAL, GRIDDIM>& v1, const std::array<REAL, GRIDDIM>& v2);
REAL len(const std::array<REAL, GRIDDIM>& v);
std::array<REAL, GRIDDIM> scale(const std::array<REAL, GRIDDIM>& v, const REAL s);
std::array<REAL, GRIDDIM> add(const std::array<REAL, GRIDDIM>& v1, const std::array<REAL, GRIDDIM>& v2);
std::array<REAL, GRIDDIM> sub(const std::array<REAL, GRIDDIM>& v1, const std::array<REAL, GRIDDIM>& v2);

#endif

#endif