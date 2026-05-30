#include "catch2/catch.hpp"
#include "helpers.H"
#include "Mesh.H"

#include <cmath>
#include <filesystem>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
// Geometry: cell centres, dx, index-to-position inverse
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Geometry: getDx matches (hi-lo)/res", "[geometry]")
{
    std::array<REAL, GRIDDIM> lo = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
    std::array<REAL, GRIDDIM> hi = {GRIDDIM_DECL(1.0, 2.0, 3.0)};
    std::array<int, GRIDDIM>  res = {GRIDDIM_DECL(10, 20, 30)};
    Geometry geom(lo, hi, res);

    auto dx = geom.getDx();
    REQUIRE(dx[0] == Approx(0.1).epsilon(1e-12));
    #if GRIDDIM >= 2
        REQUIRE(dx[1] == Approx(0.1).epsilon(1e-12));
    #endif
    #if GRIDDIM >= 3
        REQUIRE(dx[2] == Approx(0.1).epsilon(1e-12));
    #endif
}

TEST_CASE("Geometry: getPos gives cell centres at lo + (i+0.5)*dx", "[geometry]")
{
    std::array<REAL, GRIDDIM> lo  = {GRIDDIM_DECL(1.0, 2.0, 3.0)};
    std::array<REAL, GRIDDIM> hi  = {GRIDDIM_DECL(2.0, 4.0, 6.0)};
    std::array<int, GRIDDIM>  res = {GRIDDIM_DECL(4, 4, 4)};
    Geometry geom(lo, hi, res);
    auto dx = geom.getDx();  // {0.25, 0.5, 0.75}

    // First cell centre
    std::array<int, GRIDDIM>  idx0 = {GRIDDIM_DECL(0, 0, 0)};
    std::array<REAL, GRIDDIM> pos0;
    geom.getPos(pos0, idx0);
    REQUIRE(pos0[0] == Approx(lo[0] + 0.5 * dx[0]).epsilon(1e-12));
    #if GRIDDIM >= 2
        REQUIRE(pos0[1] == Approx(lo[1] + 0.5 * dx[1]).epsilon(1e-12));
    #endif
    #if GRIDDIM >= 3
        REQUIRE(pos0[2] == Approx(lo[2] + 0.5 * dx[2]).epsilon(1e-12));
    #endif

    // Last cell centre
    std::array<int, GRIDDIM>  idxN = {GRIDDIM_DECL(3, 3, 3)};
    std::array<REAL, GRIDDIM> posN;
    geom.getPos(posN, idxN);
    REQUIRE(posN[0] == Approx(hi[0] - 0.5 * dx[0]).epsilon(1e-12));
    #if GRIDDIM >= 2
        REQUIRE(posN[1] == Approx(hi[1] - 0.5 * dx[1]).epsilon(1e-12));
    #endif
    #if GRIDDIM >= 3
        REQUIRE(posN[2] == Approx(hi[2] - 0.5 * dx[2]).epsilon(1e-12));
    #endif
}

TEST_CASE("Geometry: getIdx is left-inverse of getPos", "[geometry]")
{
    std::array<REAL, GRIDDIM> lo  = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
    std::array<REAL, GRIDDIM> hi  = {GRIDDIM_DECL(1.0, 1.0, 1.0)};
    std::array<int, GRIDDIM>  res = {GRIDDIM_DECL(8, 8, 8)};
    Geometry geom(lo, hi, res);

    for (int i = 0; i < res[0]; ++i)
    {
        std::array<int, GRIDDIM>  idx_in  = {GRIDDIM_DECL(i, i, i)};
        std::array<REAL, GRIDDIM> pos, idx_out;
        geom.getPos(pos, idx_in);
        geom.getIdx(idx_out, pos);
        REQUIRE(idx_out[0] == Approx((REAL)i).margin(1e-12));
        #if GRIDDIM >= 2
            REQUIRE(idx_out[1] == Approx((REAL)i).margin(1e-12));
        #endif
        #if GRIDDIM >= 3
            REQUIRE(idx_out[2] == Approx((REAL)i).margin(1e-12));
        #endif
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// DataArray: interior + ghost-cell indexing, getNGhost
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("DataArray: interior indexing and independent storage", "[dataarray]")
{
    constexpr int NV = 3;
    std::array<int, GRIDDIM> res = {GRIDDIM_DECL(4, 4, 4)};
    DataArray<NV> arr(res, 0);

    // Write distinguishable values and read them back
    arr({GRIDDIM_DECL(0, 0, 0)})[0] = 1.0;
    arr({GRIDDIM_DECL(1, 2, 3)})[1] = 2.0;
    arr({GRIDDIM_DECL(3, 3, 3)})[2] = 3.0;

    REQUIRE(arr({GRIDDIM_DECL(0, 0, 0)})[0] == 1.0);
    REQUIRE(arr({GRIDDIM_DECL(1, 2, 3)})[1] == 2.0);
    REQUIRE(arr({GRIDDIM_DECL(3, 3, 3)})[2] == 3.0);

    // Verify that writing one cell does not corrupt a neighbour
    REQUIRE(arr({GRIDDIM_DECL(0, 0, 0)})[1] == 0.0);
    REQUIRE(arr({GRIDDIM_DECL(0, 0, 0)})[2] == 0.0);
}

TEST_CASE("DataArray: ghost-cell access with negative / overflow indices", "[dataarray]")
{
    constexpr int NV = 2;
    std::array<int, GRIDDIM> res    = {GRIDDIM_DECL(4, 4, 4)};
    const int nGhost = 2;
    DataArray<NV> arr(res, nGhost);

    // Write to ghost cells at -1 and -2 in the x-direction
    arr({GRIDDIM_DECL(-1, 0, 0)})[0] = 99.0;
    arr({GRIDDIM_DECL(-2, 0, 0)})[0] = 88.0;
    // Write to ghost cells at res[0] and res[0]+1
    arr({GRIDDIM_DECL(4, 0, 0)})[1] = 77.0;
    arr({GRIDDIM_DECL(5, 0, 0)})[1] = 66.0;

    REQUIRE(arr({GRIDDIM_DECL(-1, 0, 0)})[0] == 99.0);
    REQUIRE(arr({GRIDDIM_DECL(-2, 0, 0)})[0] == 88.0);
    REQUIRE(arr({GRIDDIM_DECL(4,  0, 0)})[1] == 77.0);
    REQUIRE(arr({GRIDDIM_DECL(5,  0, 0)})[1] == 66.0);

    REQUIRE(arr.getNGhost() == nGhost);
}


// ─────────────────────────────────────────────────────────────────────────────
// Mesh / fillGhost: transmissive, periodic, reflective boundary conditions
// ─────────────────────────────────────────────────────────────────────────────
namespace
{
    // Set up a 4×4×4 mesh with distinguishable interior values and return it.
    Mesh<Euler::NVARS> makeMesh(int nGhost = 2)
    {
        std::array<REAL, GRIDDIM> lo  = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
        std::array<REAL, GRIDDIM> hi  = {GRIDDIM_DECL(4.0, 4.0, 4.0)};
        std::array<int, GRIDDIM>  res = {GRIDDIM_DECL(4, 4, 4)};
        Geometry geom(lo, hi, res);
        Mesh<Euler::NVARS> mesh(geom, nGhost);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    for (int v = 0; v < Euler::NVARS; ++v)
                        mesh({GRIDDIM_DECL(i, j, k)})[v] =
                            100.0 * (i + 1) + 10.0 * (j + 1) + (k + 1) + 0.1 * v;
        return mesh;
    }

    // Momentum index vector used by fillGhost to know which components to flip.
    static const std::vector<std::array<int, GRIDDIM>> g_vecIdx = {
        {GRIDDIM_DECL(Euler::MOM[0], Euler::MOM[1], Euler::MOM[2])}
    };
}

TEST_CASE("Mesh fillGhost: periodic BC wraps both ends in x", "[mesh][bc]")
{
    auto mesh = makeMesh();
    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for (int s = 0; s < 2; ++s)
        for (int d = 0; d < GRIDDIM; ++d)
            bc[s][d] = PERIODIC;

    mesh.fillGhost(bc, g_vecIdx);

    // Ghost at x = -1  should equal interior at x = 3 (last cell)
    for (int v = 0; v < Euler::NVARS; ++v)
        REQUIRE(mesh({GRIDDIM_DECL(-1, 0, 0)})[v] ==
                Approx(mesh({GRIDDIM_DECL(3, 0, 0)})[v]).margin(1e-15));

    // Ghost at x = 4  should equal interior at x = 0 (first cell)
    for (int v = 0; v < Euler::NVARS; ++v)
        REQUIRE(mesh({GRIDDIM_DECL(4, 0, 0)})[v] ==
                Approx(mesh({GRIDDIM_DECL(0, 0, 0)})[v]).margin(1e-15));

    // Second ghost layer (depth 2): -2 ↔ 2, 5 ↔ 1
    for (int v = 0; v < Euler::NVARS; ++v)
        REQUIRE(mesh({GRIDDIM_DECL(-2, 0, 0)})[v] ==
                Approx(mesh({GRIDDIM_DECL(2, 0, 0)})[v]).margin(1e-15));
    for (int v = 0; v < Euler::NVARS; ++v)
        REQUIRE(mesh({GRIDDIM_DECL(5, 0, 0)})[v] ==
                Approx(mesh({GRIDDIM_DECL(1, 0, 0)})[v]).margin(1e-15));
}

TEST_CASE("Mesh fillGhost: transmissive BC copies the boundary cell", "[mesh][bc]")
{
    auto mesh = makeMesh();
    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for (int s = 0; s < 2; ++s)
        for (int d = 0; d < GRIDDIM; ++d)
            bc[s][d] = TRANSMISSIVE;

    mesh.fillGhost(bc, g_vecIdx);

    // Ghost at x = -1 (lo side) → same as first interior cell x = 0
    for (int v = 0; v < Euler::NVARS; ++v)
        REQUIRE(mesh({GRIDDIM_DECL(-1, 0, 0)})[v] ==
                Approx(mesh({GRIDDIM_DECL(0, 0, 0)})[v]).margin(1e-15));

    // Second ghost layer (x = -2) also copies from x = 0
    for (int v = 0; v < Euler::NVARS; ++v)
        REQUIRE(mesh({GRIDDIM_DECL(-2, 0, 0)})[v] ==
                Approx(mesh({GRIDDIM_DECL(0, 0, 0)})[v]).margin(1e-15));

    // Ghost at x = 4 (hi side) → same as last interior cell x = 3
    for (int v = 0; v < Euler::NVARS; ++v)
        REQUIRE(mesh({GRIDDIM_DECL(4, 0, 0)})[v] ==
                Approx(mesh({GRIDDIM_DECL(3, 0, 0)})[v]).margin(1e-15));
}

TEST_CASE("Mesh fillGhost: reflective BC mirrors position and flips normal momentum",
          "[mesh][bc]")
{
    // Build a fresh mesh where every interior cell has a known x-momentum of 0.5.
    const IdealGas eos(1.4);
    std::array<REAL, GRIDDIM> lo  = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
    std::array<REAL, GRIDDIM> hi  = {GRIDDIM_DECL(4.0, 4.0, 4.0)};
    std::array<int, GRIDDIM>  res = {GRIDDIM_DECL(4, 4, 4)};
    Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 2);

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
            {
                auto U = makeState(eos, 1.0, 0.5, 0.2, 0.0, 1.0);
                for (int v = 0; v < Euler::NVARS; ++v)
                    mesh({GRIDDIM_DECL(i, j, k)})[v] = U[v];
            }

    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    for (int s = 0; s < 2; ++s)
        for (int d = 0; d < GRIDDIM; ++d)
            bc[s][d] = TRANSMISSIVE;
    bc[1][0] = REFLECTIVE;   // only the hi-x face is reflective

    mesh.fillGhost(bc, g_vecIdx);

    // Ghost at x=4 mirrors cell x=3: same density and energy, but x-momentum flipped.
    const REAL rho = mesh({GRIDDIM_DECL(3, 0, 0)})[Euler::RHO];
    const REAL mom_x_interior = mesh({GRIDDIM_DECL(3, 0, 0)})[Euler::MOM[0]];

    REQUIRE(mesh({GRIDDIM_DECL(4, 0, 0)})[Euler::RHO]    ==  Approx(rho).margin(1e-15));
    REQUIRE(mesh({GRIDDIM_DECL(4, 0, 0)})[Euler::MOM[0]] == -Approx(mom_x_interior).margin(1e-15));
    // Transverse momentum (y) is NOT flipped for a reflective x-face
    REQUIRE(mesh({GRIDDIM_DECL(4, 0, 0)})[Euler::MOM[1]] ==
            Approx(mesh({GRIDDIM_DECL(3, 0, 0)})[Euler::MOM[1]]).margin(1e-15));
}


// ─────────────────────────────────────────────────────────────────────────────
// Mesh::save / createFromFile  (binary round-trip)
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Mesh: binary save/load round-trip preserves all cell data", "[mesh][io]")
{
    // Build a 5×5×5 mesh with cells filled by a distinctive pattern.
    std::array<REAL, GRIDDIM> lo  = {GRIDDIM_DECL(0.0, 0.0, 0.0)};
    std::array<REAL, GRIDDIM> hi  = {GRIDDIM_DECL(1.0, 1.0, 1.0)};
    std::array<int, GRIDDIM>  res = {GRIDDIM_DECL(5, 5, 5)};
    Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> src(geom, 0);

    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            for (int k = 0; k < 5; ++k)
                for (int v = 0; v < Euler::NVARS; ++v)
                    src({GRIDDIM_DECL(i, j, k)})[v] =
                        static_cast<REAL>(i + 10*j + 100*k + 1000*v) + 0.5;

    const std::string hdr = "/tmp/ilmatar_unit_roundtrip.txt";
    const int   saved_step = 77;
    const REAL  saved_time = 2.718281828459045;
    src.save(hdr, saved_step, saved_time);

    int  loaded_step; REAL loaded_time;
    auto dst = Mesh<Euler::NVARS>::createFromFile(hdr, loaded_step, loaded_time, 0);

    REQUIRE(loaded_step == saved_step);
    REQUIRE(loaded_time == Approx(saved_time).epsilon(1e-12));

    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            for (int k = 0; k < 5; ++k)
                for (int v = 0; v < Euler::NVARS; ++v)
                {
                    INFO("i=" << i << " j=" << j << " k=" << k << " v=" << v);
                    REQUIRE(dst({GRIDDIM_DECL(i, j, k)})[v] ==
                            Approx(src({GRIDDIM_DECL(i, j, k)})[v]).epsilon(1e-15));
                }

    // Clean up temp files
    std::filesystem::remove(hdr);
    std::filesystem::remove(getDataFilename(hdr));
}

TEST_CASE("Mesh: header stores correct domain metadata", "[mesh][io]")
{
    std::array<REAL, GRIDDIM> lo  = {GRIDDIM_DECL(-1.0, -2.0, -3.0)};
    std::array<REAL, GRIDDIM> hi  = {GRIDDIM_DECL(1.0,   2.0,  3.0)};
    std::array<int, GRIDDIM>  res = {GRIDDIM_DECL(8, 16, 24)};
    Geometry geom(lo, hi, res);
    Mesh<Euler::NVARS> mesh(geom, 0);

    const std::string hdr = "/tmp/ilmatar_unit_metadata.txt";
    mesh.save(hdr, 5, 1.23456789012345);

    int step; REAL time;
    auto loaded = Mesh<Euler::NVARS>::createFromFile(hdr, step, time, 0);

    REQUIRE(step == 5);
    REQUIRE(time == Approx(1.23456789012345).epsilon(1e-12));

    const auto& lo_out  = loaded.getGeometry().getLo();
    const auto& hi_out  = loaded.getGeometry().getHi();
    const auto& res_out = loaded.getGeometry().getRes();
    for (int d = 0; d < GRIDDIM; ++d)
    {
        REQUIRE(lo_out[d]   == Approx(lo[d]).margin(1e-15));
        REQUIRE(hi_out[d]   == Approx(hi[d]).margin(1e-15));
        REQUIRE(res_out[d]  == res[d]);
    }

    std::filesystem::remove(hdr);
    std::filesystem::remove(getDataFilename(hdr));
}
