#include "b_spline_surface.h"
#include "KnotInsertion.h"
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <set>

using namespace Eigen;

BSplineSurface::BSplineSurface() : BSplineSurface({2,2}, {1,1}) {}

BSplineSurface::BSplineSurface(const std::array<int,2>& deg_in, const std::array<vector<double>,2>& kv_in) : err_max(0.)
{
    set_bspline_basis(deg_in, kv_in);
}

BSplineSurface::BSplineSurface(const std::array<int,2>& deg_in, const std::array<int,2>& nel) : err_max(0.)
{
    std::array<vector<double>,2> kv;
    std::array<double,2> dt { 1. / nel[0], 1. / nel[1] };
    for (int i=0; i<2; i++)
    {
        int j = 0;
        for (j=0; j<deg_in[i]; j++)
        {
            kv[i].push_back(0.);
        }
        for (j=0; j<nel[i]; j++)
        {
            kv[i].push_back(j * dt[i]);
        }
        for (j=0; j<=deg_in[i]; j++)
        {
            kv[i].push_back(1.);
        }
    }
    set_bspline_basis(deg_in, kv);
}

void BSplineSurface::set_bspline_basis(const std::array<int,2>& deg_in, const std::array<vector<double>,2>& kv_in)
{
    for (int i=0; i<2; i++)
    {
        bsp[i].Set(deg_in[i], kv_in[i]);
    }

    cp.clear();
    cp.resize(bsp[0].nbf * bsp[1].nbf);
}

void BSplineSurface::set_bspline_basis(int deg_in, const vector<double> &kv_in, int dir)
{
    if (dir != 0 && dir != 1)
    {
        cerr << "Invalid direction indicator!\n";
        return;
    }

    bsp[dir].Set(deg_in, kv_in);
    cp.clear();
    cp.resize(bsp[0].nbf * bsp[1].nbf);
}

void BSplineSurface::set_bspline_basis(int nel, int dir)
{
    if (dir != 0 && dir != 1)
    {
        cerr << "Invalid direction indicator!\n";
        return;
    }

    vector<double> kv;
    const double dt = 1. / nel;
    const int& deg = bsp[dir].p;
    int j = 0;
    for (j=0; j<deg; j++)
    {
        kv.push_back(0.);
    }
    for (j=0; j<nel; j++)
    {
        kv.push_back(j * dt);
    }
    for (j=0; j<=deg; j++)
    {
        kv.push_back(1.);
    }

    bsp[dir].Set(bsp[dir].p, kv);
    cp.clear();
    cp.resize(bsp[0].nbf * bsp[1].nbf);
}

void BSplineSurface::set_bspline_basis(const std::array<int,2>& nel)
{
    const std::array<int,2> deg = { bsp[0].p, bsp[1].p };
    std::array<vector<double>,2> kv;
    std::array<double,2> dt { 1. / nel[0], 1. / nel[1] };
    for (int i=0; i<2; i++)
    {
        int j = 0;
        for (j=0; j<deg[i]; j++)
        {
            kv[i].push_back(0.);
        }
        for (j=0; j<nel[i]; j++)
        {
            kv[i].push_back(j * dt[i]);
        }
        for (j=0; j<=deg[i]; j++)
        {
            kv[i].push_back(1.);
        }
    }
    set_bspline_basis(deg, kv);
}

void BSplineSurface::set_degree(int deg_in)
{
    bsp[0].p = deg_in;
    bsp[1].p = deg_in;
    set_bspline_basis({1,1});
}

std::array<int,2> BSplineSurface::get_cp_num() const
{
    return {bsp[0].nbf, bsp[1].nbf};
}

std::array<int,2> BSplineSurface::get_degree() const
{
    return {bsp[0].p, bsp[1].p};
}

const vector<double>& BSplineSurface::get_knot_vector(int dir) const
{
    return bsp[dir].kv;
}

const vector<std::array<double,3>>& BSplineSurface::get_control_points() const
{
    return cp;
}

void BSplineSurface::insert_C0_knots(const std::array<set<double>, 2> &knots)
{
    const double eps = 1.e-10;
    std::array<vector<double>, 2> kv;
    for (int i=0; i<2; i++)
    {
        const auto& k0 = bsp[i].kv;
        const auto& p = bsp[i].p;
        const auto& od = bsp[i].order;
        const auto nsp = k0.size() - 2 * od + 1;
        for (auto j=0; j<od; j++)
        {
            kv[i].push_back(k0.front());
        }
        for (auto j=p; j<nsp; j++)
        {
            bool is_repeat = false;
            for (auto kt : knots[i])
            {
                if (kt > k0[j] + eps && kt < k0[j+1] - eps)
                {
                    for (int k=0; k<p; k++)
                        kv[i].push_back(kt);
                }
                else if (kt >= k0[j+1] - eps && kt <= k0[j+1] + eps)
                {
                    is_repeat = true;
                    break;
                }
                else if (kt > k0[j+1] + eps)
                {
                    break;
                }
            }
            if (j < nsp - 1)
            {
                if (is_repeat)
                {
                    for (int k=0; k<p; k++)
                        kv[i].push_back(k0[j+1]);
                }
                else
                {
                    kv[i].push_back(k0[j+1]);
                }
            }
        }
        for (auto j=0; j<od; j++)
        {
            kv[i].push_back(k0.back());
        }
    }

    set_bspline_basis({bsp[0].p, bsp[1].p}, kv);
}

void BSplineSurface::find_IEN(const std::array<double, 2> &u, vector<int> &IEN) const
{
    IEN.clear();
    IEN.resize(bsp[0].nbf * bsp[1].nbf);

    std::array<int,2> ist { -1, -1 };
    for(int i=0; i<2; i++)
    {
        ist[i] = bsp[i].FindSpan(u[i]);
    }

    int pid = 0;
    for(int i=0; i<bsp[1].order; i++)
    {
        for(int j=0; j<bsp[0].order; j++)
        {
            IEN[pid++] = bsp[0].nbf * (ist[1] - bsp[1].p + i) + ist[0] - bsp[0].p + j;
        }
    }
}

void BSplineSurface::evaluate_basis(const std::array<double, 2> &u, vector<double> &Nx) const
{
    Nx.clear();
    Nx.resize(bsp[0].nbf * bsp[1].nbf);

    std::array<int,2> ist { -1, -1 };
    std::array<vector<double>,2> Nt;
    for(int i=0; i<2; i++)
    {
        ist[i] = bsp[i].FindSpan(u[i]);
        bsp[i].BasisFunctionVal(ist[i], u[i], Nt[i]);
    }

    int pid = 0;
    for(int i=0; i<bsp[1].order; i++)
    {
        for(int j=0; j<bsp[0].order; j++)
        {
            Nx[pid++] = Nt[0][j] * Nt[1][i];
        }
    }
}

bool BSplineSurface::fit(const vector<Point> &pts)
{
    const double eps = 1.e-10; // threshold to check uv coordinates
    const auto& ncpu = bsp[0].nbf;
    const auto& ncpv = bsp[1].nbf;
    const auto ncp = ncpu * ncpv; // == cp.size()
    const auto& ku = bsp[0].kv;
    const auto& kv = bsp[1].kv;
    const auto& pu = bsp[0].p;
    const auto& pv = bsp[1].p;

    vector<int> IDBC(ncp, 0);
    vector<std::array<double,3>> gh(ncp,{0.,0.,0.});

    const std::array<int, 4> cnid { 0, ncpu-1, (ncpv-1)*ncpu, ncpu*ncpv-1 };
    const std::array<std::array<double,2>,4> cncoor { {{0.,0.},{1.,0.},{0.,1.},{1.,1.}} };

    auto same_uv = [&eps](const std::array<double,2>& a, const std::array<double,2>& b)
    {
        return fabs(a[0]-b[0]) < eps && fabs(a[1]-b[1]) < eps;
    };

    for (auto i=0; i<4; i++)
    {
        bool is_found = false;
        for (auto& pt : pts)
        {
            if (pt.is_bnd && same_uv(cncoor[i], pt.u))
            {
                is_found = true;
                gh[cnid[i]] = pt.x;
                break;
            }
        }
        if (is_found)
            IDBC[cnid[i]] = -1;
        else
            cerr << "Can't find right corners when enforcing corner constraints!\n";
    }

//    for (auto i=pu+1; i<ku.size()-pu-1; i++) // boundary v == 0 and v == 1
//    {
//        bool is_dup = true;
//        for (auto j=0; j<pu-1; j++)
//        {
//            if (fabs(ku[i+j+1] - ku[i]) > eps)
//            {
//                is_dup = false;
//                break;
//            }
//        }
//        if (is_dup)
//        {
//            auto icp0 = i-1;
//            auto icp1 = (ncpv-1)*ncpu+i-1;
//            for (auto& pt : pts)
//            {
//                if (pt.is_bnd && same_uv({ku[i], 0.}, pt.u))
//                {
//                    IDBC[icp0] = -1;
//                    gh[icp0] = pt.x;
//                    break;
//                }
//            }
//            for (auto& pt : pts)
//            {
//                if (pt.is_bnd && same_uv({ku[i], 1.}, pt.u))
//                {
//                    IDBC[icp1] = -1;
//                    gh[icp1] = pt.x;
//                    break;
//                }
//            }
//        }
//    }
//
//    for (auto i=pv+1; i<kv.size()-pv-1; i++) // boundary u == 0 and u == 1
//    {
//        bool is_dup = true;
//        for (auto j=0; j<pv-1; j++)
//        {
//            if (fabs(kv[i+j+1] - kv[i]) > eps)
//            {
//                is_dup = false;
//                break;
//            }
//        }
//        if (is_dup)
//        {
//            auto icp0 = (i-1)*ncpu;
//            auto icp1 = i*ncpu-1;
//            for (auto& pt : pts)
//            {
//                if (pt.is_bnd && same_uv({0., kv[i]}, pt.u))
//                {
//                    IDBC[icp0] = -1;
//                    gh[icp0] = pt.x;
//                    break;
//                }
//            }
//            for (auto& pt : pts)
//            {
//                if (pt.is_bnd && same_uv({1., kv[i]}, pt.u))
//                {
//                    IDBC[icp1] = -1;
//                    gh[icp1] = pt.x;
//                    break;
//                }
//            }
//        }
//    }


    int neq = 0;
    for (auto& id : IDBC)
    {
        if (id != -1)
            id = neq++;
    }

    vector<Triplet<double>> trilist;
    SparseMatrix<double> mat(neq, neq);
    std::array<VectorXd, 3> rhs;
    VectorXd sol(neq);
    mat.setZero();
    sol.setZero();
    for (int i=0; i<3; i++)
    {
        rhs[i].resize(neq);
        rhs[i].setZero();
    }

    // initialization
    const int nloc = bsp[0].order * bsp[1].order;
    vector<int> IEN (nloc, -1);
    vector<double> Nx (nloc, 0.);
    for (const auto& pt : pts)
    {
        find_IEN(pt.u, IEN);
        for (int i=0; i< nloc; i++)
        {
            for (int j=0; j<nloc; j++)
            {
                if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
                {
                    trilist.emplace_back(IDBC[IEN[i]], IDBC[IEN[j]], 0.);
                }
            }
        }
    }
    mat.setFromTriplets(trilist.begin(), trilist.end());
    mat.makeCompressed();

    // assembly
    for (const auto& pt : pts)
    {
        find_IEN(pt.u, IEN);
        evaluate_basis(pt.u, Nx);
        double mij = 0.;
        for (int i=0; i< nloc; i++)
        {
            if (IDBC[IEN[i]] != -1)
            {
                for (int j=0; j<nloc; j++)
                {
                    mij = Nx[i] * Nx[j];
                    if (IDBC[IEN[j]] != -1)
                    {
                        mat.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += mij;
                    }
                    else
                    {
                        for (int k=0; k<3; k++)
                        {
                            rhs[k](IDBC[IEN[i]]) -= mij * gh[IEN[j]][k];
                        }
                    }
                }
                for (int k=0; k<3; k++)
                {
                    rhs[k](IDBC[IEN[i]]) += Nx[i] * pt.x[k];
                }
            }
        }
    }

    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(mat);
    if (solver.info() != Success)
    {
        cerr << "The matrix in fitting can't be solved!\nIncrease the density of the input triangle mesh.\n";
        return false;
    }

    for (int i=0; i<3; i++)
    {
        sol = solver.solve(rhs[i]);
        for (auto j=0; j<cp.size(); j++)
        {
            if (IDBC[j] != -1)
            {
                cp[j][i] = sol(IDBC[j]);
            }
            else
            {
                cp[j][i] = gh[j][i];
            }
        }
    }

    compute_error(pts);

    return true;
}

bool BSplineSurface::fit_grid(const vector<Point> &pts, const int nu, const int nv)
{
    for (int i=0; i<4; i++)
    {
        auto is_success = fit_grid_boundary(pts, nu, nv, i);
        if (!is_success)
        {
            cerr << "Fail to fit the patch boundary!\nIncrease the quad mesh resolution and retry!\n";
            return false;
        }
    }

    auto is_success = fit_grid_interior(pts, nu, nv);
    if (!is_success)
    {
        cerr << "Fail to the patch interior!\nIncrease the quad mesh resolution and retry!\n";
        return false;
    }

    compute_error(pts);

    return true;
}

bool BSplineSurface::fit_grid_boundary(const vector<Point> &pts, const int nu, const int nv, const int bnd_id)
{
    const int dir = bnd_id % 2 == 0 ? 0 : 1;
    const int nsp = dir == 0 ? nu : nv;
    const int& ncp = bsp[dir].nbf;
    const int& ncpu = bsp[0].nbf;
    const int& ncpv = bsp[1].nbf;

    int col_id = -1;
    int row_id = -1;
    if (bnd_id == 0)
        row_id = 0;
    else if (bnd_id == 2)
        row_id = nv-1;

    if (bnd_id == 1)
        col_id = nu-1;
    else if (bnd_id == 3)
        col_id = 0;

    vector<std::array<double,3>> x(ncp);
    x[0] = dir == 0 ? pts[row_id*nu].x : pts[col_id].x;
    x.back() = dir == 0 ? pts[row_id*nu+nu-1].x : pts[col_id+(nv-1)*nu].x;

    vector<Point> crv_pts;
    crv_pts.reserve(nsp);
    if (dir == 0)
    {
        for (int i=0; i<nu; i++)
        {
            crv_pts.push_back(pts[row_id*nu+i]);
        }
    }
    else
    {
        for (int i=0; i<nv; i++)
        {
            crv_pts.push_back(pts[i*nu+col_id]);
        }
    }

    // check if there is degeneration
    bool is_degen = false;
    if (nsp == ncp)
    {
        const double eps = 1.e-8;
        for (auto i=0; i<crv_pts.size()-1; i++)
        {
            if (fabs(crv_pts[i+1].x[0]-crv_pts[i].x[0]) < eps &&
                    fabs(crv_pts[i+1].x[1]-crv_pts[i].x[1]) < eps &&
                fabs(crv_pts[i+1].x[2]-crv_pts[i].x[2]) < eps)
            {
                is_degen = true;
                break;
            }
        }
    }

    if (nsp < ncp || is_degen) // too few of sample points, linear approximation
    {
        vector<double> grev;
        bsp[dir].get_greville_points(grev);
        for (int i=1; i<ncp-1; i++)
        {
            for (int j=0; j<3; j++)
                x[i][j] = (1-grev[i])*x[0][j] + grev[i]*x.back()[j];
        }

        if (bnd_id == 0)
        {
            for (int i=0; i<ncp; i++)
            {
                cp[i] = x[i];
            }
        }
        else if (bnd_id == 1)
        {
            for (int i=0; i<ncp; i++)
            {
                cp[i*ncpu+ncpu-1] = x[i];
            }
        }
        else if (bnd_id == 2)
        {
            for (int i=0; i<ncp; i++)
            {
                cp[(ncpv-1)*ncpu+i] = x[i];
            }
        }
        else // if (bnd_id == 3)
        {
            for (int i=0; i<ncp; i++)
            {
                cp[i*ncpu] = x[i];
            }
        }

        return true;
    }

    vector<int> IDBC(ncp, 0);
    vector<std::array<double,3>> gh(ncp,{0.,0.,0.});
    IDBC[0] = -1;
    int count = 0;
    for (int i=1; i<ncp-1; i++)
        IDBC[i] = count++;
    IDBC.back() = -1;
    gh[0] = crv_pts[0].x;
    gh.back() = crv_pts.back().x;

    vector<Triplet<double>> trilist;
    SparseMatrix<double> mat(count, count);
    std::array<VectorXd, 3> rhs;
    VectorXd sol(count);
    mat.setZero();
    sol.setZero();
    for (int i=0; i<3; i++)
    {
        rhs[i].resize(count);
        rhs[i].setZero();
    }

    // initialization
    const int nloc = bsp[dir].order;
    vector<int> IEN (nloc, -1);
    vector<double> Nx (nloc, 0.);
    int ist = -1;
    for (const auto& pt : crv_pts)
    {
        ist = bsp[dir].FindSpan(pt.u[dir]);
        for (int i=0; i<nloc; i++)
            IEN[i] = ist - nloc + 1 + i;
        for (int i=0; i< nloc; i++)
        {
            for (int j=0; j<nloc; j++)
            {
                if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
                {
                    trilist.emplace_back(IDBC[IEN[i]], IDBC[IEN[j]], 0.);
                }
            }
        }
    }
    mat.setFromTriplets(trilist.begin(), trilist.end());
    mat.makeCompressed();

    // assembly
    for (const auto& pt : crv_pts)
    {
        ist = bsp[dir].FindSpan(pt.u[dir]);
        for (int i=0; i<nloc; i++)
            IEN[i] = ist - nloc + 1 + i;
        bsp[dir].BasisFunctionVal(ist, pt.u[dir], Nx);
        double mij = 0.;
        for (int i=0; i< nloc; i++)
        {
            if (IDBC[IEN[i]] != -1)
            {
                for (int j=0; j<nloc; j++)
                {
                    mij = Nx[i] * Nx[j];
                    if (IDBC[IEN[j]] != -1)
                    {
                        mat.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += mij;
                    }
                    else
                    {
                        for (int k=0; k<3; k++)
                        {
                            rhs[k](IDBC[IEN[i]]) -= mij * gh[IEN[j]][k];
                        }
                    }
                }
                for (int k=0; k<3; k++)
                {
                    rhs[k](IDBC[IEN[i]]) += Nx[i] * pt.x[k];
                }
            }
        }
    }

    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(mat);
    if (solver.info() != Success)
    {
        cerr << "The matrix in fitting can't be solved!\nIncrease the density of the input triangle mesh.\n";
        return false;
    }

    for (int i=0; i<3; i++)
    {
        sol = solver.solve(rhs[i]);
        for (int j=0; j<ncp; j++)
        {
            if (IDBC[j] != -1)
                x[j][i] = sol(IDBC[j]);
        }
    }

    if (bnd_id == 0)
    {
        for (int i=0; i<ncp; i++)
        {
            cp[i] = x[i];
        }
    }
    else if (bnd_id == 1)
    {
        for (int i=0; i<ncp; i++)
        {
            cp[i*ncpu+ncpu-1] = x[i];
        }
    }
    else if (bnd_id == 2)
    {
        for (int i=0; i<ncp; i++)
        {
            cp[(ncpv-1)*ncpu+i] = x[i];
        }
    }
    else // if (bnd_id == 3)
    {
        for (int i=0; i<ncp; i++)
        {
            cp[i*ncpu] = x[i];
        }
    }

    return true;
}

bool BSplineSurface::fit_grid_interior(const vector<Point> &pts, const int nu, const int nv)
{
    const double eps = 1.e-10; // threshold to check uv coordinates
    const auto& ncpu = bsp[0].nbf;
    const auto& ncpv = bsp[1].nbf;
    const auto ncp = ncpu * ncpv; // == cp.size()
    const auto& ku = bsp[0].kv;
    const auto& kv = bsp[1].kv;
    const auto& pu = bsp[0].p;
    const auto& pv = bsp[1].p;

    if (nu < ncpu)
    {
        vector<double> grev;
        bsp[0].get_greville_points(grev);
        for (int r=1; r<ncpv-1; r++)
        {
            int st = r * ncpu;
            for (int c=1; c<ncpu-1; c++)
            {
                for (int i=0; i<3; i++)
                    cp[st+c][i] = (1 - grev[c]) * cp[st][i] + grev[c] * cp[st+ncpu-1][i];
            }
        }
        return true;
    }

    if (nv < ncpv)
    {
        vector<double> grev;
        bsp[1].get_greville_points(grev);
        for (int c=1; c<ncpu-1; c++)
        {
            for (int r=1; r<ncpv-1; r++)
            {
                for (int i=0; i<3; i++)
                    cp[r*ncpu+c][i] = (1 - grev[r]) * cp[c][i] + grev[r] * cp[(ncpv-1)*ncpu+c][i];
            }
        }
        return true;
    }

    vector<int> IDBC(ncp, -1);
    int neq = 0;
    for (int i=1; i<ncpv-1; i++)
    {
        for (int j=1; j<ncpu-1; j++)
        {
            IDBC[i*ncpu+j] = neq++;
        }
    }

    vector<Triplet<double>> trilist;
    SparseMatrix<double> mat(neq, neq);
    std::array<VectorXd, 3> rhs;
    VectorXd sol(neq);
    mat.setZero();
    sol.setZero();
    for (int i=0; i<3; i++)
    {
        rhs[i].resize(neq);
        rhs[i].setZero();
    }

    // initialization
    const int nloc = bsp[0].order * bsp[1].order;
    vector<int> IEN (nloc, -1);
    vector<double> Nx (nloc, 0.);
    for (const auto& pt : pts)
    {
        find_IEN(pt.u, IEN);
        for (int i=0; i< nloc; i++)
        {
            for (int j=0; j<nloc; j++)
            {
                if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
                {
                    trilist.emplace_back(IDBC[IEN[i]], IDBC[IEN[j]], 0.);
                }
            }
        }
    }
    mat.setFromTriplets(trilist.begin(), trilist.end());
    mat.makeCompressed();

    // assembly
    for (const auto& pt : pts)
    {
        find_IEN(pt.u, IEN);
        evaluate_basis(pt.u, Nx);
        double mij = 0.;
        for (int i=0; i< nloc; i++)
        {
            if (IDBC[IEN[i]] != -1)
            {
                for (int j=0; j<nloc; j++)
                {
                    mij = Nx[i] * Nx[j];
                    if (IDBC[IEN[j]] != -1)
                    {
                        mat.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += mij;
                    }
                    else
                    {
                        for (int k=0; k<3; k++)
                        {
                            rhs[k](IDBC[IEN[i]]) -= mij * cp[IEN[j]][k];
                        }
                    }
                }
                for (int k=0; k<3; k++)
                {
                    rhs[k](IDBC[IEN[i]]) += Nx[i] * pt.x[k];
                }
            }
        }
    }

    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(mat);
    if (solver.info() != Success)
    {
        cerr << "The matrix in fitting can't be solved!\nIncrease the density of the input triangle mesh.\n";
        return false;
    }

    for (int i=0; i<3; i++)
    {
        sol = solver.solve(rhs[i]);
        for (auto j=0; j<cp.size(); j++)
        {
            if (IDBC[j] != -1)
            {
                cp[j][i] = sol(IDBC[j]);
            }
        }
    }

    return true;
}

void BSplineSurface::compute_error(const vector<Point> &pts)
{
    err.clear();
    err.resize(pts.size(), 0.);
    err_max = 0.;

    const int nloc = bsp[0].order * bsp[1].order;
    vector<int> IEN (nloc, -1);
    vector<double> Nx (nloc, 0.);
    Vector3d xh, x, diff;
    int loc = 0;
    for (auto& pt : pts)
    {
        find_IEN(pt.u, IEN);
        evaluate_basis(pt.u, Nx);
        xh.setZero();
        x << pt.x[0], pt.x[1], pt.x[2];
        for (int i=0; i< nloc; i++)
        {
            for (int k=0; k<3; k++)
            {
                xh(k) += Nx[i] * cp[IEN[i]][k];
            }
        }
        diff = xh - x;
        err[loc] = diff.norm();
        if (err[loc] > err_max)
            err_max = err[loc];
        loc++;
    }
}

double BSplineSurface::get_max_error() const
{
    return err_max;
}

void BSplineSurface::identify(const vector<Point> &pts, std::array<vector<double>,2> &knots) const
{
    knots[0].clear();
    knots[1].clear();
    const double tol = 0.7 * err_max;
    std::array<set<int>, 2> ks;
    int i = 0;
    for (auto& pt : pts)
    {
        if (err[i] > tol)
        {
            auto iu = bsp[0].FindSpan(pt.u[0]);
            auto iv = bsp[1].FindSpan(pt.u[1]);
            ks[0].insert(iu);
            ks[1].insert(iv);
        }
        i++;
    }
    double kt = 0.;
    for (i=0; i<2; i++)
    {
        for (auto is : ks[i])
        {
            kt = 0.5 * (bsp[i].kv[is] + bsp[i].kv[is+1]);
            knots[i].push_back(kt);
        }
        sort(knots[i].begin(), knots[i].end());
    }
}

void BSplineSurface::refine(const vector<double>& knots, int dir)
{
    const vector<double> kv0 (bsp[dir].kv);
    MatrixXd cp0(cp.size(),3);
    for(auto i=0; i<cp.size(); i++)
    {
        for(int j=0; j<3; j++)
            cp0(i,j) = cp[i][j];
    }
    const int nbf0 = bsp[dir].nbf;
    const int deg = bsp[dir].p;

    vector<double> kv1;
    int i1 = 0;
    for (auto i=0; i<kv0.size()-1; i++)
    {
        kv1.push_back(kv0[i]);
        while (i1 < knots.size() && knots[i1] >= kv0[i] && knots[i1] < kv0[i+1])
        {
            kv1.push_back(knots[i1]);
            i1++;
        }
    }
    kv1.push_back(kv0.back());
    const int nbf1 = kv1.size() - bsp[dir].order;

    set_bspline_basis(deg, kv1, dir);

    SparseMatrix<double> submat(nbf1, nbf0);
    TMatrix(kv0, kv1, deg, submat);

    int id;
    if(dir == 0)
    {
        MatrixXd cp1;
        for(int i=0; i<bsp[1].nbf; i++) // nbf[1]
        {
            cp1 = submat * cp0.block(i*nbf0,0,nbf0,3);
            for(int j=0; j<bsp[dir].nbf; j++) // nbf[dir]
            {
                id = i * bsp[0].nbf + j; // nbf[0]
                cp[id][0] = cp1(j,0);
                cp[id][1] = cp1(j,1);
                cp[id][2] = cp1(j,2);
            }
        }
    }
    else if(dir == 1)
    {
        MatrixXd cp0b(nbf0,3);
        MatrixXd cp1;
        for(int i=0; i<bsp[0].nbf; i++) // nbf[0]
        {
            for(int j=0; j<nbf0; j++)
            {
                id = j * bsp[0].nbf + i; // nbf[0]
                cp0b(j,0) = cp0(id,0);
                cp0b(j,1) = cp0(id,1);
                cp0b(j,2) = cp0(id,2);
            }
            cp1 = submat * cp0b;
            for(int j=0; j<bsp[dir].nbf; j++) // nbf[dir]
            {
                id = j * bsp[0].nbf + i; // nbf[0]
                cp[id][0] = cp1(j,0);
                cp[id][1] = cp1(j,1);
                cp[id][2] = cp1(j,2);
            }
        }
    }
}

void BSplineSurface::refine(const std::array<vector<double>, 2> &knots)
{
    refine(knots[0], 0);
    refine(knots[1], 1);
}

void BSplineSurface::refine()
{
    const std::array<vector<double>, 2> kv0 { bsp[0].kv, bsp[1].kv };
    std::array<vector<double>, 2> kv;
    for (int i=0; i<2; i++)
    {
        kv[i].reserve(2*kv0[i].size());
        for (auto j=0; j<kv0[i].size()-1; j++)
        {
            kv[i].push_back(kv0[i][j]);
            if (kv0[i][j] < kv0[i][j+1])
            {
                kv[i].push_back(0.5*(kv0[i][j]+kv0[i][j+1]));
            }
        }
        kv[i].push_back(kv0[i].back());
    }

    set_bspline_basis({bsp[0].p,bsp[1].p}, kv);
}

void BSplineSurface::refine(int dir)
{
    const auto kv0 = bsp[dir].kv;
    vector<double> kv;
    kv.reserve(2*kv0.size());
    for (auto i=0; i<kv0.size()-1; i++)
    {
        kv.push_back(kv0[i]);
        if (kv0[i] < kv0[i+1])
        {
            kv.push_back(0.5*(kv0[i]+kv0[i+1]));
        }
    }
    kv.push_back(kv0.back());

    set_bspline_basis(bsp[dir].p, kv, dir);
}

void BSplineSurface::geom_map(const std::array<double, 2> &u, std::array<double, 3> &pt) const
{
    vector<int> IEN;
    vector<double> Nx;
    find_IEN(u, IEN);
    evaluate_basis(u, Nx);

    pt = { 0., 0., 0. };
    for (auto i=0; i<IEN.size(); i++)
    {
        for (int j=0; j<3; j++)
        {
            pt[j] += cp[IEN[i]][j] * Nx[i];
        }
    }
}

void BSplineSurface::visualize_vtk(const std::string &fn) const
{
    //sampling points
    const int nsp = 7; // # interior sampling points per direction
    std::array<int, 2> nel { (int) bsp[0].ks.size(), (int) bsp[1].ks.size() };
    std::array<vector<double>,2> spu;
    spu[0].resize(nsp * nel[0] + nel[0] + 1);
    spu[1].resize(nsp * nel[1] + nel[1] + 1);
    double tmp = 0.;
    for(int i=0; i<2; i++)
    {
        int loc = 0;
        auto& ksi = i == 0 ? bsp[0].ks : bsp[1].ks;
        auto& kv = i == 0 ? bsp[0].kv : bsp[1].kv;
        for(int j=0; j<nel[i]; j++)
        {
            spu[i][loc] = kv[ksi[j]];
            loc++;
            tmp = (kv[ksi[j]+1] - kv[ksi[j]]) / (nsp+1);
            for(int k=0; k<nsp; k++)
            {
                spu[i][loc] = kv[ksi[j]] + (k + 1) * tmp;
                loc++;
            }
        }
        spu[i][loc] = kv.back();
    }
    vector<std::array<double,3>> spx(spu[0].size() * spu[1].size());
    int loc = 0;
    for(int i=0; i<spu[1].size(); i++)
    {
        for(int j=0; j<spu[0].size(); j++)
        {
            geom_map({spu[0][j], spu[1][i]}, spx[loc]);
            loc++;
        }
    }

    string fname = fn + "_bsp.vtk";
    ofstream fout;
    fout.open(fname);
    if(fout.is_open())
    {
        fout<<"# vtk DataFile Version 2.0\nSurface test\nASCII\n";
        fout<<"DATASET STRUCTURED_GRID\n";
        fout<<"DIMENSIONS " << spu[0].size() << " " << spu[1].size() <<" 1\n";
        fout<<"POINTS "<< spx.size() <<" float\n";
        for(auto & pt : spx)
        {
            fout << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
        }
        fout.close();
    }
    else
    {
        cout << "Cannot open " << fname << "!\n";
    }

    // physical element boundary
    int ned[2] = {(nel[1]+1)*nel[0]*(nsp+1), (nel[0]+1)*nel[1]*(nsp+1)};
    int nedpt[2] = {ned[0]+nel[1]+1, ned[1]+nel[0]+1};
    vector<std::array<int,2>> ed(ned[0]+ned[1]);
    vector<int> edpid(nedpt[0]+nedpt[1]);
    loc = 0;
    for(int i=0; i<nel[1]+1; i++)
    {
        for(int j=0; j<spu[0].size(); j++)
        {
            edpid[loc] = i * (nsp + 1) * spu[0].size() + j;
            loc++;
        }
    }
    for(int i=0; i<nel[0]+1; i++)
    {
        for(int j=0; j<spu[1].size(); j++)
        {
            edpid[loc] = j * spu[0].size() + i * (nsp + 1);
            loc++;
        }
    }
    loc = 0;
    for(int i=0; i<nel[1]+1; i++)
    {
        for(int j=0; j<spu[0].size()-1; j++)
        {
            ed[loc][0] = i*spu[0].size()+j;
            ed[loc][1] = i*spu[0].size()+j+1;
            loc++;
        }
    }
    for(int i=0; i<nel[0]+1; i++)
    {
        for(int j=0; j<spu[1].size()-1; j++)
        {
            ed[loc][0] = nedpt[0]+i*spu[1].size()+j;
            ed[loc][1] = nedpt[0]+i*spu[1].size()+j+1;
            loc++;
        }
    }
    string fname1 = fn + "_bsp_lines.vtk";
    fout.open(fname1);
    if(fout.is_open())
    {
        fout<<"# vtk DataFile Version 2.0\nLine ends\nASCII\n";
        fout<<"DATASET UNSTRUCTURED_GRID\n";
        fout<<"POINTS "<< edpid.size() <<" float\n";
        for(int i : edpid)
        {
            fout<<spx[i][0]<<" "<<spx[i][1]<<" "<<spx[i][2]<<"\n";
        }
        fout<<"\nCELLS " << ed.size() << " " << 3*ed.size() << "\n";
        for(auto & i : ed)
            fout<<"2 "<<i[0]<<" "<<i[1]<<"\n";
        fout<<"\nCELL_TYPES " << ed.size() <<"\n";
        for(int i=0; i<ed.size(); i++)
            fout << "3\n";
        fout.close();
    }
    else
    {
        cout<<"Cannot open "<<fname1<<"!\n";
    }
}

void BSplineSurface::visualize_control_mesh(const std::string &fn) const
{
    string fname = fn + "_cm.vtk";
    ofstream fout;
    fout.open(fname);
    if (!fout.is_open())
    {
        cerr << "Can't open " << fname << "!\n";
        return;
    }

    fout<<"# vtk DataFile Version 2.0\nTriangle mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    fout<<"POINTS "<<cp.size()<<" float\n";
    for(const auto& pt : cp)
    {
        fout << pt[0] << " " << pt[1] << " " << pt[2] << endl;
    }

    const auto& nu = bsp[0].nbf;
    const auto& nv = bsp[1].nbf;
    auto nfc = (nu-1) * (nv-1);
    fout << "\nCELLS " << nfc << " " << 5 * nfc << endl;
    for (int i=0; i<nv-1; i++)
    {
        for (int j=0; j<nu-1; j++)
        {
            fout << "4 " << i*nu+j << " " << i*nu+j+1 << " " << (i+1)*nu+j+1 << " " << (i+1)*nu+j << "\n";
        }
    }

    fout << "\nCELL_TYPES " << nfc << endl;
    for(int i=0; i<nfc; i++)
    {
        fout<<"9\n";
    }
}


