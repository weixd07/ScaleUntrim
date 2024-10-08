#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include "surface_parameterization.h"

SurfaceParameterization::SurfaceParameterization() : SurfaceParameterization(120.)
{}

SurfaceParameterization::SurfaceParameterization(double tol_angle_in) :
tol_angle(tol_angle_in), cnid({-1,-1,-1,-1}) {}

void SurfaceParameterization::read_tri_mesh_obj(const std::string &fname)
{
    ifstream fin;
    fin.open(fname);
    if(!fin.is_open())
    {
        cerr << "Can't open " << fname << "!\n";
        return;
    }

    clear();
    char type;
    Point pt;
    TriFace face;
    while (fin >> type)
    {
        if (type == 'v')
        {
            fin >> pt.x[0] >> pt.x[1] >> pt.x[2];
            pts.push_back(pt);
        }
        else if (type == 'f')
        {
            fin >> face.pid[0] >> face.pid[1] >> face.pid[2];
            for (int i=0; i<3; i++)
                face.pid[i]--;
            faces.push_back(face);
        }
    }

    fin.close();
}

void SurfaceParameterization::read_corner_IDs(const std::string& fname)
{
    if(fname.empty())
        return;

    ifstream fin;
    fin.open(fname);
    if(!fin.is_open())
    {
        cerr << "Can't open " << fname << "!\n";
        return;
    }

    int tmp = -1;
    int i = 0;
    while (fin >> tmp)
    {
        cnid[i++] = tmp;
        if (i == 4)
            break;
    }

    fin.close();
}

void SurfaceParameterization::retrieve_from_quad_mesh(const vector<std::array<double, 3>> &quad_pts,
                                                      const vector<vector<int>> &grid_ids)
{
    clear();
    const int nv = grid_ids.size();
    const int nu = grid_ids[0].size();
    auto npt = nu * nv;
    auto nfc = 2 * (nu-1) * (nv-1);
    pts.reserve(npt);
    faces.reserve(nfc);

    for (auto& row : grid_ids)
    {
        for (auto& i : row)
        {
            pts.emplace_back(quad_pts[i]);
        }
    }

    for (auto i=0; i<nv-1; i++)
    {
        for (auto j=0; j<nu-1; j++)
        {
            std::array<int,3> fc1 = { i*nu+j, i*nu+j+1, (i+1)*nu+j+1 };
            std::array<int,3> fc2 = { i*nu+j, (i+1)*nu+j+1, (i+1)*nu+j };
            faces.emplace_back(fc1);
            faces.emplace_back(fc2);
        }
    }

    std::array<int,4> cn = { 0, nu-1, nv*nu-1, (nv-1)*nu };
    set_corners(cn);
}

void SurfaceParameterization::set_corners(const std::array<int, 4>& cn_in)
{
    cnid = cn_in;
}

void SurfaceParameterization::build_edges()
{
    edges.clear();
    for(auto& fc : faces)
    {
        for(int i=0; i<3; i++)
        {
            Edge ed({fc.pid[i], fc.pid[(i+1)%3]});
            auto it = find(edges.begin(), edges.end(), ed);
            fc.edid[i] = it - edges.begin();
            if (it == edges.end())
                edges.push_back(ed);
        }
    }
}

void SurfaceParameterization::build_connectivity()
{
    p2e.clear();
    p2f.clear();
    e2f.clear();
    p2e.resize(pts.size());
    p2f.resize(pts.size());
    e2f.resize(edges.size());

    for(auto i=0; i<edges.size(); i++)
    {
        p2e[edges[i].pid[0]].push_back({i,0});
        p2e[edges[i].pid[1]].push_back({i,1});
    }

    for(auto i=0; i<faces.size(); i++)
    {
        p2f[faces[i].pid[0]].push_back({i,0});
        p2f[faces[i].pid[1]].push_back({i,1});
        p2f[faces[i].pid[2]].push_back({i,2});

        e2f[faces[i].edid[0]].push_back({i,0});
        e2f[faces[i].edid[1]].push_back({i,1});
        e2f[faces[i].edid[2]].push_back({i,2});
    }
}

void SurfaceParameterization::set_boundary_flags()
{
    for(auto i=0; i<e2f.size(); i++)
    {
        if(e2f[i].size() == 1)
        {
            edges[i].is_bnd = true;
            for(const auto & pid : edges[i].pid)
            {
                pts[pid].is_bnd = true;
                for(const auto& fcid : p2f[pid])
                {
                    faces[fcid[0]].is_bnd = true;
                }
            }
        }
    }
}

bool SurfaceParameterization::check_mesh_orientation()
{
    for (const auto & fc : e2f)
    {
        if (fc.size() == 2)
        {
            if (faces[fc[0][0]].pid[fc[0][1]] == faces[fc[1][0]].pid[fc[1][1]])
                return false;
        }
    }
    return true;
}

void SurfaceParameterization::find_boundary_points()
{
    pid_bnd.clear();
    int pid_cur = -1;
    for(auto i=0; i<pts.size(); i++)
    {
        if (pts[i].is_bnd)
        {
            pid_cur = i;
            break;
        }
    }
    if (pid_cur == -1)
    {
        cerr << "Can't find the initial boundary point!\n";
        abort();
    }

    const auto pid_end = pid_cur;
    while (true)
    {
        pid_bnd.push_back(pid_cur);
        int pid_next = -1;
        for (const auto& fcloc : p2f[pid_cur])
        {
            auto& fcid = fcloc[0];
            auto& loc = fcloc[1];
            if (!faces[fcid].is_bnd)
                continue;

            auto ed_loc = faces[fcid].edid[loc];
            if (!edges[ed_loc].is_bnd)
                continue;

            pid_next = edges[ed_loc].pid[0] == pid_cur ? edges[ed_loc].pid[1] : edges[ed_loc].pid[0];
            break;
        }
        if (pid_next == -1)
        {
            cerr << "Can't find the next boundary point!\n";
            abort();
        }

        if (pid_next == pid_end)
            break;

        pid_cur = pid_next;
    }
}

void SurfaceParameterization::find_sharp_points()
{
    sharp_loc.clear();
    sharp_sort.clear();
    vector<pair<int, double>> sh;
    double dot_product = 0.;
    for(auto i=0; i<pid_bnd.size(); i++)
    {
        if (is_sharp_point(i, dot_product))
        {
            sharp_loc.push_back(i);
            sh.emplace_back(i, dot_product);
        }
    }

    auto cmp = [](const pair<int,double>& p1, const pair<int,double>& p2)
    {
        return p1.second > p2.second;
    };

    sort(sh.begin(), sh.end(), cmp);
    sharp_sort.resize(sh.size());
    int i = 0;
    for (const auto& iloc : sh)
    {
        sharp_sort[i++] = iloc.first;
    }
}

void SurfaceParameterization::find_corners()
{
    if (cnid[0] != -1) // have been manually set
    {
        remove_corners_from_sharp_points();
        return;
    }
//    else
//    {
//        cerr << "Automatically finding corners works not well.\n";
//        cerr << "Try to specify the four corners manually.\n";
//        abort();
//    }

    vector<double> arc_len (pid_bnd.size()+1, 0.);
    for(auto i=0; i<pid_bnd.size(); i++)
    {
        const auto& x0 = pts[pid_bnd[i]].x;
        const auto& x1 = pts[pid_bnd[(i+1) % pid_bnd.size()]].x;
        arc_len[i+1] = arc_len[i] + get_distance(x0, x1);
    }

    const double bnd_len = arc_len.back();
    for (auto& len : arc_len)
    {
        len /= bnd_len;
    }

    if (sharp_loc.empty()) // subdivide into 4 equal pieces w.r.t. arc length
    {
        cnid[0] = pid_bnd[0];
        split(0, 0, arc_len);
        return;
    }

    vector<bool> is_concv (pid_bnd.size(), false);
    for (auto i : sharp_sort)
    {
        if (is_concave(i))
            is_concv[pid_bnd[i]] = true;
    }

    std::array<int,4> cnloc { -1, -1, -1, -1 }; // corner local ID in pid_bnd
    for (auto i : sharp_sort)
    {
        if (!is_concv[pid_bnd[i]])
        {
            cnloc[0] = i;
            break;
        }
    }
    if (cnloc[0] == -1)
        cnloc[0] = sharp_sort.front();

    cnid[0] = pid_bnd[cnloc[0]];
    if (cnloc[0] != 0)
    {
        const double offset = arc_len[cnloc[0]];
        for (auto& len : arc_len)
        {
            len -= offset;
            if (len < 0.)
                len += 1.;
        }
    }

    const double min_side_len = 0.125;
    const double max_side_len = 0.375;
    auto it = find(sharp_loc.begin(), sharp_loc.end(), cnloc[0]);
    auto ist = it - sharp_loc.begin();
    int k = 1;
    for (auto i=1; i<sharp_loc.size(); i++)
    {
        auto j = (i + ist) % sharp_loc.size();
        auto loc = sharp_loc[j];
        auto len = arc_len[loc] - arc_len[cnloc[k-1]];
        if (len < 0.)
            len += 1.;
        if (!is_concv[pid_bnd[loc]] && len > min_side_len && len < max_side_len)
        {
            cnloc[k] = loc;
            cnid[k++] = pid_bnd[loc];
            if (k == 4)
                break;
        }
    }

    auto it1 = find(cnid.begin(), cnid.end(), -1);
    if (it1 != cnid.end())
    {
        ist = it1 - cnid.begin(); // may be 1,2,3
        split(ist-1, cnloc[ist-1], arc_len);
    }

    remove_corners_from_sharp_points();
}

bool SurfaceParameterization::is_sharp_point(int bnd_loc, double& dot_product) const
{
    const double PI = 3.1415926535897932;
    const double tol = cos(tol_angle / 180 * PI);

    Vector3d vp, vn;
    get_directions(bnd_loc, vp, vn);
    dot_product = vp.dot(vn);

    if (dot_product > tol)
        return true;
    else
        return false;
}

bool SurfaceParameterization::is_concave(int bnd_loc) const
{
    const double eps = 1.e-8;
    Vector3d vp, x;
    get_directions(bnd_loc, vp, x);
    auto xc = vp.dot(x);
    if (fabs(xc-1.) < eps)
        return true;
    else if (fabs(xc+1.) < eps)
        return false;

    auto nm = x.cross(vp);
    nm.normalize();
    auto y = nm.cross(x);
    if (vp.dot(y) < 0.)
        return true;
    else
        return false;
}

void SurfaceParameterization::get_directions(int bnd_loc, Vector3d &v_prev, Vector3d &v_next) const
{
    v_prev.setZero();
    v_next.setZero();
    const int ip_loc = bnd_loc == 0 ? pid_bnd.size() - 1 : bnd_loc - 1;
    const int in_loc = (bnd_loc + 1) % pid_bnd.size();
    const int pid = pid_bnd[bnd_loc];
    const int ip = pid_bnd[ip_loc];
    const int in = pid_bnd[in_loc];
    v_prev << pts[ip].x[0]-pts[pid].x[0], pts[ip].x[1]-pts[pid].x[1], pts[ip].x[2]-pts[pid].x[2];
    v_next << pts[in].x[0]-pts[pid].x[0], pts[in].x[1]-pts[pid].x[1], pts[in].x[2]-pts[pid].x[2];
    v_prev.normalize();
    v_next.normalize();
}

void SurfaceParameterization::split(int ic, int ib, const vector<double>& arc_len)
{
    const int n_seg = 4 - ic;
    cnid[ic] = pid_bnd[ib];
    double ds = (1 - arc_len[ib]) / n_seg;
    for (int j=ic+1; j<4; j++)
    {
        auto t = arc_len[ib] + (j-ic) * ds;
        for (auto i=0; i<pid_bnd.size(); i++)
        {
            auto i0 = i == 0 ? pid_bnd.size()-1 : i-1;
            if (arc_len[i0] < t && arc_len[i] > t)
            {
                cnid[j] = pid_bnd[i];
                break;
            }
        }

    }

    set<int> cn_uniq (cnid.begin(), cnid.end());
    if (cn_uniq.size() != 4)
    {
        cerr << "Some corners are found to be the same!\nUse a denser input triangle mesh.\n";
        abort();
    }
    for (auto i : cnid)
    {
        if (i == -1)
        {
            cerr << "Can't find certain corners!\nUse a denser input triangle mesh.\n";
            abort();
        }
    }
}

void SurfaceParameterization::remove_corners_from_sharp_points()
{
    if (cnid[0] == -1)
        return;

    auto sharp_old = sharp_loc;
    sharp_loc.clear();
    for(auto i : sharp_old)
    {
        auto it = find(cnid.begin(), cnid.end(), pid_bnd[i]);
        if (it == cnid.end())
        {
            sharp_loc.push_back(i);
        }
    }
}

void SurfaceParameterization::set_C0()
{
    for (auto& pt : pts)
        pt.is_c0 = false;

    for (auto i : cnid)
        pts[i].is_c0 = true;

    for (auto i : sharp_loc)
        pts[pid_bnd[i]].is_c0 = true;
}

void SurfaceParameterization::reorder_boundary_points()
{
    if (cnid[0] == -1 || cnid[0] == pid_bnd[0])
        return;

    auto it = find(pid_bnd.begin(), pid_bnd.end(), cnid[0]);
    if (it == pid_bnd.end())
    {
        cerr << "Can't find the first corner point " << cnid[0] << " in the boundary point list!\n";
        abort();
    }
    auto ist = it - pid_bnd.begin();

    auto cn_old = cnid;
    auto pid_old = pid_bnd;
    pid_bnd.clear();
    pid_bnd.resize(pid_old.size());
    int count = 1;
    for(auto i=0; i<pid_old.size(); i++)
    {
        pid_bnd[i] = pid_old[(ist+i)%pid_old.size()];
        for (int j=1; j<4; j++)
        {
            if (pid_bnd[i] == cn_old[j])
            {
                cnid[count++] = cn_old[j];
                break;
            }
        }
    }

    for (auto& loc : sharp_loc)
    {
        loc -= ist;
        if (loc < 0)
            loc += pid_old.size();
    }

    for (auto& loc : sharp_sort)
    {
        loc -= ist;
        if (loc < 0)
            loc += pid_old.size();
    }
}

void SurfaceParameterization::parameterize_boundary()
{
    std::array<double, 4> bnd_len { 0., 0., 0., 0. };
    int bnd_id = 0;
    for (auto i=0; i<pid_bnd.size(); i++)
    {
        int pid0 = pid_bnd[i];
        int pid1 = pid_bnd[(i + 1) % pid_bnd.size()];
        const auto& x0 = pts[pid0].x;
        const auto& x1 = pts[pid1].x;
        bnd_len[bnd_id] += get_distance(x0, x1);

        if (bnd_id == 0)
            pts[pid1].u = { bnd_len[bnd_id], 0. };
        else if (bnd_id == 1)
            pts[pid1].u = { 1., bnd_len[bnd_id] };
        else if (bnd_id == 2)
            pts[pid1].u = { bnd_len[bnd_id], 1. };
        else if (bnd_id == 3)
            pts[pid1].u = { 0., bnd_len[bnd_id] };

        if (pid1 == cnid[(bnd_id+1)%4])
            bnd_id++;
    }

    bnd_id = 0;
    for (auto i=0; i<pid_bnd.size(); i++)
    {
        int pid1 = pid_bnd[(i + 1) % pid_bnd.size()];
        if (bnd_id == 0)
            pts[pid1].u[0] /= bnd_len[bnd_id];
        else if (bnd_id == 1)
            pts[pid1].u[1] /= bnd_len[bnd_id];
        else if (bnd_id == 2)
        {
            pts[pid1].u[0] /= bnd_len[bnd_id];
            pts[pid1].u[0] = 1 - pts[pid1].u[0];
        }
        else if (bnd_id == 3)
        {
            pts[pid1].u[1] /= bnd_len[bnd_id];
            pts[pid1].u[1] = 1 - pts[pid1].u[1];
        }

        if (pid1 == cnid[(bnd_id+1)%4])
            bnd_id++;
    }

    pts[cnid[0]].u = { 0., 0. };
    pts[cnid[1]].u = { 1., 0. };
    pts[cnid[2]].u = { 1., 1. };
    pts[cnid[3]].u = { 0., 1. };
}

void SurfaceParameterization::parameterize_interior()
{
    vector<int> IDBC(pts.size(), 0);
    int neq = 0;
    for (auto i=0; i<pts.size(); i++)
    {
        if (pts[i].is_bnd)
            IDBC[i] = -1;
        else
            IDBC[i] = neq++;
    }
    if (neq == 0) // all points are boundary points
        return;

    vector<Triplet<double>> trilist;
    SparseMatrix<double> mat(neq, neq);
    VectorXd b0(neq);
    VectorXd b1(neq);
    mat.setZero();
    b0.setZero();
    b1.setZero();
    for (auto i=0; i<pts.size(); i++)
    {
        if (pts[i].is_bnd)
            continue;

        trilist.emplace_back(IDBC[i], IDBC[i], 1.);

        vector<double> w (p2e.size(), 0.);
        double sum = 0.;
        int k = 0;
        for (const auto& edloc : p2e[i])
        {
            w[k] = compute_mean_value_coordinates(i, edloc[0]);
            sum += w[k];
            k++;
        }

        k = 0;
        for (const auto& edloc : p2e[i])
        {
            auto& edid = edloc[0];
            auto& loc = edloc[1];
            int j = edges[edid].pid[1-loc];
            double lambda = w[k++] / sum;
            if (pts[j].is_bnd)
            {
                b0(IDBC[i]) += lambda * pts[j].u[0];
                b1(IDBC[i]) += lambda * pts[j].u[1];
            }
            else
                trilist.emplace_back(IDBC[i], IDBC[j], -lambda);
        }
    }

    mat.setFromTriplets(trilist.begin(), trilist.end());
    mat.makeCompressed();

    SparseLU<SparseMatrix<double> > solver;
    solver.compute(mat);
    VectorXd u0 = solver.solve(b0);
    VectorXd u1 = solver.solve(b1);

    for (auto i=0; i<pts.size(); i++)
    {
        if (!pts[i].is_bnd)
        {
            pts[i].u[0] = u0[IDBC[i]];
            pts[i].u[1] = u1[IDBC[i]];
        }
    }
}

double SurfaceParameterization::compute_mean_value_coordinates(int pid, int edid)
{
    if (edges[edid].is_bnd)
    {
        cerr << "Computing lambda for boundary points!\n";
        abort();
    }

    auto len = get_distance(pts[edges[edid].pid[0]].x, pts[edges[edid].pid[1]].x);
    double w = 0.;
    for (const auto& fcloc : e2f[edid])
    {
        auto& fcid = fcloc[0];
        auto& loc = fcloc[1];
        if (faces[fcid].pid[loc] == pid) // alpha_i
        {
            auto i = faces[fcid].edid[(loc+1)%3];
            auto a = get_distance(pts[edges[i].pid[0]].x, pts[edges[i].pid[1]].x);
            i = faces[fcid].edid[(loc+2)%3];
            auto b = get_distance(pts[edges[i].pid[0]].x, pts[edges[i].pid[1]].x);
            auto cos_alpha = (b * b + len * len - a * a) / (2. * b * len);
            auto sin_alpha = sqrt(1. - cos_alpha * cos_alpha);
            w += (1. - cos_alpha) / sin_alpha;
        }
        else // alpha_i-1
        {
            auto i = faces[fcid].edid[(loc+2)%3];
            auto a = get_distance(pts[edges[i].pid[0]].x, pts[edges[i].pid[1]].x);
            i = faces[fcid].edid[(loc+1)%3];
            auto b = get_distance(pts[edges[i].pid[0]].x, pts[edges[i].pid[1]].x);
            auto cos_alpha = (b * b + len * len - a * a) / (2. * b * len);
            auto sin_alpha = sqrt(1. - cos_alpha * cos_alpha);
            w += (1. - cos_alpha) / sin_alpha;
        }
    }
    w /= len;

    return w;
}

double SurfaceParameterization::get_distance(const std::array<double, 3> &x0, const std::array<double, 3> &x1) const
{
    std::array<double, 3> diff { x0[0] - x1[0], x0[1] - x1[1], x0[2] - x1[2] };
    double len = sqrt( diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2] );
    return len;
}

double SurfaceParameterization::get_dimension_scale() const
{
    auto mx = numeric_limits<double>::max();
    auto mn = numeric_limits<double>::lowest();
    std::array<double, 3> xmin { mx, mx, mx };
    std::array<double, 3> xmax { mn, mn, mn };

    for (auto & pt : pts)
    {
        for (int i=0; i<3; i++)
        {
            if (xmin[i] > pt.x[i])
                xmin[i] = pt.x[i];
            if (xmax[i] < pt.x[i])
                xmax[i] = pt.x[i];
        }
    }

    return get_distance(xmin, xmax);
}

const vector<Point>& SurfaceParameterization::get_points() const
{
    return pts;
}

void SurfaceParameterization::get_sharp_point_param_coor(std::array<set<double>, 2> &ush) const
{
    ush[0].clear();
    ush[1].clear();
    const double eps = 1.e-10;

//    auto my_equal = [&eps](double a, double b)
//    {
//        auto diff = a - b;
//        return  diff > -eps && diff < eps;
//    };
//
//    auto my_insert = [my_equal](set<double>& vec, double item)
//    {
//        bool is_in = false;
//        for (auto val : vec)
//        {
//            if (my_equal(val, item))
//            {
//                is_in = true;
//                break;
//            }
//        }
//        if (!is_in)
//            vec.insert(item);
//    };

//    array<set<double>,2> utmp;
    for (auto i : sharp_loc)
    {
        const auto& u = pts[pid_bnd[i]].u;
        if (fabs(u[1]) < eps || fabs(u[1]-1) < eps)
        {
            ush[0].insert(u[0]);
//            my_insert(ush[0], u[0]);
        }
        if (fabs(u[0]) < eps || fabs(u[0]-1) < eps)
        {
            ush[1].insert(u[1]);
//            my_insert(ush[1], u[1]);
        }
    }
}

void SurfaceParameterization::visualize_tri_mesh(const std::string &fn) const
{
    string fname = fn + "_tri.vtk";
    ofstream fout;
    fout.open(fname);
    if (!fout.is_open())
    {
        cerr << "Can't open " << fname << "!\n";
        return;
    }

    fout<<"# vtk DataFile Version 2.0\nTriangle mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    fout<<"POINTS "<<pts.size()<<" float\n";
    for(const auto& pt : pts)
    {
        fout << pt.x[0] << " " << pt.x[1] << " " << pt.x[2] << endl;
    }

    fout << "\nCELLS " << faces.size() << " " << 4 * faces.size() << endl;
    for(const auto& fc : faces)
    {
        fout<< "3 " << fc.pid[0] << " " << fc.pid[1] << " " << fc.pid[2] << endl;
    }

    fout << "\nCELL_TYPES " << faces.size() << endl;
    for(const auto& fc : faces)
    {
        fout<<"5\n";
    }

//    fout << "POINT_DATA " << pts.size() << "\nVECTORS disp float\n";
//    for (const auto& pt : pts)
//    {
//    	fout << pt.u[0] << " " << pt.u[1] << " 0\n";
//    }

    vector<int> flag (pts.size(), 0);
    for (auto i : sharp_loc)
    {
        flag[pid_bnd[i]] = 1;
    }
    for (auto i : cnid)
    {
        flag[i] = 1;
    }
    fout << "POINT_DATA " << pts.size() << "\nSCALARS sharp float 1\nLOOKUP_TABLE default\n";
    for(auto i : flag)
    {
        fout << i << endl;
    }
//    for(auto& pt : pts)
//    {
//        fout << pt.is_c0 << endl;
//    }

    //fout<<"\nCELL_DATA "<<(ns-1)*(ns-1)*errL2.size()<<"\nSCALARS err float 1\nLOOKUP_TABLE default\n";
    //for(i=0;i<errL2.size();i++)
    //{
    //	for(int j=0; j<(ns-1)*(ns-1); j++)
    //	fout<<errL2[i]<<"\n";
    //}

    fout.close();
}

void SurfaceParameterization::visualize_param_mesh(const std::string &fn) const
{
    string fname = fn + "_param.vtk";
    ofstream fout;
    fout.open(fname);
    if (!fout.is_open())
    {
        cerr << "Can't open " << fname << "!\n";
        return;
    }

    fout<<"# vtk DataFile Version 2.0\nTriangle mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    fout<<"POINTS "<<pts.size()<<" float\n";
    for(const auto& pt : pts)
    {
        fout << pt.u[0] << " " << pt.u[1] << " 0\n";
    }

    fout << "\nCELLS " << faces.size() << " " << 4 * faces.size() << endl;
    for(const auto& fc : faces)
    {
        fout<< "3 " << fc.pid[0] << " " << fc.pid[1] << " " << fc.pid[2] << endl;
    }

    fout << "\nCELL_TYPES " << faces.size() << endl;
    for(const auto& fc : faces)
    {
        fout<<"5\n";
    }

    fout.close();
}

void SurfaceParameterization::clear()
{
    pts.clear();
    edges.clear();
    faces.clear();
    p2e.clear();
    p2f.clear();
    e2f.clear();
    pid_bnd.clear();
    sharp_loc.clear();
    sharp_sort.clear();

    for(auto& i : cnid)
        i = -1;
}

void SurfaceParameterization::run()
{
    build_edges();
    build_connectivity();
    set_boundary_flags();

    find_boundary_points();
    find_sharp_points();
    find_corners();
    set_C0();
    reorder_boundary_points();

    parameterize_boundary();
    parameterize_interior();
}