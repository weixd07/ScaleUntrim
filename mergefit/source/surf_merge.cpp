#include "surf_merge.h"
#include <iostream>
#include <fstream>

SurfMerge::SurfMerge() : SurfMerge(1.e-3, 150., 10)
{}

SurfMerge::SurfMerge(double fit_tol_in, double sharp_tol_in, int itr_max_in) :
fit_tol (fit_tol_in), sharp_tol (sharp_tol_in), itr_max (itr_max_in), surf({2,2}, {4,4})
{}

void SurfMerge::set_surf_info(const std::array<int, 2> &deg, const std::array<vector<double>, 2> &kv)
{
    surf.set_bspline_basis(deg, kv);
}

void SurfMerge::set_fit_tol(double tol)
{
    fit_tol = tol;
}

void SurfMerge::set_sharp_tol(double tol)
{
    sharp_tol = tol;
}

void SurfMerge::set_max_iteration(int itr_max_in)
{
    itr_max = itr_max_in;
}

void SurfMerge::read(const std::string& tri_mesh, const std::string corners)
{
    mesh.read_tri_mesh_obj(tri_mesh);
    mesh.read_corner_IDs(corners);
}

void SurfMerge::run(const string& fn)
{
    mesh.run();
    const auto& tri_pts = mesh.get_points();

//    array<set<double>,2> sharp_knots;
//    mesh.get_sharp_point_param_coor(sharp_knots);
//    surf.insert_C0_knots(sharp_knots);

    bool is_success = surf.fit(tri_pts);
    if (!is_success)
        return;

    int itr = 0;
    double err = surf.get_max_error();
    const auto len_scl = mesh.get_dimension_scale();

    surf.visualize_vtk(fn+"_" + to_string(itr));
    cout << "Iteration: " << itr++ << endl;
    cout << "Fitting error: " << err / len_scl << endl;

    std::array<vector<double>, 2> knots;
    while (err > fit_tol * len_scl && itr < itr_max)
    {
        surf.identify(tri_pts, knots);
        surf.refine(knots);
        is_success = surf.fit(tri_pts);
        if (!is_success)
            break;
        err = surf.get_max_error();

        surf.visualize_vtk(fn+"_" + to_string(itr));
        cout << "Iteration: " << itr++ << endl;
        cout << "Fitting error: " << err / len_scl << endl;
    }

    cout << "Surface fitting done!\n";
}

void SurfMerge::write(const std::string& fn)
{
    const string fname = fn + ".pat";
    ofstream fout;
    fout.open(fname);
    if (!fout.is_open())
    {
        cerr << "Can't open " << fname << "!\n";
        return;
    }

    fout << "#Patch 1\n";
    fout << "Patch 0\n";

    fout << "Degree\n";
    auto deg = surf.get_degree();
    fout << deg[0] << " " << deg[1] << '\n';

    fout << "KU\n";
    auto& ku = surf.get_knot_vector(0);
    int count = 0;
    for (auto& kt : ku)
    {
        fout << kt << " ";
        if (++count == 16)
        {
            fout << '\n';
            count = 0;
        }
    }
    if (count != 0)
        fout << "\n";

    fout << "KV\n";
    auto& kv = surf.get_knot_vector(1);
    count = 0;
    for (auto& kt : kv)
    {
        fout << kt << " ";
        if (++count == 16)
        {
            fout << '\n';
            count = 0;
        }
    }
    if (count != 0)
        fout << "\n";

    fout << "CP\n";
    auto& cps = surf.get_control_points();
    auto ncp = surf.get_cp_num();
    double w = 1.; // all weights are 1
    int pid = 0;
    for (int iu=0; iu<ncp[0]; iu++)
    {
        for (int iv=0; iv<ncp[1]; iv++)
        {
            pid = iv * ncp[0] + iu;
            fout << cps[pid][0] << " " << cps[pid][1] << " " <<cps[pid][2] << " " << w << '\n';
        }
    }

//        for (auto& pt : cps)
//            fout << pt[0] << " " << pt[1] << " " << pt[2] << '\n';

    fout.close();
}