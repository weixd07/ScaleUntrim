#include <iostream>
#include <fstream>
#include <algorithm>
#include "multi_patch_fit.h"

using namespace std;

MultiPatchFit::MultiPatchFit() : MultiPatchFit(1.e-3, 10, 2)
{}

MultiPatchFit::MultiPatchFit(double fit_tol_in, int itr_max_in, int deg_in) : fit_tol(fit_tol_in), itr_max(itr_max_in), degree(deg_in)
{}

void MultiPatchFit::set_fit_tol(double tol)
{
    fit_tol = tol;
}

void MultiPatchFit::set_max_iteration(int itr_max_in)
{
    itr_max = itr_max_in;
}

void MultiPatchFit::set_degree(int deg_in)
{
    degree = deg_in;
}

void MultiPatchFit::read(const std::string &fn)
{
    cout << "reading file " << fn << "...\n";
    ifstream fin;
    fin.open(fn);
    if(!fin.is_open())
    {
        cerr << "Can't open " << fn << "!\n";
        return;
    }

    string tmp;
    int npt = 0;
    int npatch = 0;
    int itmp = 0;
    fin >> tmp >> npt;
    spts.resize(npt);
    for (int i=0; i<npt; i++)
        fin >> spts[i][0] >> spts[i][1] >> spts[i][2];

    fin >> tmp >> npatch;

    patches.resize(npatch);
    for (int i=0; i<npatch; i++)
    {
        fin >> tmp >> itmp >> tmp >> patches[i].n_spt[1] >> tmp >> patches[i].n_spt[0];
        patches[i].spt_id.resize(patches[i].n_spt[1], vector<int>(patches[i].n_spt[0]));
        for (int j=0; j<patches[i].n_spt[1]; j++)
        {
            for (int k=0; k<patches[i].n_spt[0]; k++)
            {
                fin >> patches[i].spt_id[j][k];
            }
        }
    }

    fin.close();
}

void MultiPatchFit::build_connectivity()
{
    cout << "Building patch connectivity info...\n";
    if(spts.empty())
    {
        cerr << "No sample points have been read!\n";
        return;
    }

    edges.clear();
    for (auto& pat : patches)
    {
        pat.cn[0] = pat.spt_id[0][0];
        pat.cn[1] = pat.spt_id[0][pat.n_spt[0]-1];
        pat.cn[2] = pat.spt_id[pat.n_spt[1]-1][pat.n_spt[0]-1];
        pat.cn[3] = pat.spt_id[pat.n_spt[1]-1][0];
        for (int i=0; i<4; i++)
        {
            Edge ed({pat.cn[i], pat.cn[(i+1)%4]});
            auto it = find(edges.begin(), edges.end(), ed);
            pat.edid[i] = it - edges.begin();
            if (it == edges.end())
                edges.push_back(ed);
        }
    }

    e2p.clear();
    e2p.resize(edges.size());
    for (auto i=0; i<patches.size(); i++)
    {
        for (int j=0; j<4; j++)
        {
            e2p[patches[i].edid[j]].push_back({i,j});
        }
    }

    for (auto i=0; i<edges.size(); i++)
    {
        if (e2p[i].size() == 1)
            edges[i].is_bnd = true;
    }

    for (auto i=0; i<patches.size(); i++)
    {
        auto& pat = patches[i];
        for (int j=0; j<4; j++)
        {
            if (edges[pat.edid[j]].is_bnd)
                continue;

            if (e2p[pat.edid[j]][0][0] == i)
            {
                pat.ednb[j] =  e2p[pat.edid[j]][1][0];
                pat.edloc[j] = e2p[pat.edid[j]][1][1];
            }
            else
            {
                pat.ednb[j] =  e2p[pat.edid[j]][0][0];
                pat.edloc[j] = e2p[pat.edid[j]][0][1];
            }
        }
    }
}

void MultiPatchFit::initialize()
{
    cout << "Initializing...\n";
    if (patches.empty())
    {
        cerr << "No patches found!\n";
        abort();
    }

    meshes.resize(patches.size());
    surfs.resize(patches.size()); // by default: degree (2,2), element # (1,1)
    for (auto& surf : surfs)
    {
        surf.set_degree(degree);
    }

    const auto deg = surfs[0].get_degree();
    const std::array<int,2> nel2 = { 2, 2 };
    const std::array<int,2> nel4 = { 4, 4 };
    const std::array<int,2> ncp2 = { nel2[0]+2*deg[0]+1, nel2[1]+2*deg[1]+1 };
    const std::array<int,2> ncp4 = { nel4[0]+2*deg[0]+1, nel4[1]+2*deg[1]+1 };
    for (auto i=0; i<patches.size(); i++)
    {
        const auto& nsp = patches[i].n_spt;
        std::array<int,2> nel = { 1, 1 };
        for (int j=0; j<2; j++)
        {
            if (ncp4[j] < nsp[j])
                nel[j] = nel4[j];
            else if (ncp2[j] < nsp[j])
                nel[j] = nel2[j];
        }
        surfs[i].set_bspline_basis(nel);
    }
}

bool MultiPatchFit::refine()
{
    bool is_refine = false;
    const auto deg = surfs[0].get_degree();
    for (auto i=0; i<patches.size(); i++)
    {
        const auto& nsp = patches[i].n_spt;
        const auto ncp0 = surfs[i].get_cp_num();
        const std::array<int,2> ncp1 = { 2*ncp0[0]-deg[0], 2*ncp0[1]-deg[1] };
        if (ncp1[0] < nsp[0]-1 && ncp1[1] < nsp[1]-1)
        {
            surfs[i].refine();
            is_refine = true;
        }
        else if (ncp1[0] < nsp[0]-1)
        {
            surfs[i].refine(0);
            is_refine = true;
        }
        else if (ncp1[1] < nsp[1]-1)
        {
            surfs[i].refine(1);
            is_refine = true;
        }
    }
    return is_refine;
}

bool MultiPatchFit::fit()
{
    cout << "Perform surface fitting...\n";
    bool is_success_all = true;
    for (auto i=0; i<patches.size(); i++)
    {
        cerr << "patch " << i << endl;
        prepair_tri_mesh(i);
        auto& sample_pts = meshes[i].get_points();
        bool is_success = surfs[i].fit_grid(sample_pts, patches[i].n_spt[0], patches[i].n_spt[1]);
        if (!is_success)
        {
            is_success_all = false;
            break;
        }
    }
    return is_success_all;
}

void MultiPatchFit::prepair_tri_mesh(int patch_id)
{
    meshes[patch_id].retrieve_from_quad_mesh(spts, patches[patch_id].spt_id);
    meshes[patch_id].run();
}

double MultiPatchFit::get_max_error() const
{
    double error_max = 0.;
    for (const auto& surf : surfs)
    {
        double err = surf.get_max_error();
        error_max = max(error_max, err);
    }
    return error_max;
}

double MultiPatchFit::get_length_scale() const
{
    auto mx = numeric_limits<double>::max();
    auto mn = numeric_limits<double>::lowest();
    std::array<double, 3> xmin { mx, mx, mx };
    std::array<double, 3> xmax { mn, mn, mn };
    for (auto & pt : spts)
    {
        for (int i=0; i<3; i++)
        {
            if (xmin[i] > pt[i])
                xmin[i] = pt[i];
            if (xmax[i] < pt[i])
                xmax[i] = pt[i];
        }
    }
    std::array<double,3> diff = { xmax[0]-xmin[0], xmax[1]-xmin[1], xmax[2]-xmin[2] };
    double len = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);

    return len;
}

void MultiPatchFit::visualize_surface_vtk(const std::string &fn) const
{
    int i = 0;
    for (auto& surf : surfs)
    {
        surf.visualize_vtk(fn + "_" + to_string(i++));
    }
}

void MultiPatchFit::run(const string& fn)
{
    initialize();
    build_connectivity();

    bool is_success = fit();
    if (!is_success)
        return;

    int itr = 0;
    double err = get_max_error();
    const auto len_scl = get_length_scale();

//    visualize_surface_vtk(fn + "_" + to_string(itr));
    cout << "Iteration: " << itr++ << endl;
    cout << "Fitting error: " << err / len_scl << endl;

    while (err > fit_tol * len_scl && itr < itr_max)
    {
        bool is_refine = refine();
        if (!is_refine)
        {
            cout << "Cannot be refined further!\nNeed denser quad mesh to proceed!\n";
            break;
        }

        is_success = fit();
        if (!is_success)
            break;

        err = get_max_error();
//        visualize_surface_vtk(fn + "_" + to_string(itr));
        cout << "Iteration: " << itr++ << endl;
        cout << "Fitting error: " << err / len_scl << endl;
    }

    cout << "Surface fitting done!\n";
}

void MultiPatchFit::write(const std::string &fn) const
{
    const string fname = fn + ".pat";
    ofstream fout;
    fout.open(fname);
    if (!fout.is_open())
    {
        cerr << "Can't open " << fname << "!\n";
        return;
    }

    fout << "#Patch " << surfs.size() << '\n';
    int i = 0;
    double w = 1.; // all weights are 1
    for (auto& surf : surfs)
    {
        fout << "Patch " << i++ << '\n';

//        fout << "Degree\n";
        auto deg = surf.get_degree();
        fout << "Degree " <<  deg[0] << " " << deg[1] << '\n';

//        fout << "KU\n";
        auto& ku = surf.get_knot_vector(0);
        fout << "KU " << ku.size() << '\n';
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

//        fout << "KV\n";
        auto& kv = surf.get_knot_vector(1);
        fout << "KV " << kv.size() << '\n';
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

//        fout << "CP\n";
        auto& cps = surf.get_control_points();
        auto ncp = surf.get_cp_num();
        fout << "CP " << cps.size() << '\n';
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
    }

    fout.close();
}
