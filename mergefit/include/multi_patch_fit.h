#ifndef OPENSURF_MULTI_PATCH_FIT_H
#define OPENSURF_MULTI_PATCH_FIT_H

#include "surface_parameterization.h"
#include "b_spline_surface.h"

class MultiPatchFit
{
public:
    MultiPatchFit();
    MultiPatchFit(double fit_tol_in, int itr_max_in, int deg_in);

    void set_fit_tol(double tol);
    void set_max_iteration(int itr_max_in);
    void set_degree(int deg_in);

    void read(const string& fn);
    void build_connectivity();

    void initialize();
    bool refine();
    bool fit();
    void prepair_tri_mesh(int patch_id);
    double get_max_error() const;
    double get_length_scale() const;
    void visualize_surface_vtk(const string& fn) const;
    void run(const string& fn);

    void write(const string& fn) const;

private:
    vector<SurfaceParameterization> meshes;
    vector<BSplineSurface> surfs;

    vector<std::array<double,3>> spts; // sample points
    vector<Patch> patches;
    vector<Edge> edges;
    vector<vector<std::array<int,2>>> e2p; // edge to patch list, e2p[][][0] patch id, e2p[][][1] edge local id

    double fit_tol;
    int itr_max;
    int degree;
};

#endif //OPENSURF_MULTI_PATCH_FIT_H
