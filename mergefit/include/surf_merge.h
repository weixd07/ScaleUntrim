#ifndef OPENSURF_SURF_MERGE_H
#define OPENSURF_SURF_MERGE_H

#include "b_spline_surface.h"
#include "surface_parameterization.h"

class SurfMerge
{
public:
    SurfMerge();
    SurfMerge(double fit_tol_in, double sharp_tol_in, int itr_max_in);

    void set_surf_info(const std::array<int,2>& deg, const std::array<vector<double>,2>& kv);
    void set_fit_tol(double tol);
    void set_sharp_tol(double tol);
    void set_max_iteration(int itr_max_in);

    void read(const string& tri_mesh, const string corners = "");
    void run(const string& fn);
    void write(const string& fn);

private:
    SurfaceParameterization mesh;
    BSplineSurface surf;
    double fit_tol;
    double sharp_tol; // in angle
    int itr_max;
};

#endif //OPENSURF_SURF_MERGE_H
