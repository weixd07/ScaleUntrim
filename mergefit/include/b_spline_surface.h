#ifndef OPENSURF_B_SPLINE_SURFACE_H
#define OPENSURF_B_SPLINE_SURFACE_H

#include <array>
#include <vector>
#include <set>
#include <string>
#include "basic_data_structure.h"
#include "BSplineBasis.h"

class BSplineSurface
{
public:
    BSplineSurface();
    BSplineSurface(const std::array<int,2>& deg_in, const std::array<std::vector<double>,2>& kv_in);
    BSplineSurface(const std::array<int,2>& deg_in, const std::array<int,2>& nel);

    void set_bspline_basis(const std::array<int,2>& deg_in, const std::array<std::vector<double>,2>& kv_in);
    void set_bspline_basis(int deg_in, const std::vector<double>& kv_in, int dir); //dir == 0 for u, dir==1 for v
    void set_bspline_basis(int nel, int dir);
    void set_bspline_basis(const std::array<int,2>& nel);

    void set_degree(int deg_in);

    std::array<int,2> get_cp_num() const;
    std::array<int,2> get_degree() const;
    const std::vector<double>& get_knot_vector(int dir) const;
    const std::vector<std::array<double,3>>& get_control_points() const;
    void insert_C0_knots(const std::array<set<double>,2>& knots);
    void find_IEN(const std::array<double,2>& u, std::vector<int>& IEN) const;
    void evaluate_basis(const std::array<double,2>& u, std::vector<double>& Nx) const;
    bool fit(const std::vector<Point>& pts);
    bool fit_grid(const std::vector<Point>& pts, const int nu, const int nv);
    bool fit_grid_boundary(const std::vector<Point>& pts, const int nu, const int nv, const int bnd_id);
    bool fit_grid_interior(const std::vector<Point>& pts, const int nu, const int nv);
    void compute_error(const std::vector<Point>& pts);
    double get_max_error() const;
    void identify(const std::vector<Point> &pts, std::array<std::vector<double>,2>& knots) const;
    void refine(const std::vector<double>& knots, int dir); //dir == 0 for u, dir==1 for v
    void refine(const std::array<std::vector<double>,2>& knots);
    void refine();
    void refine(int dir);
    void geom_map(const std::array<double,2>& u, std::array<double,3>& pt) const;
    void visualize_vtk(const std::string& fn) const;
    void visualize_control_mesh(const std::string& fn) const;

private:
    std::vector<std::array<double,3>> cp;
    std::array<BSplineBasis,2> bsp;
    std::vector<double> err;
    double err_max;
};

#endif //OPENSURF_B_SPLINE_SURFACE_H
