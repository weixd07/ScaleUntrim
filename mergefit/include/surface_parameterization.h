#ifndef OPENSURF_SURFACE_PARAMETERIZATION_H
#define OPENSURF_SURFACE_PARAMETERIZATION_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "basic_data_structure.h"
#include "b_spline_surface.h"

using namespace std;
using namespace Eigen;

class SurfaceParameterization
{
public:
    SurfaceParameterization();
    SurfaceParameterization(double tol_angle_in);

    void read_tri_mesh_obj(const string& fname);
    void read_corner_IDs(const string& fname);

    void retrieve_from_quad_mesh(const vector<std::array<double,3>>& quad_pts, const vector<vector<int>>& grid_ids);
    void set_corners(const std::array<int,4>& cn_in);

    void build_edges();
    void build_connectivity();
    void set_boundary_flags();
    bool check_mesh_orientation();

    void find_boundary_points();
    void find_sharp_points();
    void find_corners();
    bool is_sharp_point(int bnd_loc, double& dot_product) const; // tol_angle is in degree
    bool is_concave(int bnd_loc) const;
    void get_directions(int bnd_loc, Vector3d& v_prev, Vector3d& v_next) const; // normalized vectors
    void split(int ic, int ist, const vector<double>& arc_len);
    void remove_corners_from_sharp_points();
    void set_C0();
    void reorder_boundary_points();

    void parameterize_boundary(); // chord length
    void parameterize_interior(); // Floater's mean value coordinates
    double compute_mean_value_coordinates(int pid, int edid);

    double get_distance(const std::array<double, 3>& pt1, const std::array<double, 3>& pt2) const;
    double get_dimension_scale() const;
    const vector<Point>& get_points() const;
    void get_sharp_point_param_coor(std::array<set<double>,2>& ush) const;

    void visualize_tri_mesh(const string& fname) const;
    void visualize_param_mesh(const string& fname) const;

    void clear();

    void run();

private:
    vector<Point> pts;
    vector<Edge> edges;
    vector<TriFace> faces;
    std::array<int, 4> cnid; // corner IDs, counter-clockwise

    vector<vector<std::array<int,2>>> p2e; // edges connecting to a point
    vector<vector<std::array<int,2>>> p2f; // faces sharing a point
    vector<vector<std::array<int,2>>> e2f; // faces sharing an edge

    vector<int> pid_bnd; // boundary point IDs, counter-clockwise
    vector<int> sharp_loc; // sharp point local IDs in pid_bnd, counter-clockwise
    vector<int> sharp_sort; // sharp point local IDs in pid_bnd, ordered according to dot product

    const double tol_angle; // angle to detect sharp points
};

#endif //OPENSURF_SURFACE_PARAMETERIZATION_H
