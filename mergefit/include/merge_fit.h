#ifndef OPENSURF_MERGEFIT_H
#define OPENSURF_MERGEFIT_H

#include "cad2tri.h"
#include "fixmesh.h"
#include "multi_patch_fit.h"
#include "GenerateShell.h"
#include <string>

class MergeFit
{
public:
    MergeFit();

    void set_temporary_directory(const string& dir);
    void set_run_from(int run_id);
    void set_tri_mesh_tol(double tol);
    void set_quad_element_num(int num);
    void set_magnitude_factor(double mag);
    void set_angle_tol(double tol);
    void set_fitting_tol(double tol);
    void set_area_tol(double tol);
    void set_mirror(const std::string& is_m);
    void set_local_layer(int n_layer);
    void set_smooth_iteration(int n_itr);
    void set_surface_degree(int deg_in);

    void read_config_file(const std::string& fn);
    bool file_exist(const std::string& fn) const;
    int generate_tri_mesh_from_CAD(const std::string& cad_in, const std::string& tri);
    int fix_tri_mesh(const std::string& tri_in, const std::string& tri_fix);
    int generate_quad_mesh_from_tri_mesh(const std::string& tri_fix, const std::string& patch);
    int fit_NURBS_from_quad_patch(const std::string& patch, const std::string& surf);
    int write_CAD_to_step(const std::string& surf, const std::string& cad_out);
    int run(const std::string& cad_file_in, const std::string& cad_file_out, const std::string& config);
    int run(const std::string& cad_file_in, const std::string& cad_file_out);

private:
    string temp_dir;
    double tri_tol;

    int local_layer;
    int smooth_iter;

    int preserve_sharp;
    int preserve_boundary;
    int minimum_cost;
    int adaptive_scale;
    double angle_tol;
    int quad_num;
    double magnitude_factor;

    double fit_tol;
    double area_tol;

    int degree;

    bool is_mirror;

    int run_from;

    double area_in;
    double area_out;
    double tri_quad_diff;
    double fit_err;
};

#endif //OPENSURF_MERGEFIT_H
