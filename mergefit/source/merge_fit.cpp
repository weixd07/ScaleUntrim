#include "merge_fit.h"
#include "multi_patch_fit.h"
#include "GenerateShell.h"
#include <iostream>
#include <fstream>
#include <quadgeneration.h>
//#include <boost/filesystem.hpp>
//#include <boost/algorithm/string.hpp>

MergeFit::MergeFit()
{
    temp_dir = "";
    tri_tol = 1.e-4;

    local_layer = -1;
    smooth_iter = -1;

    preserve_sharp = 1;
    preserve_boundary = 1;
    minimum_cost = 1;
    adaptive_scale = 1;
    quad_num = -1;
    angle_tol = 60.;
    magnitude_factor = -1.;

    fit_tol = 1.e-2;
    area_tol = 1.e-2;

    degree = 2;

    is_mirror = false;

    run_from = 0;
}

void MergeFit::set_temporary_directory(const std::string &dir)
{
    temp_dir = dir;
}

void MergeFit::set_run_from(int run_id)
{
    run_from = run_id;
}

void MergeFit::set_tri_mesh_tol(double tol)
{
    tri_tol = tol;
}

void MergeFit::set_quad_element_num(int num)
{
    quad_num = num;
}

void MergeFit::set_magnitude_factor(double mag)
{
    magnitude_factor = mag;
}

void MergeFit::set_angle_tol(double tol)
{
    angle_tol = tol;
}

void MergeFit::set_fitting_tol(double tol)
{
    fit_tol = tol;
}

void MergeFit::set_area_tol(double tol)
{
    area_tol = tol;
}

void MergeFit::set_mirror(const std::string &is_m)
{
    if (is_m == "yes")
      is_mirror = true;
    else
      is_mirror = false;
}

void MergeFit::set_local_layer(int n_layer)
{
    local_layer = n_layer;
}

void MergeFit::set_smooth_iteration(int n_itr)
{
    smooth_iter = n_itr;
}

void MergeFit::set_surface_degree(int deg_in)
{
    degree = deg_in;
}

void MergeFit::read_config_file(const std::string &fn)
{
    ifstream fin;
    fin.open(fn);
    if (!fin.is_open())
    {
        std::cout << "==============================================\n";
        cout << "Can't open the configuration file because it may not be provided.\n";
//        cout << "All parameters go default!\n";
    }
    else
    {
        std::string parameter_type;
        std::string temp;
        while (fin >> parameter_type)
        {
            if (parameter_type == "temp_dir:")
                fin >> temp_dir;
            else if (parameter_type == "tri_mesh_tolerance:")
                fin >> tri_tol;
            else if (parameter_type == "quad_face_num:")
                fin >> quad_num;
            else if (parameter_type == "magnitude_factor:")
                fin >> magnitude_factor;
            else if (parameter_type == "angle_for_sharp:")
                fin >> angle_tol;
            else if (parameter_type == "fit_tolerance:")
                fin >> fit_tol;
            else if (parameter_type == "perform_mirror:")
            {
                fin >> temp;
                if (temp == "yes")
                    is_mirror = false;
                else
                    is_mirror =false;
            }
            else if (parameter_type == "local_layer:")
                fin >> local_layer;
            else if (parameter_type == "smooth_iteration:")
                fin >> smooth_iter;
            else if (parameter_type == "run_from:")
                fin >> run_from;
            else if (parameter_type == "area_tolerance:")
                fin >> area_tol;
            else
                fin >> temp;
        }

        fin.close();
    }

    if (temp_dir.empty())
    {
        std::abort();
    }

//    if (temp_dir.empty())
//    {
//        boost::filesystem::path full_path(boost::filesystem::current_path() / "temp");
//        boost::filesystem::create_directory(full_path);
//        temp_dir = full_path.string();
//    }
//
//    if (temp_dir.back() != '\\' && temp_dir.back() != '/')
//    {
//#ifdef _WIN32
//        temp_dir += "\\";
//#else
//        temp_dir += "/";
//#endif
//    }
}

bool MergeFit::file_exist(const std::string &fn) const
{
    bool state = false;
    std::ifstream fin;
    fin.open(fn);
    if (fin.is_open())
    {
        state = true;
        fin.close();
    }
    return state;
}

int MergeFit::generate_tri_mesh_from_CAD(const std::string &cad, const std::string &tri)
{
    std::cout << "==============================================\n";
    std::cout << "Generating triangle mesh from the CAD file...\n";
    const std::string tri_stl = temp_dir + "tri.stl";
    triangular(cad.c_str(), tri_stl.c_str(), tri_tol);
    convert(tri_stl.c_str(), tri.c_str());
    area_in = compute_area(cad.c_str());
    return 0;
}

int MergeFit::fix_tri_mesh(const std::string &tri_obj, const std::string &tri_fix)
{
    std::cout << "==============================================\n";
    std::cout << "Fixing the triangle mesh...\n";
    fIxmesh fix_mesh;
    load load_mesh;
    load_mesh.Loader(tri_obj);
    load_mesh.initialize();
    fix_mesh.Fixmesh(load_mesh, local_layer, smooth_iter);
    load_mesh.Outer(tri_fix.c_str(), 1);
    return 0;
}

int MergeFit::generate_quad_mesh_from_tri_mesh(const std::string &tri_fix, const std::string &patch_info)
{
    std::cout << "==============================================\n";
    std::cout << "Generating quad mesh from the triangle mesh...\n";
    const std::string quad_vtk = temp_dir + "quad.vtk";
    const std::string patch_vtk = temp_dir + "patch";
    qflow::quadgeneration quad_gen;
    quad_gen.QMG(tri_fix, quad_vtk, patch_vtk, patch_info,
                 magnitude_factor, preserve_sharp, preserve_boundary, minimum_cost, adaptive_scale, angle_tol);

    // get mesh info if needed
    std::vector<Eigen::Vector3d> vertex;
    std::vector<Eigen::Vector4i> face;
    quad_gen.meshinformation(vertex, face);

    return 0;
}

int MergeFit::fit_NURBS_from_quad_patch(const std::string &patch_info, const std::string &surf)
{
    std::cout << "==============================================\n";
    std::cout << "Fitting NURBS...\n";
    MultiPatchFit fit;
    fit.set_fit_tol(fit_tol);
    fit.set_degree(degree);
    fit.read(patch_info);
    fit.run(surf);
    fit.write(surf);
    return 0;
}

int MergeFit::write_CAD_to_step(const std::string &surf_pat, const std::string &cad_out)
{
    std::cout << "==============================================\n";
    std::cout << "Writing CAD file in step format...\n";
    GenerateShell shell;
    shell.Read(surf_pat.c_str());
    shell.Shell(cad_out.c_str());
    area_out = shell.getArea(shell.sewing.SewedShape());
    if (is_mirror)
    {
        const std::string surf_mirror = temp_dir + "surf_mirror.step";
        shell.MirrorTransformation(surf_mirror.c_str());
    }
    const double area_diff = std::abs(area_in - area_out) / area_in;
    std::cout << "Area difference: " << area_diff << endl;
    std::cout << "Area tolerance : " << area_tol << endl;
    return 0;
}

int MergeFit::run(const std::string& cad_file_in, const std::string& cad_file_out, const std::string& config_file)
{
    read_config_file(config_file);

    // temporary files
    const std::string tri_obj = temp_dir + "tri.obj";
    const std::string tri_fix = temp_dir + "tri_fix.obj";
    const std::string patch_info = temp_dir + "patch.txt";
    const std::string surf = temp_dir + "surf";
    const std::string surf_pat = temp_dir + "surf.pat";

    if (run_from == 0)
    {
        generate_tri_mesh_from_CAD(cad_file_in, tri_obj);
        fix_tri_mesh(tri_obj, tri_fix);
        generate_quad_mesh_from_tri_mesh(tri_fix, patch_info);
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 1)
    {
        if (!file_exist(tri_obj))
        {
            std::cerr << "The triangle mesh does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        fix_tri_mesh(tri_obj, tri_fix);
        generate_quad_mesh_from_tri_mesh(tri_fix, patch_info);
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 2)
    {
        if (!file_exist(tri_fix))
        {
            std::cerr << "The fixed triangle mesh does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        generate_quad_mesh_from_tri_mesh(tri_fix, patch_info);
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 3)
    {
        if (!file_exist(patch_info))
        {
            std::cerr << "The patch info of the quad mesh does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 4)
    {
        if (!file_exist(patch_info))
        {
            std::cerr << "The patch info of the fitted surface does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else
    {
        std::cerr << "Unknown starting point to run the program: " << run_from << endl;
        std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
        return 1;
    }

    std::cout << "Successfully done merging and fitting the surface!\n";
    return 0;
}

int MergeFit::run(const std::string& cad_file_in, const std::string& cad_file_out)
{
    // temporary files
    const std::string tri_obj = temp_dir + "tri.obj";
    const std::string tri_fix = temp_dir + "tri_fix.obj";
    const std::string patch_info = temp_dir + "patch.txt";
    const std::string surf = temp_dir + "surf";
    const std::string surf_pat = temp_dir + "surf.pat";

    if (run_from == 0)
    {
        generate_tri_mesh_from_CAD(cad_file_in, tri_obj);
        fix_tri_mesh(tri_obj, tri_fix);
        generate_quad_mesh_from_tri_mesh(tri_fix, patch_info);
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 1)
    {
        if (!file_exist(tri_obj))
        {
            std::cerr << "The triangle mesh does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        fix_tri_mesh(tri_obj, tri_fix);
        generate_quad_mesh_from_tri_mesh(tri_fix, patch_info);
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 2)
    {
        if (!file_exist(tri_fix))
        {
            std::cerr << "The fixed triangle mesh does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        generate_quad_mesh_from_tri_mesh(tri_fix, patch_info);
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 3)
    {
        if (!file_exist(patch_info))
        {
            std::cerr << "The patch info of the quad mesh does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        fit_NURBS_from_quad_patch(patch_info, surf);
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else if (run_from == 4)
    {
        if (!file_exist(patch_info))
        {
            std::cerr << "The patch info of the fitted surface does not exist!\n";
            std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
            return 1;
        }
        write_CAD_to_step(surf_pat, cad_file_out);
    }
    else
    {
        std::cerr << "Unknown starting point to run the program: " << run_from << endl;
        std::cerr << "Set the parameter <run_from> as 0 and rerun the program.\n";
        return 1;
    }

    std::cout << "Successfully done merging and fitting the surface!\n";
    return 0;
}

