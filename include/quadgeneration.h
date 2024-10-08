#pragma execution_character_set("utf-8")



#ifndef _QUADGENERATION_H_
#define _QUADGENERATION_H_
#include <atomic>
#include <condition_variable>
#ifdef WITH_TBB
#include <tbb/tbb.h>
#endif

#include <Eigen/Core>
#include <Eigen/Dense>
#include <list>
#include <map>
#include <set>
#include <unordered_set>
#include "adjacent-matrix.hpp"
#include "disajoint-tree.hpp"
#include "field-math.hpp"
#include "hierarchy.hpp"
#include "post-solver.hpp"
#include "serialize.hpp"
#pragma once



namespace qflow {

class quadgeneration {
   public:
    void QMG(std::string input_tri, std::string output_quad, std::string output_patch_file,
             std::string output_patch_information, double magnitude_factor, int preserve_sharp,
             int preserve_boundary, int minimum_cost, int adaptive_scale,double angle);
    void meshinformation(std::vector<Vector3d> &V_, std::vector<Vector4i> &F_);


    std::vector<Vector3d> Ver;
    std::vector<Vector4i> Fer;
};

}
#endif
