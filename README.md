# Scalable Field-Aligned Reparameterization for Trimmed NURBS

## Overview
ScaleUntrim is an open-source software package to reconstruct trimmed CAD models as a watertight spline representation (i.e., multi-patch NURBS). The copyright of the code belongs to Prof. Xiaodong Wei's group at Shanghai Jiao Tong University, China.

## Software

The software supports cmake build for Linux/Mac/Windows systems, which converts a trimmed CAD modal  `(.stp/.igs)` into an untrimmed NURBS `(.step)`.

we recommend Visual Studio 2022 in Windows for compiling.

## Install

```
git clone 
mkdir build
cd build
cmake .. 
make
```

## External Dependencies

- Boost
- Eigen
- OpenCascade

## Parameter 

All parameters can be adjusted in the `setting.config` file."

```text
temp_dir: D:\\mergefit_v1\\io\\vs17\\temp1\\
run_from: 0
tri_mesh_tolerance: 1.e-4
quad_face_num: -1
magnitude_factor: -1
angle_for_sharp: 60
fit_tolerance: 1.e-2
```

`temp_dir` is to store intermediate files generated during the pipeline, such as triangular meshes, quadrilateral meshes, and so on.

`run_from` indicates from which step in the pipeline to start running

`tri_mesh_tolerance` controls the density of triangle mesh

`magnitude_factor`controls the density of the quadrilateral mesh

`angle_for _sharp`recognizes sharp feature

`fit_tolerance`indicates fitting error

## Licenses

## Authors

- <u>Zheng Wei</u>
- <u>Haoyu Su</u>
- <u>Shizhou Luo</u>
- <u>Xiaodong Wei</u>
