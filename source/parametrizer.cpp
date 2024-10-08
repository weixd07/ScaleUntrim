









#include "parametrizer.hpp"
#include "config.hpp"
#include "dedge.hpp"
#include "field-math.hpp"
#include "flow.hpp"
#include "localsat.hpp"
#include "optimizer.hpp"
#include "subdivide.hpp"

#include "dset.hpp"

#include <Eigen/Sparse>
#include <fstream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
using namespace std;



namespace qflow {

void Parametrizer::ComputeIndexMap(int with_scale) {
    // build edge info
    auto& V = hierarchy.mV[0];
    auto& F = hierarchy.mF;
    auto& Q = hierarchy.mQ[0];
    auto& N = hierarchy.mN[0];
    auto& O = hierarchy.mO[0];
    auto& S = hierarchy.mS[0];
    auto& holeP1 = hierarchy.mholeP1;
    auto& holeposition = hierarchy.mholeposition;

    
   
    // ComputeOrientationSingularities();

    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (Boundary_edges[i] == 1) {
            sharp_edges[i] = 0;
        }
    }

    BuildEdgeInfo();

    if (flag_preserve_sharp) {
        //        ComputeSharpO();
    }
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 1 || Boundary_edges[i] == 1) {
            int e = face_edgeIds[i / 3][i % 3];
            if (edge_diff[e][0] * edge_diff[e][1] != 0) {
                Vector3d d = O.col(edge_values[e].y) - O.col(edge_values[e].x);
                Vector3d q = Q.col(edge_values[e].x);
                Vector3d n = N.col(edge_values[e].x);
                Vector3d qy = n.cross(q);
                if (abs(q.dot(d)) > qy.dot(d))
                    edge_diff[e][1] = 0;
                else
                    edge_diff[e][0] = 0;
            }
        }
    }
    std::map<int, std::pair<Vector3d, Vector3d>> Boundary_constraints;
    std::set<int> Boundaryvert;
    for (int i = 0; i < Boundary_edges.size(); ++i) {
        if (Boundary_edges[i]) {
            Boundaryvert.insert(F(i % 3, i / 3));
            Boundaryvert.insert(F((i + 1) % 3, i / 3));
        }
    }
    std::map<int, std::pair<Vector3d, Vector3d>> sharp_constraints;
    std::set<int> sharpvert;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i]) {
            sharpvert.insert(F(i % 3, i / 3));
            sharpvert.insert(F((i + 1) % 3, i / 3));
        }
    }

    allow_changes.resize(edge_diff.size() * 2, 1);
    for (int i = 0; i < sharp_edges.size(); ++i) {
        int e = face_edgeIds[i / 3][i % 3];
        if (sharpvert.count(edge_values[e].x) && sharpvert.count(edge_values[e].y)) {
            if (sharp_edges[i] != 0) {
                for (int k = 0; k < 2; ++k) {
                    if (edge_diff[e][k] == 0) {
                        allow_changes[e * 2 + k] = 0;
                    }
                }
            }
        }
        if (Boundaryvert.count(edge_values[e].x) && Boundaryvert.count(edge_values[e].y)) {
            if (Boundary_edges[i] != 0) {
                for (int k = 0; k < 2; ++k) {
                    if (edge_diff[e][k] == 0) {
                        allow_changes[e * 2 + k] = 0;
                    }
                }
            }
        }
    }
#ifdef LOG_OUTPUT
    printf("Build Integer Constraints...\n");
#endif
    BuildIntegerConstraints();

    ComputeMaxFlow();
    // potential bug
#ifdef LOG_OUTPUT
    printf("subdivide...\n");
#endif
    subdivide_edgeDiff(F, V, N, Q, O, &hierarchy.mS[0], V2E, hierarchy.mE2E, boundary, nonManifold,
                       edge_diff, edge_values, face_edgeOrients, face_edgeIds, sharp_edges,
                       singularities, Boundary_edges);

    allow_changes.clear();
    allow_changes.resize(edge_diff.size() * 2, 1);
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 0 || Boundary_edges[i] == 0) continue;
        int e = face_edgeIds[i / 3][i % 3];
        for (int k = 0; k < 2; ++k) {
            if (edge_diff[e][k] == 0) allow_changes[e * 2 + k] = 0;
        }
    }

#ifdef LOG_OUTPUT
    printf("Fix flip advance...\n");
    int t1 = GetCurrentTime64();
#endif
    FixFlipHierarchy();
    subdivide_edgeDiff(F, V, N, Q, O, &hierarchy.mS[0], V2E, hierarchy.mE2E, boundary, nonManifold,
                       edge_diff, edge_values, face_edgeOrients, face_edgeIds, sharp_edges,
                       singularities, Boundary_edges);

    FixFlipSat();
    auto Vtest = V;
    for (int i = 0; i < V.cols(); ++i) {
        double a1 = V.col(i)[0];
        double a2 = V.col(i)[1];
        double a3 = V.col(i)[2];
        holeposition.push_back(a1);
        holeposition.push_back(a2);
        holeposition.push_back(a3);
    }

   

#ifdef LOG_OUTPUT
    int t2 = GetCurrentTime64();
    printf("Flip use %lf\n", (t2 - t1) * 1e-3);
    printf("Post Linear Solver...\n");
#endif
    std::set<int> sharp_vertices;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 1) {
            sharp_vertices.insert(F(i % 3, i / 3));
            sharp_vertices.insert(F((i + 1) % 3, i / 3));
        }
    }
    std::set<int> Boundary_vertices;
    for (int i = 0; i < Boundary_edges.size(); ++i) {
        if (Boundary_edges[i] == 1) {
            Boundary_vertices.insert(F(i % 3, i / 3));
            Boundary_vertices.insert(F((i + 1) % 3, i / 3));
        }
    }

    

    Optimizer::optimize_positions_sharp(hierarchy, edge_values, edge_diff, sharp_edges,
                                        sharp_vertices, sharp_constraints, Boundary_edges,
                                        Boundary_vertices, Boundary_constraints, with_scale);

    Optimizer::optimize_positions_fixed(hierarchy, edge_values, edge_diff, sharp_vertices,
                                        sharp_constraints, sharp_edges, Boundary_vertices,
                                        Boundary_constraints, Boundary_edges, flag_adaptive_scale);

    AdvancedExtractQuad();

    FixValence();

    std::vector<int> sharp_o(O_compact.size(), 0);
    std::map<int, std::pair<Vector3d, Vector3d>> compact_sharp_constraints;
    for (int i = 0; i < Vset.size(); ++i) {
        int sharpv = -1;
        for (auto& p : Vset[i]) {
            if (sharp_constraints.count(p)) {
                sharpv = p;
                sharp_o[i] = 1;
                if (compact_sharp_constraints.count(i) == 0 ||
                    compact_sharp_constraints[i].second != Vector3d::Zero()) {
                    compact_sharp_constraints[i] = sharp_constraints[sharpv];
                    O_compact[i] = O.col(sharpv);
                    compact_sharp_constraints[i].first = O_compact[i];
                }
            }
        }
    }
    std::vector<int> boundary_o(O_compact.size(), 0);
    std::map<int, std::pair<Vector3d, Vector3d>> compact_boundary_constraints;
    for (int i = 0; i < Vset.size(); ++i) {
        int boundaryv = -1;
        for (auto& p : Vset[i]) {
            if (Boundary_constraints.count(p)) {
                boundaryv = p;
                boundary_o[i] = 1;
                if (compact_boundary_constraints.count(i) == 0 ||
                    compact_boundary_constraints[i].second != Vector3d::Zero()) {
                    compact_boundary_constraints[i] = Boundary_constraints[boundaryv];
                    O_compact[i] = O.col(boundaryv);
                    compact_boundary_constraints[i].first = O_compact[i];
                }
            }
        }
    }

    std::map<std::pair<int, int>, int> o2e;
    for (int i = 0; i < F_compact.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            int v1 = F_compact[i][j];
            int v2 = F_compact[i][(j + 1) % 4];
            o2e[std::make_pair(v1, v2)] = i * 4 + j;
        }
    }
    std::vector<std::vector<int>> v2o(V.cols());
    for (int i = 0; i < Vset.size(); ++i) {
        for (auto v : Vset[i]) {
            v2o[v].push_back(i);
        }
    }
    std::vector<Vector3d> diffs(F_compact.size() * 4, Vector3d(0, 0, 0));
    std::vector<int> diff_count(F_compact.size() * 4, 0);
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(j, i);
            int v2 = F((j + 1) % 3, i);
            if (v1 != edge_values[face_edgeIds[i][j]].x) continue;
            if (edge_diff[face_edgeIds[i][j]].array().abs().sum() != 1) continue;
            if (v2o[v1].size() > 1 || v2o[v2].size() > 1) continue;
            for (auto o1 : v2o[v1]) {
                for (auto o2 : v2o[v2]) {
                    auto key = std::make_pair(o1, o2);
                    if (o2e.count(key)) {
                        int dedge = o2e[key];
                        Vector3d q_1 = Q.col(v1);
                        Vector3d q_2 = Q.col(v2);
                        Vector3d n_1 = N.col(v1);
                        Vector3d n_2 = N.col(v2);
                        Vector3d q_1_y = n_1.cross(q_1);
                        Vector3d q_2_y = n_2.cross(q_2);
                        auto index = compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
                        double s_x1 = S(0, v1), s_y1 = S(1, v1);
                        double s_x2 = S(0, v2), s_y2 = S(1, v2);
                        int rank_diff = (index.second + 4 - index.first) % 4;
                        if (rank_diff % 2 == 1) std::swap(s_x2, s_y2);
                        Vector3d qd_x = 0.5 * (rotate90_by(q_2, n_2, rank_diff) + q_1);
                        Vector3d qd_y = 0.5 * (rotate90_by(q_2_y, n_2, rank_diff) + q_1_y);
                        double scale_x = (with_scale ? 0.5 * (s_x1 + s_x2) : 1) * hierarchy.mScale;
                        double scale_y = (with_scale ? 0.5 * (s_y1 + s_y2) : 1) * hierarchy.mScale;
                        Vector2i diff = edge_diff[face_edgeIds[i][j]];
                        Vector3d C = diff[0] * scale_x * qd_x + diff[1] * scale_y * qd_y;

                        diff_count[dedge] += 1;
                        diffs[dedge] += C;
                        auto key = std::make_pair(o2, o1);
                        if (o2e.count(key)) {
                            int dedge = o2e[key];
                            diff_count[dedge] += 1;
                            diffs[dedge] -= C;
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < F.cols(); ++i) {
        Vector2i d1 = rshift90(edge_diff[face_edgeIds[i][0]], face_edgeOrients[i][0]);
        Vector2i d2 = rshift90(edge_diff[face_edgeIds[i][1]], face_edgeOrients[i][1]);
        if (d1[0] * d2[1] - d1[1] * d2[0] < 0) {
            for (int j = 0; j < 3; ++j) {
                int v1 = F(j, i);
                int v2 = F((j + 1) % 3, i);
                for (auto o1 : v2o[v1]) {
                    for (auto o2 : v2o[v2]) {
                        auto key = std::make_pair(o1, o2);
                        if (o2e.count(key)) {
                            int dedge = o2e[key];
                            diff_count[dedge] = 0;
                            diffs[dedge] = Vector3d(0, 0, 0);
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < diff_count.size(); ++i) {
        if (diff_count[i] != 0) {
            diffs[i] /= diff_count[i];
            diff_count[i] = 1;
        }
    }
    auto unchangB = hierarchy.munchangeB;
    vector<Vector3d> O_temp = O_compact;
    int need = 1;
    vector<int> quad_singularity;
    if (need) {
        vector<int> point_number(O_compact.size());
        for (int i = 0; i < F_compact.size(); ++i) {
            int v1, v2, v3, v4;
            v1 = F_compact[i][0];
            v2 = F_compact[i][1];
            v3 = F_compact[i][2];
            v4 = F_compact[i][3];
            point_number[v1] += 1;
            point_number[v2] += 1;
            point_number[v3] += 1;
            point_number[v4] += 1;
        }
        for (int i = 0; i < O_compact.size(); ++i) {
            if (boundary_o[i] == 0) continue;
            if (point_number[i] > 2) {
                quad_singularity.push_back(i);
            }
        }
        vector<int> boundaryvertex;
        for (int i = 0; i < 4 * F_compact.size(); ++i) {
            if (E2E_compact[i] == -1) {
                int u0 = F_compact[i / 4][i % 4];
                int u1 = F_compact[i / 4][(i + 1) % 4];
                boundaryvertex.push_back(u0);
                boundaryvertex.push_back(u1);
            }
        }
        for (int i = 0; i < boundaryvertex.size(); ++i) {
            if (boundaryvertex[i] == -1) {
                continue;
            }
            int point1, point2, point3, posi;
            std::vector<int> loop{};
            point1 = boundaryvertex[i];
            boundaryvertex[i] = -1;
            if (i % 2 == 0) {
                point2 = boundaryvertex[i + 1];
                boundaryvertex[i + 1] = -1;
                posi = i + 1;
            } else {
                point2 = boundaryvertex[i - 1];
                boundaryvertex[i - 1] = -1;
                posi = i - 1;
            }
            loop.push_back(point1);
            loop.push_back(point2);
            do {
                for (int j = 0; j < boundaryvertex.size(); ++j) {
                    if (boundaryvertex[j] == -1) continue;
                    if (boundaryvertex[j] == point2 && j != posi) {
                        boundaryvertex[j] = -1;
                        if (j % 2 == 0) {
                            point3 = boundaryvertex[j + 1];
                            boundaryvertex[j + 1] = -1;
                            posi = j + 1;
                        } else {
                            point3 = boundaryvertex[j - 1];
                            boundaryvertex[j - 1] = -1;
                            posi = j - 1;
                        }
                    }
                }
                loop.push_back(point2);
                loop.push_back(point3);
                point2 = point3;
            } while (point3 != point1);
            if (loop.size() < boundaryvertex.size() / 5) continue;
            for (int k = 0; k < loop.size(); ++k) {
                if (k % 2 == 1) continue;
                int i1 = loop[k];
                int i2 = loop[k + 1];
                if (find(quad_singularity.begin(), quad_singularity.end(), i2) !=
                    quad_singularity.end())
                    continue;
                if (find(quad_singularity.begin(), quad_singularity.end(), i1) !=
                    quad_singularity.end())
                    continue;
                if (find(unchangB.begin(), unchangB.end(), i2) != unchangB.end()) continue;
                if (k == loop.size() - 2) {
                    if (find(quad_singularity.begin(), quad_singularity.end(), loop[1]) !=
                        quad_singularity.end())
                        continue;
                } else {
                    if (find(quad_singularity.begin(), quad_singularity.end(), loop[k + 3]) !=
                        quad_singularity.end())
                        continue;
                }

                int i3;
                for (int m = 0; m < F_compact.size(); ++m) {
                    if (i1 == F_compact[m][0] && i2 == F_compact[m][1]) {
                        i3 = F_compact[m][2];
                        break;
                    }
                    if (i1 == F_compact[m][1] && i2 == F_compact[m][2]) {
                        i3 = F_compact[m][3];
                        break;
                    }
                    if (i1 == F_compact[m][2] && i2 == F_compact[m][3]) {
                        i3 = F_compact[m][0];
                        break;
                    }
                    if (i1 == F_compact[m][3] && i2 == F_compact[m][0]) {
                        i3 = F_compact[m][1];
                        break;
                    }
                    if (i1 == F_compact[m][1] && i2 == F_compact[m][0]) {
                        i3 = F_compact[m][3];
                        break;
                    }
                    if (i1 == F_compact[m][2] && i2 == F_compact[m][1]) {
                        i3 = F_compact[m][0];
                        break;
                    }
                    if (i1 == F_compact[m][3] && i2 == F_compact[m][2]) {
                        i3 = F_compact[m][1];
                        break;
                    }
                    if (i1 == F_compact[m][0] && i2 == F_compact[m][3]) {
                        i3 = F_compact[m][2];
                        break;
                    }
                }
                Vector3d e1 = O_compact[i1] - O_compact[i2];
                Vector3d e2 = O_compact[i3] - O_compact[i2];
                double cos = e1.dot(e2) / (sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]) *
                                           sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]));
                int i4;
                if (cos < 0.3 && cos > -0.3) continue;
                if (k == loop.size() - 2) {
                    i4 = loop[2];
                } else {
                    i4 = loop[k + 3];
                }
                Vector3d e3 = O_compact[i4] - O_compact[i2];
                double cos1 = e1.dot(e3) / (sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]) *
                                            sqrt(e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]));
                /*if (cos1 > -0.7) continue;*/
                Vector3d edge3 = e1.normalized() * e1.dot(e2) /
                                 sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);

                Vector3d line = edge3 - (O_compact[i3] - O_compact[i2]);

                O_temp[i2] = O_compact[i3] + line;
                /*if (cos > 0) {

                } else {
                    Vector3d edge3 = e3.normalized() * e3.dot(e2) /
                                     sqrt(e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]);
                    O_compact[i2] += edge3;
                }*/
            }
        }
    }
    O_compact = O_temp;
    std::vector<std::vector<int>> alignP2(O_compact.size());
    Vector3d posision9;
    for (int i = 0; i < O_compact.size(); ++i) {
        if (boundary_o[i] == 0) continue;
        posision9 = O_compact[i];
        double testnumber9 = 10000;
        int testpoint9 = 0;
        Vector3d position8;
        for (auto& v : Boundary_vertices) {
            position8[0] = holeposition[3 * v];
            position8[1] = holeposition[3 * v + 1];
            position8[2] = holeposition[3 * v + 2];
            double lengthh9 = (position8 - posision9).norm();
            if (lengthh9 < testnumber9) {
                testnumber9 = lengthh9;
                testpoint9 = v;
            }
        }
        alignP2[i].push_back(testpoint9);
    }
    Vector3d posision61;
    for (int i = 0; i < O_compact.size(); ++i) {
        if (boundary_o[i] == 0) continue;
        posision61 = O_compact[i];
        double testnumber61 = 10000;
        int testpoint61 = 0;
        for (auto& v : Boundary_vertices) {
            if (alignP2[i][0] == v) continue;
            Vector3d position71;
            position71[0] = holeposition[3 * v];
            position71[1] = holeposition[3 * v + 1];
            position71[2] = holeposition[3 * v + 2];
            double lengthh61 = (position71 - posision61).norm();
            if (lengthh61 < testnumber61) {
                testnumber61 = lengthh61;
                testpoint61 = v;
            }
        }
        alignP2[i].push_back(testpoint61);
    }
    for (int i = 0; i < O_compact.size(); ++i) {
        if (boundary_o[i] == 0) continue;
        if (find(unchangB.begin(), unchangB.end(), i) != unchangB.end()) continue;
        if (find(quad_singularity.begin(), quad_singularity.end(), i) != quad_singularity.end())
            continue;
        /*if (linetype[i]) continue;*/
        int v1 = alignP2[i][0];
        int v2 = alignP2[i][1];
        Vector3d pointv1, pointv2, pointv3;
        pointv1[0] = holeposition[3 * v1];
        pointv1[1] = holeposition[3 * v1 + 1];
        pointv1[2] = holeposition[3 * v1 + 2];
        pointv2[0] = holeposition[3 * v2];
        pointv2[1] = holeposition[3 * v2 + 1];
        pointv2[2] = holeposition[3 * v2 + 2];
        pointv3 = O_compact[i];
        Vector3d edge1 = pointv2 - pointv1;
        Vector3d edge2 = pointv3 - pointv1;
        Vector3d edge3 = edge1.normalized() * edge1.dot(edge2) /
                         sqrt(edge1[0] * edge1[0] + edge1[1] * edge1[1] + edge1[2] * edge1[2]);

        O_compact[i] = pointv1 + edge3;
    }
   

    Optimizer::optimize_positions_dynamic(
        hierarchy, F, V, N, Q, Vset, O_compact, F_compact, V2E_compact, E2E_compact,
        sqrt(surface_area / F_compact.size()), diffs, diff_count, o2e, sharp_o,
        compact_sharp_constraints, boundary_o, compact_boundary_constraints);
    vector<double> weizhi(3 * V.cols(), 0);
    vector<int> chongfu;
    for (int i = 0; i < Vset.size(); ++i) {
        for (int j = 0; j < Vset[i].size(); ++j) {
            weizhi[3 * Vset[i][j]] = O_compact[i][0];
            weizhi[3 * Vset[i][j] + 1] = O_compact[i][1];
            weizhi[3 * Vset[i][j] + 2] = O_compact[i][2];
            /*if (Vset[i][j] == 6565 || Vset[i][j] == 2406) {
                chongfu.push_back(i);
            }
            if (find(duoyu.begin(), duoyu.end(), Vset[i][j]) != duoyu.end()) {
                chongfu.push_back(i);
            }
            }*/
        }
    }
    std::map<int, Vector2i> new_sing;
    vector<int> qiyidian;
    for (int f = 0; f < F.cols(); ++f) {
        Vector2i index = Vector2i::Zero();
        uint32_t i0 = F(0, f), i1 = F(1, f), i2 = F(2, f);

        Vector3d q[3] = {Q.col(i0).normalized(), Q.col(i1).normalized(), Q.col(i2).normalized()};
        Vector3d n[3] = {N.col(i0), N.col(i1), N.col(i2)};
        Vector3d v[3] = {Vtest.col(i0), Vtest.col(i1), Vtest.col(i2)};
        Vector3d v1, v2, v3;
        v1[0] = weizhi[3 * i0];
        v1[1] = weizhi[3 * i0 + 1];
        v1[2] = weizhi[3 * i0 + 2];
        v2[0] = weizhi[3 * i1];
        v2[1] = weizhi[3 * i1 + 1];
        v2[2] = weizhi[3 * i1 + 2];
        v3[0] = weizhi[3 * i2];
        v3[1] = weizhi[3 * i2 + 1];
        v3[2] = weizhi[3 * i2 + 2];
        if (v1[0] == 0 && v1[1] == 0 && v1[2] == 0) continue;
        if (v2[0] == 0 && v2[1] == 0 && v2[2] == 0) continue;
        if (v3[0] == 0 && v3[1] == 0 && v3[2] == 0) continue;
        Vector3d o[3] = {v1, v2, v3};
        int best[3];
        double best_dp = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < 4; ++i) {
            Vector3d v0 = rotate90_by(q[0], n[0], i);
            for (int j = 0; j < 4; ++j) {
                Vector3d v1 = rotate90_by(q[1], n[1], j);
                for (int k = 0; k < 4; ++k) {
                    Vector3d v2 = rotate90_by(q[2], n[2], k);
                    double dp = std::min(std::min(v0.dot(v1), v1.dot(v2)), v2.dot(v0));
                    if (dp > best_dp) {
                        best_dp = dp;
                        best[0] = i;
                        best[1] = j;
                        best[2] = k;
                    }
                }
            }
        }
        for (int k = 0; k < 3; ++k) q[k] = rotate90_by(q[k], n[k], best[k]);

        for (int k = 0; k < 3; ++k) {
            int kn = k == 2 ? 0 : (k + 1);
            double scale_x = hierarchy.mScale, scale_y = hierarchy.mScale,
                   scale_x_1 = hierarchy.mScale, scale_y_1 = hierarchy.mScale;
            if (flag_adaptive_scale) {
                scale_x *= hierarchy.mS[0](0, F(k, f));
                scale_y *= hierarchy.mS[0](1, F(k, f));
                scale_x_1 *= hierarchy.mS[0](0, F(kn, f));
                scale_y_1 *= hierarchy.mS[0](1, F(kn, f));
                if (best[k] % 2 != 0) std::swap(scale_x, scale_y);
                if (best[kn] % 2 != 0) std::swap(scale_x_1, scale_y_1);
            }
            double inv_scale_x = 1.0 / scale_x, inv_scale_y = 1.0 / scale_y,
                   inv_scale_x_1 = 1.0 / scale_x_1, inv_scale_y_1 = 1.0 / scale_y_1;
            std::pair<Vector2i, Vector2i> value = compat_position_extrinsic_index_4(
                v[k], n[k], q[k], o[k], v[kn], n[kn], q[kn], o[kn], scale_x, scale_y, inv_scale_x,
                inv_scale_y, scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1, nullptr);
            auto diff = value.first - value.second;
            index += diff;
        }
        if (index != Vector2i::Zero()) {
            new_sing[f] = rshift90(index, best[0]);
            qiyidian.push_back(i0);
            qiyidian.push_back(i1);
            qiyidian.push_back(i2);
        }
    }
    vector<int> dian{};
    for (int i = 0; i < Vset.size(); ++i) {
        for (int j = 0; j < Vset[i].size(); ++j) {
            if (find(qiyidian.begin(), qiyidian.end(), Vset[i][j]) != qiyidian.end()) {
                if (find(dian.begin(), dian.end(), i) == dian.end())
                dian.push_back(i);
            }
        }
    }

    
    for (int i = 0; i < sharp_o.size(); ++i) {
        if (boundary_o[i]) {
            sharp_o[i] = 0;
        }
    }

    F_new = F;
    V_new = V;
    N_new = N;
    Q_new = Q;
    Vset_new = Vset;
    O_compact_new = O_compact;
    F_compact_new = F_compact;
    V2E_compact_new = V2E_compact;
    E2E_compact_new = E2E_compact;
    diffs_new = diffs;
    double mScale = sqrt(surface_area / F_compact.size());
    mScale_new = mScale;
    diff_count_new = diff_count;
    sharp_o_new = sharp_o;
    compact_sharp_constraints_new = compact_sharp_constraints;
    boundary_o_new = boundary_o;
    compact_boundary_constraints_new = compact_boundary_constraints;





    /*Optimizer::extract_patch(
        hierarchy, F, V, N, Q, Vset, O_compact, F_compact, V2E_compact, E2E_compact,
        sqrt(surface_area / F_compact.size()), diffs, diff_count, patch_compact, sharp_o,
        compact_sharp_constraints, boundary_o, compact_boundary_constraints);*/
    /*  optimize_quad_positions(O_compact, N_compact, Q_compact, F_compact, V2E_compact,
       E2E_compact,V, N, Q, O, F, V2E, hierarchy.mE2E, disajoint_tree, hierarchy.mScale, false);*/
}


}  // namespace qflow

