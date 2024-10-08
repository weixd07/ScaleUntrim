





#include "optimizer.hpp"

#include <Eigen/Sparse>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <unordered_map>
#include "config.hpp"
#include "field-math.hpp"
#include "flow.hpp"
#include "parametrizer.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
using namespace std;




namespace qflow {

#ifdef WITH_CUDA
#    include <cuda_runtime.h>
#endif

#ifndef EIGEN_MPL2_ONLY
template<class T>
using LinearSolver = Eigen::SimplicialLLT<T>;
#else
template<class T>
using LinearSolver = Eigen::SparseLU<T>;
#endif

Optimizer::Optimizer() {}

void Optimizer::optimize_orientations(Hierarchy& mRes) {
#ifdef WITH_CUDA
    optimize_orientations_cuda(mRes);
    printf("%s\n", cudaGetErrorString(cudaDeviceSynchronize()));
    cudaMemcpy(mRes.mQ[0].data(), mRes.cudaQ[0], sizeof(glm::dvec3) * mRes.mQ[0].cols(),
               cudaMemcpyDeviceToHost);

#else

    int levelIterations = 6;
    for (int level = mRes.mN.size() - 1; level >= 0; --level) {
        AdjacentMatrix& adj = mRes.mAdj[level];
        const MatrixXd& N = mRes.mN[level];
        const MatrixXd& CQ = mRes.mCQ[level];
        const VectorXd& CQw = mRes.mCQw[level];
        MatrixXd& Q = mRes.mQ[level];
        auto& phases = mRes.mPhases[level];
        
        for (int iter = 0; iter < levelIterations; ++iter) {
            for (int phase = 0; phase < phases.size(); ++phase) {
                auto& p = phases[phase];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
                vector<int> inner;
                for (int pi = 0; pi < p.size(); ++pi) {
                    int i = p[pi];
                    /*if (CQw[i] == 0) {
                        inner.push_back(i);
                        continue;
                    }*/
                    const Vector3d n_i = N.col(i);
                    double weight_sum = 0.0f;
                    Vector3d sum = Q.col(i);
                    for (auto& link : adj[i]) {
                        const int j = link.id;
                        /*if (CQw[j] == 0) continue;*/
                        const double weight = link.weight;
                        if (weight == 0) continue;
                        const Vector3d n_j = N.col(j);
                        Vector3d q_j = Q.col(j);
                        std::pair<Vector3d, Vector3d> value =
                            compat_orientation_extrinsic_4(sum, n_i, q_j, n_j);
                        sum = value.first * weight_sum + value.second * weight;
                        sum -= n_i * n_i.dot(sum);
                        weight_sum += weight;
                        double norm = sum.norm();
                        if (norm > RCPOVERFLOW) sum /= norm;
                    }

                    if (CQw.size() > 0) {
                        float cw = CQw[i];
                        if (cw != 0) {
                            std::pair<Vector3d, Vector3d> value =
                                compat_orientation_extrinsic_41(sum, n_i, CQ.col(i), n_i);
                            sum = value.first * (1 - cw) + value.second * cw;
                            /*sum = CQ.col(i);*/
                            sum -= n_i * n_i.dot(sum);

                            float norm = sum.norm();
                            if (norm > RCPOVERFLOW) sum /= norm;
                        }
                    }

                    if (weight_sum > 0) {
                        Q.col(i) = sum;
                    }
                }
                for (int pi = 0; pi < inner.size(); ++pi) {
                    int i = inner[pi];
                    if (CQw[i] != 0) continue;
                    const Vector3d n_i = N.col(i);
                    double weight_sum = 0.0f;
                    Vector3d sum = Q.col(i);
                    for (auto& link : adj[i]) {
                        const int j = link.id;
                        const double weight = link.weight;
                        if (weight == 0) continue;
                        const Vector3d n_j = N.col(j);
                        Vector3d q_j = Q.col(j);
                        std::pair<Vector3d, Vector3d> value =
                            compat_orientation_extrinsic_4(sum, n_i, q_j, n_j);
                        sum = value.first * weight_sum + value.second * weight;
                        sum -= n_i * n_i.dot(sum);
                        weight_sum += weight;
                        double norm = sum.norm();
                        if (norm > RCPOVERFLOW) sum /= norm;
                    }

                    if (CQw.size() > 0) {
                        float cw = CQw[i];
                        if (cw != 0) {
                            std::pair<Vector3d, Vector3d> value =
                                compat_orientation_extrinsic_41(sum, n_i, CQ.col(i), n_i);
                            sum = value.first * (1 - cw) + value.second * cw;
                            /*sum = CQ.col(i);*/
                            sum -= n_i * n_i.dot(sum);

                            float norm = sum.norm();
                            if (norm > RCPOVERFLOW) sum /= norm;
                        }
                    }

                    if (weight_sum > 0) {
                        Q.col(i) = sum;
                    }
                }
            }
        }
        if (level > 0) {
            const MatrixXd& srcField = mRes.mQ[level];
            const MatrixXi& toUpper = mRes.mToUpper[level - 1];
            MatrixXd& destField = mRes.mQ[level - 1];
            const MatrixXd& N = mRes.mN[level - 1];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < srcField.cols(); ++i) {
                for (int k = 0; k < 2; ++k) {
                    int dest = toUpper(k, i);
                    if (dest == -1) continue;
                    Vector3d q = srcField.col(i), n = N.col(dest);
                    destField.col(dest) = q - n * n.dot(q);
                }
            }
        }
    }
    std::vector<int> switching(mRes.mV[0].cols(), 1);

    std::queue<int> queue;
    if (mRes.mstart.size() != 0) {
        queue.push(mRes.mstart[0]);
    }
    while (!queue.empty()) {
        int vertex = queue.front();
        queue.pop();
        Vector3d crossfield = mRes.mQ[0].col(vertex);
        for (auto v : mRes.mAdj[0][vertex]) {
            int v1 = v.id;
            if (switching[v1] == 0) continue;
            switching[v1] = 0;
            queue.push(v1);
            vector<Vector3d> cross(4);
            cross[0] = mRes.mQ[0].col(v1);
            Vector3d nv = mRes.mN[0].col(v1);
            cross[1] = cross[0].cross(nv);
            cross[2] = -1 * cross[0];
            cross[3] = -1 * cross[1];
            double best_score = -10;
            int best_a = 0;
            if (std::find(mRes.mboundaryvertex.begin(), mRes.mboundaryvertex.end(), v1) == mRes.mboundaryvertex.end() ||
                std::find(mRes.msharpvertex.begin(), mRes.msharpvertex.end(), v1) == mRes.msharpvertex.end()) {
                for (int i = 0; i < 4; ++i) {
                    double score = crossfield.dot(cross[i]);
                    if (score > best_score) {
                        best_a = i;
                        best_score = score;
                    }
                }
                mRes.mQ[0].col(v1) = cross[best_a];

            } else {
                for (int i = 0; i < 2; ++i) {
                    double score = crossfield.dot(cross[2 * i]);
                    if (score > best_score) {
                        best_a = 2 * i;
                        best_score = score;
                    }
                }
                mRes.mQ[0].col(v1) = cross[best_a]; 
            }
        }
    }





     

    for (int l = 0; l < mRes.mN.size() - 1; ++l) {
        const MatrixXd& N = mRes.mN[l];
        const MatrixXd& N_next = mRes.mN[l + 1];
        const MatrixXd& Q = mRes.mQ[l];
        MatrixXd& Q_next = mRes.mQ[l + 1];
        auto& toUpper = mRes.mToUpper[l];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < toUpper.cols(); ++i) {
            Vector2i upper = toUpper.col(i);
            Vector3d q0 = Q.col(upper[0]);
            Vector3d n0 = N.col(upper[0]);
            Vector3d q;

            if (upper[1] != -1) {
                Vector3d q1 = Q.col(upper[1]);
                Vector3d n1 = N.col(upper[1]);
                auto result = compat_orientation_extrinsic_4(q0, n0, q1, n1);
                q = result.first + result.second;
            } else {
                q = q0;
            }
            Vector3d n = N_next.col(i);
            q -= n.dot(q) * n;
            if (q.squaredNorm() > RCPOVERFLOW) q.normalize();

            Q_next.col(i) = q;
        }
    }




#endif
     
     
}

void Optimizer::optimize_scale(Hierarchy& mRes, VectorXd& rho, int adaptive) {
    const MatrixXd& N = mRes.mN[0];
    MatrixXd& Q = mRes.mQ[0];
    MatrixXd& V = mRes.mV[0];
    MatrixXd& S = mRes.mS[0];
    MatrixXd& K = mRes.mK[0];
    MatrixXi& F = mRes.mF;

    if (adaptive) {
        std::vector<Eigen::Triplet<double>> lhsTriplets;

        lhsTriplets.reserve(F.cols() * 6);
        for (int i = 0; i < V.cols(); ++i) {
            for (int j = 0; j < 2; ++j) {
                S(j, i) = 1.0;
                double sc1 = std::max(0.75 * S(j, i), rho[i] * 1.0 / mRes.mScale);
                S(j, i) = std::min(S(j, i), sc1);
            }
        }

        std::vector<std::map<int, double>> entries(V.cols() * 2);
        double lambda = 1;
        for (int i = 0; i < entries.size(); ++i) {
            entries[i][i] = lambda;
        }
        for (int i = 0; i < F.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int v1 = F(j, i);
                int v2 = F((j + 1) % 3, i);
                Vector3d diff = V.col(v2) - V.col(v1);
                Vector3d q_1 = Q.col(v1);
                Vector3d q_2 = Q.col(v2);
                Vector3d n_1 = N.col(v1);
                Vector3d n_2 = N.col(v2);
                Vector3d q_1_y = n_1.cross(q_1);
                auto index = compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
                int v1_x = v1 * 2, v1_y = v1 * 2 + 1, v2_x = v2 * 2, v2_y = v2 * 2 + 1;

                double dx = diff.dot(q_1);
                double dy = diff.dot(q_1_y);

                double kx_g = K(0, v1);
                double ky_g = K(1, v1);

                if (index.first % 2 != index.second % 2) {
                    std::swap(v2_x, v2_y);
                }
                double scale_x = (fmin(fmax(1 + kx_g * dy, 0.3), 3));
                double scale_y = (fmin(fmax(1 + ky_g * dx, 0.3), 3));
                //                (v2_x - scale_x * v1_x)^2 = 0
                // x^2 - 2s xy + s^2 y^2
                entries[v2_x][v2_x] += 1;
                entries[v1_x][v1_x] += scale_x * scale_x;
                entries[v2_y][v2_y] += 1;
                entries[v1_y][v1_y] += scale_y * scale_y;
                auto it = entries[v1_x].find(v2_x);
                if (it == entries[v1_x].end()) {
                    entries[v1_x][v2_x] = -scale_x;
                    entries[v2_x][v1_x] = -scale_x;
                    entries[v1_y][v2_y] = -scale_y;
                    entries[v2_y][v1_y] = -scale_y;
                } else {
                    it->second -= scale_x;
                    entries[v2_x][v1_x] -= scale_x;
                    entries[v1_y][v2_y] -= scale_y;
                    entries[v2_y][v1_y] -= scale_y;
                }
            }
        }

        Eigen::SparseMatrix<double> A(V.cols() * 2, V.cols() * 2);
        VectorXd rhs(V.cols() * 2);
        rhs.setZero();
        for (int i = 0; i < entries.size(); ++i) {
            rhs(i) = lambda * S(i % 2, i / 2);
            for (auto& rec : entries[i]) {
                lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, rec.second));
            }
        }
        A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());
        LinearSolver<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);

        solver.factorize(A);

        VectorXd result = solver.solve(rhs);

        double total_area = 0;
        for (int i = 0; i < V.cols(); ++i) {
            S(0, i) = (result(i * 2));
            S(1, i) = (result(i * 2 + 1));
            total_area += S(0, i) * S(1, i);
        }
        total_area = sqrt(V.cols() / total_area);
        for (int i = 0; i < V.cols(); ++i) {
            //            S(0, i) *= total_area;
            //            S(1, i) *= total_area;
        }
    } else {
        for (int i = 0; i < V.cols(); ++i) {
            S(0, i) = 1;
            S(1, i) = 1;
        }
    }

    for (int l = 0; l < mRes.mS.size() - 1; ++l) {
        const MatrixXd& S = mRes.mS[l];
        MatrixXd& S_next = mRes.mS[l + 1];
        auto& toUpper = mRes.mToUpper[l];
        for (int i = 0; i < toUpper.cols(); ++i) {
            Vector2i upper = toUpper.col(i);
            Vector2d q0 = S.col(upper[0]);

            if (upper[1] != -1) {
                q0 = (q0 + S.col(upper[1])) * 0.5;
            }
            S_next.col(i) = q0;
        }
    }
}

void Optimizer::optimize_positions(Hierarchy& mRes, int with_scale) {
    int levelIterations = 6;
#ifdef WITH_CUDA
    optimize_positions_cuda(mRes);
    cudaMemcpy(mRes.mO[0].data(), mRes.cudaO[0], sizeof(glm::dvec3) * mRes.mO[0].cols(),
               cudaMemcpyDeviceToHost);
#else
    vector<Vector3d> point;

    for (int level = mRes.mAdj.size() - 1; level >= 0; --level) {
        for (int iter = 0; iter < levelIterations; ++iter) {
            AdjacentMatrix& adj = mRes.mAdj[level];
            const MatrixXd &N = mRes.mN[level], &Q = mRes.mQ[level], &V = mRes.mV[level];
            const MatrixXd& CQ = mRes.mCQ[level];
            const MatrixXd& CO = mRes.mCO[level];
            const VectorXd& COw = mRes.mCOw[level];
            MatrixXd& O = mRes.mO[level];
            MatrixXd& S = mRes.mS[level];
            auto& phases = mRes.mPhases[level];
            for (int phase = 0; phase < phases.size(); ++phase) {
                auto& p = phases[phase];
#ifdef WITH_OMP
#pragma omp parallel for
#endif          
                vector<int> inner;
                for (int pi = 0; pi < p.size(); ++pi) {
                    int i = p[pi];
                    /*if (COw[i] == 0) {
                        inner.push_back(i);
                        continue;
                    }*/
                    double scale_x = mRes.mScale;
                    double scale_y = mRes.mScale;
                    if (with_scale) {
                        scale_x *= S(0, i);
                        scale_y *= S(1, i);
                    }
                    double inv_scale_x = 1.0f / scale_x;
                    double inv_scale_y = 1.0f / scale_y;
                    const Vector3d n_i = N.col(i), v_i = V.col(i);
                    Vector3d q_i = Q.col(i);

                    Vector3d sum = O.col(i);
                    double weight_sum = 0.0f;

                    q_i.normalize();
                    for (auto& link : adj[i]) {
                        const int j = link.id;
                        //if (COw[j] == 0 /*|| COw[i] != 0*/) continue;
                        const double weight = link.weight;
                        if (weight == 0) continue;
                        double scale_x_1 = mRes.mScale;
                        double scale_y_1 = mRes.mScale;
                        if (with_scale) {
                            scale_x_1 *= S(0, j);
                            scale_y_1 *= S(1, j);
                        }
                        double inv_scale_x_1 = 1.0f / scale_x_1;
                        double inv_scale_y_1 = 1.0f / scale_y_1;

                        const Vector3d n_j = N.col(j), v_j = V.col(j);
                        Vector3d q_j = Q.col(j), o_j = O.col(j);

                        q_j.normalize();

                        std::pair<Vector3d, Vector3d> value = compat_position_extrinsic_4(
                            v_i, n_i, q_i, sum, v_j, n_j, q_j, o_j, scale_x, scale_y, inv_scale_x,
                            inv_scale_y, scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1);

                        sum = value.first * weight_sum + value.second * weight;
                        weight_sum += weight;
                        if (weight_sum > RCPOVERFLOW) sum /= weight_sum;
                        sum -= n_i.dot(sum - v_i) * n_i;
                    }

                    if (COw.size() > 0) {
                        float cw = COw[i];
                        if (cw != 0) {
                            Vector3d co = CO.col(i), cq = CQ.col(i);
                            Vector3d d = co - sum;
                            d -= cq.dot(d) * cq;
                            sum += cw * d;
                            sum -= n_i.dot(sum - v_i) * n_i;
                            /* sum = co;*/
                        }
                    }

                    if (weight_sum > 0) {
                        O.col(i) = position_round_4(sum, q_i, n_i, v_i, scale_x, scale_y,
                                                    inv_scale_x, inv_scale_y);
                    }
                }

                for (int pi = 0; pi < inner.size(); ++pi) {
                    int i = inner[pi];
                    if (COw[i] != 0) {
                        continue;
                    }
                    double scale_x = mRes.mScale;
                    double scale_y = mRes.mScale;
                    if (with_scale) {
                        scale_x *= S(0, i);
                        scale_y *= S(1, i);
                    }
                    double inv_scale_x = 1.0f / scale_x;
                    double inv_scale_y = 1.0f / scale_y;
                    const Vector3d n_i = N.col(i), v_i = V.col(i);
                    Vector3d q_i = Q.col(i);

                    Vector3d sum = O.col(i);
                    double weight_sum = 0.0f;

                    q_i.normalize();
                    for (auto& link : adj[i]) {
                        const int j = link.id;
                        const double weight = link.weight;
                        if (weight == 0) continue;
                        double scale_x_1 = mRes.mScale;
                        double scale_y_1 = mRes.mScale;
                        if (with_scale) {
                            scale_x_1 *= S(0, j);
                            scale_y_1 *= S(1, j);
                        }
                        double inv_scale_x_1 = 1.0f / scale_x_1;
                        double inv_scale_y_1 = 1.0f / scale_y_1;

                        const Vector3d n_j = N.col(j), v_j = V.col(j);
                        Vector3d q_j = Q.col(j), o_j = O.col(j);

                        q_j.normalize();

                        std::pair<Vector3d, Vector3d> value = compat_position_extrinsic_4(
                            v_i, n_i, q_i, sum, v_j, n_j, q_j, o_j, scale_x, scale_y, inv_scale_x,
                            inv_scale_y, scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1);

                        sum = value.first * weight_sum + value.second * weight;
                        weight_sum += weight;
                        if (weight_sum > RCPOVERFLOW) sum /= weight_sum;
                        sum -= n_i.dot(sum - v_i) * n_i;
                    }

                    if (COw.size() > 0) {
                        float cw = COw[i];
                        if (cw != 0) {
                            Vector3d co = CO.col(i), cq = CQ.col(i);
                            Vector3d d = co - sum;
                            d -= cq.dot(d) * cq;
                            sum += cw * d;
                            sum -= n_i.dot(sum - v_i) * n_i;
                            /* sum = co;*/
                        }
                    }
                    if (weight_sum > 0) {
                        O.col(i) = position_round_4(sum, q_i, n_i, v_i, scale_x, scale_y,
                                                    inv_scale_x, inv_scale_y);
                    }
                }
            }
        }
        if (level > 0) {
            const MatrixXd& srcField = mRes.mO[level];
            const MatrixXi& toUpper = mRes.mToUpper[level - 1];
            MatrixXd& destField = mRes.mO[level - 1];
            const MatrixXd& N = mRes.mN[level - 1];
            const MatrixXd& V = mRes.mV[level - 1];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < srcField.cols(); ++i) {
                for (int k = 0; k < 2; ++k) {
                    int dest = toUpper(k, i);
                    if (dest == -1) continue;
                    Vector3d o = srcField.col(i), n = N.col(dest), v = V.col(dest);
                    o -= n * n.dot(o - v);
                    destField.col(dest) = o;
                }
            }
        }
    }
    



#endif
    }
void Optimizer::extract_patch(
    Hierarchy& mRes, MatrixXi& F, MatrixXd& V, MatrixXd& N, MatrixXd& Q,
    std::vector<std::vector<int>>& Vset, std::vector<Vector3d>& O_compact,
    std::vector<Vector4i>& F_compact, std::vector<int>& V2E_compact, std::vector<int>& E2E_compact,
    double mScale, std::vector<Vector3d>& diffs, std::vector<int>& diff_count,
    std::vector<std::vector<std::vector<int>>>& patch_compact, std::vector<int>& sharp_o,
    std::map<int, std::pair<Vector3d, Vector3d>>& compact_sharp_constraints,
    std::vector<int>& boundary_o,
    std::map<int, std::pair<Vector3d, Vector3d>>& compact_boundary_constraints) {
    auto& cornerQ = mRes.mcornerQ;

    std::vector<Vector4i> F_compact1;
    for (int i = 0; i < F_compact.size(); ++i) {
        if (F_compact[i][0] == -1 || F_compact[i][1] == -1 || F_compact[i][2] == -1 ||
            F_compact[i][3] == -1)
            continue;
        F_compact1.push_back(F_compact[i]);
    }
    F_compact = F_compact1;

    std::vector<std::vector<int>> VtoF(O_compact.size());
    std::vector<std::vector<int>> adj_compact(O_compact.size());
    for (int i = 0; i < F_compact.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j == 3) {
                int v1 = F_compact[i][3];
                VtoF[v1].push_back(i);
                if (boundary_o[v1]) continue;
                int v2 = F_compact[i][0];
                adj_compact[v1].push_back(v2);
            } else {
                int v1 = F_compact[i][j];
                VtoF[v1].push_back(i);
                if (boundary_o[v1]) continue;
                int v2 = F_compact[i][j + 1];
                adj_compact[v1].push_back(v2);
            }
        }
    }

    for (int i = 0; i < O_compact.size(); ++i) {
        if (boundary_o[i] == 0) continue;
        for (int j = 0; j < F_compact.size(); ++j) {
            if (F_compact[j][0] == i) {
                int v1 = F_compact[j][1];
                int v2 = F_compact[j][3];
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v1) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v1);
                }
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v2) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v2);
                }
            }
            if (F_compact[j][1] == i) {
                int v1 = F_compact[j][2];
                int v2 = F_compact[j][0];
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v1) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v1);
                }
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v2) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v2);
                }
            }
            if (F_compact[j][2] == i) {
                int v1 = F_compact[j][3];
                int v2 = F_compact[j][1];
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v1) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v1);
                }
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v2) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v2);
                }
            }
            if (F_compact[j][3] == i) {
                int v1 = F_compact[j][0];
                int v2 = F_compact[j][2];
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v1) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v1);
                }
                if (find(adj_compact[i].begin(), adj_compact[i].end(), v2) ==
                    adj_compact[i].end()) {
                    adj_compact[i].push_back(v2);
                }
            }
        }
    }

    std::vector<int> singularity{};
    std::vector<int> sharp_singularity{};
    std::vector<int> boundary_singularity{};
    std::vector<std::vector<int>> adj_boundary_singularity;
    std::vector<std::vector<int>> adj_sharp_singularity;
    std::vector<std::vector<int>> adj_singularity;

    for (int i = 0; i < adj_compact.size(); ++i) {
        if (boundary_o[i]) {
            if (adj_compact[i].size() > 3) {
                boundary_singularity.push_back(i);
                adj_boundary_singularity.push_back(adj_compact[i]);
            }

        } else {
            if (adj_compact[i].size() == 3 || adj_compact[i].size() > 4) {
                if (sharp_o[i]) {
                    sharp_singularity.push_back(i);
                    adj_sharp_singularity.push_back(adj_compact[i]);
                } else {
                    singularity.push_back(i);
                    adj_singularity.push_back(adj_compact[i]);
                }
            }
        }
    }

    std::vector<std::vector<int>> adj_sharp;
    adj_sharp = adj_compact;
    std::vector<int> boundary{};
    std::vector<int> sharp{};
    std::vector<int> overlap{};
    std::vector<std::vector<int>> patch_line{};
    std::vector<int> assist(O_compact.size(), 0);
    for (int i = 0; i < assist.size(); ++i) {
        if (boundary_o[i] == 1) continue;
        if (adj_compact[i].size() == 4) {
            assist[i] = 1;
        } else {
            assist[i] = adj_compact[i].size();
        }
    }

    for (int i = 0; i < singularity.size(); ++i) {
        if (assist[singularity[i]] == 0) continue;
        for (int j = 0; j < adj_singularity[i].size(); ++j) {
            if (adj_singularity[i][j] == -1) continue;
            assist[singularity[i]] -= 1;
            std::vector<int> line{};
            line.push_back(singularity[i]);
            int v = adj_singularity[i][j];
            if (std::find(singularity.begin(), singularity.end(), v) != singularity.end()) {
                auto it = std::find(singularity.begin(), singularity.end(), v);
                int position = std::distance(singularity.begin(), it);
                assist[singularity[position]] -= 1;
                auto itt = std::find(adj_singularity[position].begin(),
                                     adj_singularity[position].end(), singularity[i]);
                int position1 = std::distance(adj_singularity[position].begin(), itt);
                adj_singularity[position][position1] = -1;
                line.push_back(v);
                patch_line.push_back(line);
                continue;
            }
            if (sharp_o[v]) {
                line.push_back(v);
                patch_line.push_back(line);
                adj_singularity[i][j] = -1;
                if (find(sharp_singularity.begin(), sharp_singularity.end(), v) !=
                    sharp_singularity.end()) {
                    auto it = std::find(sharp_singularity.begin(), sharp_singularity.end(), v);
                    int position = std::distance(sharp_singularity.begin(), it);
                    auto itt = std::find(adj_sharp_singularity[position].begin(),
                                         adj_sharp_singularity[position].end(), singularity[i]);
                    int position1 = std::distance(adj_sharp_singularity[position].begin(), itt);
                    adj_sharp_singularity[position][position1] = -1;
                } else {
                    sharp.push_back(v);
                }
                for (int k = 0; k < adj_sharp[v].size(); ++k) {
                    if (adj_sharp[v][k] == -1) continue;
                    if (adj_sharp[v][k] == singularity[i]) {
                        adj_sharp[v][k] = -1;
                    }
                }
                continue;
            }
            if (boundary_o[v]) {
                boundary.push_back(v);
                line.push_back(v);
                patch_line.push_back(line);
                adj_singularity[i][j] = -1;
                if (find(boundary_singularity.begin(), boundary_singularity.end(), v) !=
                    boundary_singularity.end()) {
                    auto it =
                        std::find(boundary_singularity.begin(), boundary_singularity.end(), v);
                    int position = std::distance(boundary_singularity.begin(), it);
                    auto itt = std::find(adj_boundary_singularity[position].begin(),
                                         adj_boundary_singularity[position].end(), singularity[i]);
                    int position1 = std::distance(adj_boundary_singularity[position].begin(), itt);
                    adj_boundary_singularity[position][position1] = -1;
                }
                continue;
            }
            line.push_back(v);
            if (assist[v]) {
                assist[v] -= 1;
            } else {
                overlap.push_back(v);
            }
            int v1 = singularity[i];
            int v2;
            do {
                int next_v;
                for (int j = 0; j < VtoF[v].size(); ++j) {
                    int x = VtoF[v][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == v1 && F_compact[x][(k + 1) % 4] == v) {
                            v2 = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                for (int j = 0; j < VtoF[v].size(); ++j) {
                    int x = VtoF[v][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == v2 && F_compact[x][(k + 1) % 4] == v) {
                            next_v = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                v1 = v;
                line.push_back(next_v);
                if (sharp_o[next_v]) {
                    for (int k = 0; k < adj_sharp[next_v].size(); ++k) {
                        if (adj_sharp[next_v][k] == -1) continue;
                        if (adj_sharp[next_v][k] == v) {
                            adj_sharp[next_v][k] = -1;
                        }
                    }
                    if (find(sharp_singularity.begin(), sharp_singularity.end(), next_v) !=
                        sharp_singularity.end()) {
                        auto it = find(sharp_singularity.begin(), sharp_singularity.end(), next_v);
                        int position = std::distance(sharp_singularity.begin(), it);
                        auto itt = find(adj_sharp_singularity[position].begin(),
                                        adj_sharp_singularity[position].end(), v);
                        int position1 =
                            std::distance(adj_sharp_singularity[position].begin(), itt);
                        adj_sharp_singularity[position][position1] = -1;

                    } else {
                        sharp.push_back(next_v);
                    }
                    v = next_v;
                    continue;
                }
                if (boundary_o[next_v]) {
                    boundary.push_back(next_v);
                    if (find(boundary_singularity.begin(), boundary_singularity.end(), next_v) !=
                        boundary_singularity.end()) {
                        auto it = std::find(boundary_singularity.begin(),
                                            boundary_singularity.end(), next_v);
                        int position = std::distance(boundary_singularity.begin(), it);
                        auto itt = std::find(adj_boundary_singularity[position].begin(),
                                             adj_boundary_singularity[position].end(), v);
                        int position1 =
                            std::distance(adj_boundary_singularity[position].begin(), itt);
                        adj_boundary_singularity[position][position1] = -1;
                    }
                    v = next_v;
                    continue;
                }
                if (std::find(singularity.begin(), singularity.end(), next_v) !=
                    singularity.end()) {
                    auto ta = std::find(singularity.begin(), singularity.end(), next_v);
                    int tasposition = std::distance(singularity.begin(), ta);
                    assist[singularity[tasposition]] -= 1;
                    auto taa = std::find(adj_singularity[tasposition].begin(),
                                         adj_singularity[tasposition].end(), v);
                    int taasposition1 = std::distance(adj_singularity[tasposition].begin(), taa);
                    adj_singularity[tasposition][taasposition1] = -1;
                    v = next_v;
                    continue;
                }
                if (assist[next_v]) {
                    assist[next_v] -= 1;
                } else {
                    overlap.push_back(next_v);
                }
                v = next_v;
            } while (boundary_o[v] == 0 && sharp_o[v] == 0 &&
                     std::find(singularity.begin(), singularity.end(), v) == singularity.end());
            patch_line.push_back(line);
        }
    }
    std::vector<int> boundaryvertex;
    for (int i = 0; i < 4 * F_compact.size(); ++i) {
        if (E2E_compact[i] == -1) {
            int u0 = F_compact[i / 4][i % 4];
            int u1 = F_compact[i / 4][(i + 1) % 4];
            boundaryvertex.push_back(u0);
            boundaryvertex.push_back(u1);
        }
    }
    std::unordered_map<int, std::vector<int>> positions;
    std::vector<int> duplicates;

    for (int i = 0; i < boundaryvertex.size(); ++i) {
        positions[boundaryvertex[i]].push_back(i);
    }

    for (int i = 0; i < boundary_singularity.size(); ++i) {
        for (int j = 0; j < adj_boundary_singularity[i].size(); ++j) {
            if (adj_boundary_singularity[i][j] == -1) continue;
            if (boundary_o[adj_boundary_singularity[i][j]]) {
                int b = boundary_singularity[i];
                int vv = adj_boundary_singularity[i][j];
                int p1 = positions.at(b)[0];
                int p2 = positions.at(b)[1];
                int adj1, adj2;
                if (p1 % 2 == 0) {
                    adj1 = boundaryvertex[p1 + 1];
                } else {
                    adj1 = boundaryvertex[p1 - 1];
                }
                if (p2 % 2 == 0) {
                    adj2 = boundaryvertex[p2 + 1];
                } else {
                    adj2 = boundaryvertex[p2 - 1];
                }
                if (vv == adj1 || vv == adj2) {
                    adj_boundary_singularity[i][j] = -1;
                    continue;
                } else {
                    std::vector<int> line{};
                    line.push_back(b);
                    line.push_back(vv);
                    patch_line.push_back(line);
                    boundary.push_back(vv);
                    continue;
                }
            }
            std::vector<int> line{};
            line.push_back(boundary_singularity[i]);
            int v = adj_boundary_singularity[i][j];
            if (sharp_o[v]) {
                line.push_back(v);
                patch_line.push_back(line);
                adj_boundary_singularity[i][j] = -1;
                for (int k = 0; k < adj_sharp[v].size(); ++k) {
                    if (adj_sharp[v][k] == -1) continue;
                    if (adj_sharp[v][k] == boundary_singularity[i]) {
                        adj_sharp[v][k] = -1;
                    }
                }
                if (find(sharp_singularity.begin(), sharp_singularity.end(), v) !=
                    sharp_singularity.end()) {
                    auto it = std::find(sharp_singularity.begin(), sharp_singularity.end(), v);
                    int position = std::distance(sharp_singularity.begin(), it);
                    auto itt =
                        std::find(adj_sharp_singularity[position].begin(),
                                  adj_sharp_singularity[position].end(), boundary_singularity[i]);
                    int position1 = std::distance(adj_sharp_singularity[position].begin(), itt);
                    adj_sharp_singularity[position][position1] = -1;
                } else {
                    sharp.push_back(v);
                }
                continue;
            }
            line.push_back(v);
            if (assist[v]) {
                assist[v] -= 1;
            } else {
                overlap.push_back(v);
            }
            int v1 = boundary_singularity[i];
            int v2;
            do {
                int next_v;
                for (int j = 0; j < VtoF[v].size(); ++j) {
                    int x = VtoF[v][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == v1 && F_compact[x][(k + 1) % 4] == v) {
                            v2 = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                for (int j = 0; j < VtoF[v].size(); ++j) {
                    int x = VtoF[v][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == v2 && F_compact[x][(k + 1) % 4] == v) {
                            next_v = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                v1 = v;
                line.push_back(next_v);
                if (sharp_o[next_v]) {
                    for (int k = 0; k < adj_sharp[next_v].size(); ++k) {
                        if (adj_sharp[next_v][k] == -1) continue;
                        if (adj_sharp[next_v][k] == v) {
                            adj_sharp[next_v][k] = -1;
                        }
                    }
                    if (find(sharp_singularity.begin(), sharp_singularity.end(), next_v) !=
                        sharp_singularity.end()) {
                        auto it =
                            std::find(sharp_singularity.begin(), sharp_singularity.end(), next_v);
                        int position = std::distance(sharp_singularity.begin(), it);
                        auto itt = find(adj_sharp_singularity[position].begin(),
                                        adj_sharp_singularity[position].end(), v);
                        int position1 =
                            std::distance(adj_sharp_singularity[position].begin(), itt);
                        adj_sharp_singularity[position][position1] = -1;
                    } else {
                        sharp.push_back(next_v);
                    }
                    v = next_v;
                    continue;
                }
                if (boundary_o[next_v]) {
                    boundary.push_back(next_v);
                    if (find(boundary_singularity.begin(), boundary_singularity.end(), next_v) !=
                        boundary_singularity.end()) {
                        auto it = std::find(boundary_singularity.begin(),
                                            boundary_singularity.end(), next_v);
                        int position = std::distance(boundary_singularity.begin(), it);
                        auto itt = std::find(adj_boundary_singularity[position].begin(),
                                             adj_boundary_singularity[position].end(), v);
                        int position1 =
                            std::distance(adj_boundary_singularity[position].begin(), itt);
                        adj_boundary_singularity[position][position1] = -1;
                    }
                    v = next_v;
                    continue;
                }
                if (assist[next_v]) {
                    assist[next_v] -= 1;
                } else {
                    overlap.push_back(next_v);
                }
                v = next_v;
            } while (boundary_o[v] == 0 && sharp_o[v] == 0 &&
                     std::find(singularity.begin(), singularity.end(), v) == singularity.end());
            patch_line.push_back(line);
        }
    }
    for (int i = 0; i < sharp_singularity.size(); ++i) {
        for (int j = 0; j < adj_sharp_singularity[i].size(); ++j) {
            if (adj_sharp_singularity[i][j] == -1) continue;
            std::vector<int> line{};
            line.push_back(sharp_singularity[i]);
            int v = adj_sharp_singularity[i][j];
            adj_sharp_singularity[i][j] = -1;
            if (boundary_o[v]) {
                boundary.push_back(v);
                line.push_back(v);
                patch_line.push_back(line);
                continue;
            }
            if (sharp_o[v]) {
                for (int k = 0; k < adj_sharp[v].size(); ++k) {
                    if (adj_sharp[v][k] == -1) continue;
                    if (adj_sharp[v][k] == sharp_singularity[i]) {
                        adj_sharp[v][k] = -1;
                    }
                }
                if (find(sharp_singularity.begin(), sharp_singularity.end(), v) !=
                    sharp_singularity.end()) {
                    auto it = std::find(sharp_singularity.begin(), sharp_singularity.end(), v);
                    int position = std::distance(sharp_singularity.begin(), it);
                    auto itt = find(adj_sharp_singularity[position].begin(),
                                    adj_sharp_singularity[position].end(), sharp_singularity[i]);
                    int position1 = std::distance(adj_sharp_singularity[position].begin(), itt);
                    adj_sharp_singularity[position][position1] = -1;
                    line.push_back(v);
                    patch_line.push_back(line);
                    continue;
                }
                if (find(sharp.begin(), sharp.end(), v) != sharp.end()) {
                    line.push_back(v);
                    patch_line.push_back(line);
                    continue;
                }
            }
            if (sharp_o[v] == 0 && boundary_o[v] == 0) {
                if (assist[v] == 0) {
                    overlap.push_back(v);
                } else {
                    assist[v] -= 1;
                }
            }
            int v1 = sharp_singularity[i];
            int v2;
            int next_v = v;
            v = v1;
            bool goon = true;
            do {
                line.push_back(next_v);
                if (sharp_o[next_v] && sharp_o[v] &&
                    find(sharp.begin(), sharp.end(), next_v) == sharp.end()) {
                    for (int k = 0; k < adj_sharp[next_v].size(); ++k) {
                        if (adj_sharp[next_v][k] == -1) continue;
                        if (adj_sharp[next_v][k] == v) {
                            adj_sharp[next_v][k] = -1;
                        }
                    }
                    if (boundary_o[next_v]) {
                        goon = false;
                        continue;
                    }
                    for (int j = 0; j < VtoF[next_v].size(); ++j) {
                        int x = VtoF[next_v][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v && F_compact[x][(k + 1) % 4] == next_v) {
                                v1 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    for (int j = 0; j < VtoF[next_v].size(); ++j) {
                        int x = VtoF[next_v][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v1 && F_compact[x][(k + 1) % 4] == next_v) {
                                v2 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    v = next_v;
                    next_v = v2;
                } else {
                    if (sharp_o[next_v]) {
                        for (int k = 0; k < adj_sharp[next_v].size(); ++k) {
                            if (adj_sharp[next_v][k] == -1) continue;
                            if (adj_sharp[next_v][k] == v) {
                                adj_sharp[next_v][k] = -1;
                            }
                        }
                        if (find(sharp.begin(), sharp.end(), next_v) == sharp.end()) {
                            sharp.push_back(next_v);
                        }
                        goon = false;
                        continue;
                    }
                    if (boundary_o[next_v]) {
                        boundary.push_back(next_v);
                        goon = false;
                        continue;
                    }
                    for (int j = 0; j < VtoF[next_v].size(); ++j) {
                        int x = VtoF[next_v][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v && F_compact[x][(k + 1) % 4] == next_v) {
                                v1 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    for (int j = 0; j < VtoF[next_v].size(); ++j) {
                        int x = VtoF[next_v][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v1 && F_compact[x][(k + 1) % 4] == next_v) {
                                v2 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    v = next_v;
                    next_v = v2;
                    for (int k = 0; k < adj_sharp[v].size(); ++k) {
                        if (adj_sharp[v][k] == -1) continue;
                        if (adj_sharp[v][k] == next_v) {
                            if (assist[adj_sharp[v][k]] == 0) {
                                if (boundary_o[adj_sharp[v][k]] == 0) {
                                    overlap.push_back(adj_sharp[v][k]);
                                }
                            } else {
                                assist[adj_sharp[v][k]] -= 1;
                            }
                        }
                    }
                }
            } while (goon);
            patch_line.push_back(line);
        }
    }
    std::queue<int> queue;
    for (int i = 0; i < sharp.size(); ++i) {
        queue.push(sharp[i]);
    }

    while (!queue.empty()) {
        int sharp_v = queue.front();
        queue.pop();
        for (int i = 0; i < adj_sharp[sharp_v].size(); ++i) {
            if (adj_sharp[sharp_v][i] == -1) continue;
            int v = sharp_v;
            std::vector<int> line{};
            line.push_back(v);
            int v_next = adj_sharp[v][i];
            if (sharp_o[v_next] == 0 && boundary_o[v_next] == 0) {
                if (assist[v_next] == 0) {
                    overlap.push_back(v_next);
                } else {
                    assist[v_next] -= 1;
                }
            }
            bool loop_is = true;
            while (loop_is) {
                if (sharp_o[v_next] && sharp_o[v] &&
                    find(sharp.begin(), sharp.end(), v_next) == sharp.end()) {
                    line.push_back(v_next);
                    for (int k = 0; k < adj_sharp[v_next].size(); ++k) {
                        if (adj_sharp[v_next][k] == -1) continue;
                        if (adj_sharp[v_next][k] == v) {
                            adj_sharp[v_next][k] = -1;
                        }
                    }
                    int v2, v3;
                    for (int j = 0; j < VtoF[v_next].size(); ++j) {
                        int x = VtoF[v_next][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v && F_compact[x][(k + 1) % 4] == v_next) {
                                v2 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    for (int j = 0; j < VtoF[v_next].size(); ++j) {
                        int x = VtoF[v_next][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v2 && F_compact[x][(k + 1) % 4] == v_next) {
                                v3 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    v = v_next;
                    v_next = v3;
                } else {
                    if (sharp_o[v_next]) {
                        line.push_back(v_next);
                        if (find(sharp.begin(), sharp.end(), v_next) == sharp.end()) {
                            sharp.push_back(v_next);
                            queue.push(v_next);
                        }
                        for (int k = 0; k < adj_sharp[v_next].size(); ++k) {
                            if (adj_sharp[v_next][k] == -1) continue;
                            if (adj_sharp[v_next][k] == v) {
                                adj_sharp[v_next][k] = -1;
                            }
                        }
                        loop_is = false;
                        continue;
                    }
                    if (boundary_o[v_next]) {
                        line.push_back(v_next);
                        boundary.push_back(v_next);
                        if (find(boundary_singularity.begin(), boundary_singularity.end(),
                                 v_next) != boundary_singularity.end()) {
                            auto it = std::find(boundary_singularity.begin(),
                                                boundary_singularity.end(), v_next);
                            int position = std::distance(boundary_singularity.begin(), it);
                            auto itt = std::find(adj_boundary_singularity[position].begin(),
                                                 adj_boundary_singularity[position].end(), v);
                            int position1 =
                                std::distance(adj_boundary_singularity[position].begin(), itt);
                            adj_boundary_singularity[position][position1] = -1;
                        }
                        loop_is = false;
                        continue;
                    }
                    line.push_back(v_next);
                    int v2, v3;
                    for (int j = 0; j < VtoF[v_next].size(); ++j) {
                        int x = VtoF[v_next][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v && F_compact[x][(k + 1) % 4] == v_next) {
                                v2 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    for (int j = 0; j < VtoF[v_next].size(); ++j) {
                        int x = VtoF[v_next][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v2 && F_compact[x][(k + 1) % 4] == v_next) {
                                v3 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    for (int k = 0; k < adj_sharp[v_next].size(); ++k) {
                        if (adj_sharp[v_next][k] == -1) continue;
                        if (adj_sharp[v_next][k] == v3) {
                            if (assist[adj_sharp[v_next][k]] == 0) {
                                if (boundary_o[adj_sharp[v_next][k]] == 0) {
                                    overlap.push_back(adj_sharp[v_next][k]);
                                }
                            } else {
                                assist[adj_sharp[v_next][k]] -= 1;
                            }
                        }
                    }
                    v = v_next;
                    v_next = v3;
                }
            }
            patch_line.push_back(line);
        }
    }

    for (int i = 0; i < boundary_singularity.size(); ++i) {
        if (find(boundary.begin(), boundary.end(), boundary_singularity[i]) != boundary.end())
            continue;
        boundary.push_back(boundary_singularity[i]);
    }

    std::vector<std::vector<int>> patch_line_inner;
    for (int i = 0; i < patch_line.size(); ++i) {
        std::vector<int> line;
        for (int j = 0; j < patch_line[i].size(); ++j) {
            if (j == patch_line[i].size() - 1) {
                int v = patch_line[i][j];
                line.push_back(v);
                patch_line_inner.push_back(line);
            } else {
                int v = patch_line[i][j];
                line.push_back(v);
                if (std::find(overlap.begin(), overlap.end(), v) != overlap.end()) {
                    patch_line_inner.push_back(line);
                    std::vector<int>().swap(line);
                    line.push_back(v);
                }
            }
        }
    }
    std::vector<std::vector<int>> patch_line_boundary;

    std::vector<int> corner;
    std::vector<int> number_F(O_compact.size());
    for (int i = 0; i < F_compact.size(); ++i) {
        number_F[F_compact[i][0]] += 1;
        number_F[F_compact[i][1]] += 1;
        number_F[F_compact[i][2]] += 1;
        number_F[F_compact[i][3]] += 1;
    }
    for (int i = 0; i < number_F.size(); ++i) {
        if (number_F[i] == 1) {
            corner.push_back(i);
        }
    }

    std::vector<int> db;
    for (int i = 0; i < boundaryvertex.size(); ++i) {
        if (find(db.begin(), db.end(), boundaryvertex[i]) != db.end()) continue;
        db.push_back(boundaryvertex[i]);
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
        std::vector<int> line{};
        std::vector<int> temp;
        std::vector<int> temp1;
        int a = 0;
        for (int j = 0; j < corner.size(); ++j) {
            if (find(boundary.begin(), boundary.end(), corner[j]) == boundary.end()) {
                boundary.push_back(corner[j]);
            }
        }
        for (int j = 0; j < loop.size(); ++j) {
            int v = loop[j];
            if (j == 0 && find(boundary.begin(), boundary.end(), v) != boundary.end()) {
                temp = loop;
                break;
            }
            if (find(boundary.begin(), boundary.end(), v) != boundary.end()) {
                a = j;
                temp1.push_back(v);
                break;
            }
            temp1.push_back(v);
        }
        for (int j = 0; j < loop.size(); ++j) {
            if (temp.size() == loop.size()) break;
            if (j > a) {
                temp.push_back(loop[j]);
            }
        }
        for (int j = 0; j < temp1.size(); ++j) {
            temp.push_back(temp1[j]);
        }

        for (int j = 0; j < temp.size(); ++j) {
            if (j % 2 == 1 && j != temp.size() - 1) continue;
            if (j == temp.size() - 1) {
                int v = temp[j];
                line.push_back(v);
                patch_line_boundary.push_back(line);
            } else {
                int v = temp[j];
                line.push_back(v);
                if (find(boundary.begin(), boundary.end(), v) != boundary.end() &&
                    line.size() != 1) {
                    patch_line_boundary.push_back(line);
                    std::vector<int>().swap(line);
                    line.push_back(v);
                }
            }
        }
    }
    std::vector<int> start;
    std::vector<std::vector<std::vector<int>>> patch;
    std::vector<int> number_line_inner(patch_line_inner.size(), 2);
    std::vector<int> number_line_boundary(patch_line_boundary.size(), 1);

    for (int i = 0; i < singularity.size(); ++i) {
        start.push_back(singularity[i]);
    }
    for (int i = 0; i < boundary_singularity.size(); ++i) {
        start.push_back(boundary_singularity[i]);
    }
    for (int i = 0; i < sharp_singularity.size(); ++i) {
        start.push_back(sharp_singularity[i]);
    }
    for (int i = 0; i < overlap.size(); ++i) {
        start.push_back(overlap[i]);
    }
    for (int i = 0; i < sharp.size(); ++i) {
        start.push_back(sharp[i]);
    }

    if (singularity.size() == 0 && boundary_singularity.size() == 0 &&
        sharp_singularity.size() == 0) {
        int v = patch_line_boundary[0][0];
        Vector4i one_quad;
        for (int i = 0; i < F_compact.size(); ++i) {
            if (F_compact[i][0] == v || F_compact[i][1] == v || F_compact[i][2] == v ||
                F_compact[i][3] == v) {
                one_quad = F_compact[i];
            }
        }
        int assist_point, start_point;
        if (one_quad[0] == v) {
            start_point = one_quad[1];
            assist_point = one_quad[2];
        }
        if (one_quad[1] == v) {
            start_point = one_quad[2];
            assist_point = one_quad[3];
        }
        if (one_quad[2] == v) {
            start_point = one_quad[3];
            assist_point = one_quad[0];
        }
        if (one_quad[3] == v) {
            start_point = one_quad[0];
            assist_point = one_quad[1];
        }
        std::vector<int> assist_line;
        do {
            assist_line.push_back(assist_point);
            int v2, v3;
            for (int j = 0; j < VtoF[assist_point].size(); ++j) {
                int x = VtoF[assist_point][j];
                for (int k = 0; k < F_compact[x].size(); ++k) {
                    if (F_compact[x][k] == start_point &&
                        F_compact[x][(k + 1) % 4] == assist_point) {
                        v2 = F_compact[x][(k + 2) % 4];
                    }
                }
            }
            for (int j = 0; j < VtoF[assist_point].size(); ++j) {
                int x = VtoF[assist_point][j];
                for (int k = 0; k < F_compact[x].size(); ++k) {
                    if (F_compact[x][k] == v2 && F_compact[x][(k + 1) % 4] == assist_point) {
                        v3 = F_compact[x][(k + 2) % 4];
                    }
                }
            }
            start_point = assist_point;
            assist_point = v3;
        } while (boundary_o[assist_point] == 0);
        std::vector<std::vector<int>> matrix;
        matrix.push_back(patch_line_boundary[0]);
        std::vector<int> counter_line;
        for (int i = 0; i < patch_line_boundary[3].size() - 1; ++i) {
            if (i == 0) continue;
            counter_line.push_back(patch_line_boundary[3][patch_line_boundary[3].size() - 1 - i]);
        }
        for (int i = 0; i < counter_line.size(); ++i) {
            std::vector<int> patchline;
            int vertex;
            vertex = assist_line[i];
            patchline.push_back(counter_line[i]);
            patchline.push_back(vertex);
            int v1, v2, v3;
            v1 = counter_line[i];
            do {
                for (int j = 0; j < VtoF[vertex].size(); ++j) {
                    int x = VtoF[vertex][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == v1 && F_compact[x][(k + 1) % 4] == vertex) {
                            v2 = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                for (int j = 0; j < VtoF[vertex].size(); ++j) {
                    int x = VtoF[vertex][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == v2 && F_compact[x][(k + 1) % 4] == vertex) {
                            v3 = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                v1 = vertex;
                vertex = v3;
                patchline.push_back(vertex);
            } while (boundary_o[vertex] == 0);
            matrix.push_back(patchline);
        }
        std::vector<int> asline;
        for (int i = 0; i < patch_line_boundary[2].size(); ++i) {
            asline.push_back(patch_line_boundary[2][patch_line_boundary[2].size() - 1 - i]);
        }
        matrix.push_back(asline);
        patch.push_back(matrix);
        patch_compact = patch;

    } else {
        std::vector<std::vector<Vector4i>> F_inner(start.size());
        std::vector<Vector4i> Inner_quad;
        std::vector<int> point;
        for (int i = 0; i < start.size(); ++i) {
            int v = start[i];
            for (int j = 0; j < F_compact.size(); ++j) {
                if (F_compact[j][0] == v || F_compact[j][1] == v || F_compact[j][2] == v ||
                    F_compact[j][3] == v) {
                    F_inner[i].push_back(F_compact[j]);
                }
            }
        }
        for (int i = 0; i < start.size(); ++i) {
            int v = start[i];

            for (int j = 0; j < F_inner[i].size(); ++j) {
                int v_adj1, v_adj2;
                int end1, end2;
                if (F_inner[i][j][0] == v) {
                    v_adj1 = F_inner[i][j][3];
                    v_adj2 = F_inner[i][j][1];
                }
                if (F_inner[i][j][1] == v) {
                    v_adj1 = F_inner[i][j][0];
                    v_adj2 = F_inner[i][j][2];
                }
                if (F_inner[i][j][2] == v) {
                    v_adj1 = F_inner[i][j][1];
                    v_adj2 = F_inner[i][j][3];
                }
                if (F_inner[i][j][3] == v) {
                    v_adj1 = F_inner[i][j][2];
                    v_adj2 = F_inner[i][j][0];
                }
                std::vector<std::vector<int>> single_patch;
                int a1, a2;
                for (int k = 0; k < patch_line_inner.size(); ++k) {
                    if (number_line_inner[k] == 0) continue;
                    if (patch_line_inner[k][0] == v && patch_line_inner[k][1] == v_adj1) {
                        single_patch.push_back(patch_line_inner[k]);
                        a1 = k;
                        end1 = patch_line_inner[k][patch_line_inner[k].size() - 1];
                    }
                    if (patch_line_inner[k][patch_line_inner[k].size() - 1] == v &&
                        patch_line_inner[k][patch_line_inner[k].size() - 2] == v_adj1) {
                        single_patch.push_back(patch_line_inner[k]);
                        a1 = k;
                        end1 = patch_line_inner[k][0];
                    }
                    if (patch_line_inner[k][0] == v && patch_line_inner[k][1] == v_adj2) {
                        single_patch.push_back(patch_line_inner[k]);
                        a2 = k;
                        end2 = patch_line_inner[k][patch_line_inner[k].size() - 1];
                    }
                    if (patch_line_inner[k][patch_line_inner[k].size() - 1] == v &&
                        patch_line_inner[k][patch_line_inner[k].size() - 2] == v_adj2) {
                        single_patch.push_back(patch_line_inner[k]);
                        a2 = k;
                        end2 = patch_line_inner[k][0];
                    }
                }
                if (single_patch.size() != 2) continue;
                number_line_inner[a1] -= 1;
                number_line_inner[a2] -= 1;
                if (boundary_o[end1] && boundary_o[end2]) {
                    for (int k = 0; k < patch_line_boundary.size(); ++k) {
                        if (number_line_boundary[k] == 0) continue;
                        int beginning = -1;
                        bool loop = true;
                        if (patch_line_boundary[k][0] == end2 &&
                            patch_line_boundary[k][patch_line_boundary[k].size() - 1] != v) {
                            beginning = patch_line_boundary[k][patch_line_boundary[k].size() - 1];
                        }
                        if (patch_line_boundary[k][patch_line_boundary[k].size() - 1] == end2 &&
                            patch_line_boundary[k][0] != v) {
                            beginning = patch_line_boundary[k][0];
                        }
                        if (beginning == -1) continue;
                        for (int m = 0; m < patch_line_boundary.size(); ++m) {
                            if (number_line_boundary[m] == 0) continue;
                            if (patch_line_boundary[m][0] == end1 &&
                                patch_line_boundary[m][patch_line_boundary[m].size() - 1] ==
                                    beginning) {
                                single_patch.push_back(patch_line_boundary[k]);
                                single_patch.push_back(patch_line_boundary[m]);
                                number_line_boundary[k] -= 1;
                                number_line_boundary[m] -= 1;
                                loop = false;
                                break;
                            }
                            if (patch_line_boundary[m][0] == beginning &&
                                patch_line_boundary[m][patch_line_boundary[m].size() - 1] ==
                                    end1) {
                                single_patch.push_back(patch_line_boundary[k]);
                                single_patch.push_back(patch_line_boundary[m]);
                                number_line_boundary[k] -= 1;
                                number_line_boundary[m] -= 1;
                                loop = false;
                                break;
                            }
                        }
                        if (loop) continue;
                        break;
                    }
                    if (single_patch.size() != 4) {
                        number_line_inner[a1] += 1;
                        number_line_inner[a2] += 1;
                        continue;
                    }
                }
                if (boundary_o[end1] == 0 && boundary_o[end2]) {
                    for (int k = 0; k < patch_line_boundary.size(); ++k) {
                        if (number_line_boundary[k] == 0) continue;
                        int beginning = -1;
                        bool loop = true;
                        if (patch_line_boundary[k][0] == end2 &&
                            patch_line_boundary[k][patch_line_boundary[k].size() - 1] != v) {
                            beginning = patch_line_boundary[k][patch_line_boundary[k].size() - 1];
                        }
                        if (patch_line_boundary[k][patch_line_boundary[k].size() - 1] == end2 &&
                            patch_line_boundary[k][0] != v) {
                            beginning = patch_line_boundary[k][0];
                        }
                        if (beginning == -1) continue;
                        for (int m = 0; m < patch_line_inner.size(); ++m) {
                            if (number_line_inner[m] == 0) continue;
                            if (patch_line_inner[m][0] == end1 &&
                                patch_line_inner[m][patch_line_inner[m].size() - 1] == beginning) {
                                single_patch.push_back(patch_line_boundary[k]);
                                single_patch.push_back(patch_line_inner[m]);
                                number_line_boundary[k] -= 1;
                                number_line_inner[m] -= 1;
                                loop = false;
                                break;
                            }
                            if (patch_line_inner[m][0] == beginning &&
                                patch_line_inner[m][patch_line_inner[m].size() - 1] == end1) {
                                single_patch.push_back(patch_line_boundary[k]);
                                single_patch.push_back(patch_line_inner[m]);
                                number_line_boundary[k] -= 1;
                                number_line_inner[m] -= 1;
                                loop = false;
                                break;
                            }
                        }
                        if (loop) continue;
                        break;
                    }
                    if (single_patch.size() != 4) {
                        number_line_inner[a1] += 1;
                        number_line_inner[a2] += 1;
                        continue;
                    }
                }
                if (boundary_o[end1] && boundary_o[end2] == 0) {
                    for (int k = 0; k < patch_line_inner.size(); ++k) {
                        if (number_line_inner[k] == 0) continue;
                        int beginning = -1;
                        bool loop = true;
                        if (patch_line_inner[k][0] == end2 &&
                            patch_line_inner[k][patch_line_inner[k].size() - 1] != v) {
                            beginning = patch_line_inner[k][patch_line_inner[k].size() - 1];
                        }
                        if (patch_line_inner[k][patch_line_inner[k].size() - 1] == end2 &&
                            patch_line_inner[k][0] != v) {
                            beginning = patch_line_inner[k][0];
                        }
                        if (beginning == -1) continue;
                        for (int m = 0; m < patch_line_boundary.size(); ++m) {
                            if (number_line_boundary[m] == 0) continue;
                            if (patch_line_boundary[m][0] == end1 &&
                                patch_line_boundary[m][patch_line_boundary[m].size() - 1] ==
                                    beginning) {
                                single_patch.push_back(patch_line_inner[k]);
                                single_patch.push_back(patch_line_boundary[m]);
                                number_line_inner[k] -= 1;
                                number_line_boundary[m] -= 1;
                                loop = false;
                                break;
                            }
                            if (patch_line_boundary[m][0] == beginning &&
                                patch_line_boundary[m][patch_line_boundary[m].size() - 1] ==
                                    end1) {
                                single_patch.push_back(patch_line_inner[k]);
                                single_patch.push_back(patch_line_boundary[m]);
                                number_line_inner[k] -= 1;
                                number_line_boundary[m] -= 1;
                                loop = false;
                                break;
                            }
                        }
                        if (loop) continue;
                        break;
                    }
                    if (single_patch.size() != 4) {
                        number_line_inner[a1] += 1;
                        number_line_inner[a2] += 1;
                        continue;
                    }
                }
                if (boundary_o[end1] == 0 && boundary_o[end2] == 0) {
                    for (int k = 0; k < patch_line_inner.size(); ++k) {
                        if (number_line_inner[k] == 0) continue;
                        int beginning = -1;
                        bool loop = true;
                        if (patch_line_inner[k][0] == end2 &&
                            patch_line_inner[k][patch_line_inner[k].size() - 1] != v) {
                            beginning = patch_line_inner[k][patch_line_inner[k].size() - 1];
                        }
                        if (patch_line_inner[k][patch_line_inner[k].size() - 1] == end2 &&
                            patch_line_inner[k][0] != v) {
                            beginning = patch_line_inner[k][0];
                        }
                        if (beginning == -1) continue;
                        for (int m = 0; m < patch_line_inner.size(); ++m) {
                            if (number_line_inner[m] == 0) continue;
                            if (patch_line_inner[m][0] == end1 &&
                                patch_line_inner[m][patch_line_inner[m].size() - 1] == beginning) {
                                single_patch.push_back(patch_line_inner[k]);
                                single_patch.push_back(patch_line_inner[m]);
                                number_line_inner[k] -= 1;
                                number_line_inner[m] -= 1;
                                loop = false;
                                break;
                            }
                            if (patch_line_inner[m][0] == beginning &&
                                patch_line_inner[m][patch_line_inner[m].size() - 1] == end1) {
                                single_patch.push_back(patch_line_inner[k]);
                                single_patch.push_back(patch_line_inner[m]);
                                number_line_inner[k] -= 1;
                                number_line_inner[m] -= 1;
                                loop = false;
                                break;
                            }
                        }
                        if (loop) continue;
                        break;
                    }
                    if (single_patch.size() != 4) {
                        number_line_inner[a1] += 1;
                        number_line_inner[a2] += 1;
                        continue;
                    }
                }
                patch.push_back(single_patch);
                Inner_quad.push_back(F_inner[i][j]);
                point.push_back(v);
            }
        }

        int need = 0;
        for (int i = 0; i < number_line_inner.size(); ++i) {
            if (number_line_inner[i] != 0) {
                need = 1;
            }
        }
        std::vector<std::vector<int>> patchline;
        std::vector<int> number_patchline;
        for (int i = 0; i < number_line_inner.size(); ++i) {
            if (number_line_inner[i] != 0) {
                for (int j = 0; j < number_line_inner[i]; ++j) {
                    patchline.push_back(patch_line_inner[i]);
                    number_patchline.push_back(1);
                }
            }
        }
        for (int i = 0; i < number_line_boundary.size(); ++i) {
            if (number_line_boundary[i] != 0) {
                for (int j = 0; j < number_line_boundary[i]; ++j) {
                    patchline.push_back(patch_line_boundary[i]);
                    number_patchline.push_back(1);
                }
            }
        }
        if (need) {
            for (int i = 0; i < patchline.size(); ++i) {
                if (number_patchline[i] == 0) continue;
                std::vector<std::vector<int>> single_patch;
                int a, b, c, d;
                int v1, v2, v3, v4;
                v1 = patchline[i][0];
                v2 = patchline[i][patchline[i].size() - 1];
                single_patch.push_back(patchline[i]);
                for (int j = 0; j < patchline.size(); ++j) {
                    if (j == i) continue;
                    if (number_patchline[j] == 0) continue;
                    if (patchline[j][0] == v2 && patchline[j][patchline[j].size() - 1] != v1) {
                        v3 = patchline[j][patchline[j].size() - 1];
                        single_patch.push_back(patchline[j]);
                        a = j;
                    }
                    if (patchline[j][patchline[j].size() - 1] == v2 && patchline[j][0] != v1) {
                        v3 = patchline[j][0];
                        std::vector<int> line;
                        for (int k = 0; k < patchline[j].size(); ++k) {
                            line.push_back(patchline[j][patchline[j].size() - 1 - k]);
                        }
                        single_patch.push_back(line);
                        a = j;
                    }
                }
                for (int j = 0; j < patchline.size(); ++j) {
                    if (j == i) continue;
                    if (j == a) continue;
                    if (number_patchline[j] == 0) continue;
                    if (patchline[j][0] == v3 && patchline[j][patchline[j].size() - 1] != v2) {
                        v4 = patchline[j][patchline[j].size() - 1];
                        single_patch.push_back(patchline[j]);
                        b = j;
                    }
                    if (patchline[j][patchline[j].size() - 1] == v3 && patchline[j][0] != v2) {
                        v4 = patchline[j][0];
                        std::vector<int> line;
                        for (int k = 0; k < patchline[j].size(); ++k) {
                            line.push_back(patchline[j][patchline[j].size() - 1 - k]);
                        }
                        single_patch.push_back(line);
                        b = j;
                    }
                }
                for (int j = 0; j < patchline.size(); ++j) {
                    if (j == i) continue;
                    if (j == a) continue;
                    if (j == b) continue;
                    if (number_patchline[j] == 0) continue;
                    if (patchline[j][0] == v4 && patchline[j][patchline[j].size() - 1] != v3) {
                        d = patchline[j][patchline[j].size() - 1];
                        single_patch.push_back(patchline[j]);
                        c = j;
                    }
                    if (patchline[j][patchline[j].size() - 1] == v4 && patchline[j][0] != v3) {
                        d = patchline[j][0];
                        std::vector<int> line;
                        for (int k = 0; k < patchline[j].size(); ++k) {
                            line.push_back(patchline[j][patchline[j].size() - 1 - k]);
                        }
                        single_patch.push_back(line);
                        c = j;
                    }
                }
                if (d == v1 && single_patch.size() == 4) {
                    int p, adj1, adj2, r;
                    /*Vector4i Qu;*/
                    if (find(start.begin(), start.end(), v1) != start.end()) {
                        r = 0;
                        p = v1;
                        adj1 = single_patch[0][1];
                        adj2 = single_patch[3][single_patch[3].size() - 2];
                    }
                    if (find(start.begin(), start.end(), v2) != start.end()) {
                        r = 1;
                        p = v2;
                        adj1 = single_patch[1][1];
                        adj2 = single_patch[0][single_patch[0].size() - 2];
                    }
                    if (find(start.begin(), start.end(), v3) != start.end()) {
                        r = 2;
                        p = v3;
                        adj1 = single_patch[2][1];
                        adj2 = single_patch[1][single_patch[1].size() - 2];
                    }
                    if (find(start.begin(), start.end(), v4) != start.end()) {
                        r = 3;
                        p = v4;
                        adj1 = single_patch[3][1];
                        adj2 = single_patch[2][single_patch[2].size() - 2];
                    }
                    for (int j = 0; j < start.size(); ++j) {
                        if (start[j] != p) continue;
                        for (int k = 0; k < F_inner[j].size(); ++k) {
                            std::vector<int> F;
                            F.push_back(F_inner[j][k][0]);
                            F.push_back(F_inner[j][k][1]);
                            F.push_back(F_inner[j][k][2]);
                            F.push_back(F_inner[j][k][3]);
                            if (find(F.begin(), F.end(), adj1) != F.end() &&
                                find(F.begin(), F.end(), adj2) != F.end() &&
                                find(F.begin(), F.end(), p) != F.end()) {
                                Inner_quad.push_back(F_inner[j][k]);
                                /*Qu = F_inner[j][k];*/
                            }
                        }
                    }
                    std::vector<std::vector<int>> Qkx;
                    if (r == 0) {
                        Qkx.push_back(single_patch[0]);
                        std::vector<int> line;
                        for (int m = 0; m < single_patch[3].size(); ++m) {
                            line.push_back(single_patch[3][single_patch[3].size() - 1 - m]);
                        }
                        Qkx.push_back(line);
                        Qkx.push_back(single_patch[1]);
                        Qkx.push_back(single_patch[2]);
                    }
                    if (r == 1) {
                        Qkx.push_back(single_patch[1]);
                        std::vector<int> line;
                        for (int m = 0; m < single_patch[0].size(); ++m) {
                            line.push_back(single_patch[0][single_patch[0].size() - 1 - m]);
                        }
                        Qkx.push_back(line);
                        Qkx.push_back(single_patch[2]);
                        Qkx.push_back(single_patch[3]);
                    }
                    if (r == 2) {
                        Qkx.push_back(single_patch[2]);
                        std::vector<int> line;
                        for (int m = 0; m < single_patch[1].size(); ++m) {
                            line.push_back(single_patch[1][single_patch[1].size() - 1 - m]);
                        }
                        Qkx.push_back(line);
                        Qkx.push_back(single_patch[3]);
                        Qkx.push_back(single_patch[0]);
                    }
                    if (r == 3) {
                        Qkx.push_back(single_patch[3]);
                        std::vector<int> line;
                        for (int m = 0; m < single_patch[2].size(); ++m) {
                            line.push_back(single_patch[2][single_patch[2].size() - 1 - m]);
                        }
                        Qkx.push_back(line);
                        Qkx.push_back(single_patch[0]);
                        Qkx.push_back(single_patch[1]);
                    }
                    patch.push_back(Qkx);
                    point.push_back(p);
                    number_patchline[i] = 0;
                    number_patchline[a] = 0;
                    number_patchline[b] = 0;
                    number_patchline[c] = 0;
                }
            }
        }
        std::vector<std::vector<std::vector<int>>> collapse_patch(patch.size());
        std::vector<int> need_collapse;
        for (int i = 0; i < patch.size(); ++i) {
            int v = point[i];
            int v_adj1, ender, starter;
            if (Inner_quad[i][0] == v) {
                v_adj1 = Inner_quad[i][1];
            }
            if (Inner_quad[i][1] == v) {
                v_adj1 = Inner_quad[i][2];
            }
            if (Inner_quad[i][2] == v) {
                v_adj1 = Inner_quad[i][3];
            }
            if (Inner_quad[i][3] == v) {
                v_adj1 = Inner_quad[i][0];
            }
            starter = v;
            for (int j = 0; j < patch[i].size(); ++j) {
                if (patch[i][j][0] == v && patch[i][j][1] == v_adj1) {
                    ender = patch[i][j][patch[i][j].size() - 1];
                    collapse_patch[i].push_back(patch[i][j]);
                }
                if (patch[i][j][patch[i][j].size() - 1] == v &&
                    patch[i][j][patch[i][j].size() - 2] == v_adj1) {
                    ender = patch[i][j][0];
                    std::vector<int> line;
                    for (int k = 0; k < patch[i][j].size(); ++k) {
                        line.push_back(patch[i][j][patch[i][j].size() - 1 - k]);
                    }
                    collapse_patch[i].push_back(line);
                }
            }
            do {
                for (int j = 0; j < patch[i].size(); ++j) {
                    if (patch[i][j][0] == ender &&
                        patch[i][j][patch[i][j].size() - 1] != starter) {
                        starter = ender;
                        ender = patch[i][j][patch[i][j].size() - 1];
                        collapse_patch[i].push_back(patch[i][j]);
                        if (ender == v) break;
                    }
                    if (patch[i][j][patch[i][j].size() - 1] == ender &&
                        patch[i][j][0] != starter) {
                        starter = ender;
                        ender = patch[i][j][0];
                        std::vector<int> line;
                        for (int k = 0; k < patch[i][j].size(); ++k) {
                            line.push_back(patch[i][j][patch[i][j].size() - 1 - k]);
                        }
                        collapse_patch[i].push_back(line);
                        if (ender == v) break;
                    }
                }
            } while (ender != v);
        }

        std::vector<int> all_singularity;
        for (int i = 0; i < singularity.size(); ++i) {
            all_singularity.push_back(singularity[i]);
        }
        for (int i = 0; i < boundary_singularity.size(); ++i) {
            all_singularity.push_back(boundary_singularity[i]);
        }
        for (int i = 0; i < sharp_singularity.size(); ++i) {
            all_singularity.push_back(sharp_singularity[i]);
        }
        std::vector<std::vector<int>> collapse_patch_edge;
        for (int i = 0; i < collapse_patch.size(); ++i) {
            std::vector<int> line;
            line.push_back(1);
            line.push_back(1);
            line.push_back(1);
            line.push_back(1);
            collapse_patch_edge.push_back(line);
        }
        std::vector<std::vector<int>> layer;

        std::vector<int> thr_for_sin_patch;
        for (int i = 0; i < collapse_patch.size(); ++i) {
            int a = 0;
            if (collapse_patch[i][0].size() == 2 && collapse_patch[i][1].size() == 2) {
                for (int k = 0; k < 4; ++k) {
                    if (find(all_singularity.begin(), all_singularity.end(),
                             collapse_patch[i][k][0]) != all_singularity.end()) {
                        a += 1;
                    }
                }
                if (a > 2) {
                    thr_for_sin_patch.push_back(i);
                }
            }
        }
        for (int i = 0; i < collapse_patch.size(); ++i) {
            std::vector<int> one_layer;
            int rev = -1;
            for (int j = 0; j < collapse_patch[i].size(); ++j) {
                if (collapse_patch_edge[i][j] == -1) continue;
                if (collapse_patch[i][j].size() > 2) {
                    collapse_patch_edge[i][j] = -1;
                    continue;
                }
                if (one_layer.size() > 1 && rev != collapse_patch[i].size() - 1 && rev != -1) {
                    if (j == rev + 1) continue;
                }
                rev = j;
                collapse_patch_edge[i][j] = -1;
                if (find(thr_for_sin_patch.begin(), thr_for_sin_patch.end(), i) !=
                    thr_for_sin_patch.end()) {
                    continue;
                }
                if (find(one_layer.begin(), one_layer.end(), i) == one_layer.end()) {
                    one_layer.push_back(i);
                }

                int v1, v2;
                v1 = collapse_patch[i][j][0];
                v2 = collapse_patch[i][j][1];
                bool loop = true;
                do {
                    int v3 = -1;
                    int v4 = -1;
                    for (int k = 0; k < collapse_patch.size(); ++k) {
                        if (k == i) continue;
                        if (find(thr_for_sin_patch.begin(), thr_for_sin_patch.end(), k) !=
                            thr_for_sin_patch.end()) {
                            continue;
                        }
                        for (int p = 0; p < collapse_patch[k].size(); ++p) {
                            if (collapse_patch_edge[k][p] == -1) continue;
                            if (collapse_patch[k][p].size() > 2) continue;
                            if (collapse_patch[k][p][0] == v1 && collapse_patch[k][p][1] == v2) {
                                collapse_patch_edge[k][p] = -1;
                                collapse_patch_edge[k][(p + 2) % 4] = -1;
                                v3 = collapse_patch[k][(p + 2) % 4][0];
                                v4 = collapse_patch[k][(p + 2) % 4][1];
                                v1 = v3;
                                v2 = v4;
                                if (j == 0 || j == 1) {
                                    one_layer.insert(one_layer.begin(), k);
                                } else {
                                    one_layer.push_back(k);
                                }

                            } else if (collapse_patch[k][p][0] == v2 &&
                                       collapse_patch[k][p][1] == v1) {
                                collapse_patch_edge[k][p] = -1;
                                collapse_patch_edge[k][(p + 2) % 4] = -1;
                                v3 = collapse_patch[k][(p + 2) % 4][0];
                                v4 = collapse_patch[k][(p + 2) % 4][1];
                                v1 = v3;
                                v2 = v4;
                                if (j == 0 || j == 1) {
                                    one_layer.insert(one_layer.begin(), k);
                                } else {
                                    one_layer.push_back(k);
                                }
                            }
                        }
                    }
                    if (v3 == -1 || v4 == -1) {
                        loop = false;
                    }

                } while (loop);
            }
            if (one_layer.size() != 0) {
                layer.push_back(one_layer);
            }
        }

        std::vector<std::vector<std::vector<int>>> relative_point;
        for (int i = 0; i < layer.size(); ++i) {
            int v = layer[i][0];
            int v1, v2, v3, v4;
            int start1, start2;
            std::vector<std::vector<int>> patch_point;
            std::vector<int> single_point;
            v1 = collapse_patch[v][0][0];
            v2 = collapse_patch[v][1][0];
            v3 = collapse_patch[v][2][0];
            v4 = collapse_patch[v][3][0];
            if (layer[i].size() == 1) {
                single_point.push_back(v1);
                single_point.push_back(v2);
                single_point.push_back(v3);
                single_point.push_back(v4);
                patch_point.push_back(single_point);
                relative_point.push_back(patch_point);
                continue;
            }
            std::vector<int> clearlove;
            clearlove.push_back(collapse_patch[layer[i][1]][0][0]);
            clearlove.push_back(collapse_patch[layer[i][1]][1][0]);
            clearlove.push_back(collapse_patch[layer[i][1]][2][0]);
            clearlove.push_back(collapse_patch[layer[i][1]][3][0]);
            for (int j = 0; j < 4; ++j) {
                int o1, o2, o3, o4;
                o1 = collapse_patch[v][j][0];
                o2 = collapse_patch[v][(j + 1) % 4][0];
                o3 = collapse_patch[v][(j + 2) % 4][0];
                o4 = collapse_patch[v][(j + 3) % 4][0];

                if (find(clearlove.begin(), clearlove.end(), o1) != clearlove.end() &&
                    find(clearlove.begin(), clearlove.end(), o2) != clearlove.end()) {
                    single_point.push_back(o3);
                    single_point.push_back(o4);
                    single_point.push_back(o1);
                    single_point.push_back(o2);
                    start1 = o1;
                    start2 = o2;
                    break;
                }
            }
            patch_point.push_back(single_point);
            for (int j = 1; j < layer[i].size(); ++j) {
                std::vector<int> single;
                v = layer[i][j];
                for (int k = 0; k < 4; ++k) {
                    v1 = collapse_patch[v][k][0];
                    v2 = collapse_patch[v][(k + 1) % 4][0];
                    v3 = collapse_patch[v][(k + 2) % 4][0];
                    v4 = collapse_patch[v][(k + 3) % 4][0];
                    if (start1 == v1 && start2 == v2) {
                        single.push_back(v1);
                        single.push_back(v2);
                        single.push_back(v3);
                        single.push_back(v4);
                        start1 = v3;
                        start2 = v4;
                        break;
                    } else if (start1 == v2 && start2 == v1) {
                        single.push_back(v1);
                        single.push_back(v2);
                        single.push_back(v3);
                        single.push_back(v4);
                        start1 = v3;
                        start2 = v4;
                        break;
                    }
                }
                patch_point.push_back(single);
            }
            relative_point.push_back(patch_point);
        }
        patch_compact.resize(patch.size());
        for (int i = 0; i < patch.size(); ++i) {
            int v = point[i];
            int v_adj1, v_adj2, next;
            if (Inner_quad[i][0] == v) {
                v_adj1 = Inner_quad[i][3];
                v_adj2 = Inner_quad[i][1];
                next = Inner_quad[i][2];
            }
            if (Inner_quad[i][1] == v) {
                v_adj1 = Inner_quad[i][0];
                v_adj2 = Inner_quad[i][2];
                next = Inner_quad[i][3];
            }
            if (Inner_quad[i][2] == v) {
                v_adj1 = Inner_quad[i][1];
                v_adj2 = Inner_quad[i][3];
                next = Inner_quad[i][0];
            }
            if (Inner_quad[i][3] == v) {
                v_adj1 = Inner_quad[i][2];
                v_adj2 = Inner_quad[i][0];
                next = Inner_quad[i][1];
            }
            if (patch[i][0].size() == 2 || patch[i][1].size() == 2 || patch[i][2].size() == 2 ||
                patch[i][3].size() == 2) {
                if (collapse_patch[i][0].size() > collapse_patch[i][1].size()) {
                    patch_compact[i].push_back(collapse_patch[i][0]);
                    std::vector<int> edge;
                    for (int j = 0; j < collapse_patch[i][2].size(); ++j) {
                        edge.push_back(collapse_patch[i][2][collapse_patch[i][2].size() - 1 - j]);
                    }
                    patch_compact[i].push_back(edge);
                } else {
                    patch_compact[i].push_back(collapse_patch[i][1]);
                    std::vector<int> edge;
                    for (int j = 0; j < collapse_patch[i][3].size(); ++j) {
                        edge.push_back(collapse_patch[i][3][collapse_patch[i][3].size() - 1 - j]);
                    }
                    patch_compact[i].push_back(edge);
                }
                continue;
            }

            std::vector<int> boundary_line;
            for (int j = 0; j < patch[i].size(); ++j) {
                for (int m = 0; m < patch[i][j].size(); ++m) {
                    if (find(boundary.begin(), boundary.end(), patch[i][j][m]) == boundary.end()) {
                        boundary_line.push_back(patch[i][j][m]);
                    }
                }
            }
            std::vector<int> assist_line;
            int point = v_adj2;
            int point_next = next;
            assist_line.push_back(v_adj2);
            assist_line.push_back(next);
            do {
                int v2, v3;
                for (int j = 0; j < VtoF[point_next].size(); ++j) {
                    int x = VtoF[point_next][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == point && F_compact[x][(k + 1) % 4] == point_next) {
                            v2 = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                for (int j = 0; j < VtoF[point_next].size(); ++j) {
                    int x = VtoF[point_next][j];
                    for (int k = 0; k < F_compact[x].size(); ++k) {
                        if (F_compact[x][k] == v2 && F_compact[x][(k + 1) % 4] == point_next) {
                            v3 = F_compact[x][(k + 2) % 4];
                        }
                    }
                }
                point = point_next;
                point_next = v3;
                assist_line.push_back(point_next);
            } while (find(boundary_line.begin(), boundary_line.end(), point_next) ==
                     boundary_line.end());
            std::vector<std::vector<int>> line;
            std::vector<int> line_line;
            std::vector<int> line_line1;
            int endl;
            for (int j = 0; j < patch[i].size(); ++j) {
                if (patch[i][j][0] == v && patch[i][j][1] == v_adj1) {
                    line.push_back(patch[i][j]);
                    endl = patch[i][j][patch[i][j].size() - 1];
                }
                if (patch[i][j][patch[i][j].size() - 1] == v &&
                    patch[i][j][patch[i][j].size() - 2] == v_adj1) {
                    endl = patch[i][j][0];
                    for (int k = 0; k < patch[i][j].size(); ++k) {
                        line_line.push_back(patch[i][j][patch[i][j].size() - k - 1]);
                    }
                    line.push_back(line_line);
                }
                if (patch[i][j][0] == v && patch[i][j][1] == v_adj2) {
                    patch_compact[i].push_back(patch[i][j]);
                }
                if (patch[i][j][patch[i][j].size() - 1] == v &&
                    patch[i][j][patch[i][j].size() - 2] == v_adj2) {
                    for (int k = 0; k < patch[i][j].size(); ++k) {
                        line_line1.push_back(patch[i][j][patch[i][j].size() - k - 1]);
                    }
                    patch_compact[i].push_back(line_line1);
                }
            }

            for (int j = 0; j < line[0].size(); ++j) {
                if (j == 0 || j == line[0].size() - 1) continue;
                std::vector<int> inner_trace;
                int v1 = line[0][j];
                int v2 = assist_line[j];
                inner_trace.push_back(v1);
                inner_trace.push_back(v2);
                do {
                    int v3, v4;
                    for (int j = 0; j < VtoF[v2].size(); ++j) {
                        int x = VtoF[v2][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v1 && F_compact[x][(k + 1) % 4] == v2) {
                                v3 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    for (int j = 0; j < VtoF[v2].size(); ++j) {
                        int x = VtoF[v2][j];
                        for (int k = 0; k < F_compact[x].size(); ++k) {
                            if (F_compact[x][k] == v3 && F_compact[x][(k + 1) % 4] == v2) {
                                v4 = F_compact[x][(k + 2) % 4];
                            }
                        }
                    }
                    v1 = v2;
                    v2 = v4;
                    inner_trace.push_back(v2);
                } while (find(boundary_line.begin(), boundary_line.end(), v2) ==
                         boundary_line.end());
                patch_compact[i].push_back(inner_trace);
            }
            std::vector<int> boundary_trace;
            for (int j = 0; j < patch[i].size(); ++j) {
                if (patch[i][j][0] == endl && patch[i][j][patch[i][j].size() - 1] != v) {
                    patch_compact[i].push_back(patch[i][j]);
                }
                if (patch[i][j][0] != v && patch[i][j][patch[i][j].size() - 1] == endl) {
                    for (int k = 0; k < patch[i][j].size(); ++k) {
                        boundary_trace.push_back(patch[i][j][patch[i][j].size() - k - 1]);
                    }
                    patch_compact[i].push_back(boundary_trace);
                }
            }
        }
        /*std::vector<int> delete_patch;
        std::vector<int> renew_O(O_compact.size(), -1);*/
        std::vector<int> valence(O_compact.size());
        for (int i = 0; i < valence.size(); ++i) {
            valence[i] = VtoF[i].size();
        }
        std::vector<int> two_sin_patch;
        std::vector<std::vector<int>> two_sin;
        std::vector<int> two_sin_po;
        std::vector<int> quad_v;
        for (int i = 0; i < thr_for_sin_patch.size(); ++i) {
            int n = thr_for_sin_patch[i];
            for (int i = 0; i < collapse_patch[n].size(); ++i) {
                int v = collapse_patch[n][i][0];
                if (find(all_singularity.begin(), all_singularity.end(), v) !=
                    all_singularity.end()) {
                    quad_v.push_back(v);
                }
            }
        }
        for (int i = 0; i < relative_point.size(); ++i) {
            if (relative_point.size() == 1) continue;
            for (int j = 0; j < relative_point[i].size(); ++j) {
                int v1 = relative_point[i][j][0];
                int v2 = relative_point[i][j][1];
                int a = 0;
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    a += 1;
                }
                if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                    all_singularity.end()) {
                    a += 1;
                }
                if (a == 2) {
                    if (find(quad_v.begin(), quad_v.end(), v2) != quad_v.end() &&
                        find(quad_v.begin(), quad_v.end(), v1) != quad_v.end())
                        break;
                    std::vector<int> mat;
                    mat.push_back(v1);
                    mat.push_back(v2);
                    two_sin.push_back(mat);
                    two_sin_patch.push_back(i);
                    two_sin_po.push_back(j);
                    break;
                }
                if (j == relative_point[i].size() - 1) {
                    int v3 = relative_point[i][j][2];
                    int v4 = relative_point[i][j][3];
                    int b = 0;
                    if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                        b += 1;
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end()) {
                        b += 1;
                    }
                    if (b == 2) {
                        if (find(quad_v.begin(), quad_v.end(), v3) != quad_v.end() &&
                            find(quad_v.begin(), quad_v.end(), v4) != quad_v.end())
                            break;
                        std::vector<int> mat;
                        mat.push_back(v3);
                        mat.push_back(v4);
                        two_sin.push_back(mat);
                        two_sin_patch.push_back(i);
                        two_sin_po.push_back(j);
                        break;
                    }
                }
            }
        }
        if (two_sin_patch.size() != 0) {
            for (int i = 0; i < two_sin_patch.size(); ++i) {
                int v = two_sin_patch[i];
                int five = 0, three = 0;
                int v1 = two_sin[i][0];
                int v2 = two_sin[i][1];
                if (valence[v1] == 5) {
                    five += 1;
                }
                if (valence[v1] == 3) {
                    three += 1;
                }
                if (valence[v2] == 5) {
                    five += 1;
                }
                if (valence[v2] == 3) {
                    three += 1;
                }
                if (five == 1 && three == 1) {
                    /*sin_ft(O_compact, patch_compact, valence, relative_point, layer,
                           collapse_patch, v, v1, v2, all_singularity, thr_for_sin_patch);*/
                }
            }
        }
        if (thr_for_sin_patch.size() != 0) {
            for (int i = 0; i < thr_for_sin_patch.size(); ++i) {
                int v = thr_for_sin_patch[i];
                int five = 0, three = 0;
                bool boundary = false;
                for (int j = 0; j < collapse_patch[v].size(); ++j) {
                    int v1 = collapse_patch[v][j][0];
                    if (find(all_singularity.begin(), all_singularity.end(), v1) ==
                        all_singularity.end())
                        continue;
                    if (boundary_o[v1]) {
                        if (valence[v1] == 3) {
                            five += 1;
                            boundary = true;
                        }
                    } else {
                        if (valence[v1] == 5) {
                            five += 1;
                        }
                        if (valence[v1] == 3) {
                            three += 1;
                        }
                    }
                }
                if (five == 2 && three == 1 && boundary == false) {
                    sin_fft(O_compact, patch_compact, valence, relative_point, layer,
                            collapse_patch, v, all_singularity, thr_for_sin_patch);
                }
                if (five == 1 && three == 2) {
                    sin_ftt(O_compact, patch_compact, valence, relative_point, layer,
                            collapse_patch, v, all_singularity, thr_for_sin_patch);
                }
                if (five == 2 && three == 1 && boundary) {
                    sin_fft_boundary(O_compact, boundary_o, patch_compact, valence, relative_point,
                                     layer, collapse_patch, v, all_singularity, thr_for_sin_patch);
                }
            }
        }
        /*relative_point.resize(0);*/
        collapse_boundary(O_compact, patch_compact, valence, relative_point, layer, collapse_patch,
                          all_singularity);
        std::vector<int> cannot_collapse;
        for (int i = 0; i < relative_point.size(); ++i) {
            if (relative_point[i].size() == 1) continue;
            int v1 = relative_point[i][0][0];
            int v2 = relative_point[i][0][1];
            int v3 = relative_point[i][relative_point[i].size() - 1][2];
            int v4 = relative_point[i][relative_point[i].size() - 1][3];
            for (int j = 0; j < thr_for_sin_patch.size(); ++j) {
                std::vector<int> quad;
                int v = thr_for_sin_patch[j];
                quad.push_back(collapse_patch[v][0][0]);
                quad.push_back(collapse_patch[v][1][0]);
                quad.push_back(collapse_patch[v][2][0]);
                quad.push_back(collapse_patch[v][3][0]);
                if (find(quad.begin(), quad.end(), v1) != quad.end() &&
                    find(quad.begin(), quad.end(), v2) != quad.end()) {
                    cannot_collapse.push_back(i);
                    break;
                }
                if (find(quad.begin(), quad.end(), v3) != quad.end() &&
                    find(quad.begin(), quad.end(), v4) != quad.end()) {
                    cannot_collapse.push_back(i);
                    break;
                }
            }
        }
        std::vector<std::vector<int>> need_collapse_part(relative_point.size());
        for (int i = 0; i < relative_point.size(); ++i) {
            if (relative_point[i].size() == 1) continue;
            if (find(cannot_collapse.begin(), cannot_collapse.end(), i) != cannot_collapse.end())
                continue;
            int v1 = relative_point[i][0][0];
            int v2 = relative_point[i][0][1];
            int v3 = relative_point[i][0][2];
            int v4 = relative_point[i][0][3];
            int begin = -1;
            int start = -1;
            int x = layer[i][0];
            int c = 0, d = 0, sizec = 0, sized = 0;
            if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end() ||
                find(all_singularity.begin(), all_singularity.end(), v2) !=
                    all_singularity.end()) {
                for (int j = 0; j < collapse_patch[x].size(); ++j) {
                    if (collapse_patch[x][j][0] == v2 &&
                        collapse_patch[x][j][collapse_patch[x][j].size() - 1] == v3) {
                        sizec = collapse_patch[x][j].size();
                        for (int k = 0; k < collapse_patch[x][j].size(); ++k) {
                            if (boundary_o[collapse_patch[x][j][k]]) {
                                c += 1;
                            }
                        }
                    }
                    if (collapse_patch[x][j][0] == v4 &&
                        collapse_patch[x][j][collapse_patch[x][j].size() - 1] == v1) {
                        sized = collapse_patch[x][j].size();
                        for (int k = 0; k < collapse_patch[x][j].size(); ++k) {
                            if (boundary_o[collapse_patch[x][j][k]]) {
                                d += 1;
                            }
                        }
                    }
                }
            }
            if (sizec == 2 && boundary_o[v2] && boundary_o[v3]) continue;
            if (sized == 2 && boundary_o[v4] && boundary_o[v1]) continue;
            if (c > 2) continue;
            if (d > 2) continue;
            if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                all_singularity.end()) {
                if (find(all_singularity.begin(), all_singularity.end(), v2) ==
                    all_singularity.end()) {
                    begin = 0;
                    start = 3;
                }
            }
            if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                all_singularity.end()) {
                if (find(all_singularity.begin(), all_singularity.end(), v1) ==
                    all_singularity.end()) {
                    begin = 0;
                    start = 4;
                }
            }
            if (begin != -1) {
                if (start == 3) {
                    if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end() &&
                        find(all_singularity.begin(), all_singularity.end(), v4) ==
                            all_singularity.end()) {
                        need_collapse_part[i].push_back(begin);
                        need_collapse_part[i].push_back(begin);
                        begin = 1;
                        start = 4;
                        if (find(need_collapse.begin(), need_collapse.end(), i) ==
                            need_collapse.end()) {
                            need_collapse.push_back(i);
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end() &&
                        find(all_singularity.begin(), all_singularity.end(), v3) ==
                            all_singularity.end()) {
                        begin = 1;
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end() &&
                        find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                        start = -1;
                    }
                }
                if (start == 4) {
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end() &&
                        find(all_singularity.begin(), all_singularity.end(), v3) ==
                            all_singularity.end()) {
                        need_collapse_part[i].push_back(begin);
                        need_collapse_part[i].push_back(begin);
                        begin = 1;
                        start = 3;
                        if (find(need_collapse.begin(), need_collapse.end(), i) ==
                            need_collapse.end()) {
                            need_collapse.push_back(i);
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end() &&
                        find(all_singularity.begin(), all_singularity.end(), v4) ==
                            all_singularity.end()) {
                        begin = 1;
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end() &&
                        find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                        start = -1;
                    }
                }
            }
            for (int j = 1; j < relative_point[i].size(); ++j) {
                if (start == -1) {
                    v1 = relative_point[i][j][0];
                    v2 = relative_point[i][j][1];
                    v3 = relative_point[i][j][2];
                    v4 = relative_point[i][j][3];
                    int y = layer[i][j];
                    int e = 0, f = 0, sizee = 0, sizef = 0;
                    for (int jjj = 0; jjj < collapse_patch[y].size(); ++jjj) {
                        if (collapse_patch[y][jjj][0] == v2 &&
                            collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v3) {
                            sizee = collapse_patch[y][jjj].size();
                            for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                                if (boundary_o[collapse_patch[y][jjj][k]]) {
                                    e += 1;
                                }
                            }
                        }
                        if (collapse_patch[y][jjj][0] == v4 &&
                            collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v1) {
                            sizef = collapse_patch[y][jjj].size();
                            for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                                if (boundary_o[collapse_patch[y][jjj][k]]) {
                                    f += 1;
                                }
                            }
                        }
                    }
                    if (sizee == 2 && boundary_o[v2] && boundary_o[v3]) continue;
                    if (sizef == 2 && boundary_o[v4] && boundary_o[v1]) continue;
                    if (e > 2) continue;
                    if (f > 2) continue;
                    if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v2) ==
                            all_singularity.end()) {
                            begin = j;
                            start = 3;
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v1) ==
                            all_singularity.end()) {
                            begin = j;
                            start = 4;
                        }
                    }
                    if (begin != -1) {
                        if (start == 3) {
                            if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                    all_singularity.end() &&
                                find(all_singularity.begin(), all_singularity.end(), v4) ==
                                    all_singularity.end()) {
                                need_collapse_part[i].push_back(begin);
                                need_collapse_part[i].push_back(begin);
                                begin = j + 1;
                                start = 4;
                                if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                    need_collapse.end()) {
                                    need_collapse.push_back(i);
                                }
                            }
                            if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                    all_singularity.end() &&
                                find(all_singularity.begin(), all_singularity.end(), v3) ==
                                    all_singularity.end()) {
                                begin = j + 1;
                            }
                            if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                    all_singularity.end() &&
                                find(all_singularity.begin(), all_singularity.end(), v3) !=
                                    all_singularity.end()) {
                                start = -1;
                            }
                        }
                        if (start == 4) {
                            if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                    all_singularity.end() &&
                                find(all_singularity.begin(), all_singularity.end(), v3) ==
                                    all_singularity.end()) {
                                need_collapse_part[i].push_back(begin);
                                need_collapse_part[i].push_back(begin);
                                begin = j + 1;
                                start = 3;
                                if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                    need_collapse.end()) {
                                    need_collapse.push_back(i);
                                }
                            }
                            if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                    all_singularity.end() &&
                                find(all_singularity.begin(), all_singularity.end(), v4) ==
                                    all_singularity.end()) {
                                begin = j + 1;
                            }
                            if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                    all_singularity.end() &&
                                find(all_singularity.begin(), all_singularity.end(), v3) !=
                                    all_singularity.end()) {
                                start = -1;
                            }
                        }
                    }
                } else {
                    v1 = relative_point[i][j][0];
                    v2 = relative_point[i][j][1];
                    v3 = relative_point[i][j][2];
                    v4 = relative_point[i][j][3];
                    int y = layer[i][j];
                    int e = 0, f = 0, sizee = 0, sizef = 0;
                    for (int jjj = 0; jjj < collapse_patch[y].size(); ++jjj) {
                        if (collapse_patch[y][jjj][0] == v2 &&
                            collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v3) {
                            sizee = collapse_patch[y][jjj].size();
                            for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                                if (boundary_o[collapse_patch[y][jjj][k]]) {
                                    e += 1;
                                }
                            }
                        }
                        if (collapse_patch[y][jjj][0] == v4 &&
                            collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v1) {
                            sizef = collapse_patch[y][jjj].size();
                            for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                                if (boundary_o[collapse_patch[y][jjj][k]]) {
                                    f += 1;
                                }
                            }
                        }
                    }
                    if (sizee == 2 && boundary_o[v2] && boundary_o[v3]) {
                        begin = -1;
                        start = -1;
                        continue;
                    }
                    if (sizef == 2 && boundary_o[v4] && boundary_o[v1]) {
                        begin = -1;
                        start = -1;
                        continue;
                    }
                    if (e > 2) {
                        begin = -1;
                        start = -1;
                        continue;
                    }
                    if (f > 2) {
                        begin = -1;
                        start = -1;
                        continue;
                    }
                    if (start == 3) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                            if (find(all_singularity.begin(), all_singularity.end(), v4) ==
                                all_singularity.end()) {
                                need_collapse_part[i].push_back(begin);
                                need_collapse_part[i].push_back(j);
                                begin = j + 1;
                                start = 4;
                                if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                    need_collapse.end()) {
                                    need_collapse.push_back(i);
                                }
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                            if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end()) {
                                begin = -1;
                                start = -1;
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v3) ==
                            all_singularity.end()) {
                            if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end()) {
                                begin = j + 1;
                            }
                        }
                    }
                    if (start == 4) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end()) {
                            if (find(all_singularity.begin(), all_singularity.end(), v3) ==
                                all_singularity.end()) {
                                need_collapse_part[i].push_back(begin);
                                need_collapse_part[i].push_back(j);
                                begin = j + 1;
                                start = 3;
                                if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                    need_collapse.end()) {
                                    need_collapse.push_back(i);
                                }
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end()) {
                            if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end()) {
                                begin = -1;
                                start = -1;
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) ==
                            all_singularity.end()) {
                            if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end()) {
                                begin = j + 1;
                            }
                        }
                    }
                }
            }
        }
        std::vector<std::vector<int>> double_q(layer.size());
        std::vector<int> aaaa;
        for (int i = 0; i < need_collapse.size(); ++i) {
            int n = need_collapse[i];
            int begin, end;
            begin = need_collapse_part[n][0];
            end = need_collapse_part[n][1];
            for (int j = begin; j < end + 1; ++j) {
                if (find(aaaa.begin(), aaaa.end(), layer[n][j]) == aaaa.end()) {
                    aaaa.push_back(layer[n][j]);
                } else {
                    double_q[n].push_back(layer[n][j]);
                }
            }
        }
        for (int i = 0; i < need_collapse.size(); ++i) {
            sin_mismatching(need_collapse_part, patch_compact, O_compact, relative_point, layer,
                            need_collapse, all_singularity, collapse_patch, i, double_q,
                            boundary_o);
        }
        size_three(patch_compact, O_compact, all_singularity, collapse_patch, boundary_o);
    }
}

void Optimizer::size_three(std::vector<std::vector<std::vector<int>>>& patch_compact,
                           std::vector<Vector3d>& O_compact, std::vector<int>& all_singularity,
                           std::vector<std::vector<std::vector<int>>>& collapse_patch,
                           std::vector<int>& boundary_o) {
    std::vector<std::vector<int>> collapse_patch_edge;
    std::vector<std::vector<int>> layer;
    for (int i = 0; i < collapse_patch.size(); ++i) {
        std::vector<int> line;
        line.push_back(1);
        line.push_back(1);
        line.push_back(1);
        line.push_back(1);
        collapse_patch_edge.push_back(line);
    }
    for (int i = 0; i < collapse_patch.size(); ++i) {
        std::vector<int> one_layer;
        int rev = -1;
        for (int j = 0; j < collapse_patch[i].size(); ++j) {
            if (collapse_patch_edge[i][j] == -1) continue;
            if (collapse_patch[i][j].size() > 5) {
                collapse_patch_edge[i][j] = -1;
                continue;
            }
            if (one_layer.size() > 1 && rev != collapse_patch[i].size() - 1 && rev != -1) {
                if (j == rev + 1) continue;
            }
            rev = j;
            collapse_patch_edge[i][j] = -1;
            if (find(one_layer.begin(), one_layer.end(), i) == one_layer.end()) {
                one_layer.push_back(i);
            }

            int v1, v2;
            v1 = collapse_patch[i][j][0];
            v2 = collapse_patch[i][j][collapse_patch[i][j].size() - 1];
            bool loop = true;
            do {
                int v3 = -1;
                int v4 = -1;
                for (int k = 0; k < collapse_patch.size(); ++k) {
                    if (k == i) continue;
                    for (int p = 0; p < collapse_patch[k].size(); ++p) {
                        if (collapse_patch_edge[k][p] == -1) continue;
                        if (collapse_patch[k][p].size() > 5) continue;
                        if (collapse_patch[k][p][0] == v1 &&
                            collapse_patch[k][p][collapse_patch[k][p].size() - 1] == v2) {
                            collapse_patch_edge[k][p] = -1;
                            collapse_patch_edge[k][(p + 2) % 4] = -1;
                            v3 = collapse_patch[k][(p + 2) % 4][0];
                            v4 = collapse_patch[k][(p + 2) % 4]
                                               [collapse_patch[k][(p + 2) % 4].size() - 1];
                            v1 = v3;
                            v2 = v4;
                            if (j == 0 || j == 1) {
                                one_layer.insert(one_layer.begin(), k);
                            } else {
                                one_layer.push_back(k);
                            }

                        } else if (collapse_patch[k][p][0] == v2 &&
                                   collapse_patch[k][p][collapse_patch[k][p].size() - 1] == v1) {
                            collapse_patch_edge[k][p] = -1;
                            collapse_patch_edge[k][(p + 2) % 4] = -1;
                            v3 = collapse_patch[k][(p + 2) % 4][0];
                            v4 = collapse_patch[k][(p + 2) % 4]
                                               [collapse_patch[k][(p + 2) % 4].size() - 1];
                            v1 = v3;
                            v2 = v4;
                            if (j == 0 || j == 1) {
                                one_layer.insert(one_layer.begin(), k);
                            } else {
                                one_layer.push_back(k);
                            }
                        }
                    }
                }
                if (v3 == -1 || v4 == -1) {
                    loop = false;
                }

            } while (loop);
        }
        if (one_layer.size() != 0) {
            layer.push_back(one_layer);
        }
    }

    std::vector<std::vector<std::vector<int>>> relative_point;
    for (int i = 0; i < layer.size(); ++i) {
        int v = layer[i][0];
        int v1, v2, v3, v4;
        int start1, start2;
        std::vector<std::vector<int>> patch_point;
        std::vector<int> single_point;
        v1 = collapse_patch[v][0][0];
        v2 = collapse_patch[v][1][0];
        v3 = collapse_patch[v][2][0];
        v4 = collapse_patch[v][3][0];
        if (layer[i].size() == 1) {
            single_point.push_back(v1);
            single_point.push_back(v2);
            single_point.push_back(v3);
            single_point.push_back(v4);
            patch_point.push_back(single_point);
            relative_point.push_back(patch_point);
            continue;
        }
        std::vector<int> clearlove;
        clearlove.push_back(collapse_patch[layer[i][1]][0][0]);
        clearlove.push_back(collapse_patch[layer[i][1]][1][0]);
        clearlove.push_back(collapse_patch[layer[i][1]][2][0]);
        clearlove.push_back(collapse_patch[layer[i][1]][3][0]);
        for (int j = 0; j < 4; ++j) {
            int o1, o2, o3, o4;
            o1 = collapse_patch[v][j][0];
            o2 = collapse_patch[v][(j + 1) % 4][0];
            o3 = collapse_patch[v][(j + 2) % 4][0];
            o4 = collapse_patch[v][(j + 3) % 4][0];

            if (find(clearlove.begin(), clearlove.end(), o1) != clearlove.end() &&
                find(clearlove.begin(), clearlove.end(), o2) != clearlove.end()) {
                single_point.push_back(o3);
                single_point.push_back(o4);
                single_point.push_back(o1);
                single_point.push_back(o2);
                start1 = o1;
                start2 = o2;
                break;
            }
        }
        patch_point.push_back(single_point);
        for (int j = 1; j < layer[i].size(); ++j) {
            std::vector<int> single;
            v = layer[i][j];
            for (int k = 0; k < 4; ++k) {
                v1 = collapse_patch[v][k][0];
                v2 = collapse_patch[v][(k + 1) % 4][0];
                v3 = collapse_patch[v][(k + 2) % 4][0];
                v4 = collapse_patch[v][(k + 3) % 4][0];
                if (start1 == v1 && start2 == v2) {
                    single.push_back(v1);
                    single.push_back(v2);
                    single.push_back(v3);
                    single.push_back(v4);
                    start1 = v3;
                    start2 = v4;
                    break;
                } else if (start1 == v2 && start2 == v1) {
                    single.push_back(v1);
                    single.push_back(v2);
                    single.push_back(v3);
                    single.push_back(v4);
                    start1 = v3;
                    start2 = v4;
                    break;
                }
            }
            patch_point.push_back(single);
        }
        relative_point.push_back(patch_point);
    }

    std::vector<std::vector<int>> need_collapse_part(relative_point.size());
    std::vector<int> need_collapse;
    for (int i = 0; i < relative_point.size(); ++i) {
        if (relative_point[i].size() == 1) continue;
        int v1 = relative_point[i][0][0];
        int v2 = relative_point[i][0][1];
        int v3 = relative_point[i][0][2];
        int v4 = relative_point[i][0][3];
        int begin = -1;
        int start = -1;
        int x = layer[i][0];
        int c = 0, d = 0, sizec = 0, sized = 0;
        if (find(all_singularity.begin(), all_singularity.end(), v1) != all_singularity.end() ||
            find(all_singularity.begin(), all_singularity.end(), v2) != all_singularity.end()) {
            for (int j = 0; j < collapse_patch[x].size(); ++j) {
                if (collapse_patch[x][j][0] == v2 &&
                    collapse_patch[x][j][collapse_patch[x][j].size() - 1] == v3) {
                    sizec = collapse_patch[x][j].size();
                    for (int k = 0; k < collapse_patch[x][j].size(); ++k) {
                        if (boundary_o[collapse_patch[x][j][k]]) {
                            c += 1;
                        }
                    }
                }
                if (collapse_patch[x][j][0] == v4 &&
                    collapse_patch[x][j][collapse_patch[x][j].size() - 1] == v1) {
                    sized = collapse_patch[x][j].size();
                    for (int k = 0; k < collapse_patch[x][j].size(); ++k) {
                        if (boundary_o[collapse_patch[x][j][k]]) {
                            d += 1;
                        }
                    }
                }
            }
        }
        if (sizec == 2 && boundary_o[v2] && boundary_o[v3]) continue;
        if (sized == 2 && boundary_o[v4] && boundary_o[v1]) continue;
        if (find(all_singularity.begin(), all_singularity.end(), v1) != all_singularity.end() &&
                find(all_singularity.begin(), all_singularity.end(), v2) !=
                    all_singularity.end() ||
            find(all_singularity.begin(), all_singularity.end(), v3) != all_singularity.end() &&
                find(all_singularity.begin(), all_singularity.end(), v4) != all_singularity.end())
            continue;
        if (c > 2) continue;
        if (d > 2) continue;
        if (find(all_singularity.begin(), all_singularity.end(), v1) != all_singularity.end()) {
            if (find(all_singularity.begin(), all_singularity.end(), v2) ==
                all_singularity.end()) {
                begin = 0;
                start = 3;
            }
        }
        if (find(all_singularity.begin(), all_singularity.end(), v2) != all_singularity.end()) {
            if (find(all_singularity.begin(), all_singularity.end(), v1) ==
                all_singularity.end()) {
                begin = 0;
                start = 4;
            }
        }
        if (begin != -1) {
            if (start == 3) {
                if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v4) ==
                        all_singularity.end()) {
                    need_collapse_part[i].push_back(begin);
                    need_collapse_part[i].push_back(begin);
                    begin = 1;
                    start = 4;
                    if (find(need_collapse.begin(), need_collapse.end(), i) ==
                        need_collapse.end()) {
                        need_collapse.push_back(i);
                    }
                }
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) ==
                        all_singularity.end()) {
                    begin = 1;
                }
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                    start = -1;
                }
            }
            if (start == 4) {
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) ==
                        all_singularity.end()) {
                    need_collapse_part[i].push_back(begin);
                    need_collapse_part[i].push_back(begin);
                    begin = 1;
                    start = 3;
                    if (find(need_collapse.begin(), need_collapse.end(), i) ==
                        need_collapse.end()) {
                        need_collapse.push_back(i);
                    }
                }
                if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v4) ==
                        all_singularity.end()) {
                    begin = 1;
                }
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                    start = -1;
                }
            }
        }
        for (int j = 1; j < relative_point[i].size(); ++j) {
            if (start == -1) {
                v1 = relative_point[i][j][0];
                v2 = relative_point[i][j][1];
                v3 = relative_point[i][j][2];
                v4 = relative_point[i][j][3];
                int y = layer[i][j];
                int e = 0, f = 0, sizee = 0, sizef = 0;
                for (int jjj = 0; jjj < collapse_patch[y].size(); ++jjj) {
                    if (collapse_patch[y][jjj][0] == v2 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v3) {
                        sizee = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                e += 1;
                            }
                        }
                    }
                    if (collapse_patch[y][jjj][0] == v4 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v1) {
                        sizef = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                f += 1;
                            }
                        }
                    }
                }
                if (sizee == 2 && boundary_o[v2] && boundary_o[v3]) continue;
                if (sizef == 2 && boundary_o[v4] && boundary_o[v1]) continue;
                if (e > 2) continue;
                if (f > 2) continue;
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    if (find(all_singularity.begin(), all_singularity.end(), v2) ==
                        all_singularity.end()) {
                        begin = j;
                        start = 3;
                    }
                }
                if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                    all_singularity.end()) {
                    if (find(all_singularity.begin(), all_singularity.end(), v1) ==
                        all_singularity.end()) {
                        begin = j;
                        start = 4;
                    }
                }
                if (begin != -1) {
                    if (start == 3) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v4) ==
                                all_singularity.end()) {
                            need_collapse_part[i].push_back(begin);
                            need_collapse_part[i].push_back(begin);
                            begin = j + 1;
                            start = 4;
                            if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                need_collapse.end()) {
                                need_collapse.push_back(i);
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) ==
                                all_singularity.end()) {
                            begin = j + 1;
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end()) {
                            start = -1;
                        }
                    }
                    if (start == 4) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) ==
                                all_singularity.end()) {
                            need_collapse_part[i].push_back(begin);
                            need_collapse_part[i].push_back(begin);
                            begin = j + 1;
                            start = 3;
                            if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                need_collapse.end()) {
                                need_collapse.push_back(i);
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v4) ==
                                all_singularity.end()) {
                            begin = j + 1;
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end()) {
                            start = -1;
                        }
                    }
                }
            } else {
                v1 = relative_point[i][j][0];
                v2 = relative_point[i][j][1];
                v3 = relative_point[i][j][2];
                v4 = relative_point[i][j][3];
                int y = layer[i][j];
                int e = 0, f = 0, sizee = 0, sizef = 0;
                for (int jjj = 0; jjj < collapse_patch[y].size(); ++jjj) {
                    if (collapse_patch[y][jjj][0] == v2 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v3) {
                        sizee = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                e += 1;
                            }
                        }
                    }
                    if (collapse_patch[y][jjj][0] == v4 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v1) {
                        sizef = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                f += 1;
                            }
                        }
                    }
                }
                if (sizee == 2 && boundary_o[v2] && boundary_o[v3]) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (sizef == 2 && boundary_o[v4] && boundary_o[v1]) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (e > 2) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (f > 2) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (start == 3) {
                    if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) ==
                            all_singularity.end()) {
                            need_collapse_part[i].push_back(begin);
                            need_collapse_part[i].push_back(j);
                            begin = j + 1;
                            start = 4;
                            if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                need_collapse.end()) {
                                need_collapse.push_back(i);
                            }
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end()) {
                            begin = -1;
                            start = -1;
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v3) ==
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end()) {
                            begin = j + 1;
                        }
                    }
                }
                if (start == 4) {
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) ==
                            all_singularity.end()) {
                            need_collapse_part[i].push_back(begin);
                            need_collapse_part[i].push_back(j);
                            begin = j + 1;
                            start = 3;
                            if (find(need_collapse.begin(), need_collapse.end(), i) ==
                                need_collapse.end()) {
                                need_collapse.push_back(i);
                            }
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                            begin = -1;
                            start = -1;
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) ==
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                            begin = j + 1;
                        }
                    }
                }
            }
        }
    }
    std::vector<std::vector<int>> double_q(layer.size());
    std::vector<int> aaaa;
    for (int i = 0; i < need_collapse.size(); ++i) {
        int n = need_collapse[i];
        int begin, end;
        begin = need_collapse_part[n][0];
        end = need_collapse_part[n][1];
        for (int j = begin; j < end + 1; ++j) {
            if (find(aaaa.begin(), aaaa.end(), layer[n][j]) == aaaa.end()) {
                aaaa.push_back(layer[n][j]);
            } else {
                double_q[n].push_back(layer[n][j]);
            }
        }
    }
    for (int i = 0; i < need_collapse.size(); ++i) {
        if (i == 1) continue;
        sin_mismatching(need_collapse_part, patch_compact, O_compact, relative_point, layer,
                        need_collapse, all_singularity, collapse_patch, i, double_q, boundary_o);
    }
}

void Optimizer::sin_fft_boundary(std::vector<Vector3d>& O_compact, std::vector<int> boundary_o,
                                 std::vector<std::vector<std::vector<int>>>& patch_compact,
                                 std::vector<int>& valence,
                                 std::vector<std::vector<std::vector<int>>>& relative_point,
                                 std::vector<std::vector<int>>& layer,
                                 std::vector<std::vector<std::vector<int>>>& collapse_patch, int v,
                                 std::vector<int>& all_singularity,
                                 std::vector<int>& thr_for_sin_patch) {
    std::vector<int> position;
    std::vector<int> renew(O_compact.size(), -1);
    std::vector<int> delete_patch;
    for (int i = 0; i < collapse_patch[v].size(); ++i) {
        int v1 = collapse_patch[v][i][0];
        int v2 = collapse_patch[v][i][1];
        for (int j = 0; j < relative_point.size(); ++j) {
            if (relative_point[j].size() == 1) {
                int n1 = relative_point[j][0][0];
                int n2 = relative_point[j][0][1];
                int n3 = relative_point[j][0][2];
                int n4 = relative_point[j][0][3];
                if (n1 == v2 && n2 == v1) {
                    position.push_back(j);
                    break;
                } else if (n2 == v2 && n3 == v1) {
                    position.push_back(j);
                    break;
                } else if (n3 == v2 && n4 == v1) {
                    position.push_back(j);
                    break;
                } else if (n4 == v2 && n1 == v1) {
                    position.push_back(j);
                    break;
                }
            } else {
                int n1 = relative_point[j][0][0];
                int n2 = relative_point[j][0][1];
                int n3 = relative_point[j][relative_point[j].size() - 1][2];
                int n4 = relative_point[j][relative_point[j].size() - 1][3];
                if (v1 == n1 && v2 == n2) {
                    position.push_back(j);
                    break;
                } else if (v1 == n2 && v2 == n1) {
                    position.push_back(j);
                    break;
                } else if (v1 == n3 && v2 == n4) {
                    position.push_back(j);
                    break;
                } else if (v1 == n4 && v2 == n3) {
                    position.push_back(j);
                    break;
                }
            }
        }
    }
    std::vector<int> singularity;
    std::vector<int> quad;
    int v1, v2, v3, v4;
    int x1, x2, x3, x4;
    v1 = collapse_patch[v][0][0];
    v2 = collapse_patch[v][1][0];
    v3 = collapse_patch[v][2][0];
    v4 = collapse_patch[v][3][0];
    if (boundary_o[v1]) {
        x1 = v1;
        x2 = v2;
        x3 = v3;
        x4 = v4;
    } else if (boundary_o[v2]) {
        x1 = v2;
        x2 = v1;
        x3 = v3;
        x4 = v4;
    } else if (boundary_o[v3]) {
        x1 = v3;
        x2 = v1;
        x3 = v2;
        x4 = v4;
    } else if (boundary_o[v4]) {
        x1 = v4;
        x2 = v1;
        x3 = v2;
        x4 = v3;
    }
    v1 = x1;
    v2 = x2;
    v3 = x3;
    v4 = x4;
    renew[v2] = v1;
    renew[v3] = v1;
    renew[v4] = v1;
    singularity.push_back(v2);
    singularity.push_back(v3);
    singularity.push_back(v4);
    quad.push_back(v1);
    quad.push_back(v2);
    quad.push_back(v3);
    quad.push_back(v4);
    Vector3d p1 = O_compact[v1];
    O_compact[v1] = p1;
    O_compact[v2] = p1;
    O_compact[v3] = p1;
    O_compact[v4] = p1;
    std::vector<int> twice;
    for (int i = 0; i < position.size(); ++i) {
        int n = position[i];
        for (int j = 0; j < layer[n].size(); ++j) {
            if (find(delete_patch.begin(), delete_patch.end(), layer[n][j]) ==
                delete_patch.end()) {
                delete_patch.push_back(layer[n][j]);
            } else {
                twice.push_back(layer[n][j]);
            }
        }
    }
    delete_patch.push_back(v);
    std::vector<int> boundary_layer(position.size(), 0);
    std::vector<int> need_collapse(position.size(), 0);
    for (int i = 0; i < position.size(); ++i) {
        int n = position[i];
        for (int j = 0; j < layer[n].size(); ++j) {
            if (relative_point[n][j][0] == v1 || relative_point[n][j][1] == v1 ||
                relative_point[n][j][2] == v1 || relative_point[n][j][3] == v1) {
                boundary_layer[i] = 1;
                break;
            }
        }
    }
    std::vector<std::vector<int>> corr(position.size());
    for (int i = 0; i < position.size(); ++i) {
        int n = position[i];
        if (layer[n].size() == 1) continue;
        if (boundary_layer[i]) continue;
        for (int j = 0; j < layer[n].size(); ++j) {
            for (int k = 0; k < 4; ++k) {
                int c1 = relative_point[n][j][k];
                if (find(all_singularity.begin(), all_singularity.end(), c1) !=
                        all_singularity.end() &&
                    find(quad.begin(), quad.end(), c1) == quad.end()) {
                    need_collapse[i] = 1;
                }
            }
        }
    }
    std::vector<std::vector<int>> sin_p(position.size());
    std::vector<std::vector<int>> collapse_size(position.size());
    for (int i = 0; i < position.size(); ++i) {
        if (need_collapse[i] == 0) continue;
        int n = position[i];
        int ver1 = collapse_patch[v][i][0];
        int ver2 = collapse_patch[v][i][1];
        if (relative_point[n][0][1] == ver1) {
            for (int j = 0; j < relative_point[n].size(); ++j) {
                int v1 = relative_point[n][j][2];
                int v2 = relative_point[n][j][3];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    sin_p[i].push_back(1);
                } else if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                           all_singularity.end()) {
                    sin_p[i].push_back(0);
                } else {
                    sin_p[i].push_back(2);
                }
            }
            if (sin_p[i][0] == 2) {
                int a1, a2;
                for (int j = 0; j < sin_p[i][0]; ++j) {
                    if (sin_p[i][j] == 2) continue;
                    a1 = j;
                    a2 = sin_p[i][j];
                    break;
                }
                for (int j = 0; j < a1; ++j) {
                    sin_p[i][j] = a2;
                    break;
                }
            }
            bool spicy = false;
            int a1 = -1;
            for (int j = 0; j < sin_p[i].size(); ++j) {
                if (spicy) {
                    if (sin_p[i][j] != 2) {
                        if (sin_p[i][j] == sin_p[i][a1]) {
                            sin_p[i][j] = 2;
                            spicy = false;
                        } else {
                            for (int k = a1 + 1; k < j; ++k) {
                                sin_p[i][k] = sin_p[i][j];
                            }
                            spicy = false;
                        }
                    }
                } else {
                    if (sin_p[i][j] == 2) {
                        spicy = true;
                        if (a1 == -1) {
                            a1 = j - 1;
                        }
                    }
                }
            }
            int begin = 0, end = sin_p[i].size() - 1;
            int size = 0, sign = sin_p[i][0];
            for (int j = 0; j < sin_p[i].size(); ++j) {
                int n1 = layer[n][j];
                if (sin_p[i][j] == sign) {
                    int aaa;
                    if (collapse_patch[n1][0].size() > collapse_patch[n1][1].size()) {
                        aaa = collapse_patch[n1][0].size();
                    } else {
                        aaa = collapse_patch[n1][1].size();
                    }
                    size = size + aaa - 1;
                    if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                        end = j;
                        for (int k = begin; k < end + 1; ++k) {
                            if (sin_p[i][k] == 2) {
                                collapse_size[i].push_back(2);
                            } else {
                                collapse_size[i].push_back(size);
                            }
                        }
                        size = 0;
                        begin = j;
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != 2) {
                            sign = sin_p[i][j + 1];
                            begin = j + 1;
                        }
                    }

                } else {
                    if (size == 0) {
                        begin = j;
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != 2) {
                            sign = sin_p[i][j + 1];
                        } else if (j != sin_p[i].size() - 1) {
                            collapse_size[i].push_back(2);
                        }
                    }
                }
                if (j == sin_p[i].size() - 1) {
                    if (sin_p[i][j] == sign) {
                        end = j;
                        for (int k = begin; k < end + 1; ++k) {
                            collapse_size[i].push_back(size);
                        }
                    } else if (sin_p[i][j] == 2) {
                        collapse_size[i].push_back(2);
                    }
                }
            }
        } else {
            for (int j = relative_point[n].size() - 1; j > -1; --j) {
                int v1 = relative_point[n][j][0];
                int v2 = relative_point[n][j][1];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    sin_p[i].insert(sin_p[i].begin(), 0);
                } else if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                           all_singularity.end()) {
                    sin_p[i].insert(sin_p[i].begin(), 1);
                } else {
                    sin_p[i].insert(sin_p[i].begin(), 2);
                }
            }
            if (sin_p[i][sin_p[i].size() - 1] == 2) {
                int a1, a2;
                for (int j = sin_p[i].size() - 1; j > -1; --j) {
                    if (sin_p[i][j] == 2) continue;
                    a1 = j;
                    a2 = sin_p[i][j];
                    break;
                }
                for (int j = a1; j < sin_p[i].size(); ++j) {
                    sin_p[i][j] = a2;
                    break;
                }
            }
            bool spicy = false;
            int a1 = -1;
            for (int j = sin_p[i].size() - 1; j > -1; --j) {
                if (spicy) {
                    if (sin_p[i][j] != 2) {
                        if (sin_p[i][j] == sin_p[i][a1]) {
                            sin_p[i][j] = 2;
                            spicy = false;
                        } else {
                            for (int k = a1 - 1; k > j; --k) {
                                sin_p[i][k] = sin_p[i][j];
                            }
                            spicy = false;
                        }
                    }
                } else {
                    if (sin_p[i][j] == 2) {
                        spicy = true;
                        if (a1 == -1) {
                            a1 = j + 1;
                        }
                    }
                }
            }
            int begin = sin_p[i].size() - 1, end = 0;
            int size = 0, sign = sin_p[i][sin_p[i].size() - 1];
            for (int j = sin_p[i].size() - 1; j > -1; --j) {
                int n1 = layer[n][j];
                if (sin_p[i][j] == sign) {
                    int aaa;
                    if (collapse_patch[n1][0].size() > collapse_patch[n1][1].size()) {
                        aaa = collapse_patch[n1][0].size();
                    } else {
                        aaa = collapse_patch[n1][1].size();
                    }
                    size = size + aaa - 1;
                    if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                        end = j;
                        for (int k = begin; k > end - 1; --k) {
                            if (sin_p[i][k] == 2) {
                                collapse_size[i].insert(collapse_size[i].begin(), 2);
                            } else {
                                collapse_size[i].insert(collapse_size[i].begin(), size);
                            }
                        }
                        size = 0;
                        begin = j;
                        if (j != 0 && sin_p[i][j - 1] != 2) {
                            sign = sin_p[i][j - 1];
                            begin = j - 1;
                        }
                    }
                } else {
                    if (size == 0) {
                        begin = j;
                        if (j != 0 && sin_p[i][j - 1] != 2) {
                            sign = sin_p[i][j - 1];
                        } else if (j != 0) {
                            collapse_size[i].insert(collapse_size[i].begin(), 2);
                        }
                    }
                }
                if (j == 0) {
                    if (sin_p[i][j] == sign) {
                        end = j;
                        for (int k = begin; k > end - 1; --k) {
                            collapse_size[i].insert(collapse_size[i].begin(), size);
                        }
                    } else if (sin_p[i][j] == 2) {
                        collapse_size[i].insert(collapse_size[i].begin(), 2);
                    }
                }
            }
        }
    }

    for (int i = 0; i < position.size(); ++i) {
        if (need_collapse[i] == 0) {
            if (boundary_layer[i]) {
                int n = position[i];
                int ver1 = collapse_patch[v][i][0];
                int ver2 = collapse_patch[v][i][1];
                if (relative_point[n].size() == 1) {
                    std::vector<int> line1, line2;
                    int p = layer[n][0];
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == ver2 && collapse_patch[p][k][1] == ver1) {
                            line1 = collapse_patch[p][(k + 3) % 4];
                            line2 = collapse_patch[p][(k + 1) % 4];
                            break;
                        }
                    }
                    for (int k = 0; k < line1.size(); ++k) {
                        int o1 = line1[k];
                        int o2 = line2[line2.size() - 1 - k];
                        Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                        O_compact[o1] = p2;
                        O_compact[o2] = p2;
                        if (renew[o2] == -1) {
                            renew[o2] = o1;
                        }
                    }
                } else {
                    int a;
                    for (int j = 0; j < layer[n].size(); ++j) {
                        if (j == layer[n].size() - 1 && layer[n].size() != 1) {
                            v1 = relative_point[n][j][2];
                            v2 = relative_point[n][j][3];
                            if (x1 == v1) {
                                a = 1;
                                continue;
                            } else if (x1 == v2) {
                                a = 0;
                                continue;
                            }
                        }
                        v1 = relative_point[n][j][0];
                        v2 = relative_point[n][j][1];
                        if (x1 == v1) {
                            a = 0;
                            break;
                        } else if (x1 == v2) {
                            a = 1;
                            break;
                        }
                    }
                    for (int j = 0; j < layer[n].size(); ++j) {
                        int p = layer[n][j];
                        if (a == 0) {
                            v1 = relative_point[n][j][0];
                            v2 = relative_point[n][j][1];
                        } else {
                            v1 = relative_point[n][j][1];
                            v2 = relative_point[n][j][0];
                        }
                        std::vector<int> line1, line2;
                        for (int k = 0; k < collapse_patch[p].size(); ++k) {
                            if (collapse_patch[p][k][0] == v2 &&
                                collapse_patch[p][k][collapse_patch[p][k].size() - 1] != v1) {
                                line2 = collapse_patch[p][k];
                                line1 = collapse_patch[p][(k + 2) % 4];
                                break;
                            }
                            if (collapse_patch[p][k][0] == v1 &&
                                collapse_patch[p][k][collapse_patch[p][k].size() - 1] != v2) {
                                line1 = collapse_patch[p][k];
                                line2 = collapse_patch[p][(k + 2) % 4];
                                break;
                            }
                        }
                        for (int k = 0; k < line1.size(); ++k) {
                            int o1 = line1[k];
                            int o2 = line2[line2.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                    }
                }

            } else {
                int n = position[i];
                int ver1 = collapse_patch[v][i][0];
                int ver2 = collapse_patch[v][i][1];
                if (valence[ver1] == 3 || valence[ver2] == 3) {
                    int sin, no_sin;
                    if (valence[ver1] == 3) {
                        sin = ver1;
                        no_sin = ver2;
                    }
                    if (valence[ver2] == 3) {
                        sin = ver2;
                        no_sin = ver1;
                    }
                    int a = -1;
                    if (relative_point[n].size() != 1) {
                        if (relative_point[n][0][0] == sin ||
                            relative_point[n][relative_point[n].size() - 1][3] == sin) {
                            a = 0;
                        } else {
                            a = 1;
                        }
                    }
                    if (a == -1) {
                        std::vector<int> line1, line2;
                        int p = layer[n][0];
                        for (int k = 0; k < collapse_patch[p].size(); ++k) {
                            if (collapse_patch[p][k][0] == sin &&
                                collapse_patch[p][k][1] == no_sin) {
                                line1 = collapse_patch[p][(k + 3) % 4];
                                line2 = collapse_patch[p][(k + 1) % 4];
                                break;
                            }
                            if (collapse_patch[p][k][0] == no_sin &&
                                collapse_patch[p][k][1] == sin) {
                                line1 = collapse_patch[p][(k + 1) % 4];
                                line2 = collapse_patch[p][(k + 3) % 4];
                                break;
                            }
                        }
                        for (int k = 0; k < line1.size(); ++k) {
                            int o1 = line1[k];
                            int o2 = line2[line2.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                    } else {
                        for (int j = 0; j < layer[n].size(); ++j) {
                            int p = layer[n][j];
                            v1 = relative_point[n][j][a];
                            v2 = relative_point[n][j][(a + 1) % 2];
                            std::vector<int> line1, line2;
                            for (int k = 0; k < collapse_patch[p].size(); ++k) {
                                if (collapse_patch[p][k][0] == v2 &&
                                    collapse_patch[p][k][collapse_patch[p][k].size() - 1] != v1) {
                                    if (a == 0) {
                                        line2 = collapse_patch[p][k];
                                        line1 = collapse_patch[p][(k + 2) % 4];
                                    } else {
                                        line1 = collapse_patch[p][k];
                                        line2 = collapse_patch[p][(k + 2) % 4];
                                    }
                                    break;
                                } else if (collapse_patch[p][k][0] == v1 &&
                                           collapse_patch[p][k][collapse_patch[p][k].size() - 1] !=
                                               v2) {
                                    if (a == 1) {
                                        line1 = collapse_patch[p][k];
                                        line2 = collapse_patch[p][(k + 2) % 4];
                                    } else {
                                        line2 = collapse_patch[p][k];
                                        line1 = collapse_patch[p][(k + 2) % 4];
                                    }
                                    break;
                                }
                            }
                            for (int k = 0; k < line1.size(); ++k) {
                                int o1 = line1[k];
                                int o2 = line2[line2.size() - 1 - k];
                                O_compact[o2] = O_compact[o1];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        }
                    }
                }
            }
        } else {
            int n = position[i];
            int ver1 = collapse_patch[v][i][0];
            int ver2 = collapse_patch[v][i][1];
            if (relative_point[n][0][1] == ver1) {
                bool begin = true;
                int b1 = 0;
                int b2 = 0;
                int c1 = sin_p[i][0];
                for (int j = 0; j < layer[n].size(); ++j) {
                    int p = layer[n][j];
                    v1 = relative_point[n][j][0];
                    v2 = relative_point[n][j][1];
                    std::vector<int> line1, line2;
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == v2) {
                            line2 = collapse_patch[p][k];
                            line1 = collapse_patch[p][(k + 2) % 4];
                            break;
                        }
                    }
                    if (valence[v2] == 3 && begin && sin_p[i][j]) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (valence[v1] == 3 && begin && sin_p[i][j] == 0) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o1] = O_compact[o2];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (begin && sin_p[i][j]) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 - b1 * d / 2 / size;
                            O_compact[o2] = p2 - b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                            if (sin_p[i][j + 1] != 2) {
                                c1 = sin_p[i][j + 1];
                            }
                            begin = false;
                        }
                    } else if (begin && sin_p[i][j] == 0) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;

                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 + b1 * d / 2 / size;
                            O_compact[o2] = p2 + b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                            if (sin_p[i][j + 1] != 2) {
                                c1 = sin_p[i][j + 1];
                            }
                            begin = false;
                        }
                    } else {
                        if (sin_p[i][j] != 2) {
                            c1 = sin_p[i][j];
                        }
                        if (c1 == 0 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o1] = O_compact[o2];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else if (c1 == 1 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o2] = O_compact[o1];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else {
                            if (c1 == 0) {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o2] - O_compact[o1];
                                    O_compact[o1] = O_compact[o1] + d * b1 / size;
                                    O_compact[o2] = O_compact[o1];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            } else {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o1] - O_compact[o2];
                                    O_compact[o2] = O_compact[o2] + d * b1 / size;
                                    O_compact[o1] = O_compact[o2];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            }
                        }
                    }
                }
            } else {
                /////////////////////暂时未修改///////////////////////////////////////////////////////////5944-6065
                bool begin = true;
                int b1 = 0;
                int b2 = 0;
                int c1 = sin_p[i][sin_p[i].size() - 1];
                for (int j = layer[n].size() - 1; j > -1; --j) {
                    int p = layer[n][j];
                    v1 = relative_point[n][j][2];
                    v2 = relative_point[n][j][3];
                    std::vector<int> line1, line2;
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == v2) {
                            line2 = collapse_patch[p][k];
                            line1 = collapse_patch[p][(k + 2) % 4];
                            break;
                        }
                    }
                    if (valence[v2] == 3 && begin && sin_p[i][j]) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (valence[v1] == 3 && begin && sin_p[i][j] == 0) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o1] = O_compact[o2];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (begin && sin_p[i][j]) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 + b1 * d / 2 / size;
                            O_compact[o2] = p2 + b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                            c1 = sin_p[i][j - 1];
                            begin = false;
                        }
                    } else if (begin && sin_p[i][j] == 0) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 - b1 * d / 2 / size;
                            O_compact[o2] = p2 - b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                            c1 = sin_p[i][j - 1];
                            begin = false;
                        }
                    } else {
                        if (sin_p[i][j] != 2) {
                            c1 = sin_p[i][j];
                        }
                        if (c1 == 0 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o2] = O_compact[o1];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else if (c1 == 1 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o1] = O_compact[o2];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else {
                            if (c1 == 0) {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o1] - O_compact[o2];
                                    O_compact[o2] = O_compact[o2] + d * b1 / size;
                                    O_compact[o1] = O_compact[o2];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            } else {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o2] - O_compact[o1];
                                    O_compact[o1] = O_compact[o1] + d * b1 / size;
                                    O_compact[o2] = O_compact[o1];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            }
                        }
                    }
                }
                /////////////////////暂时未修改///////////////////////////////////////////////////////////5944-6065
            }
        }
    }

    std::sort(delete_patch.begin(), delete_patch.end());
    // update information
    //  update valence
    valence[x1] = 4;
    valence[x2] = 0;
    valence[x3] = 0;
    valence[x4] = 0;
    // update singularity
    std::vector<int> temp_singularity;
    for (int i = 0; i < all_singularity.size(); ++i) {
        if (find(singularity.begin(), singularity.end(), all_singularity[i]) != singularity.end())
            continue;
        temp_singularity.push_back(all_singularity[i]);
    }
    all_singularity = temp_singularity;
    int singularity_point;
    for (int i = 0; i < twice.size(); ++i) {
        int u = twice[i];
        for (int j = 0; j < collapse_patch[u].size(); ++j) {
            if (find(all_singularity.begin(), all_singularity.end(), collapse_patch[u][j][0]) !=
                all_singularity.end()) {
                singularity_point = collapse_patch[u][j][0];
            }
        }
        for (int j = 0; j < collapse_patch[u].size(); ++j) {
            if (collapse_patch[u][j][0] == singularity_point) continue;
            renew[collapse_patch[u][j][0]] = singularity_point;
        }
    }
    for (int i = 0; i < renew.size(); ++i) {
        if (renew[i] == -1) continue;
        if (find(all_singularity.begin(), all_singularity.end(), i) != all_singularity.end()) {
            if (boundary_o[i]) {
                valence[renew[i]] = valence[i] - 1;
            } else {
                valence[renew[i]] = valence[i];
            }
            for (int j = 0; j < all_singularity.size(); ++j) {
                if (all_singularity[j] == i) {
                    all_singularity[j] = renew[i];
                }
            }
        }
        int v = renew[i];
        if (boundary_o[v]) {
            boundary_o[i] = 1;
        }
    }
    // update collapse_patch
    std::vector<std::vector<std::vector<int>>> temp_collapse_patch;
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < collapse_patch[i].size(); ++j) {
            for (int k = 0; k < collapse_patch[i][j].size(); ++k) {
                int v = collapse_patch[i][j][k];
                if (renew[v] != -1) {
                    collapse_patch[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_collapse_patch.push_back(collapse_patch[i]);
    }
    collapse_patch = temp_collapse_patch;
    // update patch_compact
    std::vector<std::vector<std::vector<int>>> temp_patch_compact;
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int v = patch_compact[i][j][k];
                if (renew[v] != -1) {
                    patch_compact[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_patch_compact.push_back(patch_compact[i]);
    }
    patch_compact = temp_patch_compact;
    // update relative_point
    std::vector<std::vector<std::vector<int>>> temp_relative_point(relative_point.size());
    std::vector<std::vector<int>> temp_layer(layer.size());
    for (int i = 0; i < relative_point.size(); ++i) {
        for (int j = 0; j < relative_point[i].size(); ++j) {
            for (int k = 0; k < relative_point[i][j].size(); ++k) {
                int v = relative_point[i][j][k];
                if (renew[v] != -1) {
                    relative_point[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < layer.size(); ++i) {
        if (find(position.begin(), position.end(), i) != position.end()) continue;
        for (int j = 0; j < layer[i].size(); ++j) {
            int v = layer[i][j];
            if (find(delete_patch.begin(), delete_patch.end(), v) != delete_patch.end()) continue;
            temp_relative_point[i].push_back(relative_point[i][j]);
            for (int k = 0; k < delete_patch.size(); ++k) {
                if (v < delete_patch[0]) {
                    temp_layer[i].push_back(v);
                    break;
                }
                if (k == delete_patch.size() - 1) {
                    if (v > delete_patch[k]) {
                        temp_layer[i].push_back(v - k - 1);
                    }
                } else {
                    if (v > delete_patch[k] && v < delete_patch[k + 1]) {
                        temp_layer[i].push_back(v - k - 1);
                        break;
                    }
                }
            }
        }
    }
    for (int i = 0; i < thr_for_sin_patch.size(); ++i) {
        int n = thr_for_sin_patch[i];
        for (int k = 0; k < delete_patch.size(); ++k) {
            if (k == delete_patch.size() - 1) {
                if (n > delete_patch[k]) {
                    thr_for_sin_patch[i] = n - k - 1;
                }
            } else {
                if (n > delete_patch[k] && n < delete_patch[k + 1]) {
                    thr_for_sin_patch[i] = n - k - 1;
                    break;
                }
            }
        }
    }
    std::vector<std::vector<std::vector<int>>> temp;
    std::vector<std::vector<int>> temp1;
    for (int i = 0; i < temp_relative_point.size(); ++i) {
        if (temp_relative_point[i].size() == 0) continue;
        temp.push_back(temp_relative_point[i]);
    }
    relative_point = temp;
    for (int i = 0; i < temp_layer.size(); ++i) {
        if (temp_layer[i].size() == 0) continue;
        temp1.push_back(temp_layer[i]);
    }
    layer = temp1;
}
void Optimizer::sin_fft(std::vector<Vector3d>& O_compact,
                        std::vector<std::vector<std::vector<int>>>& patch_compact,
                        std::vector<int>& valence,
                        std::vector<std::vector<std::vector<int>>>& relative_point,
                        std::vector<std::vector<int>>& layer,
                        std::vector<std::vector<std::vector<int>>>& collapse_patch, int v,
                        std::vector<int>& all_singularity, std::vector<int>& thr_for_sin_patch) {
    std::vector<int> position;
    std::vector<int> renew(O_compact.size(), -1);
    std::vector<int> delete_patch;
    for (int i = 0; i < collapse_patch[v].size(); ++i) {
        int v1 = collapse_patch[v][i][0];
        int v2 = collapse_patch[v][i][1];
        for (int j = 0; j < relative_point.size(); ++j) {
            if (relative_point[j].size() == 1) {
                int n1 = relative_point[j][0][0];
                int n2 = relative_point[j][0][1];
                int n3 = relative_point[j][0][2];
                int n4 = relative_point[j][0][3];
                if (n1 == v2 && n2 == v1) {
                    position.push_back(j);
                    break;
                } else if (n2 == v2 && n3 == v1) {
                    position.push_back(j);
                    break;
                } else if (n3 == v2 && n4 == v1) {
                    position.push_back(j);
                    break;
                } else if (n4 == v2 && n1 == v1) {
                    position.push_back(j);
                    break;
                }
            } else {
                int n1 = relative_point[j][0][0];
                int n2 = relative_point[j][0][1];
                int n3 = relative_point[j][relative_point[j].size() - 1][2];
                int n4 = relative_point[j][relative_point[j].size() - 1][3];
                if (v1 == n1 && v2 == n2) {
                    position.push_back(j);
                    break;
                } else if (v1 == n2 && v2 == n1) {
                    position.push_back(j);
                    break;
                } else if (v1 == n3 && v2 == n4) {
                    position.push_back(j);
                    break;
                } else if (v1 == n4 && v2 == n3) {
                    position.push_back(j);
                    break;
                }
            }
        }
    }
    std::vector<int> singularity;
    std::vector<int> quad;
    int v1, v2, v3, v4;
    int x1, x2, x3, x4;
    v1 = collapse_patch[v][0][0];
    v2 = collapse_patch[v][1][0];
    v3 = collapse_patch[v][2][0];
    v4 = collapse_patch[v][3][0];
    x1 = v1;
    x2 = v2;
    x3 = v3;
    x4 = v4;
    renew[v2] = v1;
    renew[v3] = v1;
    renew[v4] = v1;
    singularity.push_back(v2);
    singularity.push_back(v3);
    singularity.push_back(v4);
    quad.push_back(v1);
    quad.push_back(v2);
    quad.push_back(v3);
    quad.push_back(v4);
    Vector3d p1 = (O_compact[v1] + O_compact[v2] + O_compact[v3] + O_compact[v4]) / 4;
    O_compact[v1] = p1;
    O_compact[v2] = p1;
    O_compact[v3] = p1;
    O_compact[v4] = p1;
    std::vector<int> twice;
    for (int i = 0; i < position.size(); ++i) {
        int n = position[i];
        for (int j = 0; j < layer[n].size(); ++j) {
            if (find(delete_patch.begin(), delete_patch.end(), layer[n][j]) ==
                delete_patch.end()) {
                delete_patch.push_back(layer[n][j]);
            } else {
                twice.push_back(layer[n][j]);
            }
        }
    }
    delete_patch.push_back(v);
    std::vector<int> need_collapse(position.size(), 0);
    std::vector<std::vector<int>> corr(position.size());
    for (int i = 0; i < position.size(); ++i) {
        int n = position[i];
        if (layer[n].size() == 1) continue;
        for (int j = 0; j < layer[n].size(); ++j) {
            for (int k = 0; k < 4; ++k) {
                int v1 = relative_point[n][j][k];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                        all_singularity.end() &&
                    find(quad.begin(), quad.end(), v1) == quad.end()) {
                    need_collapse[i] = 1;
                }
            }
        }
    }
    std::vector<std::vector<int>> sin_p(position.size());
    std::vector<std::vector<int>> collapse_size(position.size());
    for (int i = 0; i < position.size(); ++i) {
        if (need_collapse[i] == 0) continue;
        int n = position[i];
        int ver1 = collapse_patch[v][i][0];
        int ver2 = collapse_patch[v][i][1];
        if (relative_point[n][0][1] == ver1) {
            for (int j = 0; j < relative_point[n].size(); ++j) {
                int v1 = relative_point[n][j][2];
                int v2 = relative_point[n][j][3];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    sin_p[i].push_back(1);
                } else if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                           all_singularity.end()) {
                    sin_p[i].push_back(0);
                } else {
                    sin_p[i].push_back(2);
                }
            }
            if (sin_p[i][0] == 2) {
                int a1, a2;
                for (int j = 0; j < sin_p[i][0]; ++j) {
                    if (sin_p[i][j] == 2) continue;
                    a1 = j;
                    a2 = sin_p[i][j];
                    break;
                }
                for (int j = 0; j < a1; ++j) {
                    sin_p[i][j] = a2;
                    break;
                }
            }
            bool spicy = false;
            int a1 = -1;
            for (int j = 0; j < sin_p[i].size(); ++j) {
                if (spicy) {
                    if (sin_p[i][j] != 2) {
                        if (sin_p[i][j] == sin_p[i][a1]) {
                            sin_p[i][j] = 2;
                            spicy = false;
                        } else {
                            for (int k = a1 + 1; k < j; ++k) {
                                sin_p[i][k] = sin_p[i][j];
                            }
                            spicy = false;
                        }
                    }
                } else {
                    if (sin_p[i][j] == 2) {
                        spicy = true;
                        if (a1 == -1) {
                            a1 = j - 1;
                        }
                    }
                }
            }
            int begin = 0, end = sin_p[i].size() - 1;
            int size = 0, sign = sin_p[i][0];
            for (int j = 0; j < sin_p[i].size(); ++j) {
                int n1 = layer[n][j];
                if (sin_p[i][j] == sign) {
                    int aaa;
                    if (collapse_patch[n1][0].size() > collapse_patch[n1][1].size()) {
                        aaa = collapse_patch[n1][0].size();
                    } else {
                        aaa = collapse_patch[n1][1].size();
                    }
                    size = size + aaa - 1;
                    if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                        end = j;
                        for (int k = begin; k < end + 1; ++k) {
                            if (sin_p[i][k] == 2) {
                                collapse_size[i].push_back(2);
                            } else {
                                collapse_size[i].push_back(size);
                            }
                        }
                        size = 0;
                        begin = j;
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != 2) {
                            sign = sin_p[i][j + 1];
                            begin = j + 1;
                        }
                    }

                } else {
                    if (size == 0) {
                        begin = j;
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != 2) {
                            sign = sin_p[i][j + 1];
                        } else if (j != sin_p[i].size() - 1) {
                            collapse_size[i].push_back(2);
                        }
                    }
                }
                if (j == sin_p[i].size() - 1) {
                    if (sin_p[i][j] == sign) {
                        end = j;
                        for (int k = begin; k < end + 1; ++k) {
                            collapse_size[i].push_back(size);
                        }
                    } else if (sin_p[i][j] == 2) {
                        collapse_size[i].push_back(2);
                    }
                }
            }
        } else {
            for (int j = relative_point[n].size() - 1; j > -1; --j) {
                int v1 = relative_point[n][j][0];
                int v2 = relative_point[n][j][1];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    sin_p[i].insert(sin_p[i].begin(), 0);
                } else if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                           all_singularity.end()) {
                    sin_p[i].insert(sin_p[i].begin(), 1);
                } else {
                    sin_p[i].insert(sin_p[i].begin(), 2);
                }
            }
            if (sin_p[i][sin_p[i].size() - 1] == 2) {
                int a1, a2;
                for (int j = sin_p[i].size() - 1; j > -1; --j) {
                    if (sin_p[i][j] == 2) continue;
                    a1 = j;
                    a2 = sin_p[i][j];
                    break;
                }
                for (int j = a1; j < sin_p[i].size(); ++j) {
                    sin_p[i][j] = a2;
                    break;
                }
            }
            bool spicy = false;
            int a1 = -1;
            for (int j = sin_p[i].size() - 1; j > -1; --j) {
                if (spicy) {
                    if (sin_p[i][j] != 2) {
                        if (sin_p[i][j] == sin_p[i][a1]) {
                            sin_p[i][j] = 2;
                            spicy = false;
                        } else {
                            for (int k = a1 - 1; k > j; --k) {
                                sin_p[i][k] = sin_p[i][j];
                            }
                            spicy = false;
                        }
                    }
                } else {
                    if (sin_p[i][j] == 2) {
                        spicy = true;
                        if (a1 == -1) {
                            a1 = j + 1;
                        }
                    }
                }
            }
            int begin = sin_p[i].size() - 1, end = 0;
            int size = 0, sign = sin_p[i][sin_p[i].size() - 1];
            for (int j = sin_p[i].size() - 1; j > -1; --j) {
                int n1 = layer[n][j];
                if (sin_p[i][j] == sign) {
                    int aaa;
                    if (collapse_patch[n1][0].size() > collapse_patch[n1][1].size()) {
                        aaa = collapse_patch[n1][0].size();
                    } else {
                        aaa = collapse_patch[n1][1].size();
                    }
                    size = size + aaa - 1;
                    if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                        end = j;
                        for (int k = begin; k > end - 1; --k) {
                            if (sin_p[i][k] == 2) {
                                collapse_size[i].insert(collapse_size[i].begin(), 2);
                            } else {
                                collapse_size[i].insert(collapse_size[i].begin(), size);
                            }
                        }
                        size = 0;
                        begin = j;
                        if (j != 0 && sin_p[i][j - 1] != 2) {
                            sign = sin_p[i][j - 1];
                            begin = j - 1;
                        }
                    }
                } else {
                    if (size == 0) {
                        begin = j;
                        if (j != 0 && sin_p[i][j - 1] != 2) {
                            sign = sin_p[i][j - 1];
                        } else if (j != 0) {
                            collapse_size[i].insert(collapse_size[i].begin(), 2);
                        }
                    }
                }
                if (j == 0) {
                    if (sin_p[i][j] == sign) {
                        end = j;
                        for (int k = begin; k > end - 1; --k) {
                            collapse_size[i].insert(collapse_size[i].begin(), size);
                        }
                    } else if (sin_p[i][j] == 2) {
                        collapse_size[i].insert(collapse_size[i].begin(), 2);
                    }
                }
            }
        }
    }

    for (int i = 0; i < position.size(); ++i) {
        if (need_collapse[i] == 0) {
            int n = position[i];
            int ver1 = collapse_patch[v][i][0];
            int ver2 = collapse_patch[v][i][1];
            if (valence[ver1] == 3 || valence[ver2] == 3) {
                int sin, no_sin;
                if (valence[ver1] == 3) {
                    sin = ver1;
                    no_sin = ver2;
                }
                if (valence[ver2] == 3) {
                    sin = ver2;
                    no_sin = ver1;
                }
                int a = -1;
                if (relative_point[n].size() != 1) {
                    if (relative_point[n][0][0] == sin ||
                        relative_point[n][relative_point[n].size() - 1][3] == sin) {
                        a = 0;
                    } else {
                        a = 1;
                    }
                }
                if (a == -1) {
                    std::vector<int> line1, line2;
                    int p = layer[n][0];
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == sin && collapse_patch[p][k][1] == no_sin) {
                            line1 = collapse_patch[p][(k + 3) % 4];
                            line2 = collapse_patch[p][(k + 1) % 4];
                            break;
                        }
                        if (collapse_patch[p][k][0] == no_sin && collapse_patch[p][k][1] == sin) {
                            line1 = collapse_patch[p][(k + 1) % 4];
                            line2 = collapse_patch[p][(k + 3) % 4];
                            break;
                        }
                    }
                    for (int k = 0; k < line1.size(); ++k) {
                        int o1 = line1[k];
                        int o2 = line2[line2.size() - 1 - k];
                        O_compact[o2] = O_compact[o1];
                        if (renew[o2] == -1) {
                            renew[o2] = o1;
                        }
                    }
                } else {
                    for (int j = 0; j < layer[n].size(); ++j) {
                        int p = layer[n][j];
                        v1 = relative_point[n][j][a];
                        v2 = relative_point[n][j][(a + 1) % 2];
                        std::vector<int> line1, line2;
                        for (int k = 0; k < collapse_patch[p].size(); ++k) {
                            if (collapse_patch[p][k][0] == v2 &&
                                collapse_patch[p][k][collapse_patch[p][k].size() - 1] != v1) {
                                if (a == 0) {
                                    line2 = collapse_patch[p][k];
                                    line1 = collapse_patch[p][(k + 2) % 4];
                                } else {
                                    line1 = collapse_patch[p][k];
                                    line2 = collapse_patch[p][(k + 2) % 4];
                                }
                                break;
                            } else if (collapse_patch[p][k][0] == v1 &&
                                       collapse_patch[p][k][collapse_patch[p][k].size() - 1] !=
                                           v2) {
                                if (a == 1) {
                                    line1 = collapse_patch[p][k];
                                    line2 = collapse_patch[p][(k + 2) % 4];
                                } else {
                                    line2 = collapse_patch[p][k];
                                    line1 = collapse_patch[p][(k + 2) % 4];
                                }
                                break;
                            }
                        }
                        for (int k = 0; k < line1.size(); ++k) {
                            int o1 = line1[k];
                            int o2 = line2[line2.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                    }
                }
            } else {
                if (relative_point[n].size() == 1) {
                    std::vector<int> line1, line2;
                    int p = layer[n][0];
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == ver2 && collapse_patch[p][k][1] == ver1) {
                            line1 = collapse_patch[p][(k + 3) % 4];
                            line2 = collapse_patch[p][(k + 1) % 4];
                            break;
                        }
                    }
                    for (int k = 0; k < line1.size(); ++k) {
                        int o1 = line1[k];
                        int o2 = line2[line2.size() - 1 - k];
                        Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                        O_compact[o1] = p2;
                        O_compact[o2] = p2;
                        if (renew[o2] == -1) {
                            renew[o2] = o1;
                        }
                    }
                } else {
                    for (int j = 0; j < layer[n].size(); ++j) {
                        int p = layer[n][j];
                        v1 = relative_point[n][j][0];
                        v2 = relative_point[n][j][1];
                        std::vector<int> line1, line2;
                        for (int k = 0; k < collapse_patch[p].size(); ++k) {
                            if (collapse_patch[p][k][0] == v2) {
                                line2 = collapse_patch[p][k];
                                line1 = collapse_patch[p][(k + 2) % 4];
                                break;
                            }
                        }
                        for (int k = 0; k < line1.size(); ++k) {
                            int o1 = line1[k];
                            int o2 = line2[line2.size() - 1 - k];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2;
                            O_compact[o2] = p2;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                    }
                }
            }
        } else {
            int n = position[i];
            int ver1 = collapse_patch[v][i][0];
            int ver2 = collapse_patch[v][i][1];
            if (relative_point[n][0][1] == ver1) {
                bool begin = true;
                int b1 = 0;
                int b2 = 0;
                int c1 = sin_p[i][0];
                for (int j = 0; j < layer[n].size(); ++j) {
                    int p = layer[n][j];
                    v1 = relative_point[n][j][0];
                    v2 = relative_point[n][j][1];
                    std::vector<int> line1, line2;
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == v2) {
                            line2 = collapse_patch[p][k];
                            line1 = collapse_patch[p][(k + 2) % 4];
                            break;
                        }
                    }
                    if (valence[v2] == 3 && begin && sin_p[i][j]) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (valence[v1] == 3 && begin && sin_p[i][j] == 0) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o1] = O_compact[o2];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (begin && sin_p[i][j]) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 - b1 * d / 2 / size;
                            O_compact[o2] = p2 - b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                            if (sin_p[i][j + 1] != 2) {
                                c1 = sin_p[i][j + 1];
                            }
                            begin = false;
                        }
                    } else if (begin && sin_p[i][j] == 0) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;

                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 + b1 * d / 2 / size;
                            O_compact[o2] = p2 + b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                            if (sin_p[i][j + 1] != 2) {
                                c1 = sin_p[i][j + 1];
                            }
                            begin = false;
                        }
                    } else {
                        if (sin_p[i][j] != 2) {
                            c1 = sin_p[i][j];
                        }
                        if (c1 == 0 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o1] = O_compact[o2];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else if (c1 == 1 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o2] = O_compact[o1];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else {
                            if (c1 == 0) {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o2] - O_compact[o1];
                                    O_compact[o1] = O_compact[o1] + d * b1 / size;
                                    O_compact[o2] = O_compact[o1];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            } else {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o1] - O_compact[o2];
                                    O_compact[o2] = O_compact[o2] + d * b1 / size;
                                    O_compact[o1] = O_compact[o2];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            }
                        }
                    }
                }
            } else {
                /////////////////////暂时未修改///////////////////////////////////////////////////////////5944-6065
                bool begin = true;
                int b1 = 0;
                int b2 = 0;
                int c1 = sin_p[i][sin_p[i].size() - 1];
                for (int j = layer[n].size() - 1; j > -1; --j) {
                    int p = layer[n][j];
                    v1 = relative_point[n][j][2];
                    v2 = relative_point[n][j][3];
                    std::vector<int> line1, line2;
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == v2) {
                            line2 = collapse_patch[p][k];
                            line1 = collapse_patch[p][(k + 2) % 4];
                            break;
                        }
                    }
                    if (valence[v2] == 3 && begin && sin_p[i][j]) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (valence[v1] == 3 && begin && sin_p[i][j] == 0) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o1] = O_compact[o2];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (begin && sin_p[i][j]) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 + b1 * d / 2 / size;
                            O_compact[o2] = p2 + b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                            c1 = sin_p[i][j - 1];
                            begin = false;
                        }
                    } else if (begin && sin_p[i][j] == 0) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 - b1 * d / 2 / size;
                            O_compact[o2] = p2 - b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                            c1 = sin_p[i][j - 1];
                            begin = false;
                        }
                    } else {
                        if (sin_p[i][j] != 2) {
                            c1 = sin_p[i][j];
                        }
                        if (c1 == 0 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o2] = O_compact[o1];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else if (c1 == 1 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o1] = O_compact[o2];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else {
                            if (c1 == 0) {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o1] - O_compact[o2];
                                    O_compact[o2] = O_compact[o2] + d * b1 / size;
                                    O_compact[o1] = O_compact[o2];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            } else {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o2] - O_compact[o1];
                                    O_compact[o1] = O_compact[o1] + d * b1 / size;
                                    O_compact[o2] = O_compact[o1];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            }
                        }
                    }
                }
                /////////////////////暂时未修改///////////////////////////////////////////////////////////5944-6065
            }
        }
    }

    std::sort(delete_patch.begin(), delete_patch.end());
    // update information
    //  update valence
    valence[x1] = 5;
    valence[x2] = 0;
    valence[x3] = 0;
    valence[x4] = 0;
    // update singularity
    std::vector<int> temp_singularity;
    for (int i = 0; i < all_singularity.size(); ++i) {
        if (find(singularity.begin(), singularity.end(), all_singularity[i]) != singularity.end())
            continue;
        temp_singularity.push_back(all_singularity[i]);
    }
    all_singularity = temp_singularity;
    int singularity_point;
    for (int i = 0; i < twice.size(); ++i) {
        int u = twice[i];
        for (int j = 0; j < collapse_patch[u].size(); ++j) {
            if (find(all_singularity.begin(), all_singularity.end(), collapse_patch[u][j][0]) !=
                all_singularity.end()) {
                singularity_point = collapse_patch[u][j][0];
            }
        }
        for (int j = 0; j < collapse_patch[u].size(); ++j) {
            if (collapse_patch[u][j][0] == singularity_point) continue;
            renew[collapse_patch[u][j][0]] = singularity_point;
        }
    }
    for (int i = 0; i < renew.size(); ++i) {
        if (renew[i] == -1) continue;
        if (find(all_singularity.begin(), all_singularity.end(), i) != all_singularity.end()) {
            std::cout << -1 << "\n";
            valence[renew[i]] = valence[i];
            for (int j = 0; j < all_singularity.size(); ++j) {
                if (all_singularity[j] == i) {
                    all_singularity[j] = renew[i];
                }
            }
        }
    }
    // update collapse_patch
    std::vector<std::vector<std::vector<int>>> temp_collapse_patch;
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < collapse_patch[i].size(); ++j) {
            for (int k = 0; k < collapse_patch[i][j].size(); ++k) {
                int v = collapse_patch[i][j][k];
                if (renew[v] != -1) {
                    collapse_patch[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_collapse_patch.push_back(collapse_patch[i]);
    }
    collapse_patch = temp_collapse_patch;
    // update patch_compact
    std::vector<std::vector<std::vector<int>>> temp_patch_compact;
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int v = patch_compact[i][j][k];
                if (renew[v] != -1) {
                    patch_compact[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_patch_compact.push_back(patch_compact[i]);
    }
    patch_compact = temp_patch_compact;
    // update relative_point
    std::vector<std::vector<std::vector<int>>> temp_relative_point(relative_point.size());
    std::vector<std::vector<int>> temp_layer(layer.size());
    for (int i = 0; i < relative_point.size(); ++i) {
        for (int j = 0; j < relative_point[i].size(); ++j) {
            for (int k = 0; k < relative_point[i][j].size(); ++k) {
                int v = relative_point[i][j][k];
                if (renew[v] != -1) {
                    relative_point[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < layer.size(); ++i) {
        if (find(position.begin(), position.end(), i) != position.end()) continue;
        for (int j = 0; j < layer[i].size(); ++j) {
            int v = layer[i][j];
            if (find(delete_patch.begin(), delete_patch.end(), v) != delete_patch.end()) continue;
            temp_relative_point[i].push_back(relative_point[i][j]);
            for (int k = 0; k < delete_patch.size(); ++k) {
                if (v < delete_patch[0]) {
                    temp_layer[i].push_back(v);
                    break;
                }
                if (k == delete_patch.size() - 1) {
                    if (v > delete_patch[k]) {
                        temp_layer[i].push_back(v - k - 1);
                    }
                } else {
                    if (v > delete_patch[k] && v < delete_patch[k + 1]) {
                        temp_layer[i].push_back(v - k - 1);
                        break;
                    }
                }
            }
        }
    }
    for (int i = 0; i < thr_for_sin_patch.size(); ++i) {
        int n = thr_for_sin_patch[i];
        for (int k = 0; k < delete_patch.size(); ++k) {
            if (k == delete_patch.size() - 1) {
                if (n > delete_patch[k]) {
                    thr_for_sin_patch[i] = n - k - 1;
                }
            } else {
                if (n > delete_patch[k] && n < delete_patch[k + 1]) {
                    thr_for_sin_patch[i] = n - k - 1;
                    break;
                }
            }
        }
    }
    std::vector<std::vector<std::vector<int>>> temp;
    std::vector<std::vector<int>> temp1;
    for (int i = 0; i < temp_relative_point.size(); ++i) {
        if (temp_relative_point[i].size() == 0) continue;
        temp.push_back(temp_relative_point[i]);
    }
    relative_point = temp;
    for (int i = 0; i < temp_layer.size(); ++i) {
        if (temp_layer[i].size() == 0) continue;
        temp1.push_back(temp_layer[i]);
    }
    layer = temp1;
}
void Optimizer::sin_ftt(std::vector<Vector3d>& O_compact,
                        std::vector<std::vector<std::vector<int>>>& patch_compact,
                        std::vector<int>& valence,
                        std::vector<std::vector<std::vector<int>>>& relative_point,
                        std::vector<std::vector<int>>& layer,
                        std::vector<std::vector<std::vector<int>>>& collapse_patch, int v,
                        std::vector<int>& all_singularity, std::vector<int>& thr_for_sin_patch) {
    std::vector<int> position;
    std::vector<int> renew(O_compact.size(), -1);
    std::vector<int> delete_patch;
    for (int i = 0; i < collapse_patch[v].size(); ++i) {
        int v1 = collapse_patch[v][i][0];
        int v2 = collapse_patch[v][i][1];
        for (int j = 0; j < relative_point.size(); ++j) {
            if (relative_point[j].size() == 1) {
                int n1 = relative_point[j][0][0];
                int n2 = relative_point[j][0][1];
                int n3 = relative_point[j][0][2];
                int n4 = relative_point[j][0][3];
                if (n1 == v2 && n2 == v1) {
                    position.push_back(j);
                    break;
                } else if (n2 == v2 && n3 == v1) {
                    position.push_back(j);
                    break;
                } else if (n3 == v2 && n4 == v1) {
                    position.push_back(j);
                    break;
                } else if (n4 == v2 && n1 == v1) {
                    position.push_back(j);
                    break;
                }
            } else {
                int n1 = relative_point[j][0][0];
                int n2 = relative_point[j][0][1];
                int n3 = relative_point[j][relative_point[j].size() - 1][2];
                int n4 = relative_point[j][relative_point[j].size() - 1][3];
                if (v1 == n1 && v2 == n2) {
                    position.push_back(j);
                    break;
                } else if (v1 == n2 && v2 == n1) {
                    position.push_back(j);
                    break;
                } else if (v1 == n3 && v2 == n4) {
                    position.push_back(j);
                    break;
                } else if (v1 == n4 && v2 == n3) {
                    position.push_back(j);
                    break;
                }
            }
        }
    }
    std::vector<int> singularity;
    std::vector<int> quad;
    int v1, v2, v3, v4;
    int x1, x2, x3, x4;
    v1 = collapse_patch[v][0][0];
    v2 = collapse_patch[v][1][0];
    v3 = collapse_patch[v][2][0];
    v4 = collapse_patch[v][3][0];
    x1 = v1;
    x2 = v2;
    x3 = v3;
    x4 = v4;
    renew[v2] = v1;
    renew[v3] = v1;
    renew[v4] = v1;
    singularity.push_back(v2);
    singularity.push_back(v3);
    singularity.push_back(v4);
    quad.push_back(v1);
    quad.push_back(v2);
    quad.push_back(v3);
    quad.push_back(v4);
    Vector3d p1 = (O_compact[v1] + O_compact[v2] + O_compact[v3] + O_compact[v4]) / 4;
    O_compact[v1] = p1;
    O_compact[v2] = p1;
    O_compact[v3] = p1;
    O_compact[v4] = p1;
    std::vector<int> twice;
    for (int i = 0; i < position.size(); ++i) {
        int n = position[i];
        for (int j = 0; j < layer[n].size(); ++j) {
            if (find(delete_patch.begin(), delete_patch.end(), layer[n][j]) ==
                delete_patch.end()) {
                delete_patch.push_back(layer[n][j]);
            } else {
                twice.push_back(layer[n][j]);
            }
        }
    }
    delete_patch.push_back(v);
    std::vector<int> need_collapse(position.size(), 0);
    std::vector<std::vector<int>> corr(position.size());
    for (int i = 0; i < position.size(); ++i) {
        int n = position[i];
        if (layer[n].size() == 1) continue;
        for (int j = 0; j < layer[n].size(); ++j) {
            for (int k = 0; k < 4; ++k) {
                int v1 = relative_point[n][j][k];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                        all_singularity.end() &&
                    find(quad.begin(), quad.end(), v1) == quad.end()) {
                    need_collapse[i] = 1;
                }
            }
        }
    }
    std::vector<std::vector<int>> sin_p(position.size());
    std::vector<std::vector<int>> collapse_size(position.size());
    for (int i = 0; i < position.size(); ++i) {
        if (need_collapse[i] == 0) continue;
        int n = position[i];
        int ver1 = collapse_patch[v][i][0];
        int ver2 = collapse_patch[v][i][1];
        if (relative_point[n][0][1] == ver1) {
            for (int j = 0; j < relative_point[n].size(); ++j) {
                int v1 = relative_point[n][j][2];
                int v2 = relative_point[n][j][3];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    sin_p[i].push_back(1);
                } else if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                           all_singularity.end()) {
                    sin_p[i].push_back(0);
                } else {
                    sin_p[i].push_back(2);
                }
            }
            if (sin_p[i][0] == 2) {
                int a1, a2;
                for (int j = 0; j < sin_p[i][0]; ++j) {
                    if (sin_p[i][j] == 2) continue;
                    a1 = j;
                    a2 = sin_p[i][j];
                    break;
                }
                for (int j = 0; j < a1; ++j) {
                    sin_p[i][j] = a2;
                    break;
                }
            }
            bool spicy = false;
            int a1 = -1;
            for (int j = 0; j < sin_p[i].size(); ++j) {
                if (spicy) {
                    if (sin_p[i][j] != 2) {
                        if (sin_p[i][j] == sin_p[i][a1]) {
                            sin_p[i][j] = 2;
                            spicy = false;
                        } else {
                            for (int k = a1 + 1; k < j; ++k) {
                                sin_p[i][k] = sin_p[i][j];
                            }
                            spicy = false;
                        }
                    }
                } else {
                    if (sin_p[i][j] == 2) {
                        spicy = true;
                        if (a1 == -1) {
                            a1 = j - 1;
                        }
                    }
                }
            }
            int begin = 0, end = sin_p[i].size() - 1;
            int size = 0, sign = sin_p[i][0];
            for (int j = 0; j < sin_p[i].size(); ++j) {
                int n1 = layer[n][j];
                if (sin_p[i][j] == sign) {
                    int aaa;
                    if (collapse_patch[n1][0].size() > collapse_patch[n1][1].size()) {
                        aaa = collapse_patch[n1][0].size();
                    } else {
                        aaa = collapse_patch[n1][1].size();
                    }
                    size = size + aaa - 1;
                    if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                        end = j;
                        for (int k = begin; k < end + 1; ++k) {
                            if (sin_p[i][k] == 2) {
                                collapse_size[i].push_back(2);
                            } else {
                                collapse_size[i].push_back(size);
                            }
                        }
                        size = 0;
                        begin = j;
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != 2) {
                            sign = sin_p[i][j + 1];
                            begin = j + 1;
                        }
                    }

                } else {
                    if (size == 0) {
                        begin = j;
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != 2) {
                            sign = sin_p[i][j + 1];
                        } else if (j != sin_p[i].size() - 1) {
                            collapse_size[i].push_back(2);
                        }
                    }
                }
                if (j == sin_p[i].size() - 1) {
                    if (sin_p[i][j] == sign) {
                        end = j;
                        for (int k = begin; k < end + 1; ++k) {
                            collapse_size[i].push_back(size);
                        }
                    } else if (sin_p[i][j] == 2) {
                        collapse_size[i].push_back(2);
                    }
                }
            }
        } else {
            for (int j = relative_point[n].size() - 1; j > -1; --j) {
                int v1 = relative_point[n][j][0];
                int v2 = relative_point[n][j][1];
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    sin_p[i].insert(sin_p[i].begin(), 0);
                } else if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                           all_singularity.end()) {
                    sin_p[i].insert(sin_p[i].begin(), 1);
                } else {
                    sin_p[i].insert(sin_p[i].begin(), 2);
                }
            }
            if (sin_p[i][sin_p[i].size() - 1] == 2) {
                int a1, a2;
                for (int j = sin_p[i].size() - 1; j > -1; --j) {
                    if (sin_p[i][j] == 2) continue;
                    a1 = j;
                    a2 = sin_p[i][j];
                    break;
                }
                for (int j = a1; j < sin_p[i].size(); ++j) {
                    sin_p[i][j] = a2;
                    break;
                }
            }
            bool spicy = false;
            int a1 = -1;
            for (int j = sin_p[i].size() - 1; j > -1; --j) {
                if (spicy) {
                    if (sin_p[i][j] != 2) {
                        if (sin_p[i][j] == sin_p[i][a1]) {
                            sin_p[i][j] = 2;
                            spicy = false;
                        } else {
                            for (int k = a1 - 1; k > j; --k) {
                                sin_p[i][k] = sin_p[i][j];
                            }
                            spicy = false;
                        }
                    }
                } else {
                    if (sin_p[i][j] == 2) {
                        spicy = true;
                        if (a1 == -1) {
                            a1 = j + 1;
                        }
                    }
                }
            }
            int begin = sin_p[i].size() - 1, end = 0;
            int size = 0, sign = sin_p[i][sin_p[i].size() - 1];
            for (int j = sin_p[i].size() - 1; j > -1; --j) {
                int n1 = layer[n][j];
                if (sin_p[i][j] == sign) {
                    int aaa;
                    if (collapse_patch[n1][0].size() > collapse_patch[n1][1].size()) {
                        aaa = collapse_patch[n1][0].size();
                    } else {
                        aaa = collapse_patch[n1][1].size();
                    }
                    size = size + aaa - 1;
                    if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                        end = j;
                        for (int k = begin; k > end - 1; --k) {
                            if (sin_p[i][k] == 2) {
                                collapse_size[i].insert(collapse_size[i].begin(), 2);
                            } else {
                                collapse_size[i].insert(collapse_size[i].begin(), size);
                            }
                        }
                        size = 0;
                        begin = j;
                        if (j != 0 && sin_p[i][j - 1] != 2) {
                            sign = sin_p[i][j - 1];
                            begin = j - 1;
                        }
                    }
                } else {
                    if (size == 0) {
                        begin = j;
                        if (j != 0 && sin_p[i][j - 1] != 2) {
                            sign = sin_p[i][j - 1];
                        } else if (j != 0) {
                            collapse_size[i].insert(collapse_size[i].begin(), 2);
                        }
                    }
                }
                if (j == 0) {
                    if (sin_p[i][j] == sign) {
                        end = j;
                        for (int k = begin; k > end - 1; --k) {
                            collapse_size[i].insert(collapse_size[i].begin(), size);
                        }
                    } else if (sin_p[i][j] == 2) {
                        collapse_size[i].insert(collapse_size[i].begin(), 2);
                    }
                }
            }
        }
    }

    for (int i = 0; i < position.size(); ++i) {
        if (need_collapse[i] == 0) {
            int n = position[i];
            int ver1 = collapse_patch[v][i][0];
            int ver2 = collapse_patch[v][i][1];
            if (valence[ver1] == 3 || valence[ver2] == 3) {
                int sin, no_sin;
                if (valence[ver1] == 3) {
                    sin = ver1;
                    no_sin = ver2;
                }
                if (valence[ver2] == 3) {
                    sin = ver2;
                    no_sin = ver1;
                }
                int a = -1;
                if (relative_point[n].size() != 1) {
                    if (relative_point[n][0][0] == sin ||
                        relative_point[n][relative_point[n].size() - 1][3] == sin) {
                        a = 0;
                    } else {
                        a = 1;
                    }
                }
                if (a == -1) {
                    std::vector<int> line1, line2;
                    int p = layer[n][0];
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == sin && collapse_patch[p][k][1] == no_sin) {
                            line1 = collapse_patch[p][(k + 3) % 4];
                            line2 = collapse_patch[p][(k + 1) % 4];
                            break;
                        }
                        if (collapse_patch[p][k][0] == no_sin && collapse_patch[p][k][1] == sin) {
                            line1 = collapse_patch[p][(k + 1) % 4];
                            line2 = collapse_patch[p][(k + 3) % 4];
                            break;
                        }
                    }
                    for (int k = 0; k < line1.size(); ++k) {
                        int o1 = line1[k];
                        int o2 = line2[line2.size() - 1 - k];
                        O_compact[o2] = O_compact[o1];
                        if (renew[o2] == -1) {
                            renew[o2] = o1;
                        }
                    }
                } else {
                    for (int j = 0; j < layer[n].size(); ++j) {
                        int p = layer[n][j];
                        if (a == 0) {
                            v1 = relative_point[n][j][a];
                            v2 = relative_point[n][j][(a + 1) % 2];
                        } else {
                            v2 = relative_point[n][j][a];
                            v1 = relative_point[n][j][(a + 1) % 2];
                        }

                        std::vector<int> line1, line2;
                        for (int k = 0; k < collapse_patch[p].size(); ++k) {
                            if (collapse_patch[p][k][0] == v2) {
                                if (a == 0) {
                                    line2 = collapse_patch[p][k];
                                    line1 = collapse_patch[p][(k + 2) % 4];
                                } else {
                                    line1 = collapse_patch[p][k];
                                    line2 = collapse_patch[p][(k + 2) % 4];
                                }
                                break;
                            }
                        }
                        for (int k = 0; k < line1.size(); ++k) {
                            int o1 = line1[k];
                            int o2 = line2[line2.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                    }
                }
            }
        } else {
            int n = position[i];
            int ver1 = collapse_patch[v][i][0];
            int ver2 = collapse_patch[v][i][1];
            if (relative_point[n][0][1] == ver1) {
                bool begin = true;
                int b1 = 0;
                int b2 = 0;
                int c1 = sin_p[i][0];
                for (int j = 0; j < layer[n].size(); ++j) {
                    int p = layer[n][j];
                    v1 = relative_point[n][j][0];
                    v2 = relative_point[n][j][1];
                    std::vector<int> line1, line2;
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == v2) {
                            line2 = collapse_patch[p][k];
                            line1 = collapse_patch[p][(k + 2) % 4];
                            break;
                        }
                    }
                    if (valence[v2] == 3 && begin && sin_p[i][j]) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (valence[v1] == 3 && begin && sin_p[i][j] == 0) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o1] = O_compact[o2];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (begin && sin_p[i][j]) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 - b1 * d / 2 / size;
                            O_compact[o2] = p2 - b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                            if (sin_p[i][j + 1] != 2) {
                                c1 = sin_p[i][j + 1];
                            }
                            begin = false;
                        }
                    } else if (begin && sin_p[i][j] == 0) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;

                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 + b1 * d / 2 / size;
                            O_compact[o2] = p2 + b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != sin_p[i].size() - 1 && sin_p[i][j + 1] != sin_p[i][j]) {
                            if (sin_p[i][j + 1] != 2) {
                                c1 = sin_p[i][j + 1];
                            }
                            begin = false;
                        }
                    } else {
                        if (sin_p[i][j] != 2) {
                            c1 = sin_p[i][j];
                        }
                        if (c1 == 0 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o1] = O_compact[o2];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else if (c1 == 1 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o2] = O_compact[o1];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else {
                            if (c1 == 0) {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o2] - O_compact[o1];
                                    O_compact[o1] = O_compact[o1] + d * b1 / size;
                                    O_compact[o2] = O_compact[o1];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            } else {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o1] - O_compact[o2];
                                    O_compact[o2] = O_compact[o2] + d * b1 / size;
                                    O_compact[o1] = O_compact[o2];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            }
                        }
                    }
                }
            } else {
                /////////////////////暂时未修改///////////////////////////////////////////////////////////5944-6065
                bool begin = true;
                int b1 = 0;
                int b2 = 0;
                int c1 = sin_p[i][sin_p[i].size() - 1];
                for (int j = layer[n].size() - 1; j > -1; --j) {
                    int p = layer[n][j];
                    v1 = relative_point[n][j][2];
                    v2 = relative_point[n][j][3];
                    std::vector<int> line1, line2;
                    for (int k = 0; k < collapse_patch[p].size(); ++k) {
                        if (collapse_patch[p][k][0] == v2) {
                            line2 = collapse_patch[p][k];
                            line1 = collapse_patch[p][(k + 2) % 4];
                            break;
                        }
                    }
                    if (valence[v2] == 3 && begin && sin_p[i][j]) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o2] = O_compact[o1];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (valence[v1] == 3 && begin && sin_p[i][j] == 0) {
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            O_compact[o1] = O_compact[o2];
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                        }
                        begin = false;
                        continue;
                    }
                    if (begin && sin_p[i][j]) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 + b1 * d / 2 / size;
                            O_compact[o2] = p2 + b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                            c1 = sin_p[i][j - 1];
                            begin = false;
                        }
                    } else if (begin && sin_p[i][j] == 0) {
                        int size = collapse_size[i][j];
                        for (int k = 0; k < line2.size(); ++k) {
                            if (k == 0 && b2 == 0) {
                                b1 += 1;
                                continue;
                            }
                            if (k == 0) continue;
                            int o1 = line2[k];
                            int o2 = line1[line1.size() - 1 - k];
                            Vector3d d = O_compact[o2] - O_compact[o1];
                            Vector3d p2 = (O_compact[o1] + O_compact[o2]) / 2;
                            O_compact[o1] = p2 - b1 * d / 2 / size;
                            O_compact[o2] = p2 - b1 * d / 2 / size;
                            if (renew[o2] == -1) {
                                renew[o2] = o1;
                            }
                            b1 += 1;
                        }
                        if (b2 == 0) {
                            b2 += 1;
                        }
                        if (b1 == size + 1) {
                            b1 = 0;
                            b2 = 0;
                        }
                        if (j != 0 && sin_p[i][j - 1] != sin_p[i][j]) {
                            c1 = sin_p[i][j - 1];
                            begin = false;
                        }
                    } else {
                        if (sin_p[i][j] != 2) {
                            c1 = sin_p[i][j];
                        }
                        if (c1 == 0 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o2] = O_compact[o1];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else if (c1 == 1 && sin_p[i][j] == 2) {
                            for (int k = 0; k < line2.size(); ++k) {
                                int o1 = line2[k];
                                int o2 = line1[line1.size() - 1 - k];
                                O_compact[o1] = O_compact[o2];
                                if (renew[o2] == -1) {
                                    renew[o2] = o1;
                                }
                            }
                        } else {
                            if (c1 == 0) {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o1] - O_compact[o2];
                                    O_compact[o2] = O_compact[o2] + d * b1 / size;
                                    O_compact[o1] = O_compact[o2];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            } else {
                                int size = collapse_size[i][j];
                                for (int k = 0; k < line2.size(); ++k) {
                                    if (k == 0 && b2 == 0) {
                                        b1 += 1;
                                        continue;
                                    }
                                    if (k == 0) continue;
                                    int o1 = line2[k];
                                    int o2 = line1[line1.size() - 1 - k];
                                    Vector3d d = O_compact[o2] - O_compact[o1];
                                    O_compact[o1] = O_compact[o1] + d * b1 / size;
                                    O_compact[o2] = O_compact[o1];
                                    if (renew[o2] == -1) {
                                        renew[o2] = o1;
                                    }
                                    b1 += 1;
                                }
                                if (b2 == 0) {
                                    b2 += 1;
                                }
                                if (b1 == size + 1) {
                                    b1 = 0;
                                    b2 = 0;
                                }
                            }
                        }
                    }
                }
                /////////////////////暂时未修改///////////////////////////////////////////////////////////5944-6065
            }
        }
    }

    std::sort(delete_patch.begin(), delete_patch.end());
    // update information
    //  update valence
    valence[x1] = 3;
    valence[x2] = 0;
    valence[x3] = 0;
    valence[x4] = 0;
    // update singularity
    std::vector<int> temp_singularity;
    for (int i = 0; i < all_singularity.size(); ++i) {
        if (find(singularity.begin(), singularity.end(), all_singularity[i]) != singularity.end())
            continue;
        temp_singularity.push_back(all_singularity[i]);
    }
    all_singularity = temp_singularity;
    int singularity_point;
    for (int i = 0; i < twice.size(); ++i) {
        int u = twice[i];
        for (int j = 0; j < collapse_patch[u].size(); ++j) {
            if (find(all_singularity.begin(), all_singularity.end(), collapse_patch[u][j][0]) !=
                all_singularity.end()) {
                singularity_point = collapse_patch[u][j][0];
            }
        }
        for (int j = 0; j < collapse_patch[u].size(); ++j) {
            if (collapse_patch[u][j][0] == singularity_point) continue;
            renew[collapse_patch[u][j][0]] = singularity_point;
        }
    }
    for (int i = 0; i < renew.size(); ++i) {
        if (renew[i] == -1) continue;
        if (find(all_singularity.begin(), all_singularity.end(), i) != all_singularity.end()) {
            valence[renew[i]] = valence[i];
            for (int j = 0; j < all_singularity.size(); ++j) {
                if (all_singularity[j] == i) {
                    all_singularity[j] = renew[i];
                }
            }
        }
    }
    // update collapse_patch
    std::vector<std::vector<std::vector<int>>> temp_collapse_patch;
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < collapse_patch[i].size(); ++j) {
            for (int k = 0; k < collapse_patch[i][j].size(); ++k) {
                int v = collapse_patch[i][j][k];
                if (renew[v] != -1) {
                    collapse_patch[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_collapse_patch.push_back(collapse_patch[i]);
    }
    collapse_patch = temp_collapse_patch;
    // update patch_compact
    std::vector<std::vector<std::vector<int>>> temp_patch_compact;
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int v = patch_compact[i][j][k];
                if (renew[v] != -1) {
                    patch_compact[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_patch_compact.push_back(patch_compact[i]);
    }
    patch_compact = temp_patch_compact;
    // update relative_point
    std::vector<std::vector<std::vector<int>>> temp_relative_point(relative_point.size());
    std::vector<std::vector<int>> temp_layer(layer.size());
    for (int i = 0; i < relative_point.size(); ++i) {
        for (int j = 0; j < relative_point[i].size(); ++j) {
            for (int k = 0; k < relative_point[i][j].size(); ++k) {
                int v = relative_point[i][j][k];
                if (renew[v] != -1) {
                    relative_point[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < layer.size(); ++i) {
        if (find(position.begin(), position.end(), i) != position.end()) continue;
        for (int j = 0; j < layer[i].size(); ++j) {
            int v = layer[i][j];
            if (find(delete_patch.begin(), delete_patch.end(), v) != delete_patch.end()) continue;
            temp_relative_point[i].push_back(relative_point[i][j]);
            for (int k = 0; k < delete_patch.size(); ++k) {
                if (v < delete_patch[0]) {
                    temp_layer[i].push_back(v);
                    break;
                }
                if (k == delete_patch.size() - 1) {
                    if (v > delete_patch[k]) {
                        temp_layer[i].push_back(v - k - 1);
                    }
                } else {
                    if (v > delete_patch[k] && v < delete_patch[k + 1]) {
                        temp_layer[i].push_back(v - k - 1);
                        break;
                    }
                }
            }
        }
    }
    for (int i = 0; i < thr_for_sin_patch.size(); ++i) {
        int n = thr_for_sin_patch[i];
        for (int k = 0; k < delete_patch.size(); ++k) {
            if (k == delete_patch.size() - 1) {
                if (n > delete_patch[k]) {
                    thr_for_sin_patch[i] = n - k - 1;
                }
            } else {
                if (n > delete_patch[k] && n < delete_patch[k + 1]) {
                    thr_for_sin_patch[i] = n - k - 1;
                    break;
                }
            }
        }
    }
    std::vector<std::vector<std::vector<int>>> temp;
    std::vector<std::vector<int>> temp1;
    for (int i = 0; i < temp_relative_point.size(); ++i) {
        if (temp_relative_point[i].size() == 0) continue;
        temp.push_back(temp_relative_point[i]);
    }
    relative_point = temp;
    for (int i = 0; i < temp_layer.size(); ++i) {
        if (temp_layer[i].size() == 0) continue;
        temp1.push_back(temp_layer[i]);
    }
    layer = temp1;
}

void Optimizer::sin_ft(std::vector<Vector3d>& O_compact,
                       std::vector<std::vector<std::vector<int>>>& patch_compact,
                       std::vector<int>& valence,
                       std::vector<std::vector<std::vector<int>>>& relative_point,
                       std::vector<std::vector<int>>& layer,
                       std::vector<std::vector<std::vector<int>>>& collapse_patch, int v, int v1,
                       int v2, std::vector<int>& all_singularity,
                       std::vector<int>& thr_for_sin_patch) {
    std::vector<int> renew(O_compact.size(), -1);
    std::vector<int> delete_patch;
    int coll = 0;
    for (int i = 0; i < relative_point[v].size(); ++i) {
        int p1, p2;
        p1 = relative_point[v][i][1];
        p2 = relative_point[v][i][2];
        if (p1 != v1 && p1 != v2 &&
            find(all_singularity.begin(), all_singularity.end(), p1) != all_singularity.end()) {
            coll = 1;
        }
        if (p2 != v1 && p2 != v2 &&
            find(all_singularity.begin(), all_singularity.end(), p2) != all_singularity.end()) {
            coll = 1;
        }
    }

    for (int i = 0; i < layer[v].size(); ++i) {
        int n = layer[v][i];
        delete_patch.push_back(n);
        std::vector<int> line1, line2;
        int p1, p2;
        p1 = relative_point[v][i][0];
        p2 = relative_point[v][i][3];
        if (coll == 0) {
            for (int j = 0; j < collapse_patch[n].size(); ++j) {
                if (collapse_patch[n][j][0] == p2 &&
                    collapse_patch[n][j][collapse_patch[n][j].size() - 1] == p1) {
                    line1 = collapse_patch[n][j];
                    line2 = collapse_patch[n][(j + 2) % 4];
                }
            }
        } else {
            for (int j = 0; j < collapse_patch[n].size(); ++j) {
                if (collapse_patch[n][j][0] == p2 &&
                    collapse_patch[n][j][collapse_patch[n][j].size() - 1] == p1) {
                    line2 = collapse_patch[n][j];
                    line1 = collapse_patch[n][(j + 2) % 4];
                }
            }
        }
        for (int j = 0; j < line1.size(); ++j) {
            int o1 = line1[j];
            int o2 = line2[line2.size() - 1 - j];
            O_compact[o2] = O_compact[o1];
            if (renew[o2] == -1) {
                renew[o2] = o1;
            }
        }
    }
    std::sort(delete_patch.begin(), delete_patch.end());
    // update information
    //  update valence
    valence[v1] = 4;
    valence[v2] = 4;
    // update singularity

    std::vector<int> temp_singularity;
    for (int i = 0; i < all_singularity.size(); ++i) {
        if (all_singularity[i] == v1) continue;
        if (all_singularity[i] == v2) continue;
        temp_singularity.push_back(all_singularity[i]);
    }
    all_singularity = temp_singularity;
    // update collapse_patch
    std::vector<std::vector<std::vector<int>>> temp_collapse_patch;
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < collapse_patch[i].size(); ++j) {
            for (int k = 0; k < collapse_patch[i][j].size(); ++k) {
                int n = collapse_patch[i][j][k];
                if (renew[n] != -1) {
                    collapse_patch[i][j][k] = renew[n];
                }
            }
        }
    }
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_collapse_patch.push_back(collapse_patch[i]);
    }
    collapse_patch = temp_collapse_patch;
    // update patch_compact
    std::vector<std::vector<std::vector<int>>> temp_patch_compact;
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int n = patch_compact[i][j][k];
                if (renew[n] != -1) {
                    patch_compact[i][j][k] = renew[n];
                }
            }
        }
    }
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_patch_compact.push_back(patch_compact[i]);
    }
    patch_compact = temp_patch_compact;
    // update relative_point
    std::vector<std::vector<std::vector<int>>> temp_relative_point(relative_point.size());
    std::vector<std::vector<int>> temp_layer(layer.size());
    for (int i = 0; i < relative_point.size(); ++i) {
        for (int j = 0; j < relative_point[i].size(); ++j) {
            for (int k = 0; k < relative_point[i][j].size(); ++k) {
                int n = relative_point[i][j][k];
                if (renew[n] != -1) {
                    relative_point[i][j][k] = renew[n];
                }
            }
        }
    }
    for (int i = 0; i < layer.size(); ++i) {
        if (i == v) continue;
        for (int j = 0; j < layer[i].size(); ++j) {
            int n = layer[i][j];
            if (find(delete_patch.begin(), delete_patch.end(), n) != delete_patch.end()) continue;
            temp_relative_point[i].push_back(relative_point[i][j]);
            for (int k = 0; k < delete_patch.size(); ++k) {
                if (n < delete_patch[0]) {
                    temp_layer[i].push_back(n);
                    break;
                }
                if (k == delete_patch.size() - 1) {
                    if (n > delete_patch[k]) {
                        temp_layer[i].push_back(n - k - 1);
                    }
                } else {
                    if (n > delete_patch[k] && n < delete_patch[k + 1]) {
                        temp_layer[i].push_back(n - k - 1);
                        break;
                    }
                }
            }
        }
    }
    for (int i = 0; i < thr_for_sin_patch.size(); ++i) {
        int n = thr_for_sin_patch[i];
        for (int k = 0; k < delete_patch.size(); ++k) {
            if (k == delete_patch.size() - 1) {
                if (n > delete_patch[k]) {
                    thr_for_sin_patch[i] = n - k - 1;
                }
            } else {
                if (n > delete_patch[k] && n < delete_patch[k + 1]) {
                    thr_for_sin_patch[i] = n - k - 1;
                    break;
                }
            }
        }
    }
    std::vector<std::vector<std::vector<int>>> temp;
    std::vector<std::vector<int>> temp1;
    for (int i = 0; i < temp_relative_point.size(); ++i) {
        if (temp_relative_point[i].size() == 0) continue;
        temp.push_back(temp_relative_point[i]);
    }
    relative_point = temp;
    for (int i = 0; i < temp_layer.size(); ++i) {
        if (temp_layer[i].size() == 0) continue;
        temp1.push_back(temp_layer[i]);
    }
    layer = temp1;
}

void Optimizer::match_patch(int i, int j, int v, std::vector<int>& edge_collapse,
                            std::vector<std::vector<std::vector<int>>>& collapse_patch,
                            std::vector<std::vector<std::vector<int>>>& patch_compact,
                            std::vector<std::vector<std::vector<int>>>& assist_patch,
                            std::vector<int>& can, std::vector<int>& relative_can,
                            std::vector<int>& relative_double, std::vector<int>& double_following,
                            std::vector<int>& thenum, std ::vector<int>& renew) {
    int v1, v2;
    int beg1, beg2;
    v1 = edge_collapse[j * 2];
    v2 = edge_collapse[j * 2 + 1];
    std::vector<int> no_need_edge, need_collapse_edge;
    if (collapse_patch[v][0][0] == v1 &&
            collapse_patch[v][0][collapse_patch[v][0].size() - 1] == v2 ||
        collapse_patch[v][1][0] == v1 &&
            collapse_patch[v][1][collapse_patch[v][1].size() - 1] == v2 ||
        collapse_patch[v][2][0] == v1 &&
            collapse_patch[v][2][collapse_patch[v][2].size() - 1] == v2 ||
        collapse_patch[v][3][0] == v1 &&
            collapse_patch[v][3][collapse_patch[v][3].size() - 1] == v2) {
    } else {
        int temp = v1;
        v1 = v2;
        v2 = temp;
    }
    if (collapse_patch[v][0].size() > collapse_patch[v][1].size()) {
        if (collapse_patch[v][0][0] == v1 ||
            collapse_patch[v][0][collapse_patch[v][0].size() - 1] == v1) {
            no_need_edge = collapse_patch[v][0];
            need_collapse_edge = collapse_patch[v][2];
        }
        if (collapse_patch[v][2][0] == v1 ||
            collapse_patch[v][2][collapse_patch[v][2].size() - 1] == v1) {
            no_need_edge = collapse_patch[v][2];
            need_collapse_edge = collapse_patch[v][0];
        }
    } else if (collapse_patch[v][0].size() < collapse_patch[v][1].size()) {
        if (collapse_patch[v][1][0] == v1 ||
            collapse_patch[v][1][collapse_patch[v][1].size() - 1] == v1) {
            no_need_edge = collapse_patch[v][1];
            need_collapse_edge = collapse_patch[v][3];
        }
        if (collapse_patch[v][3][0] == v1 ||
            collapse_patch[v][3][collapse_patch[v][3].size() - 1] == v1) {
            no_need_edge = collapse_patch[v][3];
            need_collapse_edge = collapse_patch[v][1];
        }
    } else {
        for (int k = 0; k < collapse_patch[v].size(); ++k) {
            if (collapse_patch[v][k][0] == v1 &&
                collapse_patch[v][k][collapse_patch[v][k].size() - 1] == v2) {
                no_need_edge = collapse_patch[v][k];
                need_collapse_edge = collapse_patch[v][(k + 2) % 4];
                break;
            }
            if (collapse_patch[v][k][0] == v2 &&
                collapse_patch[v][k][collapse_patch[v][k].size() - 1] == v1) {
                no_need_edge = collapse_patch[v][k];
                need_collapse_edge = collapse_patch[v][(k + 2) % 4];
                break;
            }
        }
    }
    for (int i = 0; i < no_need_edge.size(); ++i) {
        int o1 = no_need_edge[i];
        int o2 = need_collapse_edge[no_need_edge.size() - 1 - i];
        if (renew[o2] == -1) {
            renew[o2] = o1;
        }
    }
}

void Optimizer::update_O(std::vector<Vector3d>& O_compact,
                         std::vector<std::vector<std::vector<int>>>& patch_compact,
                         std::vector<int>& all_singularity, std::vector<int>& boundary_o) {
    std::vector<int> temp_O(O_compact.size(), -1);
    for (int i = 0; i < patch_compact.size(); ++i) {
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int v = patch_compact[i][j][k];
                if (temp_O[v] == -1) {
                    temp_O[v] = v;
                }
            }
        }
    }
    std::vector<Vector3d> temp_compact;
    std::vector<int> delete_O;
    std::vector<int> boundary;
    for (int i = 0; i < temp_O.size(); ++i) {
        if (temp_O[i] != -1) {
            temp_compact.push_back(O_compact[i]);
            boundary.push_back(boundary_o[i]);
            continue;
        }
        delete_O.push_back(i);
    }
    std::cout << O_compact.size() << " " << temp_compact.size() << "\n";
    O_compact = temp_compact;
    boundary_o = boundary;
    for (int i = 0; i < patch_compact.size(); ++i) {
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int v = patch_compact[i][j][k];
                for (int p = 0; p < delete_O.size(); ++p) {
                    if (v < delete_O[0]) break;
                    if (p == delete_O.size() - 1) {
                        if (v > delete_O[p]) {
                            patch_compact[i][j][k] = patch_compact[i][j][k] - p - 1;
                        }
                    } else {
                        if (v > delete_O[p] && v < delete_O[p + 1]) {
                            patch_compact[i][j][k] = patch_compact[i][j][k] - p - 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < all_singularity.size(); ++i) {
        int v = all_singularity[i];
        for (int p = 0; p < delete_O.size(); ++p) {
            if (v < delete_O[0]) break;
            if (p == delete_O.size() - 1) {
                if (v > delete_O[p]) {
                    all_singularity[i] = all_singularity[i] - p - 1;
                }
            } else {
                if (v > delete_O[p] && v < delete_O[p + 1]) {
                    all_singularity[i] = all_singularity[i] - p - 1;
                    break;
                }
            }
        }
    }
}
void Optimizer::sin_mismatching(std::vector<std::vector<int>>& need_collapse_part,
                                std::vector<std::vector<std::vector<int>>>& patch_compact,
                                std::vector<Vector3d>& O_compact,
                                std::vector<std::vector<std::vector<int>>>& relative_point,
                                std::vector<std::vector<int>>& layer,
                                std::vector<int>& need_collapse, std::vector<int>& all_singularity,
                                std::vector<std::vector<std::vector<int>>>& collapse_patch, int n,
                                std::vector<std::vector<int>>& double_q,
                                std::vector<int>& boundary_o) {
    std::vector<Vector3d> patch_position = O_compact;
    std::vector<std::vector<std::vector<int>>> assist_patch(patch_compact.size());
    std::vector<std::vector<int>> following_collapse(need_collapse.size());
    std::vector<int> edge_collapse;
    int v = need_collapse[n];
    std::vector<int> pre_lay;
    std::vector<int> renew(O_compact.size(), -1);
    for (int i = 0; i < need_collapse_part[v].size(); ++i) {
        if (i % 2 == 1) continue;
        int begin = need_collapse_part[v][i];
        int end = need_collapse_part[v][i + 1];
        if (i == 0) {
            for (int j = 0; j < begin; ++j) {
                following_collapse[n].push_back(layer[v][j]);
            }
        }
        if (i == need_collapse_part[v].size() - 2) {
            for (int j = end + 1; j < layer[v].size(); ++j) {
                following_collapse[n].push_back(layer[v][j]);
            }
        }
        if (i != 0) {
            for (int j = need_collapse_part[v][i - 1] + 1; j < begin; ++j) {
                following_collapse[n].push_back(layer[v][j]);
            }
        }
    }

    for (int j = 0; j < need_collapse_part[v].size(); ++j) {
        if (j % 2 == 1) continue;  ////  don't consider the size larger 2
        if (j == 0) {
            int begin = need_collapse_part[v][j];
            int end = need_collapse_part[v][j + 1];
            int start;
            int v1;
            v1 = relative_point[v][begin][0];
            if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                all_singularity.end()) {
                start = 0;
            } else {
                start = 1;
            }
            for (int k = 0; k < begin; ++k) {
                if (start == 0) {
                    edge_collapse.push_back(relative_point[v][k][3]);
                    edge_collapse.push_back(relative_point[v][k][0]);
                }
                if (start == 1) {
                    edge_collapse.push_back(relative_point[v][k][1]);
                    edge_collapse.push_back(relative_point[v][k][2]);
                }
            }
        }
        if (j == need_collapse_part[v].size() - 2) {
            int begin = need_collapse_part[v][j];
            int end = need_collapse_part[v][j + 1];
            int start;
            int v1;
            v1 = relative_point[v][end][3];
            if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                all_singularity.end()) {
                start = 0;
            } else {
                start = 1;
            }
            for (int k = end + 1; k < relative_point[v].size(); ++k) {
                if (start == 0) {
                    edge_collapse.push_back(relative_point[v][k][3]);
                    edge_collapse.push_back(relative_point[v][k][0]);
                }
                if (start == 1) {
                    edge_collapse.push_back(relative_point[v][k][1]);
                    edge_collapse.push_back(relative_point[v][k][2]);
                }
            }
        }
        if (j != 0) {
            int begin = need_collapse_part[v][j - 1] + 1;
            int end = need_collapse_part[v][j];
            int start;
            int v1;
            v1 = relative_point[v][begin][0];
            if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                all_singularity.end()) {
                start = 0;
            } else {
                start = 1;
            }
            for (int k = begin; k < end; ++k) {
                if (start == 0) {
                    edge_collapse.push_back(relative_point[v][k][3]);
                    edge_collapse.push_back(relative_point[v][k][0]);
                }
                if (start == 1) {
                    edge_collapse.push_back(relative_point[v][k][1]);
                    edge_collapse.push_back(relative_point[v][k][2]);
                }
            }
        }
    }
    int number = need_collapse[n];
    for (int i = 0; i < need_collapse_part[number].size(); ++i) {
        if (i % 2 == 1) continue;
        int begin = need_collapse_part[number][i];
        int end = need_collapse_part[number][i + 1];
        int total_size = 0;
        int b = 0;
        for (int k = begin; k < end + 1; ++k) {
            int pa = layer[number][k];
            if (collapse_patch[pa][0].size() > collapse_patch[pa][1].size()) {
                total_size += collapse_patch[pa][0].size() - 1;
            } else {
                total_size += collapse_patch[pa][1].size() - 1;
            }
        }
        int inner_vertex_ = 0;
        for (int s = begin; s < end + 1; ++s) {
            int pa = layer[number][s];
            pre_lay.push_back(pa);
        }
        for (int s = begin; s < end + 1; ++s) {
            int pa = layer[number][s];
            std::vector<int> edge1, edge2;
            int singularity_vertex;
            int v1, v2, v3, v4;
            v1 = relative_point[number][s][0];
            v2 = relative_point[number][s][1];
            v3 = relative_point[number][s][2];
            v4 = relative_point[number][s][3];
            if (s == begin) {
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    singularity_vertex = v1;
                } else {
                    singularity_vertex = v2;
                }
            } else {
                singularity_vertex = inner_vertex_;
            }
            for (int p = 0; p < collapse_patch[pa].size(); ++p) {
                if (collapse_patch[pa][p][0] == v2 &&
                    collapse_patch[pa][p][collapse_patch[pa][p].size() - 1] == v3) {
                    edge1 = collapse_patch[pa][p];
                    edge2 = collapse_patch[pa][(p + 2) % 4];
                    break;
                }
            }
            int inner_vertex;
            if (singularity_vertex == edge1[0]) {
                inner_vertex = edge1[edge1.size() - 1];
                inner_vertex_ = edge1[edge1.size() - 1];
                for (int j = 0; j < edge1.size(); ++j) {
                    if (s != begin && j == 0) continue;
                    Vector3d direction =
                        O_compact[edge2[edge2.size() - 1 - j]] - O_compact[edge1[j]];
                    Vector3d position =
                        O_compact[edge1[j]] + (double)b * direction / (double)total_size;
                    b += 1;
                    int o1 = edge1[j];
                    int o2 = edge2[edge2.size() - 1 - j];
                    if (renew[o2] == -1) {
                        renew[o2] = o1;
                    }
                    O_compact[o2] = position;
                    O_compact[o1] = position;
                }
            } else {
                singularity_vertex = edge1[0];
                inner_vertex = edge1[edge1.size() - 1];
                inner_vertex_ = edge2[0];
                for (int j = 0; j < edge2.size(); ++j) {
                    if (s != begin && j == 0) continue;
                    Vector3d direction =
                        O_compact[edge1[j]] - O_compact[edge2[edge2.size() - 1 - j]];
                    Vector3d position = O_compact[edge2[edge2.size() - 1 - j]] +
                                        (double)b * direction / (double)total_size;
                    b += 1;
                    int o1 = edge2[edge2.size() - 1 - j];
                    int o2 = edge1[j];
                    if (renew[o2] == -1) {
                        renew[o2] = o1;
                    }
                    O_compact[o2] = position;
                    O_compact[o1] = position;
                }
            }
        }
    }

    std::vector<int> double_following;
    std::vector<int> correspond;
    for (int i = 0; i < following_collapse.size(); ++i) {
        for (int j = 0; j < following_collapse[i].size(); ++j) {
            if (find(correspond.begin(), correspond.end(), following_collapse[i][j]) !=
                correspond.end()) {
                double_following.push_back(following_collapse[i][j]);
            } else {
                correspond.push_back(following_collapse[i][j]);
            }
        }
    }
    std::vector<int> thenum(double_following.size(), 1);
    std::vector<std::vector<int>> pocess(double_following.size());
    std::vector<int> can;
    std::vector<int> no_can;
    std::vector<std::vector<int>> double_need_collapse_edge(double_following.size());
    for (int i = 0; i < double_following.size(); ++i) {
        int v = double_following[i];
        for (int j = 0; j < following_collapse.size(); ++j) {
            if (following_collapse[j].size() == 0) continue;
            for (int k = 0; k < following_collapse[j].size(); ++k) {
                if (following_collapse[j][k] == v) {
                    double_need_collapse_edge[i].push_back(edge_collapse[k * 2]);
                    double_need_collapse_edge[i].push_back(edge_collapse[k * 2 + 1]);
                    pocess[i].push_back(j);
                }
            }
        }
    }
    for (int i = 0; i < double_need_collapse_edge.size(); ++i) {
        int v1, v2, v3, v4;
        v1 = double_need_collapse_edge[i][0];
        v2 = double_need_collapse_edge[i][1];
        v3 = double_need_collapse_edge[i][2];
        v4 = double_need_collapse_edge[i][3];
        int m1, m2, m3, m4;
        if (v1 == v3) {
            m1 = v2;
            m3 = v4;
            m2 = v1;
        }
        if (v1 == v4) {
            m1 = v2;
            m3 = v3;
            m2 = v1;
        }
        if (v2 == v3) {
            m1 = v1;
            m3 = v4;
            m2 = v2;
        }
        if (v2 == v4) {
            m1 = v1;
            m3 = v3;
            m2 = v2;
        }
        int v = double_following[i];
        bool is_break = false;
        std::vector<int> quad;
        for (int j = 0; j < layer.size(); ++j) {
            for (int k = 0; k < layer[j].size(); ++k) {
                if (layer[j][k] == v) {
                    is_break = true;
                    quad = relative_point[j][k];
                    break;
                }
            }
            if (is_break) break;
        }

        for (int j = 0; j < quad.size(); ++j) {
            if (quad[j] == m2) {
                m4 = quad[(j + 2) % 4];
            }
        }
        int n0 = pocess[i][0];
        for (int k = 0; k < following_collapse[n0].size(); ++k) {
            int x = following_collapse[n0][k];
            if (x == v) continue;
            for (int l = 0; l < collapse_patch[x].size(); ++l) {
                if (collapse_patch[x][l][0] == m1 &&
                    collapse_patch[x][l][collapse_patch[x][l].size() - 1] == m4) {
                    can.push_back(x);
                    break;
                }
                if (collapse_patch[x][l][0] == m4 &&
                    collapse_patch[x][l][collapse_patch[x][l].size() - 1] == m1) {
                    can.push_back(x);
                    break;
                }
            }
        }
        int n1 = pocess[i][1];
        for (int k = 0; k < following_collapse[n1].size(); ++k) {
            int x = following_collapse[n1][k];
            if (x == v) continue;
            for (int l = 0; l < collapse_patch[x].size(); ++l) {
                if (collapse_patch[x][l][0] == m3 &&
                    collapse_patch[x][l][collapse_patch[x][l].size() - 1] == m4) {
                    no_can.push_back(x);
                    break;
                }
                if (collapse_patch[x][l][0] == m4 &&
                    collapse_patch[x][l][collapse_patch[x][l].size() - 1] == m3) {
                    no_can.push_back(x);
                    break;
                }
            }
        }
    }
    std::vector<int> relative_can(can.size());
    std::vector<int> relative_double(double_following.size());

    for (int i = 0; i < following_collapse.size(); ++i) {
        if (following_collapse[i].size() == 0) continue;
        for (int j = 0; j < following_collapse[i].size(); ++j) {
            int v = following_collapse[i][j];
            if (find(no_can.begin(), no_can.end(), v) != no_can.end()) continue;
            match_patch(i, j, v, edge_collapse, collapse_patch, patch_compact, assist_patch, can,
                        relative_can, relative_double, double_following, thenum, renew);
        }
    }
    for (int i = 0; i < renew.size(); ++i) {
        if (renew[i] == -1) continue;
        int x = renew[i];
        if (renew[x] == i) {
            renew[x] = -1;
        }
    }
    for (int i = 0; i < all_singularity.size(); ++i) {
        int v = all_singularity[i];
        if (renew[v] != -1) {
            all_singularity[i] = renew[v];
        }
    }
    // update collapse_patch
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(pre_lay.begin(), pre_lay.end(), i) != pre_lay.end()) continue;
        for (int j = 0; j < collapse_patch[i].size(); ++j) {
            for (int k = 0; k < collapse_patch[i][j].size(); ++k) {
                int v = collapse_patch[i][j][k];
                if (renew[v] != -1) {
                    collapse_patch[i][j][k] = renew[v];
                }
            }
        }
    }
    // update patch_compact
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(pre_lay.begin(), pre_lay.end(), i) != pre_lay.end()) continue;
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int v = patch_compact[i][j][k];
                if (renew[v] != -1) {
                    patch_compact[i][j][k] = renew[v];
                }
            }
        }
    }
    std::vector<std::vector<std::vector<int>>> new_patch;
    for (int i = 0; i < relative_can.size(); ++i) {
        int a = relative_can[i];
        int a1 = relative_double[i];
        std::vector<std::vector<int>> matrix_patch1 = assist_patch[a];
        std::vector<std::vector<int>> matrix_patch2 = assist_patch[a1];
        std::vector<std::vector<int>> total_matrix;
        int v1, v2, v3, v4;
        v1 = assist_patch[a][0][0];
        v2 = assist_patch[a][0][assist_patch[a][0].size() - 1];
        v3 = assist_patch[a][assist_patch[a].size() - 1]
                         [assist_patch[a][assist_patch[a].size() - 1].size() - 1];
        v4 = assist_patch[a][assist_patch[a].size() - 1][0];

        int b1, b2, b3, b4;
        b1 = assist_patch[a1][0][0];
        b2 = assist_patch[a1][0][assist_patch[a1][0].size() - 1];
        b3 = assist_patch[a1][assist_patch[a1].size() - 1]
                         [assist_patch[a1][assist_patch[a1].size() - 1].size() - 1];
        b4 = assist_patch[a1][assist_patch[a1].size() - 1][0];

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (v1 == b2 && v2 == b1) {
            for (int j = matrix_patch2.size() - 1; j > 0; --j) {
                std::vector<int> line;
                for (int k = matrix_patch2[j].size() - 1; k > -1; --k) {
                    line.push_back(matrix_patch2[j][k]);
                }
                total_matrix.push_back(line);
            }
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
        }
        if (v1 == b1 && v2 == b4) {
            for (int k = matrix_patch2[0].size() - 1; k > 0; --k) {
                std::vector<int> line;
                for (int j = 0; j < matrix_patch2.size(); ++j) {
                    line.push_back(matrix_patch2[j][k]);
                    total_matrix.push_back(line);
                }
            }
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
        }
        if (v1 == b4 && v2 == b3) {
            for (int j = 0; j < matrix_patch2.size() - 1; ++j) {
                total_matrix.push_back(matrix_patch2[j]);
            }
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
        }
        if (v1 == b3 && v2 == b2) {
            for (int k = 0; k < matrix_patch2[0].size() - 1; ++k) {
                std::vector<int> line;
                for (int j = matrix_patch2.size() - 1; j > -1; --j) {
                    line.push_back(matrix_patch2[j][k]);
                    total_matrix.push_back(line);
                }
            }
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
        }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (v3 == b4 && v4 == b3) {
            for (int j = 0; j < matrix_patch1.size() - 1; ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
            for (int j = matrix_patch2.size() - 1; j > -1; --j) {
                std::vector<int> line;
                for (int k = matrix_patch2[j].size() - 1; k > -1; --k) {
                    line.push_back(matrix_patch2[j][k]);
                }
                total_matrix.push_back(line);
            }
        }
        if (v3 == b3 && v4 == b2) {
            for (int j = 0; j < matrix_patch1.size() - 1; ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
            for (int k = matrix_patch2[0].size() - 1; k > -1; --k) {
                std::vector<int> line;
                for (int j = 0; j < matrix_patch2.size(); ++j) {
                    line.push_back(matrix_patch2[j][k]);
                    total_matrix.push_back(line);
                }
            }
        }
        if (v3 == b2 && v4 == b1) {
            for (int j = 0; j < matrix_patch1.size() - 1; ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
            for (int j = 0; j < matrix_patch2.size(); ++j) {
                total_matrix.push_back(matrix_patch2[j]);
            }
        }
        if (v3 == b1 && v4 == b4) {
            for (int j = 0; j < matrix_patch1.size() - 1; ++j) {
                total_matrix.push_back(matrix_patch1[j]);
            }
            for (int k = 0; k < matrix_patch2[0].size(); ++k) {
                std::vector<int> line;
                for (int j = matrix_patch2.size() - 1; j > -1; --j) {
                    line.push_back(matrix_patch2[j][k]);
                    total_matrix.push_back(line);
                }
            }
        }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (v2 == b1 && v3 == b4) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = 0; k < matrix_patch1[j].size() - 1; ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                for (int k = 0; k < matrix_patch2[j].size(); ++k) {
                    line.push_back(matrix_patch2[j][k]);
                }
                total_matrix.push_back(line);
            }
        }
        if (v2 == b4 && v3 == b3) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = 0; k < matrix_patch1[j].size() - 1; ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                for (int k = matrix_patch2.size() - 1; k > -1; --k) {
                    line.push_back(matrix_patch2[k][j]);
                }
                total_matrix.push_back(line);
            }
        }
        if (v2 == b3 && v3 == b2) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = 0; k < matrix_patch1[j].size() - 1; ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                for (int k = matrix_patch2[j].size() - 1; k > -1; --k) {
                    line.push_back(matrix_patch2[matrix_patch2[j].size() - 1 - j][k]);
                }
                total_matrix.push_back(line);
            }
        }
        if (v2 == b2 && v3 == b1) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = 0; k < matrix_patch1[j].size() - 1; ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                for (int k = 0; k < matrix_patch2.size(); ++k) {
                    line.push_back(matrix_patch2[k][matrix_patch2[k].size() - 1 - j]);
                }
                total_matrix.push_back(line);
            }
        }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (v4 == b3 && v1 == b2) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = 0; k < matrix_patch2[j].size() - 1; ++k) {
                    line.push_back(matrix_patch2[j][k]);
                }
                for (int k = 0; k < matrix_patch1[j].size(); ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                total_matrix.push_back(line);
            }
        }
        if (v4 == b2 && v1 == b1) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = matrix_patch2.size() - 1; k > 0; ++k) {
                    line.push_back(matrix_patch2[k][j]);
                }
                for (int k = 0; k < matrix_patch1[j].size(); ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                total_matrix.push_back(line);
            }
        }
        if (v4 == b1 && v1 == b4) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = matrix_patch2[j].size() - 1; k > 0; --k) {
                    line.push_back(matrix_patch2[matrix_patch2[j].size() - 1 - j][k]);
                }
                for (int k = 0; k < matrix_patch1[j].size(); ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                total_matrix.push_back(line);
            }
        }
        if (v4 == b4 && v1 == b3) {
            for (int j = 0; j < matrix_patch1.size(); ++j) {
                std::vector<int> line;
                for (int k = 0; k < matrix_patch2.size() - 1; ++k) {
                    line.push_back(matrix_patch2[k][matrix_patch2[k].size() - 1 - j]);
                }
                for (int k = 0; k < matrix_patch1[j].size(); ++k) {
                    line.push_back(matrix_patch1[j][k]);
                }
                total_matrix.push_back(line);
            }
        }
        new_patch.push_back(total_matrix);
    }
    std::vector<int> Vect;
    for (int i = 0; i < following_collapse.size(); ++i) {
        for (int j = 0; j < following_collapse[i].size(); ++j) {
            if (find(Vect.begin(), Vect.end(), following_collapse[i][j]) != Vect.end()) continue;
            Vect.push_back(following_collapse[i][j]);
        }
    }
    for (int i = 0; i < relative_can.size(); ++i) {
        Vect.push_back(relative_can[i]);
        Vect.push_back(relative_double[i]);
    }
    int opq = need_collapse[n];
    for (int i = 0; i < need_collapse_part[opq].size(); ++i) {
        if (i % 2 == 1) continue;
        int begin = need_collapse_part[opq][i];
        int end = need_collapse_part[opq][i + 1];
        for (int j = begin; j < end + 1; ++j) {
            Vect.push_back(layer[opq][j]);
        }
    }
    /*for (int i = 0; i < delete_patch.size(); ++i) {
        Vect.push_back(delete_patch[i]);
    }*/
    for (int i = 0; i < assist_patch.size(); ++i) {
        if (assist_patch[i].size() == 0) {
            assist_patch[i] = patch_compact[i];
        }
    }
    std::vector<std::vector<std::vector<int>>> matrix;
    for (int i = 0; i < assist_patch.size(); ++i) {
        if (find(Vect.begin(), Vect.end(), i) != Vect.end()) continue;
        if (assist_patch[i].size() == 0) continue;
        matrix.push_back(assist_patch[i]);
    }
    for (int i = 0; i < new_patch.size(); ++i) {
        matrix.push_back(new_patch[i]);
    }
    patch_compact = matrix;
    update_O(O_compact, patch_compact, all_singularity, boundary_o);

    std::sort(Vect.begin(), Vect.end());
    for (int i = 0; i < need_collapse.size(); ++i) {
        int x = need_collapse[i];
        if (layer[x].size() == 0) continue;
        for (int j = 0; j < double_q[x].size(); ++j) {
            int xx = double_q[x][j];
            if (find(pre_lay.begin(), pre_lay.end(), xx) != pre_lay.end()) {
                auto it = std::find(layer[x].begin(), layer[x].end(), xx);
                int index = std::distance(layer[x].begin(), it);
                for (int k = 0; k < need_collapse_part[x].size(); ++k) {
                    if (need_collapse_part[x][k] < index) continue;
                    if (need_collapse_part[x][k] >= index) {
                        need_collapse_part[x][k] -= 1;
                    }
                }
            }
        }
    }
    std::vector<std::vector<int>> temp_layer(layer.size());
    std::vector<std::vector<int>> d_q(double_q.size());
    for (int i = 0; i < layer.size(); ++i) {
        for (int j = 0; j < layer[i].size(); ++j) {
            if (double_q[i].size() != 0 &&
                find(double_q[i].begin(), double_q[i].end(), layer[i][j]) != double_q[i].end()) {
                if (find(pre_lay.begin(), pre_lay.end(), layer[i][j]) != pre_lay.end()) continue;
                d_q[i].push_back(layer[i][j]);
            }
            int xyz = layer[i][j];
            for (int k = 0; k < Vect.size(); ++k) {
                if (xyz < Vect[0]) {
                    temp_layer[i].push_back(xyz);
                    break;
                }
                if (k == Vect.size() - 1) {
                    if (xyz > Vect[k]) {
                        temp_layer[i].push_back(xyz - k - 1);
                    }
                } else {
                    if (xyz > Vect[k] && xyz < Vect[k + 1]) {
                        temp_layer[i].push_back(xyz - k - 1);
                        break;
                    }
                }
            }
        }
    }
    for (int i = 0; i < d_q.size(); ++i) {
        if (d_q[i].size() == 0) continue;
        for (int j = 0; j < d_q[i].size(); ++j) {
            int xyz = d_q[i][j];
            for (int k = 0; k < Vect.size(); ++k) {
                if (xyz < Vect[0]) {
                    break;
                }
                if (k == Vect.size() - 1) {
                    if (xyz > Vect[k]) {
                        d_q[i][j] = xyz - k - 1;
                    }
                } else {
                    if (xyz > Vect[k] && xyz < Vect[k + 1]) {
                        d_q[i][j] = xyz - k - 1;
                        break;
                    }
                }
            }
        }
    }
    layer = temp_layer;
    double_q = d_q;
    std::vector<std::vector<std::vector<int>>> r_c_p(patch_compact.size());
    for (int i = 0; i < patch_compact.size(); ++i) {
        r_c_p[i].push_back(patch_compact[i][0]);
        std::vector<int> line1;
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            line1.push_back(patch_compact[i][j][patch_compact[i][j].size() - 1]);
        }
        r_c_p[i].push_back(line1);
        line1.resize(0);
        for (int j = patch_compact[i][patch_compact[i].size() - 1].size() - 1; j > -1; --j) {
            line1.push_back(patch_compact[i][patch_compact[i].size() - 1][j]);
        }
        r_c_p[i].push_back(line1);
        line1.resize(0);
        for (int j = patch_compact[i].size() - 1; j > -1; --j) {
            line1.push_back(patch_compact[i][j][0]);
        }
        r_c_p[i].push_back(line1);
    }
    collapse_patch = r_c_p;
    std::vector<std::vector<std::vector<int>>> r_p;
    for (int i = 0; i < need_collapse[n] + 1; ++i) {
        std::vector<std::vector<int>> zero;
        zero.resize(0);
        r_p.push_back(zero);
    }
    for (int i = need_collapse[n] + 1; i < layer.size(); ++i) {
        if (layer[i].size() == 0) {
            std::vector<std::vector<int>> zero;
            zero.resize(0);
            r_p.push_back(zero);
            continue;
        }
        int v = layer[i][0];
        int v1, v2, v3, v4;
        int start1, start2;
        std::vector<std::vector<int>> patch_point;
        std::vector<int> single_point;
        v1 = collapse_patch[v][0][0];
        v2 = collapse_patch[v][1][0];
        v3 = collapse_patch[v][2][0];
        v4 = collapse_patch[v][3][0];
        if (layer[i].size() == 1) {
            single_point.push_back(v1);
            single_point.push_back(v2);
            single_point.push_back(v3);
            single_point.push_back(v4);
            patch_point.push_back(single_point);
            r_p.push_back(patch_point);
            continue;
        }
        std::vector<int> clearlove;
        clearlove.push_back(collapse_patch[layer[i][1]][0][0]);
        clearlove.push_back(collapse_patch[layer[i][1]][1][0]);
        clearlove.push_back(collapse_patch[layer[i][1]][2][0]);
        clearlove.push_back(collapse_patch[layer[i][1]][3][0]);
        for (int j = 0; j < 4; ++j) {
            int o1, o2, o3, o4;
            o1 = collapse_patch[v][j][0];
            o2 = collapse_patch[v][(j + 1) % 4][0];
            o3 = collapse_patch[v][(j + 2) % 4][0];
            o4 = collapse_patch[v][(j + 3) % 4][0];

            if (find(clearlove.begin(), clearlove.end(), o1) != clearlove.end() &&
                find(clearlove.begin(), clearlove.end(), o2) != clearlove.end()) {
                single_point.push_back(o3);
                single_point.push_back(o4);
                single_point.push_back(o1);
                single_point.push_back(o2);
                start1 = o1;
                start2 = o2;
                break;
            }
        }
        patch_point.push_back(single_point);
        for (int j = 1; j < layer[i].size(); ++j) {
            std::vector<int> single;
            v = layer[i][j];
            for (int k = 0; k < 4; ++k) {
                int v5, v6, v7, v8;
                v5 = collapse_patch[v][k][0];
                v6 = collapse_patch[v][(k + 1) % 4][0];
                v7 = collapse_patch[v][(k + 2) % 4][0];
                v8 = collapse_patch[v][(k + 3) % 4][0];
                if (start1 == v5 && start2 == v6) {
                    single.push_back(v5);
                    single.push_back(v6);
                    single.push_back(v7);
                    single.push_back(v8);
                    start1 = v7;
                    start2 = v8;
                    break;
                } else if (start1 == v6 && start2 == v5) {
                    single.push_back(v5);
                    single.push_back(v6);
                    single.push_back(v7);
                    single.push_back(v8);
                    start1 = v7;
                    start2 = v8;
                    break;
                }
            }
            patch_point.push_back(single);
        }
        r_p.push_back(patch_point);
    }
    relative_point = r_p;
    std::vector<int> n_c;
    std::vector<std::vector<int>> n_c_p(relative_point.size());
    for (int i = need_collapse[n] + 1; i < relative_point.size(); ++i) {
        if (relative_point[i].size() == 1) continue;
        int v1 = relative_point[i][0][0];
        int v2 = relative_point[i][0][1];
        int v3 = relative_point[i][0][2];
        int v4 = relative_point[i][0][3];
        int begin = -1;
        int start = -1;
        int x = layer[i][0];
        int c = 0, d = 0, sizec = 0, sized = 0;
        if (find(all_singularity.begin(), all_singularity.end(), v1) != all_singularity.end() ||
            find(all_singularity.begin(), all_singularity.end(), v2) != all_singularity.end()) {
            for (int j = 0; j < collapse_patch[x].size(); ++j) {
                if (collapse_patch[x][j][0] == v2 &&
                    collapse_patch[x][j][collapse_patch[x][j].size() - 1] == v3) {
                    sizec = collapse_patch[x][j].size();
                    for (int k = 0; k < collapse_patch[x][j].size(); ++k) {
                        if (boundary_o[collapse_patch[x][j][k]]) {
                            c += 1;
                        }
                    }
                }
                if (collapse_patch[x][j][0] == v4 &&
                    collapse_patch[x][j][collapse_patch[x][j].size() - 1] == v1) {
                    sized = collapse_patch[x][j].size();
                    for (int k = 0; k < collapse_patch[x][j].size(); ++k) {
                        if (boundary_o[collapse_patch[x][j][k]]) {
                            d += 1;
                        }
                    }
                }
            }
        }
        if (sizec == 2 && boundary_o[v2] && boundary_o[v3]) continue;
        if (sized == 2 && boundary_o[v4] && boundary_o[v1]) continue;
        if (c > 2) continue;
        if (d > 2) continue;
        if (find(all_singularity.begin(), all_singularity.end(), v1) != all_singularity.end()) {
            if (find(all_singularity.begin(), all_singularity.end(), v2) ==
                all_singularity.end()) {
                begin = 0;
                start = 3;
            }
        }
        if (find(all_singularity.begin(), all_singularity.end(), v2) != all_singularity.end()) {
            if (find(all_singularity.begin(), all_singularity.end(), v1) ==
                all_singularity.end()) {
                begin = 0;
                start = 4;
            }
        }
        if (begin != -1) {
            if (start == 3) {
                if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v4) ==
                        all_singularity.end()) {
                    n_c_p[i].push_back(begin);
                    n_c_p[i].push_back(begin);
                    begin = 1;
                    start = 4;
                    if (find(n_c.begin(), n_c.end(), i) == n_c.end()) {
                        n_c.push_back(i);
                    }
                }
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) ==
                        all_singularity.end()) {
                    begin = 1;
                }
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                    start = -1;
                }
            }
            if (start == 4) {
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) ==
                        all_singularity.end()) {
                    n_c_p[i].push_back(begin);
                    n_c_p[i].push_back(begin);
                    begin = 1;
                    start = 3;
                    if (find(n_c.begin(), n_c.end(), i) == n_c.end()) {
                        n_c.push_back(i);
                    }
                }
                if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v4) ==
                        all_singularity.end()) {
                    begin = 1;
                }
                if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end() &&
                    find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                    start = -1;
                }
            }
        }
        for (int j = 1; j < relative_point[i].size(); ++j) {
            if (start == -1) {
                v1 = relative_point[i][j][0];
                v2 = relative_point[i][j][1];
                v3 = relative_point[i][j][2];
                v4 = relative_point[i][j][3];
                int y = layer[i][j];
                int e = 0, f = 0, sizee = 0, sizef = 0;
                for (int jjj = 0; jjj < collapse_patch[y].size(); ++jjj) {
                    if (collapse_patch[y][jjj][0] == v2 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v3) {
                        sizee = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                e += 1;
                            }
                        }
                    }
                    if (collapse_patch[y][jjj][0] == v4 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v1) {
                        sizef = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                f += 1;
                            }
                        }
                    }
                }
                if (sizee == 2 && boundary_o[v2] && boundary_o[v3]) continue;
                if (sizef == 2 && boundary_o[v4] && boundary_o[v1]) continue;
                if (e > 2) continue;
                if (f > 2) continue;
                if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                    all_singularity.end()) {
                    if (find(all_singularity.begin(), all_singularity.end(), v2) ==
                        all_singularity.end()) {
                        begin = j;
                        start = 3;
                    }
                }
                if (find(all_singularity.begin(), all_singularity.end(), v2) !=
                    all_singularity.end()) {
                    if (find(all_singularity.begin(), all_singularity.end(), v1) ==
                        all_singularity.end()) {
                        begin = j;
                        start = 4;
                    }
                }
                if (begin != -1) {
                    if (start == 3) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v4) ==
                                all_singularity.end()) {
                            n_c_p[i].push_back(begin);
                            n_c_p[i].push_back(begin);
                            begin = j + 1;
                            start = 4;
                            if (find(n_c.begin(), need_collapse.end(), i) == n_c.end()) {
                                n_c.push_back(i);
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) ==
                                all_singularity.end()) {
                            begin = j + 1;
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end()) {
                            start = -1;
                        }
                    }
                    if (start == 4) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) ==
                                all_singularity.end()) {
                            n_c_p[i].push_back(begin);
                            n_c_p[i].push_back(begin);
                            begin = j + 1;
                            start = 3;
                            if (find(n_c.begin(), n_c.end(), i) == n_c.end()) {
                                n_c.push_back(i);
                            }
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v4) ==
                                all_singularity.end()) {
                            begin = j + 1;
                        }
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                                all_singularity.end() &&
                            find(all_singularity.begin(), all_singularity.end(), v3) !=
                                all_singularity.end()) {
                            start = -1;
                        }
                    }
                }
            } else {
                v1 = relative_point[i][j][0];
                v2 = relative_point[i][j][1];
                v3 = relative_point[i][j][2];
                v4 = relative_point[i][j][3];
                int y = layer[i][j];
                int e = 0, f = 0, sizee = 0, sizef = 0;
                for (int jjj = 0; jjj < collapse_patch[y].size(); ++jjj) {
                    if (collapse_patch[y][jjj][0] == v2 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v3) {
                        sizee = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                e += 1;
                            }
                        }
                    }
                    if (collapse_patch[y][jjj][0] == v4 &&
                        collapse_patch[y][jjj][collapse_patch[y][jjj].size() - 1] == v1) {
                        sizef = collapse_patch[y][jjj].size();
                        for (int k = 0; k < collapse_patch[y][jjj].size(); ++k) {
                            if (boundary_o[collapse_patch[y][jjj][k]]) {
                                f += 1;
                            }
                        }
                    }
                }
                if (sizee == 2 && boundary_o[v2] && boundary_o[v3]) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (sizef == 2 && boundary_o[v4] && boundary_o[v1]) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (e > 2) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (f > 2) {
                    begin = -1;
                    start = -1;
                    continue;
                }
                if (start == 3) {
                    if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) ==
                            all_singularity.end()) {
                            n_c_p[i].push_back(begin);
                            n_c_p[i].push_back(j);
                            begin = j + 1;
                            start = 4;
                            if (find(n_c.begin(), n_c.end(), i) == n_c.end()) {
                                n_c.push_back(i);
                            }
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end()) {
                            begin = -1;
                            start = -1;
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v3) ==
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                            all_singularity.end()) {
                            begin = j + 1;
                        }
                    }
                }
                if (start == 4) {
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) ==
                            all_singularity.end()) {
                            n_c_p[i].push_back(begin);
                            n_c_p[i].push_back(j);
                            begin = j + 1;
                            start = 3;
                            if (find(n_c.begin(), n_c.end(), i) == n_c.end()) {
                                n_c.push_back(i);
                            }
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) !=
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                            begin = -1;
                            start = -1;
                        }
                    }
                    if (find(all_singularity.begin(), all_singularity.end(), v4) ==
                        all_singularity.end()) {
                        if (find(all_singularity.begin(), all_singularity.end(), v3) !=
                            all_singularity.end()) {
                            begin = j + 1;
                        }
                    }
                }
            }
        }
    }
    for (int i = need_collapse[n] + 1; i < need_collapse_part.size(); ++i) {
        need_collapse_part[i] = n_c_p[i];
    }
    vector<int> aaaa;
    for (int i = n + 1; i < need_collapse.size(); ++i) {
        int abc = need_collapse[i];
        int begin, end;
        for (int j = 0; j < need_collapse_part[abc].size(); ++j) {
            if (j % 2 == 1) continue;
            begin = need_collapse_part[abc][j];
            end = need_collapse_part[abc][j + 1];
            for (int j = begin; j < end + 1; ++j) {
                if (find(aaaa.begin(), aaaa.end(), layer[abc][j]) == aaaa.end()) {
                    aaaa.push_back(layer[abc][j]);
                } else {
                    if (find(double_q[abc].begin(), double_q[abc].end(), layer[abc][j]) ==
                        double_q[abc].end()) {
                        double_q[abc].push_back(layer[abc][j]);
                    }
                }
            }
        }
    }
    int a = 10;
    /*if (n != need_collapse.size() - 1) {
        for (int i = need_collapse[n + 1]; i < need_collapse_part.size(); ++i) {
            if (need_collapse_part[i].size() != 2) continue;
            int begin = need_collapse_part[i][0];
            int end = need_collapse_part[i][1];
            int v1 = relative_point[i][begin][0];
            int a = 0;
            if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                all_singularity.end()) {
                a = 0;
            } else {
                a = 1;
            }
            if (begin == end) continue;
            for (int j = begin; j < end + 1; ++j) {
                if (find(all_singularity.begin(), all_singularity.end(),
                         relative_point[i][j][a]) != all_singularity.end()) {
                    need_collapse_part[i][0] = j;
                }
            }
        }
        for (int i = need_collapse[n + 1]; i < need_collapse_part.size(); ++i) {
            if (need_collapse_part[i].size() != 2) continue;
            int begin = need_collapse_part[i][0];
            int end = need_collapse_part[i][1];
            int v1 = relative_point[i][begin][0];
            int a = 0;
            if (find(all_singularity.begin(), all_singularity.end(), v1) !=
                all_singularity.end()) {
                a = 2;
            } else {
                a = 3;
            }
            if (begin == end) continue;
            for (int j = end; j > begin - 1; --j) {
                if (find(all_singularity.begin(), all_singularity.end(),
                         relative_point[i][j][a]) != all_singularity.end()) {
                    need_collapse_part[i][1] = j;
                }
            }
        }
    }*/
}
void Optimizer::optimize_positions_dynamic(
    Hierarchy&mRes, MatrixXi& F, MatrixXd& V, MatrixXd& N, MatrixXd& Q,
    std::vector<std::vector<int>>& Vset,
    std::vector<Vector3d>& O_compact, std::vector<Vector4i>& F_compact,
    std::vector<int>& V2E_compact, std::vector<int>& E2E_compact, double mScale,
    std::vector<Vector3d>& diffs, std::vector<int>& diff_count,
    std::map<std::pair<int, int>, int>& o2e, std::vector<int>& sharp_o,
        std::map<int, std::pair<Vector3d, Vector3d>>& compact_sharp_constraints,
        std::vector<int>& boundary_o, std::map<int, std::pair<Vector3d, Vector3d>>& compact_boundary_constraints) {
    std::set<int> uncertain;
    int with_scale=1;
    auto& cornerP = mRes.mcornerP;
    auto& corner = mRes.mcorner;
    auto& holeP1 = mRes.mholeP1;
    auto& holeP = mRes.mholeP;
    auto& hole = mRes.mhole;
    auto& holeposition = mRes.mholeposition;
    auto& cornerQ = mRes.mcornerQ;
    vector<int> unchang{};
    vector<int> temphole{};
    for (int i = 0; i < corner.size(); ++i) {
        int knum;
        double lengthe = 10000;
        Vector3d poit;
        poit[0] = cornerP[i * 3];
        poit[1] = cornerP[i * 3 + 1];
        poit[2] = cornerP[i * 3 + 2];
        for (int j = 0; j < O_compact.size(); ++j) {
            Vector3d ooo;
            ooo[0] = O_compact[j][0];
            ooo[1] = O_compact[j][1];
            ooo[2] = O_compact[j][2];
            double lengthe1 = (poit - ooo).norm();
            if (lengthe1 < lengthe) {
                knum = j;
                lengthe = lengthe1;
            }
        }
        unchang.push_back(knum);
    }
    for (int i = 0; i < unchang.size(); ++i) {
        cornerQ.push_back(unchang[i]);
    }
    
    for (auto& info : o2e) {
        if (diff_count[info.second] == 0) {
            uncertain.insert(info.first.first);
            uncertain.insert(info.first.second);
        }
    }
    std::vector<int> Vind(O_compact.size(), -1);
    std::vector<std::list<int>> links(O_compact.size());
    std::vector<std::list<int>> dedges(O_compact.size());
    std::vector<std::vector<int>> adj(V.cols());
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(j, i);
            int v2 = F((j + 1) % 3, i);
            adj[v1].push_back(v2);
        }
    }
    auto FindNearest = [&]() {
        for (int i = 0; i < O_compact.size(); ++i) {
            if (Vind[i] == -1) {
                double min_dis = 1e30;
                int min_ind = -1;
                for (auto v : Vset[i]) {
                    double dis = (V.col(v) - O_compact[i]).squaredNorm();
                    if (dis < min_dis) {
                        min_dis = dis;
                        min_ind = v;
                    }
                }
                if (min_ind > -1) {
                    Vind[i] = min_ind;
                    double x = (O_compact[i] - V.col(min_ind)).dot(N.col(min_ind));
                    O_compact[i] -= x * N.col(min_ind);
                }
            } else {
                int current_v = Vind[i];
                Vector3d n = N.col(current_v);
                double current_dis = (O_compact[i] - V.col(current_v)).squaredNorm();
                while (true) {
                    int next_v = -1;
                    for (auto& v : adj[current_v]) {
                        if (N.col(v).dot(n) < cos(10.0 / 180.0 * 3.141592654)) continue;
                        double dis = (O_compact[i] - V.col(v)).squaredNorm();
                        if (dis < current_dis) {
                            current_dis = dis;
                            next_v = v;
                        }
                    }
                    if (next_v == -1) break;
                    // rotate ideal distance
                    Vector3d n1 = N.col(current_v);
                    Vector3d n2 = N.col(next_v);
                    Vector3d axis = n1.cross(n2);
                    double len = axis.norm();
                    double angle = atan2(len, n1.dot(n2));
                    axis.normalized();
                    Matrix3d m = AngleAxisd(angle, axis).toRotationMatrix();
                    for (auto e : dedges[i]) {
                        Vector3d& d = diffs[e];
                        d = m * d;
                    }
                    current_v = next_v;
                }
                Vind[i] = current_v;
            }
        }
    };

    auto BuildConnection = [&]() {
        for (int i = 0; i < links.size(); ++i) {
            int deid0 = V2E_compact[i];
            if (deid0 != -1) {
                std::list<int>& connection = links[i];
                std::list<int>& dedge = dedges[i];
                int deid = deid0;
                do {
                    connection.push_back(F_compact[deid / 4][(deid + 1) % 4]);
                    dedge.push_back(deid);
                    deid = E2E_compact[deid / 4 * 4 + (deid + 3) % 4];
                } while (deid != -1 && deid != deid0);
                if (deid == -1) {
                    deid = deid0;
                    do {
                        deid = E2E_compact[deid];
                        if (deid == -1) break;
                        deid = deid / 4 * 4 + (deid + 1) % 4;
                        connection.push_front(F_compact[deid / 4][(deid + 1) % 4]);
                        dedge.push_front(deid);
                    } while (true);
                }
            }
        }
    };

    std::vector<Vector3d> lines;
    auto ComputeDistance = [&]() {
        std::set<int> unobserved;
        for (auto& info : o2e) {
            if (diff_count[info.second] == 0) {
                unobserved.insert(info.first.first);
            }
        }
        while (true) {
            bool update = false;
            std::set<int> observed;
            for (auto& p : unobserved) {
                std::vector<int> observations, edges;
                int count = 0;
                for (auto& e : dedges[p]) {
                    edges.push_back(e);
                    if (diff_count[e]) {
                        count += 1;
                        observations.push_back(1);
                    } else {
                        observations.push_back(0);
                    }
                }
                if (count <= 1) continue;
                update = true;
                observed.insert(p);
                for (int i = 0; i < observations.size(); ++i) {
                    if (observations[i] == 1) continue;
                    int j = i;
                    std::list<int> interp;
                    while (observations[j] == 0) {
                        interp.push_front(j);
                        j -= 1;
                        if (j < 0) j = edges.size() - 1;
                    }
                    j = (i + 1) % edges.size();
                    while (observations[j] == 0) {
                        interp.push_back(j);
                        j += 1;
                        if (j == edges.size()) j = 0;
                    }
                    Vector3d dl = diffs[edges[(interp.front() + edges.size() - 1) % edges.size()]];
                    double lenl = dl.norm();
                    Vector3d dr = diffs[edges[(interp.back() + 1) % edges.size()]];
                    double lenr = dr.norm();
                    dl /= lenl;
                    dr /= lenr;
                    Vector3d n = dl.cross(dr).normalized();
                    double angle = atan2(dl.cross(dr).norm(), dl.dot(dr));
                    if (angle < 0) angle += 2 * 3.141592654;
                    Vector3d nc = N.col(Vind[p]);
                    if (n.dot(nc) < 0) {
                        n = -n;
                        angle = 2 * 3.141592654 - angle;
                    }
                    double step = (lenr - lenl) / (interp.size() + 1);
                    angle /= interp.size() + 1;
                    Vector3d dlp = nc.cross(dl).normalized();
                    int t = 0;
                    for (auto q : interp) {
                        t += 1;
                        observations[q] = 1;
                        double ad = angle * t;
                        int e = edges[q];
                        int re = E2E_compact[e];
                        diff_count[e] = 2;
                        diffs[e] = (cos(ad) * dl + sin(ad) * dlp) * (lenl + step * t);
                        if (re != -1) {
                            diff_count[re] = 2;
                            diffs[re] = -diffs[e];
                        }
                    }
                    for (int i = 0; i < edges.size(); ++i) {
                        lines.push_back(O_compact[p]);
                        lines.push_back(O_compact[p] + diffs[edges[i]]);
                    }
                }
            }
            if (!update) break;
            for (auto& p : observed) unobserved.erase(p);
        }
    };

    BuildConnection();
    int max_iter = 10;
    for (int iter = 0; iter < max_iter; ++iter) {
        FindNearest();
        ComputeDistance();

        std::vector<std::unordered_map<int, double>> entries(O_compact.size() * 2);
        std::vector<int> fixed_dim(O_compact.size() * 2, 0);
        std::vector<int> fixed_dimS(O_compact.size() * 2, 0);
        std::vector<int> fixed_dimB(O_compact.size() * 2, 0);
        for (auto& info : compact_sharp_constraints) {
            fixed_dim[info.first * 2 + 1] = 1;
            fixed_dimS[info.first * 2 + 1] = 1;
            if (info.second.second.norm() < 0.5) {
                    fixed_dim[info.first * 2] = 1;
                    fixed_dimS[info.first * 2] = 1;
            }
            
        }
        for (auto& info : compact_boundary_constraints) {
            fixed_dim[info.first * 2 + 1] = 1;
            fixed_dimB[info.first * 2 + 1] = 1;
            if (info.second.second.norm() < 0.5) {
                    fixed_dim[info.first * 2] = 1;
                    fixed_dimB[info.first * 2] = 1;
            }
        }
        std::vector<double> b(O_compact.size() * 2);
        std::vector<double> x(O_compact.size() * 2);
        std::vector<Vector3d> Q_compact(O_compact.size());
        std::vector<Vector3d> N_compact(O_compact.size());
        std::vector<Vector3d> V_compact(O_compact.size());
#ifdef WITH_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < O_compact.size(); ++i) {
            Q_compact[i] = Q.col(Vind[i]);
            N_compact[i] = N.col(Vind[i]);
            V_compact[i] = V.col(Vind[i]);
            if (fixed_dim[i * 2 + 1] && !fixed_dim[i * 2] && fixed_dimS[i*2+1]) {
                    Q_compact[i] = compact_sharp_constraints[i].second;
                    V_compact[i] = compact_sharp_constraints[i].first;
            }
            if (fixed_dim[i * 2 + 1] && !fixed_dim[i * 2] && fixed_dimB[i * 2 + 1]) {
                    Q_compact[i] = compact_boundary_constraints[i].second;
                    V_compact[i] = compact_boundary_constraints[i].first;
            }
        }
        for (int i = 0; i < O_compact.size(); ++i) {
            Vector3d q = Q_compact[i];
            Vector3d n = N_compact[i];
            Vector3d q_y = n.cross(q);
            auto Vi = V_compact[i];
            x[i * 2] = (O_compact[i] - Vi).dot(q);
            x[i * 2 + 1] = (O_compact[i] - Vi).dot(q_y);
        }
        for (int i = 0; i < O_compact.size(); ++i) {
            Vector3d qx = Q_compact[i];
            Vector3d qy = N_compact[i];
            qy = qy.cross(qx);
            auto dedge_it = dedges[i].begin();
            for (auto it = links[i].begin(); it != links[i].end(); ++it, ++dedge_it) {
                    int j = *it;
                    Vector3d qx2 = Q_compact[j];
                    Vector3d qy2 = N_compact[j];
                    qy2 = qy2.cross(qx2);

                    int de = o2e[std::make_pair(i, j)];
                    double lambda = (diff_count[de] == 1) ? 1 : 1;
                    Vector3d target_offset = diffs[de];

                    auto Vi = V_compact[i];
                    auto Vj = V_compact[j];

                    Vector3d offset = Vj - Vi;

                    //                target_offset.normalize();
                    //                target_offset *= mScale;
                    Vector3d C = target_offset - offset;
                    int vid[] = {j * 2, j * 2 + 1, i * 2, i * 2 + 1};
                    Vector3d weights[] = {qx2, qy2, -qx, -qy};
                    for (int ii = 0; ii < 4; ++ii) {
                    for (int jj = 0; jj < 4; ++jj) {
                        auto it = entries[vid[ii]].find(vid[jj]);
                        if (it == entries[vid[ii]].end()) {
                            entries[vid[ii]][vid[jj]] = lambda * weights[ii].dot(weights[jj]);
                        } else {
                            entries[vid[ii]][vid[jj]] += lambda * weights[ii].dot(weights[jj]);
                        }
                    }
                    b[vid[ii]] += lambda * weights[ii].dot(C);
                    }
            }
        }

        // fix sharp edges
        for (int i = 0; i < entries.size(); ++i) {
            if (entries[i].size() == 0) {
                    entries[i][i] = 1;
                    b[i] = x[i];
            }
            if (fixed_dim[i]) {
                    b[i] = x[i];
                    entries[i].clear();
                    entries[i][i] = 1;
            } else {
                    std::unordered_map<int, double> newmap;
                    for (auto& rec : entries[i]) {
                    if (fixed_dim[rec.first]) {
                        b[i] -= rec.second * x[rec.first];
                    } else {
                        newmap[rec.first] = rec.second;
                    }
                    }
                    std::swap(entries[i], newmap);
            }
        }
        std::vector<Eigen::Triplet<double>> lhsTriplets;
        lhsTriplets.reserve(F_compact.size() * 8);
        Eigen::SparseMatrix<double> A(O_compact.size() * 2, O_compact.size() * 2);
        VectorXd rhs(O_compact.size() * 2);
        rhs.setZero();
        for (int i = 0; i < entries.size(); ++i) {
            rhs(i) = b[i];
            for (auto& rec : entries[i]) {
                    lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, rec.second));
            }
        }

        A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());

#ifdef LOG_OUTPUT
        int t1 = GetCurrentTime64();
#endif

        // FIXME: IncompleteCholesky Preconditioner will fail here so I fallback to Diagonal one.
        // I suspected either there is a implementation bug in IncompleteCholesky Preconditioner
        // or there is a memory corruption somewhere.  However, g++'s address sanitizer does not
        // report anything useful.
        LinearSolver<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        //        Eigen::setNbThreads(1);
        //        ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;
        //        VectorXd x0 = VectorXd::Map(x.data(), x.size());
        //        solver.setMaxIterations(40);

        //        solver.compute(A);
        VectorXd x_new = solver.solve(rhs);  // solver.solveWithGuess(rhs, x0);

#ifdef LOG_OUTPUT
        // std::cout << "[LSQ] n_iteration:" << solver.iterations() << std::endl;
        // std::cout << "[LSQ] estimated error:" << solver.error() << std::endl;
        int t2 = GetCurrentTime64();
        printf("[LSQ] Linear solver uses %lf seconds.\n", (t2 - t1) * 1e-3);
#endif
        for (int i = 0; i < O_compact.size(); ++i) {
            // Vector3d q = Q.col(Vind[i]);
            Vector3d q = Q_compact[i];
            // Vector3d n = N.col(Vind[i]);
            Vector3d n = N_compact[i];
            Vector3d q_y = n.cross(q);
            auto Vi = V_compact[i];
            if (compact_boundary_constraints.count(i) == 0 &&
                compact_sharp_constraints.count(i) == 0) {
                    O_compact[i] = Vi + q * x_new[i * 2] + q_y * x_new[i * 2 + 1]; 
            }
        }

        // forgive my hack...
        if (iter + 1 == max_iter) {
            for (int iter = 0; iter < 5; ++iter) {
                    for (int i = 0; i < O_compact.size(); ++i) {
                    if (sharp_o[i] || boundary_o[i]) continue;
                    if (dedges[i].size() != 4 || uncertain.count(i)) {
                        Vector3d n(0, 0, 0), v(0, 0, 0);
                        Vector3d v0 = O_compact[i];
                        for (auto e : dedges[i]) {
                            Vector3d v1 = O_compact[F_compact[e / 4][(e + 1) % 4]];
                            Vector3d v2 = O_compact[F_compact[e / 4][(e + 3) % 4]];
                            n += (v1 - v0).cross(v2 - v0);
                            v += v1;
                        }
                        n.normalize();
                        Vector3d offset = v / dedges[i].size() - v0;
                        offset -= offset.dot(n) * n;
                        O_compact[i] += offset;
                    }
                    }
            }
        }
    }
    /*std::vector<std::vector<int>> adj(V.cols());
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(j, i);
            int v2 = F((j + 1) % 3, i);
            adj[v1].push_back(v2);
        }
    }*/
    std::vector<std::vector<int>> adj_compact(O_compact.size());
    for (int i = 0; i < F_compact.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j == 3) {
                    int v1 = F_compact[i][3];
                    int v2 = F_compact[i][0];
                    adj_compact[v1].push_back(v2);
            } else {
                    int v1 = F_compact[i][j];
                    int v2 = F_compact[i][j + 1];
                    adj_compact[v1].push_back(v2);
            }
        }
    }
    for (int iteration = 0; iteration < 100; ++iteration) {
        for (int i = 0; i < adj_compact.size(); ++i) {
            if (compact_boundary_constraints.count(i) == 1 ||
                compact_sharp_constraints.count(i) == 1)
                    continue;
            Vector3d position;
            position[0] = 0;
            position[1] = 0;
            position[2] = 0;
            Vector3d normal = N.col(Vset[i][0]);
            for (auto& v : adj_compact[i]) {
                    position += O_compact[v];
            }
            position /= adj_compact[i].size();
            Vector3d edge1 = position - O_compact[i];
            O_compact[i] += edge1 - edge1.dot(normal) * normal;
        }
    }
    for (int i = 0; i < F_compact.size(); ++i) {
        int nop = 0;
        int v1, v2, v3, v4;
        v1 = F_compact[i][0];
        v2 = F_compact[i][1];
        v3 = F_compact[i][2];
        v4 = F_compact[i][3];
        if (compact_boundary_constraints.count(v1)) {
            nop += 1;
        }
        if (compact_boundary_constraints.count(v2)) {
            nop += 1;
        }
        if (compact_boundary_constraints.count(v3)) {
            nop += 1;
        }
        if (compact_boundary_constraints.count(v4)) {
            nop += 1;
        }
        if (nop > 3) {
            F_compact[i][0] = -1;
            F_compact[i][1] = -1;
            F_compact[i][2] = -1;
            F_compact[i][3] = -1;
        }
    }
        
    for (int i = 0; i < unchang.size(); ++i) {
        O_compact[unchang[i]][0] = cornerP[i * 3];
        O_compact[unchang[i]][1] = cornerP[i * 3 + 1];
        O_compact[unchang[i]][2] = cornerP[i * 3 + 2];
    } 
            
}
    


void Optimizer::optimize_positions_sharp(
    Hierarchy& mRes, std::vector<DEdge>& edge_values, std::vector<Vector2i>& edge_diff,
    std::vector<int>& sharp_edges, std::set<int>& sharp_vertices,
    std::map<int, std::pair<Vector3d, Vector3d>>& sharp_constraints,
    std::vector<int>& Boundary_edges, std::set<int>& Boundary_vertices,
    std::map<int, std::pair<Vector3d, Vector3d>>& Boundary_constraints, int with_scale) {
    auto& V = mRes.mV[0];
    auto& F = mRes.mF;
    auto& Q = mRes.mQ[0];
    auto& N = mRes.mN[0];
    auto& O = mRes.mO[0];
    auto& S = mRes.mS[0];
    auto& corner = mRes.mcorner;
    auto& cornerP = mRes.mcornerP;
    auto& hole = mRes.mhole;
    auto& unchangeB = mRes.munchangeB;
    auto& linetype = mRes.mlinetype;

    DisajointTree tree(V.cols());
    for (int i = 0; i < edge_diff.size(); ++i) {
        if (edge_diff[i].array().abs().sum() == 0) {
            tree.Merge(edge_values[i].x, edge_values[i].y);
        }
    }
    tree.BuildCompactParent();
    std::map<int, int> compact_sharp_indices;
    std::set<DEdge> compact_sharp_edges;
    std::map<int, int> compact_boundary_indices;
    std::set<DEdge> compact_boundary_edges;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 1) {
            int v1 = tree.Index(F(i % 3, i / 3));
            int v2 = tree.Index(F((i + 1) % 3, i / 3));
            compact_sharp_edges.insert(DEdge(v1, v2));
        }
    }
    for (int i = 0; i < Boundary_edges.size(); ++i) {
        if (Boundary_edges[i] == 1) {
            int v1 = tree.Index(F(i % 3, i / 3));
            int v2 = tree.Index(F((i + 1) % 3, i / 3));
            compact_boundary_edges.insert(DEdge(v1, v2));
        }
    }
    for (auto& v : sharp_vertices) {
        int p = tree.Index(v);
        if (compact_sharp_indices.count(p) == 0) {
            int s = compact_sharp_indices.size();
            compact_sharp_indices[p] = s;
        }
    }
    for (auto& v : Boundary_vertices) {
        int p = tree.Index(v);
        if (compact_boundary_indices.count(p) == 0) {
            int s = compact_boundary_indices.size();
            compact_boundary_indices[p] = s;
        }
    }
    std::map<int, std::set<int>> sharp_vertices_links;
    std::set<DEdge> sharp_dedges;
    std::map<int, std::set<int>> boundary_vertices_links;
    std::set<DEdge> boundary_dedges;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i]) {
            int v1 = F(i % 3, i / 3);
            int v2 = F((i + 1) % 3, i / 3);
            if (sharp_vertices_links.count(v1) == 0) sharp_vertices_links[v1] = std::set<int>();
            sharp_vertices_links[v1].insert(v2);
            sharp_dedges.insert(DEdge(v1, v2));
        }
    }
    for (int i = 0; i < Boundary_edges.size(); ++i) {
        if (Boundary_edges[i]) {
            int v1 = F(i % 3, i / 3);
            int v2 = F((i + 1) % 3, i / 3);
            if (boundary_vertices_links.count(v1) == 0)
                boundary_vertices_links[v1] = std::set<int>();
            boundary_vertices_links[v1].insert(v2);
            boundary_dedges.insert(DEdge(v1, v2));
        }
    }
    std::vector<std::vector<int>> sharp_to_original_indices(compact_sharp_indices.size());
    std::vector<std::vector<int>> boundary_to_original_indices(compact_boundary_indices.size());
    for (auto& v : sharp_vertices_links) {
        if (v.second.size() == 2) continue;
        int p = tree.Index(v.first);
        sharp_to_original_indices[compact_sharp_indices[p]].push_back(v.first);
    }
    for (auto& v : sharp_vertices_links) {
        if (v.second.size() != 2) continue;
        int p = tree.Index(v.first);
        sharp_to_original_indices[compact_sharp_indices[p]].push_back(v.first);
    }
    for (auto& v : boundary_vertices_links) {
        if (v.second.size() == 2) continue;
        int p = tree.Index(v.first);
        boundary_to_original_indices[compact_boundary_indices[p]].push_back(v.first);
    }
    for (auto& v : boundary_vertices_links) {
        if (v.second.size() != 2) continue;
        int p = tree.Index(v.first);
        boundary_to_original_indices[compact_boundary_indices[p]].push_back(v.first);
    }
    for (int i = 0; i < V.cols(); ++i) {
        if (sharp_vertices.count(i)) continue;
        int p = tree.Index(i);
        if (compact_sharp_indices.count(p))
            sharp_to_original_indices[compact_sharp_indices[p]].push_back(i);
    }
    for (int i = 0; i < V.cols(); ++i) {
        if (Boundary_vertices.count(i)) continue;
        int p = tree.Index(i);
        if (compact_boundary_indices.count(p))
            boundary_to_original_indices[compact_boundary_indices[p]].push_back(i);
    }

    int num = sharp_to_original_indices.size();
    int num1 = boundary_to_original_indices.size();
    std::vector<std::set<int>> links(sharp_to_original_indices.size());
    std::vector<std::set<int>> links1(boundary_to_original_indices.size());
    for (int e = 0; e < edge_diff.size(); ++e) {
        int v1 = edge_values[e].x;
        int v2 = edge_values[e].y;
        int p1 = tree.Index(v1);
        int p2 = tree.Index(v2);
        if (p1 == p2 || compact_sharp_edges.count(DEdge(p1, p2)) == 0) continue;
        p1 = compact_sharp_indices[p1];
        p2 = compact_sharp_indices[p2];

        links[p1].insert(p2);
        links[p2].insert(p1);
    }
    for (int e = 0; e < edge_diff.size(); ++e) {
        int v1 = edge_values[e].x;
        int v2 = edge_values[e].y;
        int p1 = tree.Index(v1);
        int p2 = tree.Index(v2);
        if (p1 == p2 || compact_boundary_edges.count(DEdge(p1, p2)) == 0) continue;
        p1 = compact_boundary_indices[p1];
        p2 = compact_boundary_indices[p2];

        links1[p1].insert(p2);
        links1[p2].insert(p1);
    }
    std::vector<int> alignP{};
    Vector3d posision90;
    for (int i = 0; i < sharp_to_original_indices.size(); ++i) {
        posision90 = O.col(sharp_to_original_indices[i][0]);
        double testnumber = 10000;
        int testpoint = 0;
        for (auto& v : sharp_vertices) {
            Vector3d position80 = V.col(v);
            double lengthh = (position80 - posision90).norm();
            if (lengthh < testnumber) {
                testnumber = lengthh;
                testpoint = v;
            }
        }
        alignP.push_back(testpoint);
    }
    for (int i = 0; i < alignP.size(); ++i) {
        for (auto& v : sharp_to_original_indices[i]) {
            O.col(v) = V.col(alignP[i]);
        }
    }
    std::vector<int> hash(links.size(), 0);
    std::vector<std::vector<Vector3d>> loops;
    for (int i = 0; i < num; ++i) {
        if (hash[i] == 1) continue;
        if (links[i].size() == 2) {
            std::vector<int> q;
            q.push_back(i);
            hash[i] = 1;
            int v = i;
            int prev_v = -1;
            bool is_loop = false;
            while (links[v].size() == 2) {
                int next_v = -1;
                for (auto nv : links[v])
                    if (nv != prev_v) next_v = nv;
                if (hash[next_v]) {
                    is_loop = true;
                    break;
                }
                if (links[next_v].size() == 2) hash[next_v] = true;
                q.push_back(next_v);
                prev_v = v;
                v = next_v;
            }
            if (!is_loop && q.size() >= 2) {
                std::vector<int> q1;
                int v = i;
                int prev_v = q[1];
                while (links[v].size() == 2) {
                    int next_v = -1;
                    for (auto nv : links[v])
                        if (nv != prev_v) next_v = nv;
                    if (hash[next_v]) {
                        is_loop = true;
                        break;
                    }
                    if (links[next_v].size() == 2) hash[next_v] = true;
                    q1.push_back(next_v);
                    prev_v = v;
                    v = next_v;
                }
                std::reverse(q1.begin(), q1.end());
                q1.insert(q1.end(), q.begin(), q.end());
                std::swap(q1, q);
            }
            if (q.size() < 3) continue;
            if (is_loop) q.push_back(q.front());
            double len = 0, scale = 0;
            std::vector<Vector3d> o(q.size()), new_o(q.size());
            std::vector<double> sc(q.size());

            for (int i = 0; i < q.size() - 1; ++i) {
                int v1 = q[i];
                int v2 = q[i + 1];
                auto it = links[v1].find(v2);
                if (it == links[v1].end()) {
                    printf("Non exist!\n");
                    exit(0);
                }
            }

            for (int i = 0; i < q.size(); ++i) {
                if (sharp_to_original_indices[q[i]].size() == 0) {
                    continue;
                }
                o[i] = O.col(sharp_to_original_indices[q[i]][0]);
                Vector3d qx = Q.col(sharp_to_original_indices[q[i]][0]);
                Vector3d qy = Vector3d(N.col(sharp_to_original_indices[q[i]][0])).cross(qx);
                int fst = sharp_to_original_indices[q[1]][0];
                Vector3d dis = (i == 0) ? (Vector3d(O.col(fst)) - o[i]) : o[i] - o[i - 1];
                if (with_scale)
                    sc[i] = (abs(qx.dot(dis)) > abs(qy.dot(dis)))
                                ? S(0, sharp_to_original_indices[q[i]][0])
                                : S(1, sharp_to_original_indices[q[i]][0]);
                else
                    sc[i] = 1;
                new_o[i] = o[i];
            }

            if (is_loop) {
                for (int i = 0; i < q.size(); ++i) {
                    Vector3d dir =
                        (o[(i + 1) % q.size()] - o[(i + q.size() - 1) % q.size()]).normalized();
                    for (auto& ind : sharp_to_original_indices[q[i]]) {
                        sharp_constraints[ind] = std::make_pair(o[i], dir);
                    }
                }
            } else {
                for (int i = 0; i < q.size(); ++i) {
                    Vector3d dir(0, 0, 0);
                    if (i != 0 && i + 1 != q.size())
                        dir = (o[i + 1] - o[i - 1]).normalized();
                    else if (links[q[i]].size() == 1) {
                        if (i == 0)
                            dir = (o[i + 1] - o[i]).normalized();
                        else
                            dir = (o[i] - o[i - 1]).normalized();
                    }
                    for (auto& ind : sharp_to_original_indices[q[i]]) {
                        sharp_constraints[ind] = std::make_pair(o[i], dir);
                    }
                }
            }

            for (int i = 0; i < q.size() - 1; ++i) {
                len += (o[i + 1] - o[i]).norm();
                scale += sc[i];
            }

            int next_m = q.size() - 1;

            double left_norm = len * sc[0] / scale;
            int current_v = 0;
            double current_norm = (o[1] - o[0]).norm();
            for (int i = 1; i < next_m; ++i) {
                while (left_norm >= current_norm) {
                    left_norm -= current_norm;
                    current_v += 1;
                    current_norm = (o[current_v + 1] - o[current_v]).norm();
                }
                new_o[i] =
                    (o[current_v + 1] * left_norm + o[current_v] * (current_norm - left_norm)) /
                    current_norm;
                o[current_v] = new_o[i];
                current_norm -= left_norm;
                left_norm = len * sc[current_v] / scale;
            }

            for (int i = 0; i < q.size(); ++i) {
                for (auto v : sharp_to_original_indices[q[i]]) {
                    O.col(v) = new_o[i];
                }
            }

            loops.push_back(new_o);
        }
    }
    Vector3d Point1;
    for (int i = 0; i < corner.size(); ++i) {
        Point1[0] = cornerP[3 * i];
        Point1[1] = cornerP[3 * i + 1];
        Point1[2] = cornerP[3 * i + 2];
        double lengther = 100000;
        int mathma = 0;
        for (int j = 0; j < boundary_to_original_indices.size(); ++j) {
            Vector3d Point = O.col(boundary_to_original_indices[j][0]);
            double lengther1 = (Point1 - Point).norm();
            if (lengther1 < lengther) {
                lengther = lengther1;
                mathma = j;
            }
        }
        unchangeB.push_back(mathma);
    }
    for (int i = 0; i < unchangeB.size(); ++i) {
        double wilx1 = cornerP[3 * i];
        double wilx2 = cornerP[3 * i + 1];
        double wilx3 = cornerP[3 * i + 2];
        for (auto& wilx : boundary_to_original_indices[unchangeB[i]]) {
            O.col(wilx)[0] = wilx1;
            O.col(wilx)[1] = wilx2;
            O.col(wilx)[2] = wilx3;
        }
    }
    std::vector<int> boundaryunchangeP{};
    for (int i = 0; i < unchangeB.size(); ++i) {
        for (auto& wilx4 : boundary_to_original_indices[unchangeB[i]]) {
            boundaryunchangeP.push_back(wilx4);
        }
    }
    std::vector<int> alignB{};
    Vector3d posision250;
    for (int i = 0; i < boundary_to_original_indices.size(); ++i) {
        if (count(unchangeB.begin(), unchangeB.end(), i)) {
            alignB.push_back(-1);
            continue;
        }
        posision250 = O.col(boundary_to_original_indices[i][0]);
        double testnumber1 = 10000;
        int testpoint1 = 0;
        for (auto& v : Boundary_vertices) {
            if (count(boundaryunchangeP.begin(), boundaryunchangeP.end(), v)) {
                continue;
            }
            Vector3d position251 = V.col(v);
            double lengthh1 = (position251 - posision250).norm();
            if (lengthh1 < testnumber1) {
                testnumber1 = lengthh1;
                testpoint1 = v;
            }
        }
        alignB.push_back(testpoint1);
    }
    for (int i = 0; i < alignB.size(); ++i) {
        if (alignB[i] == -1) continue;
        for (auto& v : boundary_to_original_indices[i]) {
            O.col(v) = V.col(alignB[i]);
        }
    }
    linetype.resize(boundary_to_original_indices.size(), 0); 
    std::vector<int> hash1(links1.size(), 0);
    std::vector<std::vector<Vector3d>> loops1;
    for (int i = 0; i < num1; ++i) {
        if (hash1[i] == 1) continue;
        if (links1[i].size() == 2) {
            std::vector<int> q5;
            q5.push_back(i);
            hash1[i] = 1;
            int v1 = i;
            int prev_v1 = -1;
            bool is_loop1 = false;
            while (links1[v1].size() == 2) {
                int next_v1 = -1;
                for (auto nv : links1[v1])
                    if (nv != prev_v1) next_v1 = nv;
                if (hash1[next_v1]) {
                    is_loop1 = true;
                    break;
                }
                if (links1[next_v1].size() == 2) hash1[next_v1] = true;
                q5.push_back(next_v1);
                prev_v1 = v1;
                v1 = next_v1;
            }
            if (!is_loop1 && q5.size() >= 2) {
                std::vector<int> q6;
                int v6 = i;
                int prev_v1 = q5[1];
                while (links1[v1].size() == 2) {
                    int next_v1 = -1;
                    for (auto nv : links1[v1])
                        if (nv != prev_v1) next_v1 = nv;
                    if (hash1[next_v1]) {
                        is_loop1 = true;
                        break;
                    }
                    if (links1[next_v1].size() == 2) hash1[next_v1] = true;
                    q6.push_back(next_v1);
                    prev_v1 = v1;
                    v1 = next_v1;
                }
                std::reverse(q6.begin(), q6.end());
                q6.insert(q6.end(), q5.begin(), q5.end());
                std::swap(q6, q5);
            }
            if (q5.size() < 3) continue;
            if (q5.size() < num / 4) {
                for (int edg = 0; edg < q5.size(); ++edg) {
                    hole.push_back(q5[edg]);
                }
            }
            if (is_loop1) q5.push_back(q5.front());
            double len1 = 0, scale1 = 0;
            std::vector<Vector3d> o1(q5.size()), new_o1(q5.size());
            std::vector<double> sc1(q5.size());

            for (int i = 0; i < q5.size() - 1; ++i) {
                int v11 = q5[i];
                int v22 = q5[i + 1];
                auto it = links1[v11].find(v22);
                if (it == links1[v11].end()) {
                    printf("Non exist!\n");
                    exit(0);
                }
            }

            for (int i = 0; i < q5.size(); ++i) {
                if (boundary_to_original_indices[q5[i]].size() == 0) {
                    continue;
                }
                o1[i] = O.col(boundary_to_original_indices[q5[i]][0]);
                Vector3d qx1 = Q.col(boundary_to_original_indices[q5[i]][0]);
                Vector3d qy1 = Vector3d(N.col(boundary_to_original_indices[q5[i]][0])).cross(qx1);
                int fst1 = boundary_to_original_indices[q5[1]][0];
                Vector3d dis1 = (i == 0) ? (Vector3d(O.col(fst1)) - o1[i]) : o1[i] - o1[i - 1];
                if (with_scale)
                    sc1[i] = (abs(qx1.dot(dis1)) > abs(qy1.dot(dis1)))
                                 ? S(0, boundary_to_original_indices[q5[i]][0])
                                 : S(1, boundary_to_original_indices[q5[i]][0]);
                else
                    sc1[i] = 1;
                new_o1[i] = o1[i];
            }

            if (is_loop1) {
                for (int i = 0; i < q5.size(); ++i) {
                    Vector3d dir1, dir2, dir3, dir4, dir5;
                    if (i == 0) {
                        dir1 = (o1[i + 1] - o1[i]).normalized();
                        dir2 = (o1[q5.size() - 2] - o1[i]).normalized();
                    } 
                    if (i == q5.size() - 1) continue;
                    if (i != 0 && i != q5.size() - 1){
                        dir1 = (o1[i + 1] - o1[i]).normalized();
                        dir2 = (o1[i - 1] - o1[i]).normalized();
                    }
                    dir3 = (dir1 + dir2).normalized();
                    
                    if (dir3.norm()==0) {
                        dir4 = dir1;
                        linetype[q5[i]] = 1;
                    } else {
                        if (std::count(unchangeB.begin(), unchangeB.end(), q5[i])) {
                            dir4 = dir1;
                            linetype[q5[i]] = 1;
                        } else {
                            dir5 = N.col(boundary_to_original_indices[q5[i]][0]);
                            dir4 = dir5.cross(dir3);
                        }
                    }
                    for (auto& ind : boundary_to_original_indices[q5[i]]) {
                        Boundary_constraints[ind] = std::make_pair(o1[i], dir4);
                    }
                }
            } else {
                for (int i = 0; i < q5.size(); ++i) {
                    Vector3d dir1(0, 0, 0);
                    if (i != 0 && i + 1 != q5.size())
                        dir1 = (o1[i + 1] - o1[i - 1]).normalized();
                    else if (links1[q5[i]].size() == 1) {
                        if (i == 0)
                            dir1 = (o1[i + 1] - o1[i]).normalized();
                        else
                            dir1 = (o1[i] - o1[i - 1]).normalized();
                    }
                    for (auto& ind : boundary_to_original_indices[q5[i]]) {
                        Boundary_constraints[ind] = std::make_pair(o1[i], dir1);
                    }
                }
            }

            for (int i = 0; i < q5.size() - 1; ++i) {
                len1 += (o1[i + 1] - o1[i]).norm();
                scale1 += sc1[i];
            }

            int next_m1 = q5.size() - 1;

            double left_norm1 = len1 * sc1[0] / scale1;
            int current_v1 = 0;
            double current_norm1 = (o1[1] - o1[0]).norm();
            for (int i = 1; i < next_m1; ++i) {
                while (left_norm1 >= current_norm1) {
                    left_norm1 -= current_norm1;
                    current_v1 += 1;
                    current_norm1 = (o1[current_v1 + 1] - o1[current_v1]).norm();
                }

                if (std::count(unchangeB.begin(), unchangeB.end(), q5[current_v1])) {
                    current_norm1 -= left_norm1;
                    left_norm1 = len1 * sc1[current_v1] / scale1;

                } else {
                    new_o1[i] = (o1[current_v1 + 1] * left_norm1 +
                                 o1[current_v1] * (current_norm1 - left_norm1)) /
                                current_norm1;
                    o1[current_v1] = new_o1[i];
                    current_norm1 -= left_norm1;
                    left_norm1 = len1 * sc1[current_v1] / scale1;
                }
            }

            for (int i = 0; i < q5.size(); ++i) {
                if (std::count(unchangeB.begin(), unchangeB.end(), q5[i])) {
                    continue;
                } else {
                    for (auto v : boundary_to_original_indices[q5[i]]) {
                        Vector3d dirdir = new_o1[i] - O.col(v);
                        Vector3d dirdir1 = Boundary_constraints[v].second;
                        Vector3d dirdir3 =
                            dirdir1 * (dirdir.dot(dirdir1) /
                                       sqrt(dirdir1[0] * dirdir1[0] + dirdir1[1] * dirdir1[1] +
                                            dirdir1[2] * dirdir1[2]));
                        O.col(v) += dirdir3;
                    }
                }
            }

            loops1.push_back(new_o1);
        }
    }
    return;
        std::ofstream os("/Users/jingwei/Desktop/sharp.obj");
    for (int i = 0; i < loops.size(); ++i) {
        for (auto& v : loops[i]) {
            os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
    }
        int offset = 1;
    for (int i = 0; i < loops.size(); ++i) {
        for (int j = 0; j < loops[i].size() - 1; ++j) {
            os << "l " << offset + j << " " << offset + j + 1 << "\n";
        }
        offset += loops[i].size();
    }
    os.close();
    exit(0);
}


void Optimizer::collapse_boundary(std::vector<Vector3d>& O_compact,
                                  std::vector<std::vector<std::vector<int>>>& patch_compact,
                                  std::vector<int>& valence,
                                  std::vector<std::vector<std::vector<int>>>& relative_point,
                                  std::vector<std::vector<int>>& layer,
                                  std::vector<std::vector<std::vector<int>>>& collapse_patch,
                                  std::vector<int>& all_singularity) {
    std::vector<int> position;
    std::vector<int> renew(O_compact.size(), -1);
    std::vector<int> delete_patch;
    for (int i = 0; i < layer.size(); ++i) {
        cout << layer[i].size() << "\n";
        if (i != 2) continue;
        for (int j = 0; j < layer[i].size(); ++j) {
            delete_patch.push_back(layer[i][j]);
        }
    }
    for (int i = 0; i < relative_point.size(); ++i) {
        if (i != 2) continue;
        for (int j = 0; j < relative_point[i].size(); ++j) {
            int v1, v2;
            v1 = relative_point[i][j][3];
            v2 = relative_point[i][j][0];
            int p = layer[i][j];
            std::vector<int> line1, line2;
            for (int k = 0; k < collapse_patch[p].size(); ++k) {
                if (collapse_patch[p][k][0] == v1 &&
                    collapse_patch[p][k][collapse_patch[p][k].size() - 1] == v2) {
                    line1 = collapse_patch[p][k];
                    line2 = collapse_patch[p][(k + 2) % 4];
                    break;
                }
            }
            for (int k = 0; k < line1.size(); ++k) {
                int o1 = line1[k];
                int o2 = line2[line2.size() - 1 - k];
                O_compact[o1] = O_compact[o2];
                if (renew[o2] == -1) {
                    renew[o2] = o1;
                }
            }
        }
    }
    std::sort(delete_patch.begin(), delete_patch.end());
    // update information
    // update collapse_patch
    std::vector<std::vector<std::vector<int>>> temp_collapse_patch;
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < collapse_patch[i].size(); ++j) {
            for (int k = 0; k < collapse_patch[i][j].size(); ++k) {
                int v = collapse_patch[i][j][k];
                if (renew[v] != -1) {
                    collapse_patch[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < collapse_patch.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_collapse_patch.push_back(collapse_patch[i]);
    }
    collapse_patch = temp_collapse_patch;
    // update patch_compact
    std::vector<std::vector<std::vector<int>>> temp_patch_compact;
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        for (int j = 0; j < patch_compact[i].size(); ++j) {
            for (int k = 0; k < patch_compact[i][j].size(); ++k) {
                int v = patch_compact[i][j][k];
                if (renew[v] != -1) {
                    patch_compact[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < patch_compact.size(); ++i) {
        if (find(delete_patch.begin(), delete_patch.end(), i) != delete_patch.end()) continue;
        temp_patch_compact.push_back(patch_compact[i]);
    }
    patch_compact = temp_patch_compact;
    // update relative_point
    std::vector<std::vector<std::vector<int>>> temp_relative_point(relative_point.size());
    std::vector<std::vector<int>> temp_layer(layer.size());
    for (int i = 0; i < relative_point.size(); ++i) {
        for (int j = 0; j < relative_point[i].size(); ++j) {
            for (int k = 0; k < relative_point[i][j].size(); ++k) {
                int v = relative_point[i][j][k];
                if (renew[v] != -1) {
                    relative_point[i][j][k] = renew[v];
                }
            }
        }
    }
    for (int i = 0; i < layer.size(); ++i) {
        if (find(position.begin(), position.end(), i) != position.end()) continue;
        for (int j = 0; j < layer[i].size(); ++j) {
            int v = layer[i][j];
            if (find(delete_patch.begin(), delete_patch.end(), v) != delete_patch.end()) continue;
            temp_relative_point[i].push_back(relative_point[i][j]);
            for (int k = 0; k < delete_patch.size(); ++k) {
                if (v < delete_patch[0]) {
                    temp_layer[i].push_back(v);
                    break;
                }
                if (k == delete_patch.size() - 1) {
                    if (v > delete_patch[k]) {
                        temp_layer[i].push_back(v - k - 1);
                    }
                } else {
                    if (v > delete_patch[k] && v < delete_patch[k + 1]) {
                        temp_layer[i].push_back(v - k - 1);
                        break;
                    }
                }
            }
        }
    }
    std::cout << "finish"
              << "\n";
    std::vector<std::vector<std::vector<int>>> temp;
    std::vector<std::vector<int>> temp1;
    for (int i = 0; i < temp_relative_point.size(); ++i) {
        if (temp_relative_point[i].size() == 0) continue;
        temp.push_back(temp_relative_point[i]);
    }
    relative_point = temp;
    for (int i = 0; i < temp_layer.size(); ++i) {
        if (temp_layer[i].size() == 0) continue;
        temp1.push_back(temp_layer[i]);
    }
    layer = temp1;
}
void Optimizer::optimize_positions_fixed(
    Hierarchy& mRes, std::vector<DEdge>& edge_values, std::vector<Vector2i>& edge_diff,
    std::set<int>& sharp_vertices, std::map<int, std::pair<Vector3d, Vector3d>>& sharp_constraints,
    std::vector<int>& sharp_edges, std::set<int>& Boundary_vertices,
    std::map<int, std::pair<Vector3d, Vector3d>>& Boundary_constraints, std::vector<int>& Boundary_edges,
    int with_scale) {
    auto& V = mRes.mV[0];
    auto& F = mRes.mF;
    auto& Q = mRes.mQ[0];
    auto& N = mRes.mN[0];
    auto& O = mRes.mO[0];
    auto& S = mRes.mS[0];
    auto& corner = mRes.mcorner;
    auto& hole = mRes.mhole;
    auto& holeP = mRes.mholeP;
    auto& holeposition = mRes.mholeposition;
    auto& linetype = mRes.mlinetype;
    auto& unchangeB = mRes.munchangeB;

    
    /*std::vector<double> holeposition1{};
    for (int i = 0; i < V.cols(); ++i) {
        double a1 = V.col(i)[0];
        double a2 = V.col(i)[1];
        double a3 = V.col(i)[2];
        holeposition1.push_back(a1);
        holeposition1.push_back(a2);
        holeposition1.push_back(a3);
    }
    for (int i = 0; i < holeposition.size(); ++i) {
        if (holeposition[i] != holeposition1[i]) {
            int asasa = 1;
        }
    }*/
    

    DisajointTree tree(V.cols());
    for (int i = 0; i < edge_diff.size(); ++i) {
        if (edge_diff[i].array().abs().sum() == 0) {
            tree.Merge(edge_values[i].x, edge_values[i].y);
        }
    }
    tree.BuildCompactParent();
    int num = tree.CompactNum();

    // Find the most descriptive vertex
    std::vector<Vector3d> v_positions(num, Vector3d(0, 0, 0));
    std::vector<int> v_count(num);
    std::vector<double> v_distance(num, 1e30);
    std::vector<int> v_index(num, -1);

    for (int i = 0; i < V.cols(); ++i) {
        v_positions[tree.Index(i)] += O.col(i);
        v_count[tree.Index(i)] += 1;
    }
    for (int i = 0; i < num; ++i) {
        if (v_count[i] > 0) v_positions[i] /= v_count[i];
    }
    
    for (int i = 0; i < V.cols(); ++i) {
        int p = tree.Index(i);
        double dis = (v_positions[p] - V.col(i)).squaredNorm();
        if (dis < v_distance[p]) {
            v_distance[p] = dis;
            v_index[p] = i;
        }
    }
    std::vector<int> unchange{};
    for (int i = 0; i < corner.size(); ++i) {
        int pppp = tree.Index(corner[i]);
        unchange.push_back(pppp);
    }

    std::set<int> compact_sharp_vertices;
    for (auto& v : sharp_vertices) {
        v_positions[tree.Index(v)] = O.col(v);
        v_index[tree.Index(v)] = v;
        V.col(v) = O.col(v);
        compact_sharp_vertices.insert(tree.Index(v));
    }
    std::set<int> compact_boundary_vertices;
    for (auto& v : Boundary_vertices) {
        v_positions[tree.Index(v)] = O.col(v);
        v_index[tree.Index(v)] = v;
        V.col(v) = O.col(v);
        compact_boundary_vertices.insert(tree.Index(v));
    }

    std::map<int, int> compact_sharp_indices;
    std::set<DEdge> compact_sharp_edges;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 1) {
            int v1 = tree.Index(F(i % 3, i / 3));
            int v2 = tree.Index(F((i + 1) % 3, i / 3));
            compact_sharp_edges.insert(DEdge(v1, v2));
        }
    }
    for (auto& v : sharp_vertices) {
        int p = tree.Index(v);
        if (compact_sharp_indices.count(p) == 0) {
            int s = compact_sharp_indices.size();
            compact_sharp_indices[p] = s;
        }
    }
    std::map<int, int> compact_boundary_indices;
    std::set<DEdge> compact_boundary_edges;
    for (int i = 0; i < Boundary_edges.size(); ++i) {
        if (Boundary_edges[i] == 1) {
            int v1 = tree.Index(F(i % 3, i / 3));
            int v2 = tree.Index(F((i + 1) % 3, i / 3));
            compact_boundary_edges.insert(DEdge(v1, v2));
        }
    }
    for (auto& v : Boundary_vertices) {
        int p = tree.Index(v);
        if (compact_boundary_indices.count(p) == 0) {
            int s = compact_boundary_indices.size();
            compact_boundary_indices[p] = s;
        }
    }
    std::map<int, std::set<int>> sharp_vertices_links;
    std::set<DEdge> sharp_dedges;
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i]) {
            int v1 = F(i % 3, i / 3);
            int v2 = F((i + 1) % 3, i / 3);
            if (sharp_vertices_links.count(v1) == 0) sharp_vertices_links[v1] = std::set<int>();
            sharp_vertices_links[v1].insert(v2);
            sharp_dedges.insert(DEdge(v1, v2));
        }
    }
    std::map<int, std::set<int>> boundary_vertices_links;
    std::set<DEdge> boundary_dedges;
    for (int i = 0; i < Boundary_edges.size(); ++i) {
        if (Boundary_edges[i]) {
            int v1 = F(i % 3, i / 3);
            int v2 = F((i + 1) % 3, i / 3);
            if (boundary_vertices_links.count(v1) == 0)
                boundary_vertices_links[v1] = std::set<int>();
            boundary_vertices_links[v1].insert(v2);
            boundary_dedges.insert(DEdge(v1, v2));
        }
    }
    std::vector<std::vector<int>> sharp_to_original_indices(compact_sharp_indices.size());
    for (auto& v : sharp_vertices_links) {
        if (v.second.size() == 2) continue;
        int p = tree.Index(v.first);
        sharp_to_original_indices[compact_sharp_indices[p]].push_back(v.first);
    }
    for (auto& v : sharp_vertices_links) {
        if (v.second.size() != 2) continue;
        int p = tree.Index(v.first);
        sharp_to_original_indices[compact_sharp_indices[p]].push_back(v.first);
    }

    for (int i = 0; i < V.cols(); ++i) {
        if (sharp_vertices.count(i)) continue;
        int p = tree.Index(i);
        if (compact_sharp_indices.count(p))
            sharp_to_original_indices[compact_sharp_indices[p]].push_back(i);
    }
    std::vector<std::vector<int>> boundary_to_original_indices(compact_boundary_indices.size());
    for (auto& v : boundary_vertices_links) {
        if (v.second.size() == 2) continue;
        int p = tree.Index(v.first);
        boundary_to_original_indices[compact_boundary_indices[p]].push_back(v.first);
    }
    for (auto& v : boundary_vertices_links) {
        if (v.second.size() != 2) continue;
        int p = tree.Index(v.first);
        boundary_to_original_indices[compact_boundary_indices[p]].push_back(v.first);
    }

    for (int i = 0; i < V.cols(); ++i) {
        if (Boundary_vertices.count(i)) continue;
        int p = tree.Index(i);
        if (compact_boundary_indices.count(p))
            boundary_to_original_indices[compact_boundary_indices[p]].push_back(i);
    }
    std::vector<std::map<int, std::pair<int, Vector3d>>> ideal_distances(tree.CompactNum());
    for (int e = 0; e < edge_diff.size(); ++e) {
        int v1 = edge_values[e].x;
        int v2 = edge_values[e].y;

        int p1 = tree.Index(v1);
        int p2 = tree.Index(v2);
        int q1 = v_index[p1];
        int q2 = v_index[p2];

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
        double scale_x = (with_scale ? 0.5 * (s_x1 + s_x2) : 1) * mRes.mScale;
        double scale_y = (with_scale ? 0.5 * (s_y1 + s_y2) : 1) * mRes.mScale;
        Vector2i diff = edge_diff[e];

        Vector3d origin1 =
            /*(sharp_constraints.count(q1)) ? sharp_constraints[q1].first : */ V.col(q1);
        Vector3d origin2 =
            /*(sharp_constraints.count(q2)) ? sharp_constraints[q2].first : */ V.col(q2);
        Vector3d C = diff[0] * scale_x * qd_x + diff[1] * scale_y * qd_y + origin1 - origin2;
        auto it = ideal_distances[p1].find(p2);
        if (it == ideal_distances[p1].end()) {
            ideal_distances[p1][p2] = std::make_pair(1, C);
        } else {
            it->second.first += 1;
            it->second.second += C;
        }
    }

    std::vector<std::unordered_map<int, double>> entries(num * 2);
    std::vector<double> b(num * 2);

    for (int m = 0; m < num; ++m) {
        int v1 = v_index[m];
        for (auto& info : ideal_distances[m]) {
            int v2 = v_index[info.first];
            Vector3d q_1 = Q.col(v1);
            Vector3d q_2 = Q.col(v2);
            if (sharp_constraints.count(v1)) {
                Vector3d d = sharp_constraints[v1].second;
                if (d != Vector3d::Zero()) q_1 = d;
            }
            if (sharp_constraints.count(v2)) {
                Vector3d d = sharp_constraints[v2].second;
                if (d != Vector3d::Zero()) q_2 = d;
            }
            if (Boundary_constraints.count(v1)) {
                Vector3d d = Boundary_constraints[v1].second;
                if (d != Vector3d::Zero()) q_1 = d;
            }
            if (Boundary_constraints.count(v2)) {
                Vector3d d = Boundary_constraints[v2].second;
                if (d != Vector3d::Zero()) q_2 = d;
            }

            Vector3d n_1 = N.col(v1);
            Vector3d n_2 = N.col(v2);
            Vector3d q_1_y = n_1.cross(q_1);
            Vector3d q_2_y = n_2.cross(q_2);
            Vector3d weights[] = {q_2, q_2_y, -q_1, -q_1_y};
            int vid[] = {info.first * 2, info.first * 2 + 1, m * 2, m * 2 + 1};
            Vector3d dis = info.second.second / info.second.first;
            double lambda = 1;
            if (sharp_vertices.count(v1) && sharp_vertices.count(v2)) lambda = 1;
            if (Boundary_vertices.count(v1) && Boundary_vertices.count(v2)) lambda = 1;
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    auto it = entries[vid[i]].find(vid[j]);
                    if (it == entries[vid[i]].end()) {
                        entries[vid[i]][vid[j]] = weights[i].dot(weights[j]) * lambda;
                    } else {
                        entries[vid[i]][vid[j]] += weights[i].dot(weights[j]) * lambda;
                    }
                }
                b[vid[i]] += weights[i].dot(dis) * lambda;
            }
        }
    }

    std::vector<int> fixed_dim(num * 2, 0);
    std::vector<int> fixed_dimB(num * 2, 0);
    std::vector<double> x(num * 2);
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < num; ++i) {
        int p = v_index[i];
        Vector3d q = Q.col(p);

        if (sharp_constraints.count(p)) {
            Vector3d dir = sharp_constraints[p].second;
            fixed_dim[i * 2 + 1] = 1;
            if (dir != Vector3d::Zero()) {
                q = dir;
            } else {
                fixed_dim[i * 2] = 1;
            }
        }
        if (Boundary_constraints.count(p)) {
            Vector3d dir = Boundary_constraints[p].second;
            fixed_dim[i * 2 + 1] = 1;
            fixed_dimB[i * 2 + 1] = 1;

            if (dir != Vector3d::Zero()) {
                q = dir;
            } else {
                fixed_dim[i * 2] = 1;
                fixed_dimB[i * 2] = 1;
            }
        }
        Vector3d n = N.col(p);
        Vector3d q_y = n.cross(q);
        x[i * 2] = (v_positions[i] - V.col(p)).dot(q);
        x[i * 2 + 1] = (v_positions[i] - V.col(p)).dot(q_y);
    }

    // fix sharp edges
    for (int i = 0; i < entries.size(); ++i) {
        if (fixed_dim[i]) {
            b[i] = x[i];
            entries[i].clear();
            entries[i][i] = 1;
        } else {
            std::unordered_map<int, double> newmap;
            for (auto& rec : entries[i]) {
                if (fixed_dim[rec.first]) {
                    b[i] -= rec.second * x[rec.first];
                } else {
                    newmap[rec.first] = rec.second;
                }
            }
            std::swap(entries[i], newmap);
        }
    }
    for (int i = 0; i < entries.size(); ++i) {
        if (entries[i].size() == 0) {
            entries[i][i] = 1;
        }
    }

    std::vector<Eigen::Triplet<double>> lhsTriplets;
    lhsTriplets.reserve(F.cols() * 6);
    Eigen::SparseMatrix<double> A(num * 2, num * 2);
    VectorXd rhs(num * 2);
    rhs.setZero();
    for (int i = 0; i < entries.size(); ++i) {
        rhs(i) = b[i];
        if (std::isnan(b[i])) {
            printf("Equation has nan!\n");
            exit(0);
        }
        for (auto& rec : entries[i]) {
            lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, rec.second));
            if (std::isnan(rec.second)) {
                printf("Equation has nan!\n");
                exit(0);
            }
        }
    }
    A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());

#ifdef LOG_OUTPUT
    int t1 = GetCurrentTime64();
#endif
    /*
        Eigen::setNbThreads(1);
        ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;
        VectorXd x0 = VectorXd::Map(x.data(), x.size());
        solver.setMaxIterations(40);

        solver.compute(A);
     */
    LinearSolver<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    VectorXd x_new = solver.solve(rhs);
#ifdef LOG_OUTPUT
    // std::cout << "[LSQ] n_iteration:" << solver.iterations() << std::endl;
    // std::cout << "[LSQ] estimated error:" << solver.error() << std::endl;
    int t2 = GetCurrentTime64();
    printf("[LSQ] Linear solver uses %lf seconds.\n", (t2 - t1) * 1e-3);
#endif

    for (int i = 0; i < x.size(); ++i) {
        if (!std::isnan(x_new[i])) {
            if (!fixed_dim[i / 2 * 2 + 1]) {
                double total = 0;
                for (auto& res : entries[i]) {
                    double t = x_new[res.first];
                    if (std::isnan(t)) t = 0;
                    total += t * res.second;
                }
            }
            x[i] = x_new[i];
        }
    }

    for (int i = 0; i < O.cols(); ++i) {
        int p = tree.Index(i);
        int c = v_index[p];
        Vector3d q = Q.col(c);
        if (fixed_dim[p * 2 + 1] == 1 && sharp_constraints.count(c) == 1) {
            Vector3d dir = sharp_constraints[c].second;
            if (dir != Vector3d::Zero()) q = dir;
        }
        if (fixed_dim[p * 2 + 1] == 1 && Boundary_constraints.count(c) == 1) {
            Vector3d dir = Boundary_constraints[c].second;
            if (dir != Vector3d::Zero()) q = dir;
        }
        
        Vector3d n = N.col(c);
        Vector3d q_y = n.cross(q);
        if (std::count(unchange.begin(), unchange.end(), p)){
            continue;
        } else {
            if (Boundary_constraints.count(c)) {
                Vector3d O_new = V.col(c) + q * x[p * 2] + q_y * x[p * 2 + 1];
                Vector3d dirdir = O_new - O.col(i);
                Vector3d dirdir1 = Boundary_constraints[c].second;
                Vector3d dirdir3 = dirdir1 * (dirdir.dot(dirdir1) / sqrt(dirdir1[0] * dirdir1[0] +
                                                                         dirdir1[1] * dirdir1[1] +
                                                                         dirdir1[2] * dirdir1[2]));
                O.col(i) += dirdir3;
            } else {
                O.col(i) = V.col(c) + q * x[p * 2] + q_y * x[p * 2 + 1];
            }
            
        }
    }


    std::vector<int> boundary_edge;
    for (int i = 0; i < Boundary_edges.size(); ++i) {
        if (Boundary_edges[i] == 1) {
            boundary_edge.push_back(F(i % 3, i / 3));
            boundary_edge.push_back(F((i + 1) % 3, i / 3));
        }
    }

    std::vector<std::vector<int>> alignP2(boundary_to_original_indices.size());
    Vector3d posision9;
    for (int i = 0; i < boundary_to_original_indices.size(); ++i) {
        posision9 = O.col(boundary_to_original_indices[i][0]);
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
    for (int i = 0; i < boundary_to_original_indices.size(); ++i) {
        posision61 = O.col(boundary_to_original_indices[i][0]);
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
    for (int i = 0; i < boundary_to_original_indices.size(); ++i) {
        if (find(unchangeB.begin(), unchangeB.end(), i) != unchangeB.end()) continue;
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
        pointv3 = O.col(boundary_to_original_indices[i][0]);
        Vector3d edge1 = pointv2 - pointv1;
        Vector3d edge2 = pointv3 - pointv1;
        Vector3d edge3 = edge1.normalized() * edge1.dot(edge2) / sqrt(edge1[0] * edge1[0] + edge1[1] * edge1[1] + edge1[2] * edge1[2]);
        for (auto& v : boundary_to_original_indices[i]) {
            O.col(v) = pointv1 + edge3;
        }
    }
 
}

void Optimizer::optimize_integer_constraints(Hierarchy& mRes, std::map<int, int>& singularities,
                                             bool use_minimum_cost_flow) {
    int edge_capacity = 2;
    bool fullFlow = false;
    std::vector<std::vector<int>>& AllowChange = mRes.mAllowChanges;
    for (int level = mRes.mToUpperEdges.size(); level >= 0; --level) {
        auto& EdgeDiff = mRes.mEdgeDiff[level];
        auto& FQ = mRes.mFQ[level];
        auto& F2E = mRes.mF2E[level];
        auto& E2F = mRes.mE2F[level];

        int iter = 0;
        while (!fullFlow) {
            std::vector<Vector4i> edge_to_constraints(E2F.size() * 2, Vector4i(-1, 0, -1, 0));
            std::vector<int> initial(F2E.size() * 2, 0);
            for (int i = 0; i < F2E.size(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    int e = F2E[i][j];
                    Vector2i index = rshift90(Vector2i(e * 2 + 1, e * 2 + 2), FQ[i][j]);
                    for (int k = 0; k < 2; ++k) {
                        int l = abs(index[k]);
                        int s = index[k] / l;
                        int ind = l - 1;
                        int equationID = i * 2 + k;
                        if (edge_to_constraints[ind][0] == -1) {
                            edge_to_constraints[ind][0] = equationID;
                            edge_to_constraints[ind][1] = s;
                        } else {
                            edge_to_constraints[ind][2] = equationID;
                            edge_to_constraints[ind][3] = s;
                        }
                        initial[equationID] += s * EdgeDiff[ind / 2][ind % 2];
                    }
                }
            }
            std::vector<std::pair<Vector2i, int>> arcs;
            std::vector<int> arc_ids;
            for (int i = 0; i < edge_to_constraints.size(); ++i) {
                if (AllowChange[level][i] == 0) continue;
                if (edge_to_constraints[i][0] == -1 || edge_to_constraints[i][2] == -1) continue;
                if (edge_to_constraints[i][1] == -edge_to_constraints[i][3]) {
                    int v1 = edge_to_constraints[i][0];
                    int v2 = edge_to_constraints[i][2];
                    if (edge_to_constraints[i][1] < 0) std::swap(v1, v2);
                    int current_v = EdgeDiff[i / 2][i % 2];
                    arcs.push_back(std::make_pair(Vector2i(v1, v2), current_v));
                    if (AllowChange[level][i] == 1)
                        arc_ids.push_back(i + 1);
                    else {
                        arc_ids.push_back(-(i + 1));
                    }
                }
            }
            int supply = 0;
            int demand = 0;
            for (int i = 0; i < initial.size(); ++i) {
                int init_val = initial[i];
                if (init_val > 0) {
                    arcs.push_back(std::make_pair(Vector2i(-1, i), initial[i]));
                    supply += init_val;
                } else if (init_val < 0) {
                    demand -= init_val;
                    arcs.push_back(std::make_pair(Vector2i(i, initial.size()), -init_val));
                }
            }

            std::unique_ptr<MaxFlowHelper> solver = nullptr;
            if (use_minimum_cost_flow && level == mRes.mToUpperEdges.size()) {
                lprintf("network simplex MCF is used\n");
                solver = std::make_unique<NetworkSimplexFlowHelper>();
            } else if (supply < 20) {
                solver = std::make_unique<ECMaxFlowHelper>();
            } else {
                solver = std::make_unique<BoykovMaxFlowHelper>();
            }

#ifdef WITH_GUROBI
            if (use_minimum_cost_flow && level == mRes.mToUpperEdges.size()) {
                solver = std::make_unique<GurobiFlowHelper>();
            }
#endif
            solver->resize(initial.size() + 2, arc_ids.size());

            std::set<int> ids;
            for (int i = 0; i < arcs.size(); ++i) {
                int v1 = arcs[i].first[0] + 1;
                int v2 = arcs[i].first[1] + 1;
                int c = arcs[i].second;
                if (v1 == 0 || v2 == initial.size() + 1) {
                    solver->addEdge(v1, v2, c, 0, -1);
                } else {
                    if (arc_ids[i] > 0)
                        solver->addEdge(v1, v2, std::max(0, c + edge_capacity),
                                        std::max(0, -c + edge_capacity), arc_ids[i] - 1);
                    else {
                        if (c > 0)
                            solver->addEdge(v1, v2, std::max(0, c - 1),
                                            std::max(0, -c + edge_capacity), -1 - arc_ids[i]);
                        else
                            solver->addEdge(v1, v2, std::max(0, c + edge_capacity),
                                            std::max(0, -c - 1), -1 - arc_ids[i]);
                    }
                }
            }
            int flow_count = solver->compute();

            solver->applyTo(EdgeDiff);

            lprintf("flow_count = %d, supply = %d\n", flow_count, supply);
            if (flow_count == supply) fullFlow = true;
            if (level != 0 || fullFlow) break;
            edge_capacity += 1;
            iter++;
            if (iter == 10) {
              /* Probably won't converge. */
              break;
            }
            lprintf("Not full flow, edge_capacity += 1\n");
        }

        if (level != 0) {
            auto& nEdgeDiff = mRes.mEdgeDiff[level - 1];
            auto& toUpper = mRes.mToUpperEdges[level - 1];
            auto& toUpperOrients = mRes.mToUpperOrients[level - 1];
            for (int i = 0; i < toUpper.size(); ++i) {
                if (toUpper[i] >= 0) {
                    int orient = (4 - toUpperOrients[i]) % 4;
                    nEdgeDiff[i] = rshift90(EdgeDiff[toUpper[i]], orient);
                }
            }
        }
    }
}

#ifdef WITH_CUDA

void Optimizer::optimize_orientations_cuda(Hierarchy& mRes) {
    int levelIterations = 6;
    for (int level = mRes.mN.size() - 1; level >= 0; --level) {
        Link* adj = mRes.cudaAdj[level];
        int* adjOffset = mRes.cudaAdjOffset[level];
        glm::dvec3* N = mRes.cudaN[level];
        glm::dvec3* Q = mRes.cudaQ[level];
        auto& phases = mRes.cudaPhases[level];
        for (int iter = 0; iter < levelIterations; ++iter) {
            for (int phase = 0; phase < phases.size(); ++phase) {
                int* p = phases[phase];
                UpdateOrientation(p, mRes.mPhases[level][phase].size(), N, Q, adj, adjOffset,
                                  mRes.mAdj[level][phase].size());
            }
        }
        if (level > 0) {
            glm::dvec3* srcField = mRes.cudaQ[level];
            glm::ivec2* toUpper = mRes.cudaToUpper[level - 1];
            glm::dvec3* destField = mRes.cudaQ[level - 1];
            glm::dvec3* N = mRes.cudaN[level - 1];
            PropagateOrientationUpper(srcField, mRes.mQ[level].cols(), toUpper, N, destField);
        }
    }

    for (int l = 0; l < mRes.mN.size() - 1; ++l) {
        glm::dvec3* N = mRes.cudaN[l];
        glm::dvec3* N_next = mRes.cudaN[l + 1];
        glm::dvec3* Q = mRes.cudaQ[l];
        glm::dvec3* Q_next = mRes.cudaQ[l + 1];
        glm::ivec2* toUpper = mRes.cudaToUpper[l];

        PropagateOrientationLower(toUpper, Q, N, Q_next, N_next, mRes.mToUpper[l].cols());
    }
}

void Optimizer::optimize_positions_cuda(Hierarchy& mRes) {
    int levelIterations = 6;
    for (int level = mRes.mAdj.size() - 1; level >= 0; --level) {
        Link* adj = mRes.cudaAdj[level];
        int* adjOffset = mRes.cudaAdjOffset[level];
        glm::dvec3* N = mRes.cudaN[level];
        glm::dvec3* Q = mRes.cudaQ[level];
        glm::dvec3* V = mRes.cudaV[level];
        glm::dvec3* O = mRes.cudaO[level];
        std::vector<int*> phases = mRes.cudaPhases[level];
        for (int iter = 0; iter < levelIterations; ++iter) {
            for (int phase = 0; phase < phases.size(); ++phase) {
                int* p = phases[phase];
                UpdatePosition(p, mRes.mPhases[level][phase].size(), N, Q, adj, adjOffset,
                               mRes.mAdj[level][phase].size(), V, O, mRes.mScale);
            }
        }
        if (level > 0) {
            glm::dvec3* srcField = mRes.cudaO[level];
            glm::ivec2* toUpper = mRes.cudaToUpper[level - 1];
            glm::dvec3* destField = mRes.cudaO[level - 1];
            glm::dvec3* N = mRes.cudaN[level - 1];
            glm::dvec3* V = mRes.cudaV[level - 1];
            PropagatePositionUpper(srcField, mRes.mO[level].cols(), toUpper, N, V, destField);
        }
    }
}

#endif

} // namespace qflow
