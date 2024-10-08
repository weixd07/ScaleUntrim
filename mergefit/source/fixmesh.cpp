#include "fixmesh.h"
#include "load.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <algorithm>




using namespace std;


void fIxmesh::Fixmesh(load& mRes, int fix_hole, int local_layer, int iteration) {
    int local;
    auto& length = mRes.min_edge_length;
    auto& average = mRes.average_edge_length;
    auto& F = mRes.Face;
    auto& V = mRes.Vertice;
    auto& n = mRes.N;
    auto& new_f = mRes.new_F;
    auto& new_v = mRes.new_V;
    auto& nonManifold = mRes.nonmanifold;
    auto& EtE = mRes.E2E;
    int local_iteration;
    if (local_layer == -1) {
        local = 3;
    }
    else {
        local = local_layer;
    }
    if (iteration == -1) {
        local_iteration = 10;
    }
    else {
        local_iteration = iteration;
    }
    vector<int> data1(V.cols(), -1);
    vector<int> boundary{};
    vector<int> Bpoint;
    for (uint32_t i = 0; i < 3 * F.cols(); ++i) {
        if (EtE[i] == -1) {
            uint32_t u0 = F(i % 3, i / 3);
            uint32_t u1 = F((i + 1) % 3, i / 3);
            boundary.push_back(u0);
            boundary.push_back(u1);
            if (data1[u0] == -1) {
                Bpoint.push_back(u0);
                data1[u0] = 1;
            }
            if (data1[u1] == -1) {
                Bpoint.push_back(u1);
                data1[u1] = 1;
            }
        }
    }
    //
    //for (int i = 0; i < boundary.size(); ++i) {
    //    if (find(Bpoint.begin(), Bpoint.end(), boundary[i]) == Bpoint.end()) {
    //        Bpoint.push_back(boundary[i]);
    //    }
    //}
    sort(Bpoint.begin(), Bpoint.end());
    //ofstream ossss("C:\\Users\\Administrator\\Desktop\\1\\plate\\test2.vtk");
    //ossss << "# vtk DataFile Version 2.0"
    //    << "\n";
    //ossss << "name, Created by Gmsh 4.11.2-git-c089f96d2 "
    //    << "\n";
    //ossss << "ASCII"
    //    << "\n";
    //ossss << "DATASET UNSTRUCTURED_GRID"
    //    << "\n";
    //ossss << "POINTS"
    //    << " " << Bpoint.size() << " "
    //    << "double"
    //    << "\n";
    //for (int i = 0; i < Bpoint.size(); ++i) {
    //    double t1, t2, t3;
    //    int u = Bpoint[i];
    //    t1 = V.col(u)[0];
    //    t2 = V.col(u)[1];
    //    t3 = V.col(u)[2];

    //    ossss << t1 << " " << t2 << " " << t3 << "\n";
    //    /* oss << 1 << " "<< i << "\n";*/
    //    /*oss << 1 << "\n";*/
    //}
    //ossss << "\n";
    //ossss << "CELLS"
    //    << " " << Bpoint.size() << " " << 2 * Bpoint.size() << "\n";
    //for (int i = 0; i < Bpoint.size(); ++i) {
    //    ossss << 1 << " " << i << "\n";
    //    /*oss << 1 << "\n";*/
    //}
    //ossss << "\n";
    //ossss << "CELL_TYPES"
    //    << " " << Bpoint.size() << "\n";
    //for (int i = 0; i < Bpoint.size(); ++i) {
    //    /*oss << t1 << " " << t2 << " " << t3 << "\n";*/
    //    /* oss << 1 << " "<< i << "\n";*/
    //    ossss << 1 << "\n";
    //}
    //ossss.close();
    std::vector<std::vector<int>> adj_v(V.cols());
    for (int f = 0; f < F.cols(); ++f) {
        for (uint32_t i = 0; i < 3; ++i) {
            uint32_t v1 = F(i, f), v2 = F((i + 1) % 3, f), v3 = F((i + 2) % 3, f);
            if (std::find(adj_v[v1].begin(), adj_v[v1].end(), v2) == adj_v[v1].end()) {
                adj_v[v1].push_back(v2);
            }
            if (std::find(adj_v[v1].begin(), adj_v[v1].end(), v3) == adj_v[v1].end()) {
                adj_v[v1].push_back(v3);
            }
        }
    }
    vector<double> min_edge;
    for (int i = 0; i < adj_v.size(); ++i) {
        double min = numeric_limits<double>::infinity();
        for (auto& vi : adj_v[i]) {
            double length = (V.col(i) - V.col(vi)).norm();
            if (length < min) {
                min = length;
            }
        }
        min_edge.push_back(min);
    }
    cout << "combine vertice" << "\n";
    vector<int> merge; // 合并较近的顶点
    vector<int> test(Bpoint.size(), 1);
    for (int i = 0; i < Bpoint.size(); ++i) {
        if (test[i] == -1)continue;
        for (int j = i + 1; j < Bpoint.size(); ++j) {
            if (test[j] == -1)continue;
            double distance = (V.col(Bpoint[i]) - V.col(Bpoint[j])).norm();
            double min_edge1, min_edge2;
            min_edge1 = min_edge[Bpoint[i]];
            min_edge2 = min_edge[Bpoint[j]];
            double minNum = (min_edge1 < min_edge2) ? min_edge1 : min_edge2;
            if (distance < minNum / 5) {
                merge.push_back(Bpoint[i]);
                merge.push_back(Bpoint[j]);
                test[j] = -1;
            }
        }
    }
    vector<int> number;
    for (int i = 0; i < merge.size(); ++i) {
        if (i % 2 == 1)continue;
        int v = merge[i];
        number.push_back(v);
        int v1 = merge[i + 1];
        int cout = std::count(number.begin(), number.end(), v);
        Eigen::Vector3d position = (V.col(v) * cout + V.col(v1)) / (cout + 1);
        Eigen::Vector3d normal = n.col(v).normalized();
        Eigen::Vector3d d = position - V.col(v);
        V.col(v) += d - normal.dot(d) * normal;
        V.col(v) = position;
    }
    std::unordered_map<int, int> merge_index_map;
    vector<int> decrease;
    /*for (int i = 0; i < merge.size(); ++i) {
        if (i % 2 == 1)continue;
        decrease.push_back(merge[i + 1]);
    }
    vector<int> relative(decrease.size());
    sort(decrease.begin(), decrease.end());
    for (int i = 0; i < decrease.size(); ++i) {
        auto it = std::find(merge.begin(), merge.end(), decrease[i]);
        int index = std::distance(merge.begin(), it);
        relative[i] = merge[index - 1];
    }*/
    
    for (int i = 0; i < merge.size(); ++i) {
        merge_index_map[merge[i]] = i;
    }
    for (int i = 1; i < merge.size(); i += 2) {
        decrease.push_back(merge[i]);
    }
    std::sort(decrease.begin(), decrease.end());
    std::vector<int> relative(decrease.size());
    for (int i = 0; i < decrease.size(); ++i) {
        relative[i] = merge[merge_index_map[decrease[i]] - 1];
    }
    vector<int> table(V.cols(), 0);
    vector<int> r_table(V.cols(), -1);
    for (int i = 0; i < decrease.size(); ++i) {
        r_table[decrease[i]] = relative[i];
        if (i == decrease.size() - 1)continue;
        int v = decrease[i];
        int v1 = decrease[i + 1];
        for (int j = v + 1; j < v1; ++j) {
            table[j] = i + 1;
        }
    }
    if (decrease.size() != 0) {
        for (int i = decrease[decrease.size() - 1] + 1; i < table.size(); ++i) {
            table[i] = decrease.size();
        }
    }
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            auto& t1 = F.col(i)[j];
            if (r_table[t1] == -1) {
                t1 -= table[t1];
            }
            else {
                int v = r_table[t1];
                t1 = v - table[v];
            }

        }
    }
    Eigen::MatrixXd new_n;
    new_v.resize(3, V.cols() - decrease.size());
    new_n.resize(3, V.cols() - decrease.size());
    int a = 0;
    for (int i = 0; i < V.cols(); ++i) {
        if (r_table[i] != -1)continue;
        new_v.col(a) = V.col(i);
        new_n.col(a) = n.col(i);
        a += 1;
    }
    /*for (int i = 0; i < decrease.size(); ++i) {
        int v = decrease[i];
        for (int j = 0; j < F.cols(); ++j) {
            auto& t1 = F.col(j)[0];
            auto& t2 = F.col(j)[1];
            auto& t3 = F.col(j)[2];
            if (t1 == v) {
                t1 = relative[i];
            }
            if (t2 == v) {
                t2 = relative[i];
            }
            if (t3 == v) {
                t3 = relative[i];
            }

        }
    }
    for (int i = 0; i < decrease.size(); ++i) {
        int v = decrease[i];
        if (i == 0) {
            for (int j = 0; j < v + 1; ++j) {
                table[j] = -1;
            }
            continue;
        }
        for (int j = 0; j < F.cols(); ++j) {
            auto& t1 = F.col(j)[0];
            auto& t2 = F.col(j)[1];
            auto& t3 = F.col(j)[2];
            if (t1 < v && table[t1] != -1) {
                t1 -= i;
            }
            if (t2 < v && table[t2] != -1) {
                t2 -= i;
            }
            if (t3 < v && table[t3] != -1) {
                t3 -= i;
            }
        }
        for (int j = 0; j < v; ++j) {
            table[j] = -1;
        }
        table[v] = -1;
        if (i == decrease.size() - 1 && v != V.cols()) {
            for (int j = 0; j < F.cols(); ++j) {
                auto& t1 = F.col(j)[0];
                auto& t2 = F.col(j)[1];
                auto& t3 = F.col(j)[2];
                if (t1 > v) {
                    t1 -= decrease.size();
                }
                if (t2 > v) {
                    t2 -= decrease.size();
                }
                if (t3 > v) {
                    t3 -= decrease.size();
                }
            }

        }
    }*/
    /*Eigen::MatrixXd new_n;
    new_v.resize(3, V.cols() - decrease.size());
    new_n.resize(3, V.cols() - decrease.size());
    int a = 0;
    for (int i = 0; i < V.cols(); ++i) {
        if (find(decrease.begin(), decrease.end(), i) != decrease.end()) continue;
        new_v.col(a) = V.col(i);
        new_n.col(a) = n.col(i);
        a += 1;
    }*/
    cout << "combine vertice and edge" << "\n";
    computeboundary(new_v, F, EtE);
    std::vector<std::vector<int>> adjv(new_v.cols());
    for (int f = 0; f < F.cols(); ++f) {
        for (uint32_t i = 0; i < 3; ++i) {
            uint32_t v1 = F(i, f), v2 = F((i + 1) % 3, f), v3 = F((i + 2) % 3, f);
            if (std::find(adjv[v1].begin(), adjv[v1].end(), v2) == adjv[v1].end()) {
                adjv[v1].push_back(v2);
            }
            if (std::find(adjv[v1].begin(), adjv[v1].end(), v3) == adjv[v1].end()) {
                adjv[v1].push_back(v3);
            }
        }
    }
    vector<double> minedge;
    for (int i = 0; i < adjv.size(); ++i) {
        double min = numeric_limits<double>::infinity();
        for (auto& vi : adjv[i]) {
            double length = (new_v.col(i) - new_v.col(vi)).norm();
            if (length < min) {
                min = length;
            }
        }
        minedge.push_back(min);
    }
    vector<int> sequence;
    vector<int> boundaryedge{};
    vector<int> data2(new_v.cols(), -1);
    vector<int> boundarypoint; 
    vector<int> small_edge;
    for (uint32_t i = 0; i < 3 * F.cols(); ++i) {
        if (EtE[i] == -1) {
            sequence.push_back(i);
            uint32_t u0 = F(i % 3, i / 3);
            uint32_t u1 = F((i + 1) % 3, i / 3);
            boundaryedge.push_back(u0);
            boundaryedge.push_back(u1);
            if (data2[u0] == -1) {
                boundarypoint.push_back(u0);
                data2[u0] = 1;
            }
            if (data2[u1] == -1) {
                boundarypoint.push_back(u1);
                data2[u1] = 1;
            }
        }
    }
    /*for (uint32_t i = 0; i < 3 * F.cols(); ++i) {
        if (EtE[i] == -1) {
            sequence.push_back(i);
            uint32_t u0 = F(i % 3, i / 3);
            uint32_t u1 = F((i + 1) % 3, i / 3);
            boundaryedge.push_back(u0);
            boundaryedge.push_back(u1);
        }
    }

    vector<int> boundarypoint;
    vector<int> small_edge;
    for (int i = 0; i < boundaryedge.size(); ++i) {
        if (find(boundarypoint.begin(), boundarypoint.end(), boundaryedge[i]) == boundarypoint.end()) {
            boundarypoint.push_back(boundaryedge[i]);
        }
    }*/
    ofstream os("C:\\Users\\Administrator\\Desktop\\1\\plate\\test1.vtk");
    os << "# vtk DataFile Version 2.0"
        << "\n";
    os << "name, Created by Gmsh 4.11.2-git-c089f96d2 "
        << "\n";
    os << "ASCII"
        << "\n";
    os << "DATASET UNSTRUCTURED_GRID"
        << "\n";
    os << "POINTS"
        << " " << boundarypoint.size() << " "
        << "double"
        << "\n";
    for (int i = 0; i < boundarypoint.size(); ++i) {
        double t1, t2, t3;
        int u = boundarypoint[i];
        t1 = new_v.col(u)[0];
        t2 = new_v.col(u)[1];
        t3 = new_v.col(u)[2];

        os << t1 << " " << t2 << " " << t3 << "\n";
        /* oss << 1 << " "<< i << "\n";*/
        /*oss << 1 << "\n";*/
    }
    os << "\n";
    os << "CELLS"
        << " " << boundarypoint.size() << " " << 2 * boundarypoint.size() << "\n";
    for (int i = 0; i < boundarypoint.size(); ++i) {
        os << 1 << " " << i << "\n";
        /*oss << 1 << "\n";*/
    }
    os << "\n";
    os << "CELL_TYPES"
        << " " << boundarypoint.size() << "\n";
    for (int i = 0; i < boundarypoint.size(); ++i) {
        /*oss << t1 << " " << t2 << " " << t3 << "\n";*/
        /* oss << 1 << " "<< i << "\n";*/
        os << 1 << "\n";
    }
    os.close();
    vector<int> need_increase;
    vector<int> relative_best_j;
    for (int i = 0; i < boundarypoint.size(); ++i) {
        int u = boundarypoint[i];
        int best_j = -1;
        double little_distance = 100000;
        Eigen::Vector3d best_projection;
        for (int j = 0; j < sequence.size(); ++j) {
            int u0 = boundaryedge[2 * j];
            int u1 = boundaryedge[2 * j + 1];
            if (u0 == u || u1 == u)continue;
            Eigen::Vector3d edge1 = new_v.col(u0) - new_v.col(u);
            Eigen::Vector3d edge2 = new_v.col(u1) - new_v.col(u);
            Eigen::Vector3d edge = new_v.col(u1) - new_v.col(u0);
            if (edge1.norm() > edge.norm() || edge2.norm() > edge.norm()) continue;
            edge1.normalize();
            edge2.normalize();
            double criterian = edge1.dot(edge2);
            if (criterian < -0.5) {
                Eigen::Vector3d edge3 = new_v.col(u1) - new_v.col(u0);
                Eigen::Vector3d edge4 = new_v.col(u) - new_v.col(u0);
                Eigen::Vector3d edge5 = edge3.normalized();
                Eigen::Vector3d projection = edge4 - edge4.dot(edge5) * edge5;
                double distance = projection.norm();
                if (distance < little_distance) {
                    best_projection = projection;
                    best_j = j;
                    little_distance = distance;
                }
            }
        }
        relative_best_j.push_back(best_j);
        if (best_j == -1)continue;
        int n = sequence[best_j];
        int v0 = F(n % 3, n / 3);
        int v1 = F((n + 1) % 3, n / 3);
        int v2 = F((n + 2) % 3, n / 3);
        Eigen::Vector3d edge1 = new_v.col(v1) - new_v.col(v0);
        Eigen::Vector3d edge2 = new_v.col(v2) - new_v.col(v0);
        Eigen::Vector3d normal = edge1.cross(edge2);
        Eigen::Vector3d edge3 = new_v.col(u) - new_v.col(v0);
        normal.normalize();
        Eigen::Vector3d projection_u = new_v.col(v0) + edge3 - edge3.dot(normal) * normal;
        Eigen::Vector3d edge4 = (new_v.col(v0) - projection_u).normalized();
        Eigen::Vector3d edge5 = (new_v.col(v1) - projection_u).normalized();
        Eigen::Vector3d edge6 = (new_v.col(v2) - projection_u).normalized();
        Eigen::Vector3d edge7 = edge4.cross(edge5);
        Eigen::Vector3d edge8 = edge5.cross(edge6);
        Eigen::Vector3d edge9 = edge6.cross(edge4);
        if (edge7.dot(edge8) >= 0 && edge7.dot(edge9) >= 0 && edge8.dot(edge9) >= 0) {
            edge1.normalize();
            new_v.col(u) = new_v.col(v0) + (projection_u - new_v.col(v0)).dot(edge1) * edge1;
        }
    }
    for (int i = 0; i < relative_best_j.size(); ++i) {
        if (relative_best_j[i] == -1)continue;
        int u = boundarypoint[i];
        int n = sequence[relative_best_j[i]];
        int v0 = F(n % 3, n / 3);
        int v1 = F((n + 1) % 3, n / 3);
        int v2 = F((n + 2) % 3, n / 3);
        Eigen::Vector3d edge1 = new_v.col(v1) - new_v.col(v0);
        Eigen::Vector3d edge2 = new_v.col(u) - new_v.col(v0);
        Eigen::Vector3d edge3 = edge1.normalized();
        double little_distance = (edge2 - edge3.dot(edge2) * edge3).norm();
        double edge1_length = edge1.norm();
        double edge = minedge[u];
        double minNum = (edge1_length < edge) ? edge1_length : edge;
        if (little_distance < minNum / 5) {
            need_increase.push_back(i);
            small_edge.push_back(relative_best_j[i]);
        }
    }
    vector<int> relative_triangle;
    vector<vector<int>> same;
    for (int i = 0; i < need_increase.size(); ++i) {
        int n = sequence[small_edge[i]];
        int s = n / 3;
        relative_triangle.push_back(s);
    }
    vector<int> same_increase;
    vector<int> a_s(relative_triangle.size(), 1);
    for (int i = 0; i < relative_triangle.size(); ++i) {
        if (a_s[i] == -1)continue;
        int n = relative_triangle[i];
        a_s[i] = -1;
        if (count(relative_triangle.begin(), relative_triangle.end(), n) == 1)continue;
        vector<int> matrix;
        matrix.push_back(i);
        same_increase.push_back(i);
        for (int j = 0; j < relative_triangle.size(); ++j) {
            if (a_s[j] == -1)continue;
            if (j == i)continue;
            if (n == relative_triangle[j]) {
                matrix.push_back(j);
                a_s[j] = -1;
            }
        }
        same.push_back(matrix);
    }
    int p = 0;
    vector<int> relative_n;
    vector<int> a_switch(need_increase.size(), 1);
    for (int i = 0; i < need_increase.size(); ++i) {
        if (a_switch[i] == -1)continue;
        if (find(same_increase.begin(), same_increase.end(), i) == same_increase.end()) {
            int n = sequence[small_edge[i]];
            int u0 = F(n % 3, n / 3);
            int u1 = F((n + 1) % 3, n / 3);
            int u2 = F((n + 2) % 3, n / 3);
            new_f.push_back(u0);
            new_f.push_back(boundarypoint[need_increase[i]]);
            new_f.push_back(u2);
            new_f.push_back(boundarypoint[need_increase[i]]);
            new_f.push_back(u1);
            new_f.push_back(u2);
            a_switch[i] = -1;
            //F(n % 3, n / 3) = -1;
            //F((n + 1) % 3, n / 3) = -1;
            //F((n + 2) % 3, n / 3) = -1;
            relative_n.push_back(n / 3);
        }
        else {
            a_switch[i] = -1;
            int n = sequence[small_edge[i]];
            relative_n.push_back(n / 3);
            int u0 = F(n % 3, n / 3);
            int u1 = F((n + 1) % 3, n / 3);
            int u2 = F((n + 2) % 3, n / 3);
            vector<int> u0_u1, u1_u2, u2_u0;
            u0_u1.push_back(boundarypoint[need_increase[i]]);
            for (int j = 1; j < same[p].size(); ++j) {
                int n1 = same[p][j];
                a_switch[n1] = -1;
                int n2 = sequence[small_edge[n1]];
                int v00 = F(n2 % 3, n2 / 3);
                int v11 = F((n2 + 1) % 3, n2 / 3);
                if (v00 == u0) {
                    u0_u1.push_back(boundarypoint[need_increase[n1]]);
                }
                else if (v00 == u1) {
                    u1_u2.push_back(boundarypoint[need_increase[n1]]);
                }
                else if (v11 == u0) {
                    u2_u0.push_back(boundarypoint[need_increase[n1]]);
                }
            }
            p += 1;
            vector<double> dis0, dis1, dis2, sort0, sort1, sort2;
            for (int j = 0; j < u0_u1.size(); ++j) {
                double dis = (new_v.col(u0) - new_v.col(u0_u1[j])).norm();
                dis0.push_back(dis);
            }
            for (int j = 0; j < u1_u2.size(); ++j) {
                double dis = (new_v.col(u2) - new_v.col(u1_u2[j])).norm();
                dis1.push_back(dis);
            }
            for (int j = 0; j < u2_u0.size(); ++j) {
                double dis = (new_v.col(u0) - new_v.col(u2_u0[j])).norm();
                dis2.push_back(dis);
            }
            sort0 = dis0;
            sort1 = dis1;
            sort2 = dis2;
            std::sort(dis0.begin(), dis0.end());
            std::sort(dis1.begin(), dis1.end());
            std::sort(dis2.begin(), dis2.end());
            vector<int> relative0, relative1, relative2;
            for (int j = 0; j < dis0.size(); ++j) {
                auto it = std::find(sort0.begin(), sort0.end(), dis0[j]);
                int index = std::distance(sort0.begin(), it);
                relative0.push_back(u0_u1[index]);
            }
            for (int j = 0; j < dis1.size(); ++j) {
                auto it = std::find(sort1.begin(), sort1.end(), dis1[j]);
                int index = std::distance(sort1.begin(), it);
                relative1.push_back(u1_u2[index]);
            }
            for (int j = 0; j < dis2.size(); ++j) {
                auto it = std::find(sort2.begin(), sort2.end(), dis2[j]);
                int index = std::distance(sort2.begin(), it);
                relative2.push_back(u2_u0[index]);
            }
            if (relative2.size() == 0 && relative1.size() == 0 && relative0.size() > 1) {
                new_f.push_back(u0);
                new_f.push_back(relative0[0]);
                new_f.push_back(u2);
                for (int j = 0; j < relative0.size(); ++j) {
                    if (j == relative0.size() - 1) {
                        new_f.push_back(relative0[j]);
                        new_f.push_back(u1);
                        new_f.push_back(u2);
                        continue;
                    }
                    new_f.push_back(relative0[j]);
                    new_f.push_back(relative0[j + 1]);
                    new_f.push_back(u2);

                }
                continue;
            }
            if (relative0.size() == 0 && relative1.size() == 0 && relative2.size() > 1) {
                new_f.push_back(relative2[0]);
                new_f.push_back(u0);
                new_f.push_back(u1);
                for (int j = 0; j < relative2.size(); ++j) {
                    if (j == relative2.size() - 1) {
                        new_f.push_back(relative2[j]);
                        new_f.push_back(u1);
                        new_f.push_back(u2);
                        continue;
                    }
                    new_f.push_back(relative2[j]);
                    new_f.push_back(u1);
                    new_f.push_back(relative2[j + 1]);
                }
                continue;
            }
            if (relative2.size() == 0 && relative0.size() == 0 && relative1.size() > 1) {
                new_f.push_back(u0);
                new_f.push_back(relative1[0]);
                new_f.push_back(u2);
                for (int j = 0; j < relative1.size(); ++j) {
                    if (j == relative1.size() - 1) {
                        new_f.push_back(relative1[j]);
                        new_f.push_back(u0);
                        new_f.push_back(u1);
                        continue;
                    }
                    new_f.push_back(relative1[j]);
                    new_f.push_back(u0);
                    new_f.push_back(relative1[j + 1]);
                }
                continue;
            }
            if (relative2.size() == 0) {
                new_f.push_back(u0);
                new_f.push_back(relative0[0]);
                new_f.push_back(u2);
            }
            else {
                new_f.push_back(u0);
                new_f.push_back(relative0[0]);
                new_f.push_back(relative2[0]);
                for (int j = 0; j < relative2.size(); ++j) {
                    if (relative2.size() == 1) {
                        new_f.push_back(relative2[0]);
                        new_f.push_back(relative0[0]);
                        new_f.push_back(u2);
                        continue;
                    }
                    else if (j == relative2.size() - 1) {
                        new_f.push_back(relative2[j]);
                        new_f.push_back(relative0[0]);
                        new_f.push_back(u2);
                        continue;
                    }
                    else {
                        new_f.push_back(relative2[j]);
                        new_f.push_back(relative0[0]);
                        new_f.push_back(relative2[j + 1]);
                    }
                }
            }
            if (relative1.size() == 0) {
                if (relative0.size() == 1) {
                    new_f.push_back(u2);
                    new_f.push_back(relative0[0]);
                    new_f.push_back(u1);
                }
                else {
                    for (int j = 0; j < relative0.size(); ++j) {
                        if (j == relative0.size() - 1) {
                            new_f.push_back(u2);
                            new_f.push_back(relative0[j]);
                            new_f.push_back(u1);
                            continue;
                        }
                        else {
                            new_f.push_back(u2);
                            new_f.push_back(relative0[j]);
                            new_f.push_back(relative0[j + 1]);
                        }
                    }
                }
            }
            else {
                if (relative0.size() == 1) {
                    new_f.push_back(u2);
                    new_f.push_back(relative0[0]);
                    new_f.push_back(relative1[0]);
                    for (int j = 0; j < relative1.size(); ++j) {
                        if (relative1.size() == 1) {
                            new_f.push_back(relative1[0]);
                            new_f.push_back(relative0[0]);
                            new_f.push_back(u1);
                            continue;
                        }
                        if (j == relative1.size() - 1) {
                            new_f.push_back(relative1[j]);
                            new_f.push_back(relative0[0]);
                            new_f.push_back(u1);
                            continue;
                        }
                        new_f.push_back(relative1[j]);
                        new_f.push_back(relative0[0]);
                        new_f.push_back(relative1[j + 1]);
                    }
                }
                else {
                    for (int j = 0; j < relative1.size(); ++j) {
                        if (j == 0) {
                            new_f.push_back(u2);
                            new_f.push_back(relative0[0]);
                            new_f.push_back(relative1[0]);
                            continue;
                        }
                        if (j == relative1.size() - 1) {
                            continue;
                        }
                        new_f.push_back(relative1[j]);
                        new_f.push_back(relative0[0]);
                        new_f.push_back(relative1[j + 1]);
                    }
                    for (int j = 0; j < relative0.size(); ++j) {
                        if (relative0.size() == 1) {
                            new_f.push_back(relative0[0]);
                            new_f.push_back(relative1[relative1.size() - 1]);
                            new_f.push_back(u1);
                            continue;
                        }
                        if (j == relative0.size() - 1) {
                            new_f.push_back(relative1[relative1.size() - 1]);
                            new_f.push_back(u1);
                            new_f.push_back(relative0[j]);
                            continue;
                        }
                        new_f.push_back(relative0[j]);
                        new_f.push_back(relative1[relative1.size() - 1]);
                        new_f.push_back(relative0[j + 1]);
                    }
                }
            }
        }
    }
    Eigen::MatrixXi FF;
    vector<int> tab(F.cols(), 0);
    for (int i = 0; i < relative_n.size(); ++i) {
        int n = relative_n[i];
        tab[n] = -1;
    }
    int b = 0;
    FF.resize(3, F.cols() - relative_n.size() + new_f.size() / 3);
    for (int i = 0; i < F.cols(); ++i) {
        if (tab[i] == -1)continue;
        FF.col(b) = F.col(i);
        b += 1;
    }
    for (int i = 0; i < new_f.size() / 3; ++i) {
        int t1, t2, t3;
        t1 = new_f[3 * i];
        t2 = new_f[3 * i + 1];
        t3 = new_f[3 * i + 2];
        FF.col(b)[0] = t1;
        FF.col(b)[1] = t2;
        FF.col(b)[2] = t3;
        b += 1;
    }
    F = FF;
    /*for (int i = 0; i < relative_n.size(); ++i) {
        int n = relative_n[i];
        F(n % 3, n / 3) = -1;
        F((n + 1) % 3, n / 3) = -1;
        F((n + 2) % 3, n / 3) = -1;
    }

    vector<int> temp;
    for (int i = 0; i < F.cols(); ++i) {
        double t1, t2, t3;
        t1 = F.col(i)[0];
        t2 = F.col(i)[1];
        t3 = F.col(i)[2];
        if (t1 == -1 || t2 == -1 || t3 == -1)continue;
        temp.push_back(t1);
        temp.push_back(t2);
        temp.push_back(t3);
    }
    for (int i = 0; i < new_f.size(); ++i) {
        temp.push_back(new_f[i]);
    }
    Eigen::MatrixXi FF;
    FF.resize(3, temp.size() / 3);
    for (int i = 0; i < temp.size() / 3; ++i) {
        FF.col(i)[0] = temp[i * 3];
        FF.col(i)[1] = temp[i * 3 + 1];
        FF.col(i)[2] = temp[i * 3 + 2];
    }
    F.resize(3, FF.cols());
    for (int i = 0; i < temp.size() / 3; ++i) {
        F.col(i)[0] = temp[i * 3];
        F.col(i)[1] = temp[i * 3 + 1];
        F.col(i)[2] = temp[i * 3 + 2];
    }
    new_f.resize(0);*/
    cout << "combine end" << "\n";
    computeboundary(new_v, F, EtE);
    /*for (int i = 0; i < 1; ++i) {
        vector<int> n;
        n.push_back(5458150);
        n.push_back(10867491);
        vector<vector<int>> adj_face(F.size());
        for (int i = 0; i < EtE.size(); ++i) {
            if (EtE[i] == -1)continue;
            adj_face[i / 3].push_back(EtE[i] / 3);
        }
        queue<int> queue;
        vector<int> open(F.size(), 1);
        for (int i = 0; i < n.size(); ++i) {
            queue.push(n[i]);
            open[n[i]] = -1;
        }
        while (!queue.empty()) {
            int n = queue.front();
            queue.pop();
            int u0, u1, temp;
            u0 = F.col(n)[1];
            u1 = F.col(n)[2];
            F.col(n)[1] = u1;
            F.col(n)[2] = u0;
            for (auto& v : adj_face[n]) {
                if (open[v] == -1)continue;
                queue.push(v);
                open[v] = -1;
            }
        }
        computeboundary(V, F, EtE);
    }*/
    vector<vector<int>> adj_face(F.size());
    for (int i = 0; i < EtE.size(); ++i) {
        if (EtE[i] == -1)continue;
        adj_face[i / 3].push_back(EtE[i] / 3);
    }
    queue<int> queue;
    queue.push(0);
    vector<int> open(F.size(), 1);
    open[0] = -1;
    bool change = false;
    /*while (!queue.empty()) {
        int f = queue.front();
        queue.pop();
        int u0, u1, u2;
        u0 = F.col(f)[0];
        u1 = F.col(f)[1];
        u2 = F.col(f)[2];
        for (auto& v : adj_face[f]) {
            if (open[v] == -1)continue;
            int v0, v1, v2;
            v0 = F.col(v)[0];
            v1 = F.col(v)[1];
            v2 = F.col(v)[2];
            if (v0 == u0 && v1 == u1) {
                change = true;
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
            }
            else if (v1 == u0 && v2 == u1) {
                change = true;
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
            }
            else if (v2 == u0 && v0 == u1) {
                change = true;
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
            }
            else if (v0 == u1 && v1 == u2) {
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
                change = true;
            }
            else if (v1 == u1 && v2 == u2) {
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
                change = true;
            }
            else if (v2 == u1 && v0 == u2) {
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
                change = true;
            }
            else if (v0 == u2 && v1 == u0) {
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
                change = true;
            }
            else if (v1 == u2 && v2 == u0) {
                F.col(v)[0] = v1;
                F.col(v)[1] = v0; 
                change = true;
            }
            else if (v2 == u2 && v0 == u0) {
                F.col(v)[0] = v1;
                F.col(v)[1] = v0;
                change = true;
            }
            queue.push(v);
            open[v] = -1;
        }
    }
    if (change) {
        cout << change << "\n";
        computeboundary(new_v, F, EtE);
    }*/
    if (fix_hole) {
        vector<int>sharp_edges;
        sharp_edges.resize(F.cols() * 3, 0);
        vector<Eigen::Vector3d> face_normals(F.cols());
        for (int i = 0; i < F.cols(); ++i) {
            Eigen::Vector3d p1 = new_v.col(F(0, i));
            Eigen::Vector3d p2 = new_v.col(F(1, i));
            Eigen::Vector3d p3 = new_v.col(F(2, i));
            face_normals[i] = (p2 - p1).cross(p3 - p1).normalized();
        }
        double cos_thres = cos(60.0 / 180.0 * 3.141592654);
        for (int i = 0; i < sharp_edges.size(); ++i) {
            int e = i;
            int re = EtE[e];
            if (re == -1)continue;
            Eigen::Vector3d n1 = face_normals[e / 3];
            Eigen::Vector3d n2 = face_normals[re / 3];
            if (n1.dot(n2) < cos_thres) {
                sharp_edges[i] = 1;
            }
        }
        vector<int>sharpvertex;
        vector<int> optimize_point{};
        vector<int> boundaryvertex{};
        vector<int> real_boundary;
        vector<int>bvertex;
        vector<int> table_b(new_v.size(), -1);
        vector<int> table_s(new_v.size(), -1);
        for (uint32_t i = 0; i < 3 * F.cols(); ++i) {
            if (EtE[i] == -1) {
                uint32_t u0 = F(i % 3, i / 3);
                uint32_t u1 = F((i + 1) % 3, i / 3);
                boundaryvertex.push_back(u0);
                boundaryvertex.push_back(u1);
                if (table_b[u0] == -1) {
                    bvertex.push_back(u0);
                    table_b[u0] = 0;
                }
                if (table_b[u1] == -1) {
                    bvertex.push_back(u1);
                    table_b[u1] = 0;
                }
            }
            if (sharp_edges[i]) {
                uint32_t u0 = F(i % 3, i / 3);
                uint32_t u1 = F((i + 1) % 3, i / 3);
                if (table_s[u0] == -1) {
                    sharpvertex.push_back(u0);
                    table_s[u0] = 0;
                }
                if (table_s[u1] == -1) {
                    sharpvertex.push_back(u1);
                    table_s[u1] = 0;
                }
            }
            
        }
        /*for (uint32_t i = 0; i < boundaryvertex.size(); ++i) {
            int u0 = boundaryvertex[i];
            if (find(bvertex.begin(), bvertex.end(), u0) == bvertex.end()) {
                bvertex.push_back(u0);
            }
        }*/
        ofstream oss("C:\\Users\\Administrator\\Desktop\\1\\plate\\test2.vtk");
        oss << "# vtk DataFile Version 2.0"
            << "\n";
        oss << "name, Created by Gmsh 4.11.2-git-c089f96d2 "
            << "\n";
        oss << "ASCII"
            << "\n";
        oss << "DATASET UNSTRUCTURED_GRID"
            << "\n";
        oss << "POINTS"
            << " " << bvertex.size() << " "
            << "double"
            << "\n";
        for (int i = 0; i < bvertex.size(); ++i) {
            double t1, t2, t3;
            int u = bvertex[i];
            t1 = new_v.col(u)[0];
            t2 = new_v.col(u)[1];
            t3 = new_v.col(u)[2];

            oss << t1 << " " << t2 << " " << t3 << "\n";
            /* oss << 1 << " "<< i << "\n";*/
            /*oss << 1 << "\n";*/
        }
        oss << "\n";
        oss << "CELLS"
            << " " << bvertex.size() << " " << 2 * bvertex.size() << "\n";
        for (int i = 0; i < bvertex.size(); ++i) {
            oss << 1 << " " << i << "\n";
            /*oss << 1 << "\n";*/
        }
        oss << "\n";
        oss << "CELL_TYPES"
            << " " << bvertex.size() << "\n";
        for (int i = 0; i < bvertex.size(); ++i) {
            /*oss << t1 << " " << t2 << " " << t3 << "\n";*/
            /* oss << 1 << " "<< i << "\n";*/
            oss << 1 << "\n";
        }
        oss.close();

        vector<int> compare(new_v.size(), 0);
        vector<int> error_vertex;
        bool error = false;
        for (int i = 0; i < boundaryvertex.size(); ++i) {
            compare[boundaryvertex[i]] += 1;
        }
        for (int i = 0; i < compare.size(); ++i) {
            if (compare[i] > 2) {
                error_vertex.push_back(i);
            }
        }
        if (error_vertex.size() > 0) {
            error = true;
        }
        for (int i = 0; i < boundaryvertex.size(); ++i) {
            if (i % 2 == 1)continue;
            if (boundaryvertex[i] == -1) continue;
            int point1, point2, point3, posi;
            std::vector<int> loop;
            point1 = boundaryvertex[i];
            if (find(error_vertex.begin(), error_vertex.end(), point1) != error_vertex.end()) continue;
            boundaryvertex[i] = -1;
            point2 = boundaryvertex[i + 1];
            boundaryvertex[i + 1] = -1;
            posi = i + 1;
            loop.push_back(point1);
            loop.push_back(point2);
            vector<int> inner_b;
            if (find(error_vertex.begin(), error_vertex.end(), point2) != error_vertex.end()) {
                inner_b.push_back(point2);
            }
            do {
                if (error) {
                    for (int j = 0; j < boundaryvertex.size(); ++j) {
                        if (j % 2 == 1)continue;
                        if (boundaryvertex[j] == -1) continue;
                        if (boundaryvertex[j] == point2 && j != posi) {
                            boundaryvertex[j] = -1;
                            point3 = boundaryvertex[j + 1];
                            boundaryvertex[j + 1] = -1;
                            posi = j + 1;
                            break;
                        }
                    }
                    loop.push_back(point2);
                    loop.push_back(point3);
                    point2 = point3;
                    if (find(error_vertex.begin(), error_vertex.end(), point2) != error_vertex.end()) {
                        inner_b.push_back(point2);
                        if (count(inner_b.begin(), inner_b.end(), point2) > 1) {
                            vector<int> loop1, loop2;
                            int h = 0;
                            for (int j = 0; j < loop.size(); ++j) {
                                if (loop[j] == point2) {
                                    h += 1;
                                }
                                if (h == 1 || h == 0) {
                                    loop1.push_back(loop[j]);
                                }
                                if (h == 2 || h == 3) {
                                    loop2.push_back(loop[j]);
                                }
                            }
                            loop = loop1;

                            for (int j = 0; j < loop2.size(); ++j) {
                                if (j % 2 == 0)continue;
                                if (find(optimize_point.begin(), optimize_point.end(), loop[j]) != optimize_point.end())continue;
                                optimize_point.push_back(loop2[j]);
                            }
                            if (loop2.size() == 6) {
                                new_f.push_back(loop2[0]);
                                new_f.push_back(loop2[4]);
                                new_f.push_back(loop2[2]);
                            }
                            else if (loop2.size() == 8) {
                                new_f.push_back(loop2[0]);
                                new_f.push_back(loop2[6]);
                                new_f.push_back(loop2[2]);
                                new_f.push_back(loop2[2]);
                                new_f.push_back(loop2[6]);
                                new_f.push_back(loop2[4]);
                            }
                            else {
                                double max = -numeric_limits<double>::infinity();
                                double sub_max = -numeric_limits<double>::infinity();
                                int best_j, sub_best_j;
                                for (int j = 0; j < loop2.size(); ++j) {
                                    if (j % 2 == 0)continue;
                                    int v = loop2[j];
                                    int v1;
                                    if (j == loop2.size() - 1) {
                                        v1 = loop2[1];
                                    }
                                    else {
                                        v1 = loop2[j + 2];
                                    }
                                    int v2 = loop2[j - 1];
                                    Eigen::Vector3d edge1 = (new_v.col(v1) - new_v.col(v)).normalized();
                                    Eigen::Vector3d edge2 = (new_v.col(v2) - new_v.col(v)).normalized();
                                    double cos = edge1.dot(edge2);
                                    if (cos > max) {
                                        best_j = j;
                                        max = cos;
                                    }
                                }
                                for (int j = 0; j < loop2.size(); ++j) {
                                    if (j % 2 == 0)continue;
                                    if (j == best_j)continue;
                                    int v = loop2[j];
                                    int v1;
                                    if (j == loop2.size() - 1) {
                                        v1 = loop2[1];
                                    }
                                    else {
                                        v1 = loop2[j + 2];
                                    }
                                    int v2 = loop2[j - 1];
                                    Eigen::Vector3d edge1 = (new_v.col(v1) - new_v.col(v)).normalized();
                                    Eigen::Vector3d edge2 = (new_v.col(v2) - new_v.col(v)).normalized();
                                    double cos = edge1.dot(edge2);
                                    if (cos > sub_max) {
                                        sub_best_j = j;
                                        sub_max = cos;
                                    }
                                }
                                if (best_j > sub_best_j) {
                                    int temp = sub_best_j;
                                    sub_best_j = best_j;
                                    best_j = temp;
                                }
                                vector<int> E1, E2;
                                for (int j = best_j; j < sub_best_j; ++j) {
                                    if (j % 2 == 0)continue;
                                    E1.push_back(loop2[j]);
                                }
                                E1.push_back(loop2[sub_best_j]);
                                for (int j = best_j; j > 0; --j) {
                                    if (j % 2 == 0)continue;
                                    E2.push_back(loop2[j]);
                                }
                                for (int j = loop2.size() - 1; j > sub_best_j - 1; --j) {
                                    if (j % 2 == 0)continue;
                                    E2.push_back(loop2[j]);
                                }
                                if (E1.size() < E2.size()) {
                                    for (int j = 0; j < E1.size(); ++j) {
                                        if (j == 0 || j == E1.size() - 1)continue;
                                        new_f.push_back(E1[j]);
                                        new_f.push_back(E1[j - 1]);
                                        new_f.push_back(E2[j]);
                                        new_f.push_back(E2[j]);
                                        new_f.push_back(E2[j + 1]);
                                        new_f.push_back(E1[j]);
                                    }
                                    for (int j = E1.size() - 1; j < E2.size(); ++j) {
                                        if (j == E2.size() - 1)continue;
                                        new_f.push_back(E2[j]);
                                        new_f.push_back(E2[j + 1]);
                                        new_f.push_back(E1[E1.size() - 2]);
                                    }

                                }
                                else if (E2.size() < E1.size()) {
                                    for (int j = 0; j < E2.size(); ++j) {
                                        if (j == 0 || j == E2.size() - 1)continue;
                                        new_f.push_back(E2[j - 1]);
                                        new_f.push_back(E2[j]);
                                        new_f.push_back(E1[j]);
                                        new_f.push_back(E1[j + 1]);
                                        new_f.push_back(E1[j]);
                                        new_f.push_back(E2[j]);
                                    }
                                    for (int j = E2.size() - 1; j < E1.size(); ++j) {
                                        if (j == E1.size() - 1)continue;
                                        new_f.push_back(E1[j + 1]);
                                        new_f.push_back(E1[j]);
                                        new_f.push_back(E2[E2.size() - 2]);
                                    }
                                }
                                else {
                                    for (int j = 0; j < E1.size(); ++j) {
                                        if (j == 0 || j == E1.size() - 1)continue;
                                        new_f.push_back(E1[j]);
                                        new_f.push_back(E1[j - 1]);
                                        new_f.push_back(E2[j]);
                                        new_f.push_back(E2[j]);
                                        new_f.push_back(E2[j + 1]);
                                        new_f.push_back(E1[j]);
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    for (int j = 0; j < boundaryvertex.size(); ++j) {
                        if (j % 2 == 1)continue;
                        if (boundaryvertex[j] == -1) continue;
                        if (boundaryvertex[j] == point2 && j != posi) {
                            boundaryvertex[j] = -1;
                            point3 = boundaryvertex[j + 1];
                            boundaryvertex[j + 1] = -1;
                            posi = j + 1;
                            break;
                        }
                    }
                    loop.push_back(point2);
                    loop.push_back(point3);
                    point2 = point3;
                }

            } while (point3 != point1);
            if (loop.size() < boundaryvertex.size() / 5) {
                if (loop.size() == 6) {
                    new_f.push_back(loop[0]);
                    new_f.push_back(loop[4]);
                    new_f.push_back(loop[2]);
                }
                else if (loop.size() == 8) {
                    new_f.push_back(loop[0]);
                    new_f.push_back(loop[6]);
                    new_f.push_back(loop[2]);
                    new_f.push_back(loop[2]);
                    new_f.push_back(loop[6]);
                    new_f.push_back(loop[4]);
                }
                else {
                    double max = -numeric_limits<double>::infinity();
                    double sub_max = -numeric_limits<double>::infinity();
                    int best_j, sub_best_j;
                    for (int j = 0; j < loop.size(); ++j) {
                        if (j % 2 == 0)continue;
                        int v = loop[j];
                        int v1;
                        if (j == loop.size() - 1) {
                            v1 = loop[1];
                        }
                        else {
                            v1 = loop[j + 2];
                        }
                        int v2 = loop[j - 1];
                        Eigen::Vector3d edge1 = (new_v.col(v1) - new_v.col(v)).normalized();
                        Eigen::Vector3d edge2 = (new_v.col(v2) - new_v.col(v)).normalized();
                        double cos = edge1.dot(edge2);
                        if (cos > max) {
                            best_j = j;
                            max = cos;
                        }
                    }
                    for (int j = 0; j < loop.size(); ++j) {
                        if (j % 2 == 0)continue;
                        if (j == best_j)continue;
                        int v = loop[j];
                        int v1;
                        if (j == loop.size() - 1) {
                            v1 = loop[1];
                        }
                        else {
                            v1 = loop[j + 2];
                        }
                        int v2 = loop[j - 1];
                        Eigen::Vector3d edge1 = (new_v.col(v1) - new_v.col(v)).normalized();
                        Eigen::Vector3d edge2 = (new_v.col(v2) - new_v.col(v)).normalized();
                        double cos = edge1.dot(edge2);
                        if (cos > sub_max) {
                            sub_best_j = j;
                            sub_max = cos;
                        }
                    }
                    if (best_j > sub_best_j) {
                        int temp = sub_best_j;
                        sub_best_j = best_j;
                        best_j = temp;
                    }
                    vector<int> E1, E2;
                    for (int j = best_j; j < sub_best_j; ++j) {
                        if (j % 2 == 0)continue;
                        E1.push_back(loop[j]);
                    }
                    E1.push_back(loop[sub_best_j]);
                    for (int j = best_j; j > 0; --j) {
                        if (j % 2 == 0)continue;
                        E2.push_back(loop[j]);
                    }
                    for (int j = loop.size() - 1; j > sub_best_j - 1; --j) {
                        if (j % 2 == 0)continue;
                        E2.push_back(loop[j]);
                    }
                    if (E1.size() < E2.size()) {
                        for (int j = 0; j < E1.size(); ++j) {
                            if (j == 0 || j == E1.size() - 1)continue;
                            new_f.push_back(E1[j]);
                            new_f.push_back(E1[j - 1]);
                            new_f.push_back(E2[j]);
                            new_f.push_back(E2[j]);
                            new_f.push_back(E2[j + 1]);
                            new_f.push_back(E1[j]);
                        }
                        for (int j = E1.size() - 1; j < E2.size(); ++j) {
                            if (j == E2.size() - 1)continue;
                            new_f.push_back(E2[j]);
                            new_f.push_back(E2[j + 1]);
                            new_f.push_back(E1[E1.size() - 2]);
                        }

                    }
                    else if (E2.size() < E1.size()) {
                        for (int j = 0; j < E2.size(); ++j) {
                            if (j == 0 || j == E2.size() - 1)continue;
                            new_f.push_back(E2[j - 1]);
                            new_f.push_back(E2[j]);
                            new_f.push_back(E1[j]);
                            new_f.push_back(E1[j + 1]);
                            new_f.push_back(E1[j]);
                            new_f.push_back(E2[j]);
                        }
                        for (int j = E2.size() - 1; j < E1.size(); ++j) {
                            if (j == E1.size() - 1)continue;
                            new_f.push_back(E1[j + 1]);
                            new_f.push_back(E1[j]);
                            new_f.push_back(E2[E2.size() - 2]);
                        }
                    }
                    else {
                        for (int j = 0; j < E1.size(); ++j) {
                            if (j == 0 || j == E1.size() - 1)continue;
                            new_f.push_back(E1[j]);
                            new_f.push_back(E1[j - 1]);
                            new_f.push_back(E2[j]);
                            new_f.push_back(E2[j]);
                            new_f.push_back(E2[j + 1]);
                            new_f.push_back(E1[j]);
                        }
                    }
                    for (int j = 0; j < loop.size(); ++j) {
                        if (j % 2 == 0)continue;
                        if (find(optimize_point.begin(), optimize_point.end(), loop[j]) != optimize_point.end())continue;
                        optimize_point.push_back(loop[j]);
                    }
                }
            }
            else {
                real_boundary = loop;
            }
        }
        cout << "fix endding" << "\n";
        vector<int> temp1;
        for (int i = 0; i < F.cols(); ++i) {
            int t1, t2, t3;
            t1 = F.col(i)[0];
            t2 = F.col(i)[1];
            t3 = F.col(i)[2];
            if (t1 == -1 || t2 == -1 || t3 == -1)continue;

            temp1.push_back(t1);
            temp1.push_back(t2);
            temp1.push_back(t3);
        }
        for (int i = 0; i < new_f.size() / 3; ++i) {
            int t1, t2, t3;
            t1 = new_f[3 * i];
            t2 = new_f[3 * i + 1];
            t3 = new_f[3 * i + 2];
            if (t1 == t2 || t2 == t3 || t3 == t1) {
                continue;
            }
            temp1.push_back(t1);
            temp1.push_back(t2);
            temp1.push_back(t3);
        }
        Eigen::MatrixXi FFF;
        FFF.resize(3, temp1.size() / 3);
        for (int i = 0; i < temp1.size() / 3; ++i) {
            FFF.col(i)[0] = temp1[i * 3];
            FFF.col(i)[1] = temp1[i * 3 + 1];
            FFF.col(i)[2] = temp1[i * 3 + 2];
        }
        F.resize(3, FFF.cols());
        for (int i = 0; i < temp1.size() / 3; ++i) {
            F.col(i)[0] = temp1[i * 3];
            F.col(i)[1] = temp1[i * 3 + 1];
            F.col(i)[2] = temp1[i * 3 + 2];
        }
        new_f.resize(0);
        std::vector<std::vector<int>> adj_compact(new_v.size());
        for (int f = 0; f < F.cols(); ++f) {
            for (uint32_t i = 0; i < 3; ++i) {
                uint32_t v1 = F(i, f), v2 = F((i + 1) % 3, f);
                adj_compact[v1].push_back(v2);
            }
        }
        vector<int> optimize;
        for (int i = 0; i < optimize_point.size(); ++i) {
            optimize.push_back(optimize_point[i]);
        }
        for (int o = 0; o < local; ++o) {
            vector<int> temp = optimize;
            for (int i = 0; i < temp.size(); ++i) {
                for (auto& v : adj_compact[optimize[i]]) {
                    if (find(optimize.begin(), optimize.end(), v) != optimize.end()) continue;
                    optimize.push_back(v);
                }
            }
        }
        for (int ite = 0; ite < local_iteration; ++ite) {
            for (int i = 0; i < optimize.size(); ++i) {
                if (find(real_boundary.begin(), real_boundary.end(), optimize[i]) != real_boundary.end()) continue;
                if (find(sharpvertex.begin(), sharpvertex.end(), optimize[i]) != sharpvertex.end()) continue;
                Eigen::Vector3d position;
                position[0] = 0;
                position[1] = 0;
                position[2] = 0;
                for (auto& v : adj_compact[optimize[i]]) {
                    position += new_v.col(v);
                }
                position /= adj_compact[optimize[i]].size();
                new_v.col(optimize[i]) = position;
            }
        }
    }
}
void fIxmesh::computeboundary(Eigen::MatrixXd ver, Eigen::MatrixXi fac, Eigen::VectorXi &e2e) {
    Eigen::VectorXi V2E;
    Eigen::VectorXi manifold;
    V2E.resize(ver.cols());
    V2E.setConstant(-1);
    uint32_t deg = fac.rows();
    vector<pair<uint32_t, uint32_t>> tmp(fac.size());
    for (int f = 0; f < fac.cols(); ++f) {
        for (unsigned int i = 0; i < deg; ++i) {
            unsigned int idx_cur = fac(i, f), idx_next = fac((i + 1) % deg, f), edge_id = deg * f + i;
            if (idx_cur >= ver.cols() || idx_next >= ver.cols())
                throw runtime_error("Mesh data contains an out-of-bounds vertex reference!");
            if (idx_cur == idx_next) continue;

            tmp[edge_id] = make_pair(idx_next, -1);
            if (V2E[idx_cur] == -1)
                V2E[idx_cur] = edge_id;
            else {
                unsigned int idx = V2E[idx_cur];
                while (tmp[idx].second != -1) {
                    idx = tmp[idx].second;
                }
                tmp[idx].second = edge_id;
            }
        }
    }
    manifold.resize(ver.cols());
    manifold.setConstant(false);
    e2e.resize(fac.cols() * deg);
    e2e.setConstant(-1);
    for (int f = 0; f < fac.cols(); ++f) {
        for (uint32_t i = 0; i < deg; ++i) {
            uint32_t idx_cur = fac(i, f), idx_next = fac((i + 1) % deg, f), edge_id_cur = deg * f + i;

            if (idx_cur == idx_next) continue;

            uint32_t it = V2E[idx_next], edge_id_opp = -1;
            while (it != -1) {
                if (tmp[it].first == idx_cur) {
                    if (edge_id_opp == -1) {
                        edge_id_opp = it;
                    }
                    else {
                        manifold[idx_cur] = true;
                        manifold[idx_next] = true;
                        edge_id_opp = -1;
                        break;
                    }
                }
                it = tmp[it].second;
            }
            if (edge_id_opp != -1 && edge_id_cur < edge_id_opp) {
                e2e[edge_id_cur] = edge_id_opp;
                e2e[edge_id_opp] = edge_id_cur;
            }
        }
    }
}
void fIxmesh::change_normal(load& mRes, vector<int> f) {
    auto& V = mRes.new_V;
    auto& F = mRes.Face;
    auto& EtE = mRes.E2E;
    vector<vector<int>> adj_face(F.size());
    for (int i = 0; i < EtE.size(); ++i) {
        if (EtE[i] == -1)continue;
        adj_face[i / 3].push_back(EtE[i] / 3);
    }
    queue<int> queue;
    vector<int> open(F.size(), 1);
    for (int i = 0; i < f.size(); ++i) {
        queue.push(f[i]);
        open[f[i] - 1] = -1;
    }
    while (!queue.empty()) {
        int n = queue.front();
        queue.pop();
        int u0, u1, temp;
        u0 = F.col(n)[1];
        u1 = F.col(n)[2];
        F.col(n)[1] = u1;
        F.col(n)[2] = u0;
        for (auto& v : adj_face[n]) {
            if (open[v] == -1)continue;
            queue.push(v);
            open[v] = -1;
        }
    }
    computeboundary(V, F, EtE);
}