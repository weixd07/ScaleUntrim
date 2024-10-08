#include "load.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
//#include <cstdio>





void load::Loader(std::string filename) {
	std::ifstream is(filename);
	std::string line_str;
	std::vector<double> V{};
	std::vector<int> F{};
	size_t dotPosition = filename.find_last_of(".");
	int a = 0;
	if (dotPosition != std::string::npos) {
		std::string fileExtension = filename.substr(dotPosition + 1);


		if (fileExtension == "m") {
			double x, y, z;
			int v1, v2, v3;
			bool readingpos = false;
			bool readingtri = false;
			while (std::getline(is, line_str)) {

				if (line_str.find("msh.POS") != std::string::npos) {
					a += 1;
					if (a == 1)readingpos = true;
					continue;
				}
				else if (line_str.find("msh.TRIANGLES") != std::string::npos) {

					readingtri = true;
					readingpos = false;
					continue;
				}
				else if (line_str.find("]") != std::string::npos) {

					readingpos = false;
					readingtri = false;
				}

				if (readingpos) {

					if (sscanf_s(line_str.c_str(), "%lf %lf %lf", &x, &y, &z) == 3) {
						V.push_back(x);
						V.push_back(y);
						V.push_back(z);
					}
				}

				else if (readingtri) {

					if (sscanf_s(line_str.c_str(), "%d %d %d", &v1, &v2, &v3) == 3) {
						F.push_back(v1);
						F.push_back(v2);
						F.push_back(v3);
					}
				}
			}
		}
		else if (fileExtension == "obj") {
			bool vn = false;
			while (std::getline(is, line_str)) {
				std::istringstream line(line_str);
				std::string prefix;
				line >> prefix;
				if (prefix == "v") {
					double x, y, z;
					line >> x >> y >> z;
					V.push_back(x);
					V.push_back(y);
					V.push_back(z);
				}
				if (prefix == "vn") {
					vn = true;
				}

				if (prefix == "f") {
					if (vn == false) {
						int v1, v2, v3;
						line >> v1 >> v2 >> v3;
						F.push_back(v1);
						F.push_back(v2);
						F.push_back(v3);
					}
					else {
						std::string token;
						while (line >> token) {
							size_t pos = token.find("//");
							if (pos != std::string::npos) {
								std::string num_str = token.substr(0, pos);
								int value = std::stoi(num_str);
								F.push_back(value);
							}
						}
					}

				}

			}
		}
		else {
			std::cout << "please input .m or .obj file" << std::endl;
		}
	}
	V1.resize(V.size(), 0);
	F1.resize(F.size(), 0);
	std::vector<int> temp(V.size() / 3, -1);
	int b = 0;

	for (int i = 0; i < V.size(); ++i) {
		double t1 = V[i];
		V1[i] = t1;

	}

	for (int i = 0; i < F.size(); ++i) {
		int t = F[i] - 1;
		F1[i] = t;
	}
	/*for (int i = 0; i < F1.size(); ++i) {
		if (temp[F1[i]] == -1) {
			temp[F1[i]] = b;
			b += 1;
		}
	}
	for (int i = 0; i < F1.size(); ++i) {
		F1[i] = temp[F1[i]];
	}
	std::vector<int> temp1(V.size() / 3);
	for (int i = 0; i < temp.size(); ++i) {
		auto it = std::find(temp.begin(), temp.end(), i);
		int index = std::distance(temp.begin(), it);
		temp1[i] = index;
	}*/

	Vertice.resize(3, V1.size() / 3);
	Face.resize(3, F1.size() / 3);
	/*for (int i = 0; i < V1.size() / 3; ++i) {
		Vertice.col(i)[0] = V1[temp1[i] * 3];
		Vertice.col(i)[1] = V1[temp1[i] * 3 + 1];
		Vertice.col(i)[2] = V1[temp1[i] * 3 + 2];
	}*/
	for (int i = 0; i < V1.size() / 3; ++i) {
		Vertice.col(i)[0] = V1[i * 3];
		Vertice.col(i)[1] = V1[i * 3 + 1];
		Vertice.col(i)[2] = V1[i * 3 + 2];
	}
	for (int i = 0; i < F1.size() / 3; ++i) {
		Face.col(i)[0] = F1[i * 3];
		Face.col(i)[1] = F1[i * 3 + 1];
		Face.col(i)[2] = F1[i * 3 + 2];
	}
}

void load::Loader_errortri(string filename, int p1, int p2) {
	std::ifstream is(filename);
	std::string line_str;
	std::vector<double> V{};
	std::vector<int> F{};
	size_t dotPosition = filename.find_last_of(".");
	int a = 0;
	if (dotPosition != std::string::npos) {
		std::string fileExtension = filename.substr(dotPosition + 1);


		if (fileExtension == "m") {
			double x, y, z;
			int v1, v2, v3;
			bool readingpos = false;
			bool readingtri = false;
			while (std::getline(is, line_str)) {

				if (line_str.find("msh.POS") != std::string::npos) {
					a += 1;
					if (a == 1)readingpos = true;
					continue;
				}
				else if (line_str.find("msh.TRIANGLES") != std::string::npos) {

					readingtri = true;
					readingpos = false;
					continue;
				}
				else if (line_str.find("]") != std::string::npos) {

					readingpos = false;
					readingtri = false;
				}

				if (readingpos) {

					if (sscanf_s(line_str.c_str(), "%lf %lf %lf", &x, &y, &z) == 3) {
						V.push_back(x);
						V.push_back(y);
						V.push_back(z);
					}
				}

				else if (readingtri) {

					if (sscanf_s(line_str.c_str(), "%d %d %d", &v1, &v2, &v3) == 3) {
						F.push_back(v1);
						F.push_back(v2);
						F.push_back(v3);
					}
				}
			}
		}
		else if (fileExtension == "obj") {
			bool vn = false;
			while (std::getline(is, line_str)) {
				std::istringstream line(line_str);
				std::string prefix;
				line >> prefix;
				if (prefix == "v") {
					double x, y, z;
					line >> x >> y >> z;
					V.push_back(x);
					V.push_back(y);
					V.push_back(z);
				}
				if (prefix == "vn") {
					vn = true;
				}

				if (prefix == "f") {
					if (vn == false) {
						int v1, v2, v3;
						line >> v1 >> v2 >> v3;
						F.push_back(v1);
						F.push_back(v2);
						F.push_back(v3);
					}
					else {
						std::string token;
						while (line >> token) {
							size_t pos = token.find("//");
							if (pos != std::string::npos) {
								std::string num_str = token.substr(0, pos);
								int value = std::stoi(num_str);
								F.push_back(value);
							}
						}
					}

				}

			}
		}
		else {
			std::cout << "please input .m or .obj file" << std::endl;
		}
	}
	std::vector<int> temp(V.size() / 3, -1);
	int b = 0;

	for (int i = 0; i < V.size(); ++i) {
		if (i == (p1 - 1) * 3 || i == (p1 - 1) * 3 + 1 || i == (p1 - 1) * 3 + 2)continue;
		if (i == (p2 - 1) * 3 || i == (p2 - 1) * 3 + 1 || i == (p2 - 1) * 3 + 2)continue;
		double t1 = V[i];
		V1.push_back(t1);
	}

	for (int i = 0; i < F.size(); ++i) {
		int t = F[i] - 1;
		if (t < p1 - 1) {
			F1.push_back(t);
		}
		else if (t > p1 - 1 && t < p2 - 1) {
			F1.push_back(t - 1);
		}
		else if (t > p2 - 1) {
			F1.push_back(t - 2);
		}
		
	}
	/*for (int i = 0; i < F1.size(); ++i) {
		if (temp[F1[i]] == -1) {
			temp[F1[i]] = b;
			b += 1;
		}
	}
	for (int i = 0; i < F1.size(); ++i) {
		F1[i] = temp[F1[i]];
	}
	std::vector<int> temp1(V.size() / 3);
	for (int i = 0; i < temp.size(); ++i) {
		auto it = std::find(temp.begin(), temp.end(), i);
		int index = std::distance(temp.begin(), it);
		temp1[i] = index;
	}*/

	Vertice.resize(3, V1.size() / 3);
	Face.resize(3, F1.size() / 3);
	/*for (int i = 0; i < V1.size() / 3; ++i) {
		Vertice.col(i)[0] = V1[temp1[i] * 3];
		Vertice.col(i)[1] = V1[temp1[i] * 3 + 1];
		Vertice.col(i)[2] = V1[temp1[i] * 3 + 2];
	}*/
	for (int i = 0; i < V1.size() / 3; ++i) {
		Vertice.col(i)[0] = V1[i * 3];
		Vertice.col(i)[1] = V1[i * 3 + 1];
		Vertice.col(i)[2] = V1[i * 3 + 2];
	}
	for (int i = 0; i < F1.size() / 3; ++i) {
		Face.col(i)[0] = F1[i * 3];
		Face.col(i)[1] = F1[i * 3 + 1];
		Face.col(i)[2] = F1[i * 3 + 2];
	}
}
void load::converse(std::string filename, std::string filename1) {
	std::ifstream is(filename);
	std::string line_str;
	std::vector<double> V{};
	std::vector<int> F{};
	size_t dotPosition = filename.find_last_of(".");
	int a = 0;
	if (dotPosition != std::string::npos) {
		std::string fileExtension = filename.substr(dotPosition + 1);


		if (fileExtension == "m") {
			double x, y, z;
			int v1, v2, v3;
			bool readingpos = false;
			bool readingtri = false;
			while (std::getline(is, line_str)) {

				if (line_str.find("msh.POS") != std::string::npos) {
					a += 1;
					if (a == 1)readingpos = true;
					continue;
				}
				else if (line_str.find("msh.TRIANGLES") != std::string::npos) {

					readingtri = true;
					readingpos = false;
					continue;
				}
				else if (line_str.find("]") != std::string::npos) {

					readingpos = false;
					readingtri = false;
				}

				if (readingpos) {

					if (sscanf_s(line_str.c_str(), "%lf %lf %lf", &x, &y, &z) == 3) {
						V.push_back(x);
						V.push_back(y);
						V.push_back(z);
					}
				}

				else if (readingtri) {

					if (sscanf_s(line_str.c_str(), "%d %d %d", &v1, &v2, &v3) == 3) {
						F.push_back(v1);
						F.push_back(v2);
						F.push_back(v3);
					}
				}
			}
		}
		else if (fileExtension == "obj") {
			bool vn = false;
			while (std::getline(is, line_str)) {
				std::istringstream line(line_str);
				std::string prefix;
				line >> prefix;
				if (prefix == "v") {
					double x, y, z;
					line >> x >> y >> z;
					V.push_back(x);
					V.push_back(y);
					V.push_back(z);
				}
				if (prefix == "vn") {
					vn = true;
				}

				if (prefix == "f") {
					if (vn == false) {
						int v1, v2, v3;
						line >> v1 >> v2 >> v3;
						F.push_back(v1);
						F.push_back(v2);
						F.push_back(v3);
					}
					else {
						std::string token;
						while (line >> token) {
							size_t pos = token.find("//");
							if (pos != std::string::npos) {
								std::string num_str = token.substr(0, pos);
								int value = std::stoi(num_str);
								F.push_back(value);
							}
						}
					}

				}

			}
		}
		else {
			std::cout << "please input .m or .obj file" << std::endl;
		}
	}
	std::ifstream osos(filename1);
	std::string line_str1;
	std::vector<int> need_converse;
	while (std::getline(osos, line_str1)) {
		std::istringstream iss(line_str1);
		int num;
		while (iss >> num) { 
			if (find(need_converse.begin(), need_converse.end(), num) == need_converse.end()) {
				need_converse.push_back(num);
			}
		}
	}
	V1.resize(V.size(), 0);
	F1.resize(F.size(), 0);
	std::vector<int> temp(V.size() / 3, -1);
	int b = 0;

	for (int i = 0; i < V.size(); ++i) {
		double t1 = V[i];
		V1[i] = t1;

	}

	for (int i = 0; i < F.size(); ++i) {
		int t = F[i] - 1;
		F1[i] = t;
	}
	for (int i = 0; i < F1.size(); ++i) {
		if (temp[F1[i]] == -1) {
			temp[F1[i]] = b;
			b += 1;
		}
	}
	for (int i = 0; i < F1.size(); ++i) {
		F1[i] = temp[F1[i]];
	}
	std::vector<int> temp1(V.size() / 3);
	for (int i = 0; i < temp.size(); ++i) {
		auto it = std::find(temp.begin(), temp.end(), i);
		int index = std::distance(temp.begin(), it);
		temp1[i] = index;
	}
	for (int i = 0; i < need_converse.size(); ++i) {
		int v = need_converse[i];
		int t1 = F1[3 * v + 1];
		int t2 = F1[3 * v + 2];
		F1[3 * v + 1] = t2;
		F1[3 * v + 2] = t1;
	}
	Vertice.resize(3, V1.size() / 3);
	Face.resize(3, F1.size() / 3);
	for (int i = 0; i < V1.size() / 3; ++i) {
		Vertice.col(i)[0] = V1[i * 3];
		Vertice.col(i)[1] = V1[i * 3 + 1];
		Vertice.col(i)[2] = V1[i * 3 + 2];
	}
	for (int i = 0; i < F1.size() / 3; ++i) {
		Face.col(i)[0] = F1[i * 3];
		Face.col(i)[1] = F1[i * 3 + 1];
		Face.col(i)[2] = F1[i * 3 + 2];
	}
	for (int i = 0; i < V1.size() / 3; ++i) {
		Vertice.col(i)[0] = V1[temp1[i] * 3];
		Vertice.col(i)[1] = V1[temp1[i] * 3 + 1];
		Vertice.col(i)[2] = V1[temp1[i] * 3 + 2];
	}
	for (int i = 0; i < F1.size() / 3; ++i) {
		Face.col(i)[0] = F1[i * 3];
		Face.col(i)[1] = F1[i * 3 + 1];
		Face.col(i)[2] = F1[i * 3 + 2];
	}
}

void load::initialize() {
	ComputeMeshStatus();
	ComputeE2E();
	ComputeSmoothNormal();
}

void load::ComputeMeshStatus() {
	surface_area = 0;
	average_edge_length = 0;
	max_edge_length = 0;
	min_edge_length = std::numeric_limits<double>::infinity();
	for (int f = 0; f < Face.cols(); ++f) {
		Eigen::Vector3d v[3] = { Vertice.col(Face(0, f)), Vertice.col(Face(1, f)), Vertice.col(Face(2, f)) };
		double area = 0.5f * (v[1] - v[0]).cross(v[2] - v[0]).norm();
		surface_area += area;
		for (int i = 0; i < 3; ++i) {
			double len = (v[(i + 1) % 3] - v[i]).norm();
			average_edge_length += len;
			if (len > max_edge_length) max_edge_length = len;
			if (len < min_edge_length)min_edge_length = len;
		}
	}
	average_edge_length /= (Face.cols() * 3);
}
void load::ComputeE2E() {
	V2E.resize(Vertice.cols());
	V2E.setConstant(-1);
	uint32_t deg = Face.rows();
	vector<pair<uint32_t, uint32_t>> tmp(Face.size());
	for (int f = 0; f < Face.cols(); ++f) {
		for (unsigned int i = 0; i < deg; ++i) {
			unsigned int idx_cur = Face(i, f), idx_next = Face((i + 1) % deg, f), edge_id = deg * f + i;
			if (idx_cur >= Vertice.cols() || idx_next >= Vertice.cols())
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
	nonmanifold.resize(Vertice.cols());
	nonmanifold.setConstant(false);
	E2E.resize(Face.cols() * deg);
	E2E.setConstant(-1);
	for (int f = 0; f < Face.cols(); ++f) {
		for (uint32_t i = 0; i < deg; ++i) {
			uint32_t idx_cur = Face(i, f), idx_next = Face((i + 1) % deg, f), edge_id_cur = deg * f + i;

			if (idx_cur == idx_next) continue;

			uint32_t it = V2E[idx_next], edge_id_opp = -1;
			while (it != -1) {
				if (tmp[it].first == idx_cur) {
					if (edge_id_opp == -1) {
						edge_id_opp = it;
					}
					else {
						nonmanifold[idx_cur] = true;
						nonmanifold[idx_next] = true;
						edge_id_opp = -1;
						break;
					}
				}
				it = tmp[it].second;
			}
			if (edge_id_opp != -1 && edge_id_cur < edge_id_opp) {
				E2E[edge_id_cur] = edge_id_opp;
				E2E[edge_id_opp] = edge_id_cur;
			}	
		}
	}
	std::ofstream osss("C:\\Users\\Administrator\\Desktop\\compare\\E2E2.txt");
	for (int i = 0; i < E2E.size(); ++i) {
		osss << E2E[i] << "\n";
	}
	osss.close();
}
void load::ComputeSmoothNormal() {


	sharp_edges.resize(Face.cols() * 3, 0);
	Boundary_edges.resize(Face.cols() * 3, 0);

	for (int i = 0; i < Boundary_edges.size(); ++i) {
		int re = E2E[i];
		if (re == -1) {
			sharp_edges[i] = 1;
			Boundary_edges[i] = 1;
		}
	}

	/* Compute face normals */
	Nf.resize(3, Face.cols());

	for (int f = 0; f < Face.cols(); ++f) {
		Eigen::Vector3d v0 = Vertice.col(Face(0, f)), v1 = Vertice.col(Face(1, f)), v2 = Vertice.col(Face(2, f)),
			n = (v1 - v0).cross(v2 - v0);
		double norm = n.norm();
		n /= norm;
		Nf.col(f) = n;
	}

	N.resize(3, Vertice.cols());
	for (int i = 0; i < V2E.rows(); ++i) {
		int edge = V2E[i];


		int stop = edge;
		do {
			if (sharp_edges[edge])
				break;
			edge = E2E[edge];
			if (edge != -1)
				edge = dedge_next_3(edge);
		} while (edge != stop && edge != -1);
		if (edge == -1)
			edge = stop;
		else
			stop = edge;
		Eigen::Vector3d normal = Eigen::Vector3d::Zero();
		do {
			int idx = edge % 3;

			Eigen::Vector3d d0 = Vertice.col(Face((idx + 1) % 3, edge / 3)) - Vertice.col(i);
			Eigen::Vector3d d1 = Vertice.col(Face((idx + 2) % 3, edge / 3)) - Vertice.col(i);
			double angle = fast_acos(d0.dot(d1) / std::sqrt(d0.squaredNorm() * d1.squaredNorm()));

			/* "Computing Vertex Normals from Polygonal Facets"
			 by Grit Thuermer and Charles A. Wuethrich, JGT 1998, Vol 3 */
			if (std::isfinite(angle)) normal += Nf.col(edge / 3) * angle;

			int opp = E2E[edge];
			if (opp == -1) break;

			edge = dedge_next_3(opp);
			if (sharp_edges[edge])
				break;
		} while (edge != stop);
		double norm = normal.norm();
		N.col(i) = Eigen::Vector3d(normal / norm);
	}
}


void load::Outer(const char* filename, int need_fix) {
	std::ofstream os(filename);
	if (need_fix) {
		Vertice = new_V;
	}
	for (int i = 0; i < Vertice.cols(); ++i) {
		double t1 = Vertice.col(i)[0];
		double t2 = Vertice.col(i)[1];
		double t3 = Vertice.col(i)[2];

		os << "v " << t1 << " " << t2 << " " << t3 << "\n";
	}
	/*for (int i = 0; i < new_V.cols(); ++i) {
		double t1 = new_V.col(i)[0];
		double t2 = new_V.col(i)[1];
		double t3 = new_V.col(i)[2];

		os << "v " << t1 << " " << t2 << " " << t3 << "\n";
	}*/
	for (int i = 0; i < Face.cols(); ++i) {

		int r1 = Face.col(i)[0] + 1;
		int r2 = Face.col(i)[1] + 1;
		int r3 = Face.col(i)[2] + 1;
		if (r1 == 0 || r2 == 0 || r3 == 0)continue;
		os << "f " << r1 << " " << r2 << " " << r3 << "\n";

	}
	for (int i = 0; i < new_F.size() / 3; ++i) {
		int r1 = new_F[3 * i] + 1;
		int r2 = new_F[3 * i + 1] + 1;
		int r3 = new_F[3 * i + 2] + 1;
		if (r1 == 0 || r2 == 0 || r3 == 0)continue;
		os << "f " << r1 << " " << r2 << " " << r3 << "\n";
	}
	os.close();
	/*os << "# vtk DataFile Version 2.0"
		<< "\n";
	os << "name, Created by Gmsh 4.11.2-git-c089f96d2 "
		<< "\n";
	os << "ASCII"
		<< "\n";
	os << "DATASET UNSTRUCTURED_GRID"
		<< "\n";
	os << "POINTS"
		<< " " << Vertice.cols() << " "
		<< "double"
		<< "\n";

	for (int i = 0; i < Vertice.cols(); ++i) {
		double t1 = Vertice.col(i)[0];
	    double t2 = Vertice.col(i)[1];
	    double t3 = Vertice.col(i)[2];
		os << t1 << " " << t2 << " " << t3 << "\n";
	}
	os << "\n";
	os << "CELLS"
		<< " " << Face.cols() + new_F.size() / 3 << " " << 4 * (Face.cols() + new_F.size() / 3) << "\n";
	for (int i = 0; i < Face.cols(); ++i) {
		int r1 = Face.col(i)[0];
	    int r2 = Face.col(i)[1];
	    int r3 = Face.col(i)[2];
		if (r1 == -1 || r2 == -1 || r3 == -1)continue;
		os << 3 << " " << r1 << " " << r2 << " " << r3 <<"\n";
	}
	for (int i = 0; i < new_F.size() / 3; ++i) {
		int r1 = new_F[3 * i];
		int r2 = new_F[3 * i + 1];
		int r3 = new_F[3 * i + 2];
		if (r1 == -1 || r2 == -1 || r3 == -1)continue;
		os << 3 << " " << r1 << " " << r2 << " " << r3 << "\n";
	}
	os.close();
	os << "\n";
	os << "CELL_TYPES"
		<< " " << Face.cols() + new_F.size() / 3 << "\n";
	for (int i = 0; i < Face.cols() + new_F.size() / 3; ++i) {
		os << 5 << "\n";
	}
	os.close();*/
}

void load::Loadertri(const char* filename) {
	std::ifstream is(filename);
	std::string line_str;
	std::vector<double> V{};
	std::vector<int> F{};
	int a = 0;
	while (std::getline(is, line_str)) {
		std::istringstream line(line_str);
		bool vn = false;
		std::string prefix;
		line >> prefix;
		if (prefix == "v") {
			double x, y, z;
			line >> x >> y >> z;
			V.push_back(x);
			V.push_back(y);
			V.push_back(z);
		}
		if (prefix == "vn") {
			vn = true;
		}

		if (prefix == "f") {
			if (vn == false) {
				int v1, v2, v3;
				line >> v1 >> v2 >> v3;
				F.push_back(v1);
				F.push_back(v2);
				F.push_back(v3);
			}
			else {
				std::string token;
				while (line >> token) {
					size_t pos = token.find("//");
					if (pos != std::string::npos) {
						std::string num_str = token.substr(0, pos);
						int value = std::stoi(num_str);
						F.push_back(value);
					}
				}
			}

		}

	}
	V_tri.resize(3, V.size() / 3);
	F_tri.resize(3, F.size() / 3);
	for (int i = 0; i < V.size() / 3; ++i) {
		V_tri.col(i)[0] = V[3 * i];
		V_tri.col(i)[1] = V[3 * i + 1];
		V_tri.col(i)[2] = V[3 * i + 2];

	}
	for (int i = 0; i < F.size() / 3; ++i) {
		F_tri.col(i)[0] = F[3 * i] - 1;
		F_tri.col(i)[1] = F[3 * i + 1] - 1;
		F_tri.col(i)[2] = F[3 * i + 2] - 1;
	}

}



void load::Loaderquad(const char* filename) {
	std::ifstream is(filename);
	std::string line_str;
	std::vector<double> V{};
	std::vector<int> F{};
	int a = 0;
	while (getline(is, line_str)) {
		if (line_str.empty()) continue;
		if (line_str.find("POINTS") != std::string::npos) {
			a += 1;
			continue;
		}
		if (line_str.find("CELLS") != std::string::npos) {
			a += 1;
			continue;
		}
		if (line_str.find("CELL_TYPES") != std::string::npos) break;
		std::istringstream line(line_str);
		if (a == 1) {
			double x, y, z;
			line >> x >> y >> z;
			V.push_back(x);
			V.push_back(y);
			V.push_back(z);
		}
		if (a == 2) {
			int v1, v2, v3, v4, v5;
			line >> v1 >> v2 >> v3 >> v4 >> v5;
			F.push_back(v2);
			F.push_back(v3);
			F.push_back(v4);
			F.push_back(v5);
		}
	}
	V_quad.resize(3, V.size() / 3);
	F_quad.resize(4, F.size() / 4);
	for (int i = 0; i < V.size() / 3; ++i) {
		V_quad.col(i)[0] = V[3 * i];
		V_quad.col(i)[1] = V[3 * i + 1];
		V_quad.col(i)[2] = V[3 * i + 2];

	}
	for (int i = 0; i < F.size() / 4; ++i) {
		F_quad.col(i)[0] = F[4 * i];
		F_quad.col(i)[1] = F[4 * i + 1];
		F_quad.col(i)[2] = F[4 * i + 2];
		F_quad.col(i)[3] = F[4 * i + 3];
	}
}


void load::maxdeflection(double& absolute_error, double& relative_error) {
	minimum_edge = std::numeric_limits<double>::infinity();
	for (int f = 0; f < F_tri.cols(); ++f) {
		Eigen::Vector3d v[3] = { V_tri.col(F_tri(0, f)), V_tri.col(F_tri(1, f)), V_tri.col(F_tri(2, f)) };
		for (int i = 0; i < 3; ++i) {
			double len = (v[(i + 1) % 3] - v[i]).norm();
			average_edge_length += len;
			if (len < minimum_edge)minimum_edge = len;
		}
	}
	std::vector<int> boundarypoint;
	deflection.resize(V_quad.cols());
	std::vector<int>E;
	std::vector<int> relative_face(V_quad.cols());
	E.resize(4 * F_quad.cols(), -1);
	for (int f = 0; f < F_quad.cols(); ++f) {
		for (uint32_t i = 0; i < 4; ++i) {
			uint32_t idx_cur = F_quad(i, f), idx_next = F_quad((i + 1) % 4, f), edge_id_cur = 4 * f + i;
			if (E[edge_id_cur] != -1)continue;
			for (int k = 0; k < F_quad.cols(); ++k) {
				if (k == f)continue;
				for (uint32_t j = 0; j < 4; ++j) {
					uint32_t v1 = F_quad(j, k), v2 = F_quad((j + 1) % 4, k), circle = 4 * k + j;
					if (E[edge_id_cur] != -1)continue;
					if (v2 == idx_cur && v1 == idx_next) {
						E[edge_id_cur] = circle;
						E[circle] = edge_id_cur;
						break;
					}
				}
			}
		}
	}

	std::vector<int> boundaryvertex{};
	for (uint32_t i = 0; i < 4 * F_quad.cols(); ++i) {
		if (E[i] == -1) {
			uint32_t u0 = F_quad(i % 4, i / 4);
			uint32_t u1 = F_quad((i + 1) % 4, i / 4);
			boundaryvertex.push_back(u0);
			boundaryvertex.push_back(u1);
		}
	}
	for (int i = 0; i < V_quad.cols(); ++i) {
		if (std::find(boundaryvertex.begin(), boundaryvertex.end(), i) != boundaryvertex.end())continue;
		double little_distance = std::numeric_limits<double>::infinity();
		int best_j = -1;
		for (int j = 0; j < F_tri.cols(); ++j) {
			int v0, v1, v2;
			v0 = F_tri.col(j)[0];
			v1 = F_tri.col(j)[1];
			v2 = F_tri.col(j)[2];
			Eigen::Vector3d edge1 = V_tri.col(v1) - V_tri.col(v0);
			Eigen::Vector3d edge2 = V_tri.col(v2) - V_tri.col(v0);
			Eigen::Vector3d normal = edge1.cross(edge2);
			Eigen::Vector3d edge3 = V_quad.col(i) - V_tri.col(v0);
			normal.normalize();
			Eigen::Vector3d projection_i = V_tri.col(v0) + edge3 - edge3.dot(normal) * normal;
			Eigen::Vector3d edge4 = (V_tri.col(v0) - projection_i).normalized();
			Eigen::Vector3d edge5 = (V_tri.col(v1) - projection_i).normalized();
			Eigen::Vector3d edge6 = (V_tri.col(v2) - projection_i).normalized();
			Eigen::Vector3d edge7 = edge4.cross(edge5);
			Eigen::Vector3d edge8 = edge5.cross(edge6);
			Eigen::Vector3d edge9 = edge6.cross(edge4);
			if (edge7.dot(edge8) >= 0 && edge7.dot(edge9) >= 0 && edge8.dot(edge9) >= 0) {
				double dis = (V_quad.col(i) - projection_i).norm();
				if (dis < little_distance) {
					little_distance = dis;
					best_j = j;
				}
			}
		}
		if (little_distance == std::numeric_limits<double>::infinity()) {
			std::cout << "Vertice" << " " << i << "is wrong" << "\n";
		}
		else {
			deflection[i] = little_distance;
			relative_face[i] = best_j;
		}
	}
	double a = -std::numeric_limits<double>::infinity();
	for (int i = 0; i < deflection.size(); ++i) {
		if (deflection[i] == std::numeric_limits<double>::infinity())continue;
		if (a < deflection[i]) {
			a = deflection[i];
		}
	}

	double max_x = -std::numeric_limits<double>::infinity();
	double max_y = -std::numeric_limits<double>::infinity();
	double max_z = -std::numeric_limits<double>::infinity();
	double min_x = std::numeric_limits<double>::infinity();
	double min_y = std::numeric_limits<double>::infinity();
	double min_z = std::numeric_limits<double>::infinity();
	for (int i = 0; i < V_tri.cols(); ++i) {
		Eigen::Vector3d v = V_tri.col(i);
		if (v[0] > max_x) {
			max_x = v[0];
		}
		if (v[1] > max_y) {
			max_y = v[0];
		}
		if (v[2] > max_z) {
			max_z = v[0];
		}
		if (v[0] < min_x) {
			min_x = v[0];
		}
		if (v[1] < min_y) {
			min_y = v[1];
		}
		if (v[2] < min_z) {
			min_z = v[2];
		}
	}
	Eigen::Vector3d max_v, min_v;
	max_v[0] = max_x;
	max_v[1] = max_y;
	max_v[2] = max_z;
	min_v[0] = min_x;
	min_v[1] = min_y;
	min_v[2] = min_z;
	absolute_error = a;
	relative_error = a / (max_v - min_v).norm();
}