#pragma execution_character_set("utf-8")

#ifndef __LOAD_H
#define __LOAD_H
#include <stdio.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <list>
#include <map>
#include <set>

using namespace std;



class load {
public:
	vector<double> V1;
	vector<int> F1;
	vector<int> new_F;
	Eigen::MatrixXd new_V;
	Eigen::MatrixXd Vertice;
	Eigen::MatrixXd N;
	Eigen::MatrixXd Nf;
	Eigen::MatrixXi Face;
	Eigen::VectorXi V2E;
	Eigen::VectorXi E2E;
	Eigen::VectorXi nonmanifold;
	vector<int> sharp_edges;
	vector<int> Boundary_edges;
	double surface_area;
	double average_edge_length;
	double max_edge_length;
	double min_edge_length;
	double minimum_edge;
	Eigen::MatrixXd V_tri;
	Eigen::MatrixXi F_tri;
	Eigen::MatrixXd V_quad;
	Eigen::MatrixXi F_quad;
	vector<double> deflection;
	void Loadertri(const char* filename);
	void Loaderquad(const char* filename);
	void maxdeflection(double& absolute_error, double& relative_error);

	// function
	void Loader(string filename);
	/*void Loader1(const char* filename);*/
	void Outer(const char* filename, int need_fix);
	void initialize();
	void ComputeMeshStatus();
	void ComputeSmoothNormal();

	// inline function
	inline int dedge_prev_3(int e) { return (e % 3 == 0) ? e + 2 : e - 1; }
	inline int dedge_next_3(int e) { return (e % 3 == 2) ? e - 2 : e + 1; }
	inline double fast_acos(double x) {
		double negate = double(x < 0.0f);
		x = std::abs(x);
		double ret = -0.0187293f;
		ret *= x;
		ret = ret + 0.0742610f;
		ret *= x;
		ret = ret - 0.2121144f;
		ret *= x;
		ret = ret + 1.5707288f;
		ret = ret * std::sqrt(1.0f - x);
		ret = ret - 2.0f * negate * ret;
		return negate * (double)3.14159265358979323846 + ret;
	}

};

#endif