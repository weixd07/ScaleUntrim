#ifndef OPENSURF_BASIC_DATA_STRUCTURE_H
#define OPENSURF_BASIC_DATA_STRUCTURE_H

#include <array>
#include <vector>

using namespace std;

class Point
{
public:
    Point();
    Point(const array<double, 3>& pt_in);

    array<double, 3> x; // coordinates in physical space
    array<double, 2> u; // coordinates in parametric space
    bool is_bnd; // is boundary
    bool is_c0;
};

class Edge
{
public:
    Edge();
    Edge(const array<int, 2>& ed_in);

    bool operator==(const Edge& ed) const;

    array<int, 2> pid; // 1st point ID smaller than 2nd
    bool is_bnd;
};

class TriFace
{
public:
    TriFace();
    TriFace(const array<int, 3>& fc_in);

    array<int, 3> pid; // counter-clockwise
    array<int, 3> edid; // counter-clockwise
    bool is_bnd;
};

class Patch
{
public:
    Patch();

    array<int, 4> cn;
    array<int, 4> edid;
    array<int, 4> ednb;
    array<int, 4> edloc;

    array<int, 2> n_spt;
    vector<vector<int>> spt_id;
};

#endif //OPENSURF_BASIC_DATA_STRUCTURE_H
