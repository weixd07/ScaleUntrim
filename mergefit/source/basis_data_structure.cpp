#include "basic_data_structure.h"

Point::Point() : x({0.,0.,0.}), u({0.,0.}), is_bnd(false), is_c0(false)
{}

Point::Point(const array<double, 3>& pt) : x(pt), u({0.,0.}), is_bnd(false), is_c0(false)
{}

Edge::Edge() : pid({-1,-1}), is_bnd(false)
{}

Edge::Edge(const array<int, 2>& ed_in) : pid(ed_in), is_bnd(false)
{
    if (pid[0] > pid[1])
    {
        int tmp = pid[0];
        pid[0] = pid[1];
        pid[1] = tmp;
    }
}

bool Edge::operator==(const Edge &ed) const
{
    if(pid == ed.pid || (pid[0] == ed.pid[1] && pid[1] == ed.pid[0]))
        return true;
    else
        return false;
}

TriFace::TriFace() : pid({-1,-1,-1}), edid({-1,-1,-1}), is_bnd(false)
{}

TriFace::TriFace(const array<int, 3>& fc_in) : pid(fc_in), edid({-1,-1,-1}), is_bnd(false)
{}

Patch::Patch() : cn{-1,-1,-1,-1}, edid{-1,-1,-1,-1}, ednb{-1,-1,-1,-1}, edloc{-1,-1,-1,-1}, n_spt{0,0}
{}