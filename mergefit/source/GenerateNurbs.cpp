#include "GenerateNurbs.h"

TopoDS_Face GenerateNurbs::get()
{
	Handle(Geom_BSplineSurface) NURBSSurface
		= new Geom_BSplineSurface(
			cp.controlpoints_matrix,
			cp.weights_matrix,

			uvec.knots,
			vvec.knots,

			uvec.knotmultiplicity,
			vvec.knotmultiplicity,

			uvec.knotDegree,
			vvec.knotDegree
		);
	TopoDS_Face Nurbs_Surface = BRepBuilderAPI_MakeFace(NURBSSurface, false);
	return Nurbs_Surface;
};

