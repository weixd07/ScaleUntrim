#include "GenerateShell.h"

void GenerateShell::Shell(Standard_CString a)
{
	if (patchnumber_is==1)
	{
		STEPControl_Writer Facewriter;
		Facewriter.Transfer(tempFace, STEPControl_AsIs);
		Facewriter.Write(a);
	}
	else
	{
		TopoDS_Shell sewedshape = TopoDS::Shell(sewing.SewedShape());
		STEPControl_Writer writer;
		writer.Transfer(sewedshape, STEPControl_AsIs);
		writer.Write(a);
	}

}



void GenerateShell::Read(const char* filename)
{
	ifstream infile;
	infile.open(filename);
	string initial_string;
	int patchnumber;
	infile >> initial_string >> patchnumber;
	patchnumber_is = patchnumber;

	int filename_number = 1;

	for (int i = 0; i < patchnumber ; i++)
	{
		string temp;
		int ith_patch;
		int Degree_u, Degree_v;
		int number_of_uknots, number_of_vknots, number_of_controlpoints;
		vector<Standard_Real>knotvector_u, knotvector_v;
		Standard_Real cp_x, cp_y, cp_z, cp_w;
		int row;
		int col;

		GenerateNurbs Face;
		
		//-----------------degree-----------------------
		infile >> temp >> ith_patch;
		infile >> temp >> Degree_u >> Degree_v;
		
		Face.uvec.knotDegree = Degree_u;
		Face.vvec.knotDegree = Degree_v;

		//--------------uknots-----------------------------
		infile >> temp >> number_of_uknots;
		Standard_Real u_knots;
		for (int i = 0; i < number_of_uknots; i++)
		{
			infile >> u_knots;
			knotvector_u.push_back(u_knots);
		}

		vector<int> vector_of_multiplicity_u;
		int num_initialize_u = 0;
		for (int i = 0; i < knotvector_u.size(); i = i + num_initialize_u)
		{
			num_initialize_u = count(knotvector_u.begin() + num_initialize_u, knotvector_u.end(), knotvector_u[i]);
			vector_of_multiplicity_u.push_back(num_initialize_u);
		}

		Face.uvec.knotmultiplicity.Resize(1, vector_of_multiplicity_u.size(), true);
		for (int i = 0; i < vector_of_multiplicity_u.size(); i++)
		{
			Face.uvec.knotmultiplicity.SetValue(i + 1, vector_of_multiplicity_u[i]);
		}
		auto it_initialize_u = unique(knotvector_u.begin(), knotvector_u.end());
		knotvector_u.erase(it_initialize_u, knotvector_u.end());
		vector<Standard_Real> vector_of_uknots_opencascade;
		vector_of_uknots_opencascade = knotvector_u;

		Face.uvec.knots.Resize(1, vector_of_uknots_opencascade.size(), true);
		for (int i = 0; i < vector_of_uknots_opencascade.size(); i++)
		{
			Face.uvec.knots.SetValue(i + 1, vector_of_uknots_opencascade[i]);
		}


		//--------------vknots-----------------------------
		infile >> temp >> number_of_vknots;
		Standard_Real v_knots;
		for (int i = 0; i < number_of_vknots; i++)
		{
			infile >> v_knots;
			knotvector_v.push_back(v_knots);
		}

		vector<int> vector_of_multiplicity_v;
		int num_initialize_v = 0;
		for (int i = 0; i < knotvector_v.size(); i = i + num_initialize_v)
		{
			num_initialize_v = count(knotvector_v.begin() + num_initialize_v, knotvector_v.end(), knotvector_v[i]);
			vector_of_multiplicity_v.push_back(num_initialize_v);
		}

		Face.vvec.knotmultiplicity.Resize(1, vector_of_multiplicity_v.size(), true);
		for (int i = 0; i < vector_of_multiplicity_v.size(); i++)
		{
			Face.vvec.knotmultiplicity.SetValue(i + 1, vector_of_multiplicity_v[i]);
		}


		auto it_initialize_v = unique(knotvector_v.begin(), knotvector_v.end());
		knotvector_v.erase(it_initialize_v, knotvector_v.end());
		vector<Standard_Real> vector_of_vknots_opencascade;
		vector_of_vknots_opencascade = knotvector_v;

		Face.vvec.knots.Resize(1, vector_of_vknots_opencascade.size(), true);
		for (int i = 0; i < vector_of_vknots_opencascade.size(); i++)
		{
			Face.vvec.knots.SetValue(i + 1, vector_of_vknots_opencascade[i]);
		}


		//--------------dimension----------------------
		row = number_of_uknots - Degree_u - 1;
		col = number_of_vknots - Degree_v - 1;

		Face.cp.controlpoints_matrix.Resize(1, row, 1, col, true);
		Face.cp.weights_matrix.Resize(1, row, 1, col, true);


		//---------------control points-----------------
		infile >> temp >> number_of_controlpoints;
		vector<Standard_Real> v_cp;
		for (int i = 0; i < number_of_controlpoints; i++)
		{
			infile >> cp_x >> cp_y >> cp_z >> cp_w;
			v_cp.push_back(cp_x);
			v_cp.push_back(cp_y);
			v_cp.push_back(cp_z);
			v_cp.push_back(cp_w);
		}
		
		Standard_Real index = 0;
		for (Standard_Integer i = 1; i <= row; i++)
		{
			Standard_Integer a = i;

			for (Standard_Integer j = 1; j <= col; j++)
			{
				Standard_Integer b = j;
				Standard_Real x = v_cp[index];
				Standard_Real y = v_cp[index + 1];
				Standard_Real z = v_cp[index + 2];
				Standard_Real w = v_cp[index + 3];
				Face.cp.controlpoints_matrix.SetValue(a, b, gp_Pnt(x, y, z));
				Face.cp.weights_matrix.SetValue(a, b, w);
				index = index + 4;
			}
		}
		
		tempFace = Face.get();

		sewing.Add(Face.get());
	}
	sewing.Perform();
}

Standard_Real GenerateShell::getArea(TopoDS_Shape myshape)
{
	Standard_Real total_area = 0;
	for (TopExp_Explorer Face_Area_explorer(myshape, TopAbs_FACE); Face_Area_explorer.More(); Face_Area_explorer.Next())
	{
		TopoDS_Face face = TopoDS::Face(Face_Area_explorer.Current());
		GProp_GProps Area_Prop;
		BRepGProp::SurfaceProperties(face, Area_Prop);
		Standard_Real area = Area_Prop.Mass();
		total_area += area;
	}
	cout << "The Area of the selected shell is: " << total_area << endl;
	return total_area;
}


void GenerateShell::MirrorTransformation(const char* a)
{
	TopTools_IndexedDataMapOfShapeListOfShape M;
	TopExp::MapShapesAndAncestors(sewing.SewedShape(), TopAbs_EDGE, TopAbs_FACE, M);
	BRep_Builder wirebuilder;
	TopoDS_Wire profile_wire;
	wirebuilder.MakeWire(profile_wire);
	gp_Pnt A;
	gp_Pnt B;
	int loop = 0;
    vector<Standard_Real> coordinates_vector;
	
	for (Standard_Integer i = 1; i <= M.Extent(); i++)
	{
		if (M.FindFromIndex(i).Extent() == 1)
		{

			TopoDS_Edge iEdge = TopoDS::Edge(M.FindKey(i));

			TopoDS_Vertex FirstV = TopExp::FirstVertex(iEdge);
			TopoDS_Vertex LastV = TopExp::LastVertex(iEdge);
			if (loop <= 1)
			{
				A = BRep_Tool::Pnt(FirstV);
				coordinates_vector.push_back(A.X());
				coordinates_vector.push_back(A.Y());
				coordinates_vector.push_back(A.Z());
				B = BRep_Tool::Pnt(LastV);
				coordinates_vector.push_back(B.X());
				coordinates_vector.push_back(B.Y());
				coordinates_vector.push_back(B.Z());
			}
			loop++;
			wirebuilder.Add(profile_wire, iEdge);
		}
	}
	vector<Standard_Real> Vx, Vy, Vz;
	for (int i = 0; i < 3; i++)
	{
		Vx.push_back(coordinates_vector[i] - coordinates_vector[i + 3]);
		Vy.push_back(coordinates_vector[i + 6] - coordinates_vector[i + 9]);
	}
	Vz.push_back(Vx[1] * Vy[2] - Vy[1] * Vx[2]);
	Vz.push_back(Vx[2] * Vy[0] - Vy[2] * Vx[0]);
	Vz.push_back(Vx[0] * Vy[1] - Vy[0] * Vx[1]);


	gp_Trsf MirrorTransformation;
	gp_Pnt origin(coordinates_vector[0], coordinates_vector[1], coordinates_vector[2]);
	gp_Dir plane_normal(Vz[0], Vz[1], Vz[2]);
	gp_Dir plane_Vx(Vx[0], Vx[1], Vx[2]);
	gp_Ax2 plane_mirror(origin, plane_normal, plane_Vx);
	MirrorTransformation.SetMirror(plane_mirror);
	BRepBuilderAPI_Transform ApplyMirror(sewing.SewedShape(), MirrorTransformation);
	TopoDS_Shell mirroredShell = TopoDS::Shell(ApplyMirror.Shape());
	
	sewing_mirror.Add(sewing.SewedShape());
	sewing_mirror.Add(mirroredShell);
	sewing_mirror.Perform();
	TopoDS_Shell originshell_plus_mirroredshell = TopoDS::Shell(sewing_mirror.SewedShape());
	STEPControl_Writer mirrorwriter;
	mirrorwriter.Transfer(originshell_plus_mirroredshell, STEPControl_AsIs);
	mirrorwriter.Write(a);
}