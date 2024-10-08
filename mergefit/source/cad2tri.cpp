#include <BRepTools.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <StlAPI_Writer.hxx>
#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Poly_Triangulation.hxx>
#include <BRep_Tool.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI.hxx>
#include <Bnd_Box.hxx>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <IGESControl_Reader.hxx>

double compute_area(const char* in_file)
{
	double total_area = 0.0;
	// Create a STEPControl_Reader object

	std::string a = in_file;
	int num = a.find_last_of('.');
	std::string b = a.substr(num + 1);
	if (b == "step" || b == "stp")
	{
		STEPControl_Reader reader;
		// Read the file and check the status
		IFSelect_ReturnStatus status = reader.ReadFile(in_file);
		// Transfer the shape from the file to the process
		reader.TransferRoots();
		// Get the shape from the process
		TopoDS_Shape  shape = reader.OneShape();
		// Iterate over all faces in the shape
		for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next())
		{
			const TopoDS_Face& face = TopoDS::Face(explorer.Current());

			// Create a GProp_GProps object and pass it the surface
			GProp_GProps props;
			BRepGProp::SurfaceProperties(face, props);

			// Get the surface area of the face
			double area = props.Mass();

			// Add the area of this face to the total area
			total_area += area;
		}

		return total_area;
	}
	else if (b == "igs" || b == "iges")
	{
		IGESControl_Reader reader;
		// Read the file and check the status
		IFSelect_ReturnStatus status = reader.ReadFile(in_file);
		// Transfer the shape from the file to the process
		reader.TransferRoots();
		// Get the shape from the process
		TopoDS_Shape  shape = reader.OneShape();
		// Iterate over all faces in the shape
		for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next())
		{
			const TopoDS_Face& face = TopoDS::Face(explorer.Current());

			// Create a GProp_GProps object and pass it the surface
			GProp_GProps props;
			BRepGProp::SurfaceProperties(face, props);

			// Get the surface area of the face
			double area = props.Mass();

			// Add the area of this face to the total area
			total_area += area;
		}

		return total_area;
	}
	else
	{
		std::cout << "Unsupported file formats, please check" << std::endl;
	}

	
}



void triangular(const char* input_file, const char* output_file, double ratio = 0.0001)
{

	//std::string a = input_file;
	//int num = a.find_last_of('.');
	//std::string b = a.substr(num + 1);
	//if (b == "step" || b == "stp")
	//{
	//	// Create a STEPControl_Reader object
	//	STEPControl_Reader reader;
	//	// Read the file and check the status
	//	IFSelect_ReturnStatus status = reader.ReadFile(input_file);
	//	// Transfer the shape from the file to the process
	//	reader.TransferRoots();
	//	// Get the shape from the process
	//	TopoDS_Shape  shape = reader.OneShape();

	//	Bnd_Box box;
	//	BRepBndLib::Add(shape, box);
	//	Standard_Real xmin, ymin, zmin, xmax, ymax, zmax;
	//	box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
	//	Standard_Real dx = xmax - xmin;
	//	Standard_Real dy = ymax - ymin;
	//	Standard_Real dz = zmax - zmin;
	//	Standard_Real diagonal = sqrt(dx * dx + dy * dy + dz * dz);
	//	std::cout << "diagonal = " << diagonal << std::endl;

	//	const Standard_Real linear_deflection = diagonal * ratio; // The maximum distance between the original surface and the mesh
	//	const Standard_Real angular_deflection = 50; // The maximum angle between two adjacent triangles
	//	// Generate the mesh using BRepMesh_IncrementalMesh
	//	BRepMesh_IncrementalMesh mesher(shape, linear_deflection, Standard_False, angular_deflection, Standard_True);
	//	mesher.Perform();
	//	// Write the mesh to the stl file using StlAPI_Writer
	//	StlAPI_Writer writer;
	//	writer.Write(shape, output_file);
	//}
	//else if (b == "igs" || b == "iges")
	//{
	//	IGESControl_Reader reader;
	//	// Read the file and check the status
	//	IFSelect_ReturnStatus status = reader.ReadFile(input_file);
	//	// Transfer the shape from the file to the process
	//	reader.TransferRoots();
	//	// Get the shape from the process
	//	TopoDS_Shape  shape = reader.OneShape();

	//	Bnd_Box box;
	//	BRepBndLib::Add(shape, box);
	//	Standard_Real xmin, ymin, zmin, xmax, ymax, zmax;
	//	box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
	//	Standard_Real dx = xmax - xmin;
	//	Standard_Real dy = ymax - ymin;
	//	Standard_Real dz = zmax - zmin;
	//	Standard_Real diagonal = sqrt(dx * dx + dy * dy + dz * dz);
	//	std::cout << "diagonal = " << diagonal << std::endl;

	//	const Standard_Real linear_deflection = diagonal * ratio; // The maximum distance between the original surface and the mesh
	//	const Standard_Real angular_deflection = 50; // The maximum angle between two adjacent triangles
	//	// Generate the mesh using BRepMesh_IncrementalMesh
	//	BRepMesh_IncrementalMesh mesher(shape, linear_deflection, Standard_False, angular_deflection, Standard_True);
	//	mesher.Perform();
	//	// Write the mesh to the stl file using StlAPI_Writer
	//	StlAPI_Writer writer;
	//	writer.Write(shape, output_file);
	//}
	//else
	//{
	//	std::cout << "Unsupported file formats, please check" << std::endl;
	//}
	
	
};

void convert(const char* conv_in, const char* conv_out)
{
	std::ofstream outfile(conv_out);

	std::ifstream file(conv_in);
	std::vector<double> vertex_x;
	std::vector<double> vertex_y;
	std::vector<double> vertex_z;
	int j = 0; 
	int face = 0; 
	double tolerance = 1e-8;
	int flag_3 = 0; 
	int flag = 0;
	std::vector<std::vector<int> > group;
	std::vector<int> group_row;
	if (file.is_open())
	{
		std::string line;
		while (getline(file, line))
		{
			if (line.find("facet normal") != std::string::npos)
			{
				face = face + 1;
			}
			if (line.find("vertex") != std::string::npos)
			{
				std::stringstream ss(line);
				std::string vertex;
				double x, y, z;

				ss >> vertex >> x >> y >> z;
		

				for (int m = 0; m <= j; m++)
				{

					if (m == j || j == 0)
					{

						vertex_x.push_back(x);
						vertex_y.push_back(y);
						vertex_z.push_back(z);

						
						group_row.push_back(j + 1);
						j = j + 1;

						break;

					}
					
					if ((abs(x - vertex_x[m]) <= tolerance) && (abs(y - vertex_y[m]) <= tolerance) && (abs(z - vertex_z[m]) <= tolerance))
					{
						group_row.push_back(m + 1);
						break;
					}
				}
				flag_3++;
				if (flag_3 == 3)
				{
					flag_3 = 0;
					flag++;
					group.push_back(group_row);
					
					group_row.clear();
				}
			}
		}
		std::cout << "j = " << j << std::endl;
		std::cout << "face = " << face << std::endl;

		for (int k = 0; k < (j); k++)
		{
			outfile << "v" << " " << vertex_x[k] << " " << vertex_y[k] << " " << vertex_z[k] << std::endl;
		}

		for (int k = 0; k < flag; k++)
		{
			if ((group[k][0] != group[k][1]) && (group[k][1] != group[k][2]) && (group[k][0] != group[k][2]))
			{
				outfile << "f" << " " << group[k][0] << " " << group[k][1] << " " << group[k][2] << std::endl;
			}

		}
		file.close();
		outfile.close();
	}
	else
	{
		std::cout << "failed to open file" << std::endl;

	}
}
