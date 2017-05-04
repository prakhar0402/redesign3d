// thicken3d.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "halfedgememo.h"
#include "vertexmemo.h"

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h> 
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <armadillo>

#include <iostream>
#include <fstream>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef std::map<Polyhedron::Halfedge_const_handle, HalfedgeMemo> he_memo_map;
typedef std::map<Polyhedron::Vertex_const_handle, VertexMemo> v_memo_map;

typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef boost::graph_traits<Polyhedron> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;

const std::string FILENAME = "data/cactus.off";
const std::string IDENTIFIER = "0007";

const double T1 = 0.1;
const double T2 = 0.15;

const double K_SDF = 20.0;
const double K_S = 30.0;
const double K_D = 2.0;

const double VERTEX_MASS = 1.0;
const double MASS_INV = 1.0 / VERTEX_MASS;

const double TIME_STEP = 0.5;
const size_t STEPS = 10;

std::pair<double, double> min_max_sdf;
double nT1, nT2;
size_t nMove;

void generateOutput(
	const std::string filename,
	const Polyhedron& mesh,
	const boost::associative_property_map<v_memo_map>& Vertex_memo_map
)
{
	std::ofstream ofile(filename);
	ofile << "Input File:\n" << FILENAME << "\n";
	ofile << "Output Identifier:\n" << IDENTIFIER << "\n";
	ofile << "Strict Threshold:\n" << T1 << "\n";
	ofile << "Moderate Threshold:\n" << T2 << "\n";
	ofile << "Force Constant SDF:\n" << K_SDF << "\n";
	ofile << "Force Constant Spring:\n" << K_S << "\n";
	ofile << "Force Constant Damper:\n" << K_D << "\n";
	ofile << "Vertex Mass:\n" << VERTEX_MASS << "\n";
	ofile << "Time Step:\n" << TIME_STEP << "\n";
	ofile << "Steps:\n" << STEPS << "\n";
	ofile << "Min SDF:\n" << min_max_sdf.first << "\n";
	ofile << "Max SDF:\n" << min_max_sdf.second << "\n";
	ofile << "Number of Movable Vertices:\n" << nMove << "\n";
	ofile << "Total Number of Vertices:\n" << mesh.size_of_vertices() << "\n";
	ofile << "SDF Values:\n";
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
		ofile << Vertex_memo_map[vi].get_sdf() << "\n";
	ofile << "Normal Vectors:\n";
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
		ofile << arma::trans(Vertex_memo_map[vi].get_normal());
	ofile.close();
}

double normalize_diameter(
	const double diameter,
	const std::pair<double, double>& min_max_sdf
)
{
	return (diameter - min_max_sdf.first) / (min_max_sdf.second - min_max_sdf.first);
}

// compute element forces and Jacobians on edge elements
void populate_halfedge_memos(
	const Polyhedron& mesh,
	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	boost::associative_property_map<he_memo_map>& Halfedge_memo_map
)
{
	// loop over each halfedge to initialize halfedge memos
	for (Polyhedron::Halfedge_const_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi)
	{
		HalfedgeMemo hem(hi);
		boost::put(Halfedge_memo_map, hi, hem);
	}

	Polyhedron::Halfedge_const_handle he1;
	Polyhedron::Halfedge_const_handle he2;
	// loop over each vertex
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		he1 = vi->halfedge();
		he2 = he1;
		// compute forces and Jacobians
		do
		{
			if (Vertex_memo_map[vi].index > 0 && Halfedge_memo_map[he2].isMaster == 0) // perform computation for movable edges only
			{
				Halfedge_memo_map[he2].isMaster = 1; // set as master
				Halfedge_memo_map[he2->opposite()].isMaster = -1; // set the opposite as slave
				Halfedge_memo_map[he2].compute_force(K_S, K_D, Vertex_memo_map[vi].get_velocity(), Vertex_memo_map[he2->opposite()->vertex()].get_velocity());
			}

			he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
		} while (he2 != he1);
	}
}


size_t populate_memos(
	const Polyhedron& mesh,
	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	boost::associative_property_map<he_memo_map>& Halfedge_memo_map
)
{
	// computing normals
	typedef std::map<Polyhedron::Vertex_const_handle, Vector> vertex_normal_map;
	vertex_normal_map vnm;
	boost::associative_property_map<vertex_normal_map> normal_map(vnm);
	CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, normal_map);

	// create a property-map for faces
	typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
	Facet_double_map fdm;
	boost::associative_property_map<Facet_double_map> sdf_property_map(fdm);

	// compute SDF values
	// std::pair<double, double> min_max_sdf = CGAL::sdf_values(mesh, sdf_property_map);

	// It is possible to compute the raw SDF values and post-process them using
	// the following lines:
	const std::size_t number_of_rays = 25;  // cast 25 rays per facet
	const double cone_angle = 2.0 / 3.0 * CGAL_PI; // set cone opening-angle
	CGAL::sdf_values(mesh, sdf_property_map, cone_angle, number_of_rays, false);
	min_max_sdf = CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

	// print minimum & maximum SDF values
	std::cout << "minimum SDF: " << min_max_sdf.first	<< " maximum SDF: " << min_max_sdf.second << std::endl;

	nT1 = normalize_diameter(T1, min_max_sdf);
	nT2 = normalize_diameter(T2, min_max_sdf);

	// loop over each vertex to initialize vertex memos
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		VertexMemo vmemo(vi);
		boost::put(Vertex_memo_map, vi, vmemo);
	}

	double sdf_sum = 0.0, v_sdf = 0.0;
	int count = 0;
	size_t movable = 0;
	Polyhedron::Halfedge_const_handle he1;
	Polyhedron::Halfedge_const_handle he2;
	// loop over each vertex to comoute vertex sdf
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		he1 = vi->halfedge();
		he2 = he1;
		sdf_sum = 0.0;
		count = 0;
		// add sdf values from neighboring faces
		do
		{
			sdf_sum += sdf_property_map[he2->facet()];
			count++;
			he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
		} while (he2 != he1);
		// map vertex to average of neighboring face sdf values

		v_sdf = sdf_sum / count;

		Vertex_memo_map[vi].set_normal(normal_map[vi]);
		Vertex_memo_map[vi].set_sdf(v_sdf);
		Vertex_memo_map[vi].compute_sdf_force(K_SDF, nT1);

		// index is zero if the vertex is fixed
		if (v_sdf <= nT2)
			Vertex_memo_map[vi].index = ++movable;
		else
			Vertex_memo_map[vi].index = 0;
	}

	populate_halfedge_memos(mesh, Vertex_memo_map, Halfedge_memo_map);

	return movable;
}

// assembly
void global_assembly(
	arma::vec& FORCE,
	arma::mat& JPOS,
	arma::mat& JVEL,
	arma::vec& VELOCITY,
	const size_t& N,
	const Polyhedron& mesh,
	const boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	const boost::associative_property_map<he_memo_map>& Halfedge_memo_map
)
{
	FORCE.set_size(N * 3);
	JPOS.set_size(N * 3, N * 3);
	JVEL.set_size(N * 3, N * 3);
	VELOCITY.set_size(N * 3);

	FORCE.fill(0.0);
	JPOS.fill(0.0);
	JVEL.fill(0.0);
	VELOCITY.fill(0.0);

	int s_idx, e_idx;
	// loop over each halfedge
	for (Polyhedron::Halfedge_const_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi)
	{
		if (Halfedge_memo_map[hi].isMaster > 0)
		{
			s_idx = ((int)Vertex_memo_map[hi->vertex()].index - 1)*3;
			e_idx = ((int)Vertex_memo_map[hi->opposite()->vertex()].index - 1)*3;
			FORCE(arma::span(s_idx, s_idx + 2)) += Halfedge_memo_map[hi].get_force();
			JPOS(arma::span(s_idx, s_idx + 2), arma::span(s_idx, s_idx + 2)) += Halfedge_memo_map[hi].get_Jacobian_pos();
			JVEL(arma::span(s_idx, s_idx + 2), arma::span(s_idx, s_idx + 2)) += Halfedge_memo_map[hi].get_Jacobian_vel();

			if (e_idx >= 0) // if the opposite vertex is also movable
			{
				FORCE(arma::span(e_idx, e_idx + 2)) -= Halfedge_memo_map[hi].get_force();

				JPOS(arma::span(e_idx, e_idx + 2), arma::span(e_idx, e_idx + 2)) += Halfedge_memo_map[hi].get_Jacobian_pos();
				JPOS(arma::span(s_idx, s_idx + 2), arma::span(e_idx, e_idx + 2)) -= Halfedge_memo_map[hi].get_Jacobian_pos();
				JPOS(arma::span(e_idx, e_idx + 2), arma::span(s_idx, s_idx + 2)) -= Halfedge_memo_map[hi].get_Jacobian_pos();

				JVEL(arma::span(e_idx, e_idx + 2), arma::span(e_idx, e_idx + 2)) += Halfedge_memo_map[hi].get_Jacobian_vel();
				JVEL(arma::span(s_idx, s_idx + 2), arma::span(e_idx, e_idx + 2)) -= Halfedge_memo_map[hi].get_Jacobian_vel();
				JVEL(arma::span(e_idx, e_idx + 2), arma::span(s_idx, s_idx + 2)) -= Halfedge_memo_map[hi].get_Jacobian_vel();
			}
		}
	}

	// loop over each vertex
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		if (Vertex_memo_map[vi].index > 0)
		{
			s_idx = (Vertex_memo_map[vi].index - 1)*3;
			FORCE(arma::span(s_idx, s_idx + 2)) += Vertex_memo_map[vi].get_sdf_force();
			VELOCITY(arma::span(s_idx, s_idx + 2)) = Vertex_memo_map[vi].get_velocity(); //TODO: this is unnecessary during updation
		}
	}
}

// compute update in velocity and position for one time step
void march_one_time_step(
	arma::vec& delta_VEL,
	arma::vec& delta_POS,
	const size_t& N,
	const arma::vec& FORCE,
	const arma::mat& JPOS,
	const arma::mat& JVEL,
	const arma::vec& VELOCITY
)
{
	arma::mat A;
	arma::vec b;

	A = arma::eye(N * 3, N * 3) - TIME_STEP*MASS_INV*JVEL - TIME_STEP*TIME_STEP*MASS_INV*JPOS;
	b = TIME_STEP*MASS_INV*(FORCE + TIME_STEP*JPOS*VELOCITY);

	delta_VEL = arma::inv(A)*b;
	delta_POS = TIME_STEP*(VELOCITY + delta_VEL);
}

// update the positions, velocities, and memos
void update(
	const arma::vec& delta_POS,
	const arma::vec& delta_VEL,
	const size_t& N,
	Polyhedron& mesh,
	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	boost::associative_property_map<he_memo_map>& Halfedge_memo_map,
	arma::vec& FORCE,
	arma::mat& JPOS,
	arma::mat& JVEL,
	arma::vec& VELOCITY
)
{
	size_t idx;
	
	// First obtain the type of the property map:
	typedef boost::property_map<Polyhedron, CGAL::vertex_point_t>::type PM;

	// Then obtain the property map for P
	PM coord = get(CGAL::vertex_point, mesh);

	// update vertex coordinates
	BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
		if (Vertex_memo_map[vd].index > 0)
		{
			idx = (Vertex_memo_map[vd].index - 1) * 3;
			// Now access the property
			Point p = get(coord, vd);
			p = p + Vector(delta_POS(idx), delta_POS(idx + 1), delta_POS(idx + 2));
			put(coord, vd, p);
		}
	}

	// TODO: decide whether to update normals and sdf_force every time step or not
	// updating normals
	typedef std::map<Polyhedron::Vertex_const_handle, Vector> vertex_normal_map;
	vertex_normal_map vnm;
	boost::associative_property_map<vertex_normal_map> normal_map(vnm);
	CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, normal_map);

	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		Vertex_memo_map[vi].set_normal(normal_map[vi]);
		Vertex_memo_map[vi].compute_sdf_force(K_SDF, nT1);
		if (Vertex_memo_map[vi].index > 0)
		{
			idx = (Vertex_memo_map[vi].index - 1) * 3;
			Vertex_memo_map[vi].set_velocity(Vertex_memo_map[vi].get_velocity() + delta_VEL(arma::span(idx, idx + 2)));
		}
	}

	for (Polyhedron::Halfedge_const_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi)
	{
		Halfedge_memo_map[hi].compute_length();
		if (Halfedge_memo_map[hi].isMaster > 0)
			Halfedge_memo_map[hi].compute_force(K_S, K_D, Vertex_memo_map[hi->vertex()].get_velocity(), Vertex_memo_map[hi->opposite()->vertex()].get_velocity());
	}

	global_assembly(FORCE, JPOS, JVEL, VELOCITY, N, mesh, Vertex_memo_map, Halfedge_memo_map);
}

// march one step and then update
void march_and_update(
	const size_t& N,
	Polyhedron& mesh,
	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	boost::associative_property_map<he_memo_map>& Halfedge_memo_map,
	arma::vec& FORCE,
	arma::mat& JPOS,
	arma::mat& JVEL,
	arma::vec& VELOCITY
)
{
	arma::vec delta_VEL;
	arma::vec delta_POS;
	march_one_time_step(delta_VEL, delta_POS, N, FORCE, JPOS, JVEL, VELOCITY);
	update(delta_POS, delta_VEL, N, mesh, Vertex_memo_map, Halfedge_memo_map, FORCE, JPOS, JVEL, VELOCITY);
}


int main()
{
	// create and read Polyhedron
	Polyhedron mesh;
	std::ifstream input(FILENAME);
	if (!input || !(input >> mesh) || mesh.empty()) {
		std::cerr << "Not a valid off file." << std::endl;
		return EXIT_FAILURE;
	}

	// create a property-map for halfedges
	he_memo_map hmm;
	boost::associative_property_map<he_memo_map> Halfedge_memo_map(hmm);
	// create a property-map for vertices
	v_memo_map vmm;
	boost::associative_property_map<v_memo_map> Vertex_memo_map(vmm);

	// populate memos - compute and store normals, sdf, forces, and Jacobians
	std::cout << "Populating memos..." << std::endl;
	nMove = populate_memos(mesh, Vertex_memo_map, Halfedge_memo_map); // nMove is the number of movable vertices
	std::cout << "Populating completed!\nNumber of movable vertices = " << nMove << std::endl;
	std::cout << std::endl;

	arma::vec FORCE;
	arma::mat JPOS;
	arma::mat JVEL;
	arma::vec VELOCITY;

	// assemble
	std::cout << "Assembling..." << std::endl;
	global_assembly(FORCE, JPOS, JVEL, VELOCITY, nMove, mesh, Vertex_memo_map, Halfedge_memo_map);
	std::cout << "Assembling completed!" << std::endl;
	std::cout << std::endl;

	// numerical time integration
	std::cout << "Performing numerical integration..." << std::endl;
	for (size_t i = 0; i < STEPS; i++)
	{
		march_and_update(nMove, mesh, Vertex_memo_map, Halfedge_memo_map, FORCE, JPOS, JVEL, VELOCITY);
	}
	std::cout << "Integration completed!" << std::endl;
	std::cout << std::endl;

	// populate memos - compute and store normals, sdf, forces, and Jacobians
	std::cout << "Populating memos..." << std::endl;
	nMove = populate_memos(mesh, Vertex_memo_map, Halfedge_memo_map); // nMove is the number of movable vertices
	std::cout << "Populating completed!\nNumber of movable vertices = " << nMove << std::endl;
	std::cout << std::endl;

	std::string meshOutFile = FILENAME.substr(0, FILENAME.size() - 4) + "_" + IDENTIFIER + ".off";
	std::string outputFile = FILENAME.substr(0, FILENAME.size() - 4) + "_" + IDENTIFIER + ".out";

	std::ofstream output(meshOutFile);
	output << mesh;
	output.close();

	generateOutput(outputFile, mesh, Vertex_memo_map);

	//// Testing armadillo library
	//std::cout << "Armadillo version: " << arma::arma_version::as_string() << std::endl;

	//arma::mat A(2, 3);  // directly specify the matrix size (elements are uninitialised)

	//std::cout << "A.n_rows: " << A.n_rows << std::endl;  // .n_rows and .n_cols are read only
	//std::cout << "A.n_cols: " << A.n_cols << std::endl;

	//A.set_size(4, 5); // change the size (data is not preserved)

	//A.fill(5.0);     // set all elements to a particular value
	//A.print("A:");

	//// endr indicates "end of row"
	//A << 0.165300 << 0.454037 << 0.995795 << 0.124098 << 0.047084 << arma::endr
	//	<< 0.688782 << 0.036549 << 0.552848 << 0.937664 << 0.866401 << arma::endr
	//	<< 0.348740 << 0.479388 << 0.506228 << 0.145673 << 0.491547 << arma::endr
	//	<< 0.148678 << 0.682258 << 0.571154 << 0.874724 << 0.444632 << arma::endr
	//	<< 0.245726 << 0.595218 << 0.409327 << 0.367827 << 0.385736 << arma::endr;

	//A.print("A:");

	//// determinant
	//std::cout << "det(A): " << arma::det(A) << std::endl;

	//// inverse
	//std::cout << "inv(A): " << std::endl << arma::inv(A) << std::endl;

	//arma::vec v1, v2, v12, sv12;
	//v2.set_size(3);
	//v1 << 1 << 2 << 3;
	//v2.fill(1.5);
	//v12 = v1 - v2;
	//sv12 = v12 % v12;
	//v1.print("v1:");
	//v2.print("v2:");
	//arma::mat m = v1 * arma::trans(v2);
	//m.print("m:");
	//v12.print("v12:");
	//sv12.print("sv12:");
}