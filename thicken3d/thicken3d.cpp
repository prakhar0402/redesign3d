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

#include <armadillo>

#include <iostream>
#include <fstream>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef std::map<Polyhedron::Halfedge_const_handle, HalfedgeMemo> he_memo_map;
typedef std::map<Polyhedron::Vertex_const_handle, VertexMemo> v_memo_map;

const double t1 = 0.1;
const double t2 = 0.2;

const double K_sdf = 100.0;
const double K_s = 20.0;
const double K_d = 20.0;

double normalize_diameter(double diameter, std::pair<double, double>& min_max_sdf)
{
	return (diameter - min_max_sdf.first) / (min_max_sdf.second - min_max_sdf.first);
}

size_t populate_memos(Polyhedron& mesh, boost::associative_property_map<v_memo_map>& Vertex_memo_map, boost::associative_property_map<he_memo_map>& Halfedge_memo_map)
{
	// computing normals
	typedef boost::graph_traits<Polyhedron>::vertex_descriptor vd;
	typedef std::map<Polyhedron::Vertex_const_handle, Kernel::Vector_3> vertex_normal_map;
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
	std::pair<double, double> min_max_sdf = CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

	// print minimum & maximum SDF values
	std::cout << "minimum SDF: " << min_max_sdf.first
		<< " maximum SDF: " << min_max_sdf.second << std::endl;

	double nt1 = normalize_diameter(t1, min_max_sdf);
	double nt2 = normalize_diameter(t2, min_max_sdf);

	// loop over each vertex
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin();
		vi != mesh.vertices_end(); ++vi)
	{
		VertexMemo vmemo(vi);
		boost::put(Vertex_memo_map, vi, vmemo);
	}

	double sdf_sum = 0.0, v_sdf = 0.0;
	int count = 0;
	size_t movable = 0;
	Polyhedron::Halfedge_const_handle he1;
	Polyhedron::Halfedge_const_handle he2;
	// loop over each vertex
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin();
		vi != mesh.vertices_end(); ++vi)
	{
		he1 = vi->halfedge();
		he2 = he1;
		sdf_sum = 0.0;
		count = 0;
		// add sdf values from neighboring faces
		do
		{
			HalfedgeMemo hem(he2);
			hem.compute_force(K_s, K_d, Vertex_memo_map[vi].get_velocity(), Vertex_memo_map[he2->opposite()->vertex()].get_velocity());
			boost::put(Halfedge_memo_map, he2, hem);

			sdf_sum += sdf_property_map[he2->facet()];
			count++;
			he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
		} while (he2 != he1);
		// map vertex to average of neighboring face sdf values

		v_sdf = sdf_sum / count;

		Vertex_memo_map[vi].set_normal(normal_map[vi]);
		Vertex_memo_map[vi].set_sdf(v_sdf);
		Vertex_memo_map[vi].compute_force(K_sdf, nt1, Halfedge_memo_map);

		// index is zero if the vertex is fixed
		if (v_sdf <= nt2)
			Vertex_memo_map[vi].index = ++movable;
		else
			Vertex_memo_map[vi].index = 0;

		//std::cout << Vertex_memo_map[vi].get_Jacobian_vel() << std::endl;
	}
	std::cout << std::endl;

	return movable;
}


int main()
{
	// create and read Polyhedron
	Polyhedron mesh;
	std::ifstream input("data/cactus.off");
	if (!input || !(input >> mesh) || mesh.empty()) {
		std::cerr << "Not a valid off file." << std::endl;
		return EXIT_FAILURE;
	}

	// compute vertex normals
	//compute_vertex_normals(mesh);

	// compute sdf values at vertex by averaging neighboring face sdf values
	// TODO: weighted average
	//compute_vertex_sdf(mesh);


	// create a property-map for halfedges
	he_memo_map hmm;
	boost::associative_property_map<he_memo_map> Halfedge_memo_map(hmm);
	// create a property-map for vertices
	v_memo_map vmm;
	boost::associative_property_map<v_memo_map> Vertex_memo_map(vmm);

	size_t N = populate_memos(mesh, Vertex_memo_map, Halfedge_memo_map); // N is the number of movable vertices
	//populate_halfedge_memo(mesh, Halfedge_memo_map);

	
	//size_t N = mesh.size_of_vertices();

	arma::vec FORCE;
	arma::mat JPOS;
	arma::mat JVEL;

	FORCE.set_size(N * 3);
	JPOS.set_size(N * 3, N * 3);
	JVEL.set_size(N * 3, N * 3);

	FORCE.fill(0.0);
	JPOS.fill(0.0);
	JVEL.fill(0.0);

	// assemble
	size_t id = 0;
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin();
		vi != mesh.vertices_end(); ++vi)
	{
		id = Vertex_memo_map[vi].index;
		if (id > 0)
		{
			//FORCE(id - 1, N + id - 1, 2 * N + id - 1) = Vertex_memo_map[vi].get_force();
		}
	}

	// Testing armadillo library
	std::cout << "Armadillo version: " << arma::arma_version::as_string() << std::endl;

	arma::mat A(2, 3);  // directly specify the matrix size (elements are uninitialised)

	std::cout << "A.n_rows: " << A.n_rows << std::endl;  // .n_rows and .n_cols are read only
	std::cout << "A.n_cols: " << A.n_cols << std::endl;

	A.set_size(4, 5); // change the size (data is not preserved)

	A.fill(5.0);     // set all elements to a particular value
	A.print("A:");

	// endr indicates "end of row"
	A << 0.165300 << 0.454037 << 0.995795 << 0.124098 << 0.047084 << arma::endr
		<< 0.688782 << 0.036549 << 0.552848 << 0.937664 << 0.866401 << arma::endr
		<< 0.348740 << 0.479388 << 0.506228 << 0.145673 << 0.491547 << arma::endr
		<< 0.148678 << 0.682258 << 0.571154 << 0.874724 << 0.444632 << arma::endr
		<< 0.245726 << 0.595218 << 0.409327 << 0.367827 << 0.385736 << arma::endr;

	A.print("A:");

	// determinant
	std::cout << "det(A): " << arma::det(A) << std::endl;

	// inverse
	std::cout << "inv(A): " << std::endl << arma::inv(A) << std::endl;

	arma::vec v1, v2, v12, sv12;
	v2.set_size(3);
	v1 << 1 << 2 << 3;
	v2.fill(1.5);
	v12 = v1 - v2;
	sv12 = v12 % v12;
	v1.print("v1:");
	v2.print("v2:");
	arma::mat m = v1 * arma::trans(v2);
	m.print("m:");
	v12.print("v12:");
	sv12.print("sv12:");
}