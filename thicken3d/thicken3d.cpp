/*
MAD Lab, University at Buffalo
Copyright (C) 2018  Prakhar Jaiswal <prakharj@buffalo.edu>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "halfedgememo.h"
#include "vertexmemo.h"

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h> 
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <random>
#include <iostream>
#include <fstream>

typedef std::map<Polyhedron::Halfedge_const_handle, HalfedgeMemo> he_memo_map;
typedef std::map<Polyhedron::Vertex_const_handle, VertexMemo> v_memo_map;

typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef boost::graph_traits<Polyhedron> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;

// CONSTANTS
std::string FILENAME = "../data/cactus/cactus.off";
std::string IDENTIFIER = "test";

// SDF parameters
const size_t number_of_rays = 30;  // cast rays per facet
const double cone_angle = CGAL_PI * 2.0 / 3.0; // set cone opening-angle
std::default_random_engine re;

double Td = 0.1;

double K_S = 64.0;
double K_D = 1024.0;
double ksdf_multiplier = 1.0;
double K_SDF = ksdf_multiplier * K_S * Td;

double VERTEX_MASS = 1.0;
double MASS_INV = 1.0 / VERTEX_MASS;

double TIME_STEP = 1.0;
double TOTAL_TIME = 100.0;
size_t STEPS = (size_t)std::ceil(TOTAL_TIME / TIME_STEP);

size_t N_BATCH = 4; // number of times SDF and spring-damper system is reset during iterations (including initial)
size_t BATCH_SIZE = (size_t)std::ceil((double)STEPS / N_BATCH);

// Global Variables
std::pair<double, double> min_max_sdf;
double nTd;
size_t nMove = 0;
arma::vec max_change_vec(STEPS);
double rho = 1.0; // fabricability

// normalize between 0 and 1 based on min and max SDF values
double normalize_diameter(
	const double diameter,
	const std::pair<double, double>& min_max_sdf
)
{
	return (diameter - min_max_sdf.first) / (min_max_sdf.second - min_max_sdf.first);
}

// unnormalize from 0 and 1 based on min and max SDF values
double unnormalize_diameter(
	const double diameter,
	const std::pair<double, double>& min_max_sdf
)
{
	return diameter * (min_max_sdf.second - min_max_sdf.first) + min_max_sdf.first;
}

//// set n_ring neighborhood to movable
//void setMovable(
//	const Polyhedron& mesh,
//	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
//	Polyhedron::Vertex_const_iterator vi,
//	size_t n_ring
//)
//{
//	Vertex_memo_map[vi].isMovable = true;
//	if (n_ring <= 0) return;
//
//	Polyhedron::Halfedge_const_handle he1;
//	Polyhedron::Halfedge_const_handle he2;
//	he1 = vi->halfedge();
//	he2 = he1;
//	do
//	{
//		setMovable(mesh, Vertex_memo_map, he2->opposite()->vertex(), n_ring - 1);
//		he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
//	} while (he2 != he1);
//}

// specifies n-ring neighbor of a vertex and computes decaying forces on them for smoothness
void setSDFNeighbor(
	const Polyhedron& mesh,
	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	Polyhedron::Vertex_const_iterator vi,
	double force,
	size_t n_ring
)
{
	Vertex_memo_map[vi].isMovable = true;
	Vertex_memo_map[vi].setMaxSDFMag(force);
	if (n_ring <= 0) return;

	double f = force / 1.2;
	Polyhedron::Halfedge_const_handle he1;
	Polyhedron::Halfedge_const_handle he2;
	he1 = vi->halfedge();
	he2 = he1;
	do
	{
		setSDFNeighbor(mesh, Vertex_memo_map, he2->opposite()->vertex(), f, n_ring - 1);
		he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
	} while (he2 != he1);
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
			if (Vertex_memo_map[vi].index > 0 && Halfedge_memo_map[he2].isMaster == 0) // perform computation for movable master edges only
			{
				Halfedge_memo_map[he2].isMaster = 1; // set as master
				Halfedge_memo_map[he2->opposite()].isMaster = -1; // set the opposite as slave - avoids repeated computations
				Halfedge_memo_map[he2].compute_force(K_S, K_D, Vertex_memo_map[vi].get_velocity(), Vertex_memo_map[he2->opposite()->vertex()].get_velocity());
			}

			he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
		} while (he2 != he1);
	}
}

// compute element forces and Jacobians on vertex elements
size_t populate_memos(
	const Polyhedron& mesh,
	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	boost::associative_property_map<he_memo_map>& Halfedge_memo_map,
	bool SDFexists = false
)
{
	// computing vertex normals
	typedef std::map<Polyhedron::Vertex_const_handle, Vector> vertex_normal_map;
	vertex_normal_map vnm;
	boost::associative_property_map<vertex_normal_map> vnormal_map(vnm);
	CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, vnormal_map);

	// computing face normals
	typedef std::map<Polyhedron::Facet_const_handle, Vector> facet_normal_map;
	facet_normal_map fnm;
	boost::associative_property_map<facet_normal_map> fnormal_map(fnm);
	CGAL::Polygon_mesh_processing::compute_face_normals(mesh, fnormal_map);

	// create a property-map for faces
	typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
	Facet_double_map fdm;
	boost::associative_property_map<Facet_double_map> sdf_property_map(fdm);

	// use pre-existing face sdf values if already computed before
	std::string face_sdf_file = FILENAME.substr(0, FILENAME.size() - 4) + ".facesdf";
	std::ifstream ifile(face_sdf_file);
	if (SDFexists && ifile.good()) // if the file exists and can be read
	{
		// read face sdf values from the file
		ifile >> min_max_sdf.first;
		ifile >> min_max_sdf.second;
		double val;
		for (Polyhedron::Facet_const_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi)
		{
			ifile >> val;
			boost::put(sdf_property_map, fi, val);
		}
	}
	else // if file doesn't exist or can't be read, compute face sdf values and write to file
	{
		ifile.close();

		// compute SDF values
		// std::pair<double, double> min_max_sdf = CGAL::sdf_values(mesh, sdf_property_map);

		// It is possible to compute the raw SDF values and post-process them using
		// the following lines:
		CGAL::sdf_values(mesh, sdf_property_map, cone_angle, number_of_rays, false);
		min_max_sdf = CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

		// print minimum & maximum SDF values
		std::cout << "minimum SDF: " << min_max_sdf.first << " maximum SDF: " << min_max_sdf.second << std::endl;

		if (SDFexists)
		{
			// write face sdf values to a file
			std::ofstream ofile(face_sdf_file);
			ofile << min_max_sdf.first << "\n";
			ofile << min_max_sdf.second << "\n";
			for (Polyhedron::Facet_const_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi)
				ofile << sdf_property_map[fi] << "\n";
		}
	}

	// normalize the threshold diameter based on min and max SDF values
	nTd = normalize_diameter(Td, min_max_sdf);

	// compute face area
	Facet_double_map fam;
	boost::associative_property_map<Facet_double_map> Facet_area_map(fam);
	// loop over each face to compute face area and volume
	double area = 1.0;
	double vol = 0.0, numA = 0.0, denA = 0.0, numV = 0.0, denV = 0.0;
	for (Polyhedron::Facet_const_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi)
	{
		area = Kernel::Compute_area_3()(fi->halfedge()->vertex()->point(), fi->halfedge()->next()->vertex()->point(), fi->halfedge()->opposite()->vertex()->point());
		boost::put(Facet_area_map, fi, area);

		vol = area*unnormalize_diameter(sdf_property_map[fi], min_max_sdf); // should be divided by 2 for radius, but later 2 is cancelled in ratio
		if (sdf_property_map[fi] <= nTd)
		{
			numA += area;
			numV += vol;
		}
		denA += area;
		denV += vol;
	}
	// computes fabricability as fraction of estimated volume and area
	//double wA = 0.5, wV = 0.5;
	//rho = 1.0 - wA*(numA / denA) - wV*(numV / denV); // additive
	double wA = 1.0, wV = 0.5;
	rho = 1.0 - 0.5*pow((numA / denA), 1.0 / wA) - 0.5*pow((numV / denV), 1.0 / wV); // exponent

	// loop over each vertex to initialize vertex memos
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		VertexMemo vmemo(vi);
		boost::put(Vertex_memo_map, vi, vmemo);
	}

	double sdf_sum = 0.0, v_sdf = 0.0, proj_area = 0.0,  v_area = 0.0, area_avg = 0.0;
	int count = 0;
	size_t movable = 0;
	Polyhedron::Halfedge_const_handle he1;
	Polyhedron::Halfedge_const_handle he2;
	// loop over each vertex to compute vertex sdf
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		he1 = vi->halfedge();
		he2 = he1;
		sdf_sum = 0.0;
		count = 0;
		v_area = 0.0;
		// add sdf values from neighboring faces with projected area as weights
		do
		{
			proj_area = Facet_area_map[he2->facet()] * std::abs(Kernel::Compute_scalar_product_3()(vnormal_map[vi], fnormal_map[he2->facet()]));
			sdf_sum += sdf_property_map[he2->facet()] * proj_area;
			v_area += proj_area;
			count++;
			he2 = he2->next_on_vertex(); // loop over halfedges on the vertex
		} while (he2 != he1);
		// map vertex to weighted average of neighboring face sdf values
		v_sdf = sdf_sum / v_area;
		area_avg += v_area;

		Vertex_memo_map[vi].set_normal(vnormal_map[vi]);
		Vertex_memo_map[vi].set_sdf(v_sdf);
		Vertex_memo_map[vi].set_area_factor(v_area);
		Vertex_memo_map[vi].set_ref_point(vi->point() - unnormalize_diameter(v_sdf, min_max_sdf)*vnormal_map[vi] / 2.0);

		//if (v_sdf <= nTd)
		//	setMovable(mesh, Vertex_memo_map, vi, 3); // set 3-ring neighborhood vertices to movable
	}

	area_avg /= mesh.size_of_vertices();
	//std::cout << "Average area = " << area_avg << std::endl;

	double force = 0.0;
	// normalize area factor and set movable vertices neighborhood
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		Vertex_memo_map[vi].set_area_factor(Vertex_memo_map[vi].get_area_factor() / area_avg); // use average face area to normalize area factor
		if (Vertex_memo_map[vi].get_sdf() <= nTd)
		{
			force = K_SDF*Vertex_memo_map[vi].get_area_factor()*(1.0 - Vertex_memo_map[vi].get_sdf() / nTd);
			setSDFNeighbor(mesh, Vertex_memo_map, vi, force, 3);
		}
	}

	// compute vertex force and Jacobians after setting area factor and movable vertices
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		// index is zero if the vertex is fixed
		if (Vertex_memo_map[vi].isMovable)
			Vertex_memo_map[vi].index = ++movable;
		else
			Vertex_memo_map[vi].index = 0;

		Vertex_memo_map[vi].compute_force(K_SDF, K_S, K_D);
	}

	populate_halfedge_memos(mesh, Vertex_memo_map, Halfedge_memo_map);

	return movable;
}

// return location indicies of 3x3 submatrix starting at (row, col) in 2 row location format
arma::umat getLocationIndices(
	int row,
	int col
)
{
	arma::umat loc;
	loc << row << row + 1 << row + 2 << row << row + 1 << row + 2 << row << row + 1 << row + 2 << arma::endr
		<< col << col << col << col + 1 << col + 1 << col + 1 << col + 2 << col + 2 << col + 2 << arma::endr;
	return loc;
}

// assembly
void global_assembly(
	arma::vec& FORCE,
	arma::sp_mat& JPOS,
	arma::sp_mat& JVEL,
	arma::vec& VELOCITY,
	const size_t& N,
	const Polyhedron& mesh,
	const boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	const boost::associative_property_map<he_memo_map>& Halfedge_memo_map
)
{
	FORCE.set_size(N * 3);
	VELOCITY.set_size(N * 3);

	FORCE.fill(0.0);
	VELOCITY.fill(0.0);

	arma::umat locations;
	arma::vec JPOSvalues, JVELvalues;

	int s_idx, e_idx, count = 0;
	// compute the number of entries in sparse Jacobian matrices and fill FORCE and VELOCITY
	for (Polyhedron::Halfedge_const_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi)
	{
		if (Halfedge_memo_map[hi].isMaster > 0)
		{
			s_idx = ((int)Vertex_memo_map[hi->vertex()].index - 1) * 3;
			e_idx = ((int)Vertex_memo_map[hi->opposite()->vertex()].index - 1) * 3;
			FORCE(arma::span(s_idx, s_idx + 2)) += Halfedge_memo_map[hi].get_force();
			count += 9;

			if (e_idx >= 0) // if the opposite vertex is also movable
			{
				FORCE(arma::span(e_idx, e_idx + 2)) -= Halfedge_memo_map[hi].get_force();
				count += 27;
			}
		}
	}

	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		if (Vertex_memo_map[vi].index > 0)
		{
			s_idx = (Vertex_memo_map[vi].index - 1) * 3;
			FORCE(arma::span(s_idx, s_idx + 2)) += Vertex_memo_map[vi].get_force();
			VELOCITY(arma::span(s_idx, s_idx + 2)) = Vertex_memo_map[vi].get_velocity(); //TODO: this is unnecessary during updation
			count += 9;
		}
	}

	locations.set_size(2, count);
	JPOSvalues.set_size(count);
	JVELvalues.set_size(count);
	
	count = 0;
	// loop over each halfedge
	for (Polyhedron::Halfedge_const_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi)
	{
		if (Halfedge_memo_map[hi].isMaster > 0)
		{
			s_idx = ((int)Vertex_memo_map[hi->vertex()].index - 1) * 3;
			e_idx = ((int)Vertex_memo_map[hi->opposite()->vertex()].index - 1) * 3;

			JPOSvalues(arma::span(count, count + 8)) = arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_pos());
			JVELvalues(arma::span(count, count + 8)) = arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_vel());
			locations(arma::span::all, arma::span(count, count + 8)) = getLocationIndices(s_idx, s_idx);
			count += 9;

			if (e_idx >= 0) // if the opposite vertex is also movable
			{


				JPOSvalues(arma::span(count, count + 8)) = arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_pos());
				JVELvalues(arma::span(count, count + 8)) = arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_vel());
				locations(arma::span::all, arma::span(count, count + 8)) = getLocationIndices(e_idx, e_idx);
				count += 9;

				JPOSvalues(arma::span(count, count + 8)) = -arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_pos());
				JVELvalues(arma::span(count, count + 8)) = -arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_vel());
				locations(arma::span::all, arma::span(count, count + 8)) = getLocationIndices(s_idx, e_idx);
				count += 9;

				JPOSvalues(arma::span(count, count + 8)) = -arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_pos());
				JVELvalues(arma::span(count, count + 8)) = -arma::vectorise(Halfedge_memo_map[hi].get_Jacobian_vel());
				locations(arma::span::all, arma::span(count, count + 8)) = getLocationIndices(e_idx, s_idx);
				count += 9;
			}
		}
	}

	// loop over each vertex
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		if (Vertex_memo_map[vi].index > 0)
		{
			s_idx = (Vertex_memo_map[vi].index - 1) * 3;

			JPOSvalues(arma::span(count, count + 8)) = arma::vectorise(Vertex_memo_map[vi].get_Jacobian_pos());
			JVELvalues(arma::span(count, count + 8)) = arma::vectorise(Vertex_memo_map[vi].get_Jacobian_vel());
			locations(arma::span::all, arma::span(count, count + 8)) = getLocationIndices(s_idx, s_idx);
			count += 9;
		}
	}

	// fill sparse matrices JPOS and JVEL
	JPOS = arma::sp_mat(true, locations, JPOSvalues, N * 3, N * 3);
	JVEL = arma::sp_mat(true, locations, JVELvalues, N * 3, N * 3);
}

// compute update in velocity and position for one time step
void march_one_time_step(
	arma::vec& delta_VEL,
	arma::vec& delta_POS,
	const size_t& N,
	const arma::vec& FORCE,
	const arma::sp_mat& JPOS,
	const arma::sp_mat& JVEL,
	const arma::vec& VELOCITY
)
{
	arma::sp_mat A;
	arma::vec b;

	A = arma::speye<arma::sp_mat>(N * 3, N * 3) - TIME_STEP*MASS_INV*JVEL - TIME_STEP*TIME_STEP*MASS_INV*JPOS;
	b = TIME_STEP*MASS_INV*(FORCE + TIME_STEP*JPOS*VELOCITY);

	delta_VEL = arma::spsolve(A, b); // using SuperLU solver, TODO: use settings to compute faster, such as symmetric
	//delta_VEL = arma::solve(arma::mat(A), b);
	//delta_VEL = arma::inv_sympd(arma::mat(A))*b;
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
	arma::sp_mat& JPOS,
	arma::sp_mat& JVEL,
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
		if (Vertex_memo_map[vi].index > 0)
		{
			idx = (Vertex_memo_map[vi].index - 1) * 3;
			Vertex_memo_map[vi].set_velocity(Vertex_memo_map[vi].get_velocity() + delta_VEL(arma::span(idx, idx + 2)));
		}
		Vertex_memo_map[vi].set_normal(normal_map[vi]);
		Vertex_memo_map[vi].compute_force(K_SDF, K_S, K_D);
	}

	for (Polyhedron::Halfedge_const_iterator hi = mesh.halfedges_begin(); hi != mesh.halfedges_end(); ++hi)
	{
		if (Halfedge_memo_map[hi].isMaster > 0)
			Halfedge_memo_map[hi].compute_force(K_S, K_D, Vertex_memo_map[hi->vertex()].get_velocity(), Vertex_memo_map[hi->opposite()->vertex()].get_velocity());
	}
}

// computes the maxima of change of location of vertices
double max_change(
	arma::vec delta_POS
)
{
	size_t N = delta_POS.n_rows;
	arma::mat delta_XYZ = delta_POS;
	arma::vec dist;
	delta_XYZ = arma::square(delta_XYZ); // square each element
	delta_XYZ.reshape(3, N / 3); // reshape into 3x(N/3) matrix with squared [x y z] coordinate change
	arma::inplace_trans(delta_XYZ); // transpose
	dist = arma::sum(delta_XYZ, 1); // sum elements in each row
	dist = arma::sqrt(dist); // square root of each element
	return dist.max();
}

// march one step and then update
double march_and_update(
	const size_t& N,
	Polyhedron& mesh,
	boost::associative_property_map<v_memo_map>& Vertex_memo_map,
	boost::associative_property_map<he_memo_map>& Halfedge_memo_map,
	arma::vec& FORCE,
	arma::sp_mat& JPOS,
	arma::sp_mat& JVEL,
	arma::vec& VELOCITY
)
{
	arma::vec delta_VEL;
	arma::vec delta_POS;

	global_assembly(FORCE, JPOS, JVEL, VELOCITY, N, mesh, Vertex_memo_map, Halfedge_memo_map);
	march_one_time_step(delta_VEL, delta_POS, N, FORCE, JPOS, JVEL, VELOCITY);
	update(delta_POS, delta_VEL, N, mesh, Vertex_memo_map, Halfedge_memo_map, FORCE, JPOS, JVEL, VELOCITY);

	return max_change(delta_POS);
}

// creates an output file with all the necessary parameters and output values
void generateOutput(
	const std::string filename,
	const Polyhedron& mesh,
	const boost::associative_property_map<v_memo_map>& Vertex_memo_map
)
{
	std::ofstream ofile(filename);
	ofile << "Input File:\n" << FILENAME << "\n";
	ofile << "Output Identifier:\n" << IDENTIFIER << "\n";
	ofile << "Strict Threshold:\n" << Td << "\n";
	ofile << "Force Constant SDF:\n" << K_SDF << "\n";
	ofile << "Force Constant Spring:\n" << K_S << "\n";
	ofile << "Force Constant Damper:\n" << K_D << "\n";
	ofile << "Vertex Mass:\n" << VERTEX_MASS << "\n";
	ofile << "Time Step:\n" << TIME_STEP << "\n";
	ofile << "Total Time:\n" << TOTAL_TIME << "\n";
	ofile << "Steps:\n" << STEPS << "\n";
	ofile << "Min SDF:\n" << min_max_sdf.first << "\n";
	ofile << "Max SDF:\n" << min_max_sdf.second << "\n";
	ofile << "Number of Movable Vertices:\n" << nMove << "\n";
	ofile << "Total Number of Vertices:\n" << mesh.size_of_vertices() << "\n";
	ofile << "SDF Values:\n";
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
		ofile << Vertex_memo_map[vi].get_sdf() << "\n";
	ofile << "Vertex Movable Index:\n";
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
		ofile << Vertex_memo_map[vi].index << "\n";
	ofile << "Normal Vectors:\n";
	for (Polyhedron::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
		ofile << arma::trans(Vertex_memo_map[vi].get_normal());
	ofile << "Max Change:\n";
	ofile << max_change_vec;
	ofile.close();
}

// print usage
void print_usage()
{
	std::cout << std::endl;
	std::cout << "USAGE:" << std::endl;
	std::cout << "\tthicken3d.exe" << std::endl;
	std::cout << "\tthicken3d.exe" << " <flag1> <value1> <flag2> <value2> ..." << std::endl;
	std::cout << std::endl;
	std::cout << "FLAGS:" << std::endl;
	std::cout << "\t-f:\tinput file path FILENAME (default = \"data\\cactus\\cactus.off\")" << std::endl;
	std::cout << "\t-id:\texperiment index IDENTIFIER (default = \"test\")" << std::endl;
	std::cout << "\t-td:\tthreshold Td (default = 0.1)" << std::endl;
	std::cout << "\t-ks:\tspring constant K_S (default = 64.0)" << std::endl;
	std::cout << "\t-kd:\tdamper constant K_D (default = 1024.0)" << std::endl;
	std::cout << "\t-ksdf:\tforce constant multiplier ksdf_multiplier (default = 1.0)" << std::endl;
	std::cout << "\t-vm:\tvertex mass VERTEX_MASS (default = 1.0)" << std::endl;
	std::cout << "\t-ts:\tone time step TIME_STEP (default = 1.0)" << std::endl;
	std::cout << "\t-t:\ttotal time TOTAL_TIME (default = 100.0)" << std::endl;
	std::cout << "\t-nb:\tnumber of batches N_BATCH (default = 4)" << std::endl;
	std::cout << "\t-h, -help:\tprint usage" << std::endl;
	std::cout << std::endl;
	std::cout << "Running using default settings..." << std::endl;
	std::cout << std::endl;
}


int main(int argc, char *argv[])
{
	if (argc < 2)
		print_usage();
	else
	{
		int i = 1;
		while (i < argc)
		{
			std::string flag(argv[i]);
			if (flag == "-f")
				FILENAME = argv[++i];
			if (flag == "-id")
				IDENTIFIER = argv[++i];
			if (flag == "-td")
				Td = std::atof(argv[++i]);
			if (flag == "-ks")
				K_S = std::atof(argv[++i]);
			if (flag == "-kd")
				K_D = std::atof(argv[++i]);
			if (flag == "-ksdf")
				ksdf_multiplier = std::atof(argv[++i]);
			if (flag == "-vm")
				VERTEX_MASS = std::atof(argv[++i]);
			if (flag == "-ts")
				TIME_STEP = std::atof(argv[++i]);
			if (flag == "-t")
				TOTAL_TIME = std::atof(argv[++i]);
			if (flag == "-nb")
				N_BATCH = (size_t)std::atoi(argv[++i]);
			if (flag == "-h" || flag == "-help")
				print_usage();
			++i;
		}
	}

	// redefine dependent parameters
	K_SDF = ksdf_multiplier * K_S * Td;
	MASS_INV = 1.0 / VERTEX_MASS;
	STEPS = (size_t)std::ceil(TOTAL_TIME / TIME_STEP);
	BATCH_SIZE = (size_t)std::ceil((double)STEPS / N_BATCH);
	max_change_vec.set_size(STEPS);

	// output filenames
	std::string meshOutFile = FILENAME.substr(0, FILENAME.size() - 4) + "_" + IDENTIFIER + ".off";
	std::string outputFile = FILENAME.substr(0, FILENAME.size() - 4) + "_" + IDENTIFIER + ".out";

	// seed for random number generator
	re.seed(0);

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

	arma::vec FORCE;
	arma::sp_mat JPOS;
	arma::sp_mat JVEL;
	arma::vec VELOCITY;

	// numerical time integration
	max_change_vec.fill(0.0);
	std::cout << "Performing numerical integration..." << std::endl;
	for (size_t i = 0; i < STEPS; i++)
	{
		if (i % BATCH_SIZE == 0)
		{
			std::cout << "Step: " << i << ", Time: " << i*TIME_STEP << std::endl;
			// populate memos - compute and store normals, sdf, forces, and Jacobians
			std::cout << "Populating memos..." << std::endl;
			nMove = populate_memos(mesh, Vertex_memo_map, Halfedge_memo_map, false); // nMove is the number of movable vertices
			std::cout << "Populating completed!\nNumber of movable vertices = " << nMove << std::endl;
			std::cout << std::endl;
			std::cout << "Fabricability = " << rho << std::endl;
			std::cout << std::endl;

			if (nMove == 0) break;
		}
		max_change_vec[i] = march_and_update(nMove, mesh, Vertex_memo_map, Halfedge_memo_map, FORCE, JPOS, JVEL, VELOCITY);
	}
	std::cout << "Integration completed!" << std::endl;
	std::cout << std::endl;

	// populate memos - compute and store normals, sdf, forces, and Jacobians
	std::cout << "Populating memos..." << std::endl;
	nMove = populate_memos(mesh, Vertex_memo_map, Halfedge_memo_map, false); // nMove is the number of movable vertices
	std::cout << "Populating completed!\nNumber of movable vertices = " << nMove << std::endl;
	std::cout << std::endl;

	// write output mesh into a file
	std::ofstream output(meshOutFile);
	output << mesh;
	output.close();

	// write other output values and parameters into a result file
	generateOutput(outputFile, mesh, Vertex_memo_map);

	std::cout << "Final Fabricability = " << rho << std::endl;
	std::cout << std::endl;

	return EXIT_SUCCESS;
}