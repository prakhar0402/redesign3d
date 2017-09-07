#ifndef VERTEXMEMO_H
#define VERTEXMEMO_H

//#include "stdafx.h"

#include "halfedgememo.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/squared_distance_3.h>

#include <armadillo>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef std::map<Polyhedron::Halfedge_const_handle, HalfedgeMemo> he_memo_map;

class VertexMemo
{
public:
	VertexMemo();
	VertexMemo(Polyhedron::Vertex_const_handle vertex);
	VertexMemo(Polyhedron::Vertex_const_handle vertex, Kernel::Vector_3 vertex_normal, double dia, double area_value);
	void set_normal(Kernel::Vector_3 vertex_normal);
	void set_velocity(const arma::vec& vel);
	void set_ref_point(Kernel::Point_3 ref);
	void set_sdf(double dia);
	void set_area_factor(double value);
	double compute_length();
	void compute_sdf_force(const double& K_sdf, const double& threshold_dia);
	double get_sdf();
	double get_area_factor();
	arma::vec get_normal();
	arma::vec get_velocity();
	arma::vec get_sdf_force();
	size_t index;
	bool isMovable;
	void compute_force(const double& K_sdf, const double& K_s, const double& K_d, const double& threshold_dia);
	arma::vec get_force();
	arma::mat get_Jacobian_pos();
	arma::mat get_Jacobian_vel();

private:
	Polyhedron::Vertex_const_handle v;
	arma::vec normal; // vertex normal
	arma::vec velocity;
	Kernel::Point_3 ref_point;
	double initial_length;
	double current_length;
	double sdf;
	double area_factor; // sum of projected area of all faces around the vertex divided by average projected area for all vertices
	arma::vec sdf_force;
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif