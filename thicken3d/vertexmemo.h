#ifndef VERTEXMEMO_H
#define VERTEXMEMO_H

#include "stdafx.h"

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
	VertexMemo(Polyhedron::Vertex_const_handle vertex, Kernel::Vector_3 vertex_normal, double dia);
	void set_normal(Kernel::Vector_3 vertex_normal);
	void set_sdf(double dia);
	void compute_force(const double& K_sdf, const double& threshold_dia, boost::associative_property_map<he_memo_map>& Halfedge_memo_map);
	arma::vec get_velocity();
	arma::vec get_force();
	arma::mat get_Jacobian_pos();
	arma::mat get_Jacobian_vel();
	size_t index;

private:
	Polyhedron::Vertex_const_handle v;
	arma::vec normal;
	arma::vec velocity;
	double sdf;
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif