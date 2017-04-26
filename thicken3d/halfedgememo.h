#ifndef HALFEDGEMEMO_H
#define HALFEDGEMEMO_H

#include "stdafx.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/squared_distance_3.h>

#include <armadillo>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

class HalfedgeMemo
{
public:
	HalfedgeMemo();
	HalfedgeMemo(Polyhedron::Halfedge_const_handle halfedge);
	double compute_length();
	double get_current_length();
	void compute_force(const double& K_s, const double& K_d, const arma::vec& vel1, const arma::vec& vel2);
	arma::vec get_force();
	arma::mat get_Jacobian_pos();
	arma::mat get_Jacobian_vel();

private:
	Polyhedron::Halfedge_const_handle he;
	double initial_length;
	double current_length;
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif