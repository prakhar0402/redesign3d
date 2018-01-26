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


#ifndef HALFEDGEMEMO_H
#define HALFEDGEMEMO_H

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
	int isMaster; // 1 if master, -1 if slave, 0 if not known

private:
	Polyhedron::Halfedge_const_handle he;
	double initial_length;
	double current_length;
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif