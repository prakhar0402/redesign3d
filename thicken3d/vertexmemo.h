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


#ifndef VERTEXMEMO_H
#define VERTEXMEMO_H

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
	//void compute_sdf_force(const double& K_sdf, const double& threshold_dia);
	double get_sdf();
	double get_area_factor();
	arma::vec get_normal();
	arma::vec get_velocity();
	double get_sdf_force();
	size_t index;
	bool isMovable;
	void compute_force(const double& K_sdf, const double& K_s, const double& K_d);
	arma::vec get_force();
	arma::mat get_Jacobian_pos();
	arma::mat get_Jacobian_vel();
	void setSDFForceMag(double f);
	void setMaxSDFMag(double f);

private:
	Polyhedron::Vertex_const_handle v;
	arma::vec normal; // vertex normal
	arma::vec velocity;
	Kernel::Point_3 ref_point;
	double initial_length;
	double current_length;
	double sdf;
	double sdf_force;
	double area_factor; // sum of projected area of all faces around the vertex divided by average projected area for all vertices
	arma::vec force;
	arma::mat Jpos;
	arma::mat Jvel;
};

#endif