#include "vertexmemo.h"

VertexMemo::VertexMemo()
{
	v = NULL;
	normal.set_size(3);
	normal.fill(0.0);
	velocity.set_size(3);
	velocity.fill(0.0);
	sdf = 0.0;
}

VertexMemo::VertexMemo(Polyhedron::Vertex_const_handle vertex)
{
	v = vertex;
	normal.set_size(3);
	normal.fill(0.0);
	velocity.set_size(3);
	velocity.fill(0.0);
	sdf_force.set_size(3);
	sdf_force.fill(0.0);
	sdf = 0.0;
}

VertexMemo::VertexMemo(Polyhedron::Vertex_const_handle vertex, Kernel::Vector_3 vertex_normal, double dia)
{
	v = vertex;
	normal.set_size(3);
	normal << vertex_normal.x() << vertex_normal.y() << vertex_normal.z();
	velocity.set_size(3);
	velocity.fill(0.0);
	sdf_force.set_size(3);
	sdf_force.fill(0.0);
	sdf = dia;
}

void VertexMemo::set_normal(Kernel::Vector_3 vertex_normal)
{
	normal << vertex_normal.x() << vertex_normal.y() << vertex_normal.z();
}

void VertexMemo::set_velocity(const arma::vec& vel)
{
	velocity = vel;
}

void VertexMemo::set_sdf(double dia)
{
	sdf = dia;
}

void VertexMemo::compute_sdf_force(const double& K_sdf, const double& threshold_dia)
{
	double mag = 0.0;
	if (threshold_dia > sdf)
		 mag = K_sdf*CGAL::square(threshold_dia - sdf);

	sdf_force = mag*normal;
}

arma::vec VertexMemo::get_velocity()
{
	return velocity;
}

arma::vec VertexMemo::get_sdf_force()
{
	return sdf_force;
}