#include "vertexmemo.h"

VertexMemo::VertexMemo()
{
	v = NULL;
	normal.set_size(3);
	normal.fill(0.0);
	velocity.set_size(3);
	velocity.fill(0.0);
	force.set_size(3);
	force.fill(0.0);
	sdf = 0.0;
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
}

VertexMemo::VertexMemo(Polyhedron::Vertex_const_handle vertex)
{
	v = vertex;
	normal.set_size(3);
	normal.fill(0.0);
	velocity.set_size(3);
	velocity.fill(0.0);
	force.set_size(3);
	force.fill(0.0);
	sdf = 0.0;
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
}

VertexMemo::VertexMemo(Polyhedron::Vertex_const_handle vertex, Kernel::Vector_3 vertex_normal, double dia)
{
	v = vertex;
	normal.set_size(3);
	normal << vertex_normal.x() << vertex_normal.y() << vertex_normal.z();
	velocity.set_size(3);
	velocity.fill(0.0);
	force.set_size(3);
	force.fill(0.0);
	sdf = dia;
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
}

void VertexMemo::set_normal(Kernel::Vector_3 vertex_normal)
{
	normal << vertex_normal.x() << vertex_normal.y() << vertex_normal.z();
}

void VertexMemo::set_sdf(double dia)
{
	sdf = dia;
}

void VertexMemo::compute_force(const double& K_sdf, const double& threshold_dia, boost::associative_property_map<he_memo_map>& Halfedge_memo_map)
{
	double mag = 0.0;
	if (threshold_dia > sdf)
	{
		 mag = K_sdf*CGAL::square(threshold_dia - sdf);
	}

	force = mag*normal;

	Polyhedron::Halfedge_const_handle he1 = v->halfedge();
	Polyhedron::Halfedge_const_handle he2 = he1;

	do
	{
		force += Halfedge_memo_map[he2].get_force();
		Jpos += Halfedge_memo_map[he2].get_Jacobian_pos();
		Jvel += Halfedge_memo_map[he2].get_Jacobian_vel();
		he2 = he2->next_on_vertex();
	} while (he2 != he1);
}

arma::vec VertexMemo::get_velocity()
{
	return velocity;
}

arma::vec VertexMemo::get_force()
{
	return force;
}

arma::mat VertexMemo::get_Jacobian_pos()
{
	return Jpos;
}

arma::mat VertexMemo::get_Jacobian_vel()
{
	return Jvel;
}