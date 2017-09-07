#include "vertexmemo.h"

VertexMemo::VertexMemo()
{
	v = NULL;
	normal.set_size(3);
	normal.fill(0.0);
	velocity.set_size(3);
	velocity.fill(0.0);
	sdf = 0.0;
	area_factor = 1.0;
	ref_point = Kernel::Point_3(0.0, 0.0, 0.0);
	initial_length = 0.0;
	current_length = 0.0;
	force.set_size(3);
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	force.fill(0.0);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
	isMovable = false;
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
	area_factor = 1.0;
	ref_point = Kernel::Point_3(0.0, 0.0, 0.0);
	initial_length = 0.0;
	current_length = 0.0;
	force.set_size(3);
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	force.fill(0.0);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
	isMovable = false;
}

VertexMemo::VertexMemo(Polyhedron::Vertex_const_handle vertex, Kernel::Vector_3 vertex_normal, double dia = 0.0, double af_value = 1.0)
{
	v = vertex;
	normal.set_size(3);
	normal << vertex_normal.x() << vertex_normal.y() << vertex_normal.z();
	velocity.set_size(3);
	velocity.fill(0.0);
	sdf_force.set_size(3);
	sdf_force.fill(0.0);
	sdf = dia;
	area_factor = af_value;
	ref_point = Kernel::Point_3(0.0, 0.0, 0.0);
	initial_length = 0.0;
	current_length = 0.0;
	force.set_size(3);
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	force.fill(0.0);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
	isMovable = false;
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

void VertexMemo::set_area_factor(double value)
{
	area_factor = value;
}

void VertexMemo::set_ref_point(Kernel::Point_3 ref)
{
	ref_point = ref;
	initial_length = compute_length();
}

double VertexMemo::compute_length()
{
	if (v != NULL)
		current_length = CGAL::sqrt(CGAL::squared_distance(v->point(), ref_point));

	return current_length;
}

void VertexMemo::compute_sdf_force(const double& K_sdf, const double& threshold_dia)
{
	double mag = 0.0;
	if (threshold_dia > sdf)
		mag = K_sdf*area_factor*(1.0 - sdf / threshold_dia);
		//mag = -K_sdf*area*std::log(sdf / threshold_dia);
		//mag = K_sdf*area*(threshold_dia - sdf);
		//mag = K_sdf*area*CGAL::square(threshold_dia - sdf);

	sdf_force = mag*normal;
}

void VertexMemo::compute_force(const double& K_sdf, const double& K_s, const double& K_d, const double& threshold_dia)
{
	if (threshold_dia > sdf)
	{
		compute_length();
		compute_sdf_force(K_sdf, threshold_dia);

		arma::vec pos1, pos2, vec12;
		pos1 << v->point().x() << v->point().y() << v->point().z();
		pos2 << ref_point.x() << ref_point.y() << ref_point.z();
		vec12 = pos2 - pos1; // vector pointing from pos1 to pos2

		arma::mat mv12 = vec12 * arma::trans(vec12); // product of different elements of vec12

		double coef1 = 2.0*K_s*(1 - initial_length / current_length);
		double coef2 = 2.0*K_s*initial_length / pow(current_length, 3);

		force = coef1*vec12 + K_d*(-velocity) + sdf_force;

		Jpos = -coef2*mv12;
		Jpos.diag() -= coef1;

		Jvel.diag() -= K_d;
	}
}

double VertexMemo::get_sdf()
{
	return sdf;
}

double VertexMemo::get_area_factor()
{
	return area_factor;
}

arma::vec VertexMemo::get_normal()
{
	return normal;
}

arma::vec VertexMemo::get_velocity()
{
	return velocity;
}

arma::vec VertexMemo::get_sdf_force()
{
	return sdf_force;
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