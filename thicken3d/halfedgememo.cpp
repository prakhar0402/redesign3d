#include "halfedgememo.h"

HalfedgeMemo::HalfedgeMemo()
{
	he = NULL;
	initial_length = 0.0;
	current_length = initial_length;
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	force.fill(0.0);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
}

HalfedgeMemo::HalfedgeMemo(Polyhedron::Halfedge_const_handle halfedge)
{
	he = halfedge;
	initial_length = compute_length();
	current_length = initial_length;
	Jpos.set_size(3, 3);
	Jvel.set_size(3, 3);
	force.fill(0.0);
	Jpos.fill(0.0);
	Jvel.fill(0.0);
}

double HalfedgeMemo::compute_length()
{
	double dist = 0.0;
	if (he != NULL)
	{
		dist = CGAL::sqrt(CGAL::squared_distance(he->vertex()->point(), he->opposite()->vertex()->point()));
	}
	return dist;
}

double HalfedgeMemo::get_current_length()
{
	return current_length;
}

void HalfedgeMemo::compute_force(const double& K_s, const double& K_d, const arma::vec& vel1, const arma::vec& vel2)
{
	Polyhedron::Vertex_const_handle v1, v2;
	v1 = he->vertex();
	v2 = he->opposite()->vertex();

	arma::vec pos1, pos2, vec12;
	pos1 << v1->point().x() << v1->point().y() << v1->point().z();
	pos2 << v2->point().x() << v2->point().y() << v2->point().z();
	vec12 = pos2 - pos1; // vector pointing from pos1 to pos2

	arma::mat mv12 = vec12 * arma::trans(vec12); // product of different elements of vec12

	double coef1 = K_s*(1 - initial_length / current_length);
	double coef2 = K_s*initial_length / pow(current_length, 3);

	force = coef1*vec12 + K_d*(vel2 - vel1);

	Jpos = -coef2*mv12;
	Jpos.diag() -= coef1;

	Jvel.diag() -= K_d;

	//Jpos(1, 1) = -K_s + coef*(mv12(2, 2) + mv12(3, 3));
	//Jpos(1, 2) = -coef*v12(1)*v12(2);
	//Jpos(1, 3) = -coef*v12(1)*v12(3);

	//Jpos(2, 1) = -K_s + coef*(mv12(2, 2) + mv12(3, 3));
	//Jpos(2, 2) = -coef*v12(1)*v12(2);
	//Jpos(2, 3) = -coef*v12(1)*v12(3);

	//Jpos(3, 1) = -K_s + coef*(mv12(2, 2) + mv12(3, 3));
	//Jpos(3, 2) = -coef*v12(1)*v12(2);
	//Jpos(3, 3) = -coef*v12(1)*v12(3);
	
}

arma::vec HalfedgeMemo::get_force()
{
	return force;
}

arma::mat HalfedgeMemo::get_Jacobian_pos()
{
	return Jpos;
}

arma::mat HalfedgeMemo::get_Jacobian_vel()
{
	return Jvel;
}