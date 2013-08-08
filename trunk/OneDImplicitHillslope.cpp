//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// oneD_implicit_nonlinear_hillslope.cpp
//
// an object for solving nonlinear hillslope evolution implicitly
// Simon M. Mudd, University of Edinburgh, July 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "TNT/tnt.h"
#include "TNT/jama_lu.h"
#include "OneDImplicitHillslope.hpp"
using namespace std;
using namespace TNT;
using namespace JAMA;

#ifndef OneDImplicitHillslope_CPP
#define OneDImplicitHillslope_CPP

// operators
OneDImplicitHillslope& OneDImplicitHillslope::operator=(const OneDImplicitHillslope& rhs)
 {
  if (&rhs != this)
   {
    create(get_n_nodes(),get_ridgetop_nod() ,get_t_hat_peak(),get_U_hat_peak(),
    	   get_U_hat_width(),get_zeta_hat(),get_zeta_last_timestep(),get_zeta_intermediate(),
    	   get_f(),get_Coeff_matrix(),get_x_hat(),get_A_hat_denom(),get_B_hat_denom(),
    	   get_A_slope_denom2(),get_B_slope_denom2());
   }
  return *this;
 }


// this creates an object with default values
void OneDImplicitHillslope::create()
{
	double tp_temp = 1.0;
	double Up_temp = 1.0;
	double Uw_temp = 1.0;
	double dx = 0.2;
	create(tp_temp,Up_temp,Uw_temp, dx);
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function initializes the hillslope
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::create(double tp_temp,double Up_temp,double Uw_temp)
{
	t_hat_peak = tp_temp;
	U_hat_peak = Up_temp;
	U_hat_width = Uw_temp;

	double dx = 0.05;
	int N_nodes1 = int(0.5/dx) + 1;
	int N_nodes2 = int(2*0.4/dx);
	int N_nodes3 = int(4*0.1/dx);

	vector<double> x_hat_temp;

	// do irregular spacing, with a point exactly on the divide
	// we have 5 segments
	// from 0 to 0.5, the spacing is 0.015
	double last = 0.0;
	double latest;
	x_hat_temp.push_back(last);
	for(int i = 1; i<N_nodes1; i++)
	{
		latest = last+dx;
		x_hat_temp.push_back(latest);
		last= latest;
	}
	for(int i = N_nodes1; i<(N_nodes1+N_nodes2); i++)
	{
		latest = last+0.5*dx;
		x_hat_temp.push_back(latest);
		last= latest;
	}
	for(int i = N_nodes1+N_nodes2; i<(N_nodes1+N_nodes2+N_nodes3); i++)
	{
		latest = last+0.25*dx;
		x_hat_temp.push_back(latest);
		last= latest;
	}

	int temp_n_nodes = x_hat_temp.size();
	for (int i = temp_n_nodes-2; i>= 0; i--)
	{
		//cout << "Pushing_back: " << x_hat_temp[i] << endl;
		x_hat_temp.push_back((1-x_hat_temp[i])+1);
	}
	ridgetop_node = temp_n_nodes-1;


	temp_n_nodes = x_hat_temp.size();
	n_nodes = temp_n_nodes;
	x_hat = x_hat_temp;
	//cout << "LINE 76; N_nodes: " << n_nodes  << endl;
	//for (int i = 0; i< n_nodes; i++)
	//{
	//	cout << "x_hat["<<i<<"]: " << x_hat[i] << endl;
	//}

	Array1D<double> zeta_temp(n_nodes,0.0);
	zeta_hat = zeta_temp.copy();
	zeta_last_timestep = zeta_temp.copy();
	zeta_intermediate = zeta_temp.copy();

	vector<double> A_hat_denom_t(n_nodes,0.0);
	vector<double> B_hat_denom_t(n_nodes,0.0);
	vector<double> A_slope_denom2_t(n_nodes,0.0);
	vector<double> B_slope_denom2_t(n_nodes,0.0);

	for (int i = 1; i<n_nodes-1; i++)
	{
		A_hat_denom_t[i] = 2/( (x_hat[i+1]-x_hat[i-1])*(x_hat[i+1]-x_hat[i]) );
		B_hat_denom_t[i] = 2/( (x_hat[i+1]-x_hat[i-1])*(x_hat[i]-x_hat[i-1]) );
		A_slope_denom2_t[i] = 1/( (x_hat[i+1]-x_hat[i])*(x_hat[i+1]-x_hat[i]) );
		B_slope_denom2_t[i] = 1/( (x_hat[i]-x_hat[i-1])*(x_hat[i]-x_hat[i-1]) );
	}

	A_hat_denom = A_hat_denom_t;
	B_hat_denom = B_hat_denom_t;
	A_slope_denom2 = A_slope_denom2_t;
	B_slope_denom2 = B_slope_denom2_t;

	// set up the arrays for solving the system of equations
	Array2D<double> Temp_coeff_matx(n_nodes,n_nodes,0.0);
	Temp_coeff_matx[0][0] = 1.0;
	Temp_coeff_matx[n_nodes-1][n_nodes-1] = 1.0;
	Array1D<double> temp_f(n_nodes,0.0);
	Coeff_matrix = Temp_coeff_matx.copy();
	f = temp_f.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function initializes the hillslope
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::create(double tp_temp,double Up_temp,double Uw_temp, double dx)
{
	t_hat_peak = tp_temp;
	U_hat_peak = Up_temp;
	U_hat_width = Uw_temp;

	//cout << "dx: " << dx << " and in 0.7: " << 0.7/dx << endl;

	int N_nodes1 = int(0.8/dx+0.5) + 1;
	int N_nodes2 = int(10*0.18/dx+0.5);
	int N_nodes3 = int(20*0.02/dx+0.5);

	vector<double> x_hat_temp;

	// do irregular spacing, with a point exactly on the divide
	// we have 5 segments
	// from 0 to 0.5, the spacing is 0.015
	double last = 0.0;
	double latest;
	x_hat_temp.push_back(last);
	for(int i = 1; i<N_nodes1; i++)
	{
		latest = last+dx;
		x_hat_temp.push_back(latest);
		last= latest;
	}
	for(int i = N_nodes1; i<(N_nodes1+N_nodes2); i++)
	{
		latest = last+0.1*dx;
		x_hat_temp.push_back(latest);
		last= latest;
	}
	for(int i = N_nodes1+N_nodes2; i<(N_nodes1+N_nodes2+N_nodes3); i++)
	{
		latest = last+0.05*dx;
		x_hat_temp.push_back(latest);
		last= latest;
	}

	int temp_n_nodes = x_hat_temp.size();
	for (int i = temp_n_nodes-2; i>= 0; i--)
	{
	//	//cout << "Pushing_back: " << x_hat_temp[i] << endl;
		x_hat_temp.push_back((1-x_hat_temp[i])+1);
	}
	ridgetop_node = temp_n_nodes-1;


	temp_n_nodes = x_hat_temp.size();
	n_nodes = temp_n_nodes;
	x_hat = x_hat_temp;
	//cout << "LINE 76; N_nodes: " << n_nodes  << endl;
	//for (int i = 0; i< n_nodes; i++)
	//{
	//	cout << "x_hat["<<i<<"]: " << x_hat[i] << endl;
	//}

	Array1D<double> zeta_temp(n_nodes,0.0);
	zeta_hat = zeta_temp.copy();
	zeta_last_timestep = zeta_temp.copy();
	zeta_intermediate = zeta_temp.copy();

	vector<double> A_hat_denom_t(n_nodes,0.0);
	vector<double> B_hat_denom_t(n_nodes,0.0);
	vector<double> A_slope_denom2_t(n_nodes,0.0);
	vector<double> B_slope_denom2_t(n_nodes,0.0);

	for (int i = 1; i<n_nodes-1; i++)
	{
		A_hat_denom_t[i] = 2/( (x_hat[i+1]-x_hat[i-1])*(x_hat[i+1]-x_hat[i]) );
		B_hat_denom_t[i] = 2/( (x_hat[i+1]-x_hat[i-1])*(x_hat[i]-x_hat[i-1]) );
		A_slope_denom2_t[i] = 1/( (x_hat[i+1]-x_hat[i])*(x_hat[i+1]-x_hat[i]) );
		B_slope_denom2_t[i] = 1/( (x_hat[i]-x_hat[i-1])*(x_hat[i]-x_hat[i-1]) );
	}

	A_hat_denom = A_hat_denom_t;
	B_hat_denom = B_hat_denom_t;
	A_slope_denom2 = A_slope_denom2_t;
	B_slope_denom2 = B_slope_denom2_t;

	// set up the arrays for solving the system of equations
	Array2D<double> Temp_coeff_matx(n_nodes,n_nodes,0.0);
	Temp_coeff_matx[0][0] = 1.0;
	Temp_coeff_matx[n_nodes-1][n_nodes-1] = 1.0;
	Array1D<double> temp_f(n_nodes,0.0);
	Coeff_matrix = Temp_coeff_matx.copy();
	f = temp_f.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// this copy constructor is for use with the = operator
void OneDImplicitHillslope::create(int tn_nodes, int tridgetop_node, double tt_hat_peak, double tU_hat_peak,
					double tU_hat_width, Array1D<double> tzeta_hat, Array1D<double> tzeta_last_timestep,
					Array1D<double> tzeta_intermediate, Array1D<double> tf, Array2D<double> tCoeff_matrix,
					vector<double> tx_hat,vector<double> tA_hat_denom, vector<double> tB_hat_denom,
					vector<double> tA_slope_denom2, vector<double> tB_slope_denom2)
{
	n_nodes = n_nodes;
	ridgetop_node = tridgetop_node;
	t_hat_peak = tt_hat_peak;
	U_hat_peak = tU_hat_peak;
	U_hat_width = tU_hat_width;
	zeta_hat = tzeta_hat.copy();
	zeta_last_timestep = tzeta_last_timestep.copy();
	zeta_intermediate = tzeta_intermediate.copy();
	f = tf.copy();
	Coeff_matrix =  tCoeff_matrix.copy();
	x_hat = tx_hat;
	A_hat_denom = tA_hat_denom;
	B_hat_denom = tB_hat_denom;
	A_slope_denom2  = tA_slope_denom2;
	B_slope_denom2 = tB_slope_denom2;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this resets the 1D hillslope so that the coefficient matrices and
// x locations don't have to be recaluclated
// it starts from a flat surface
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::reset_hillslope(double tp_temp, double Up_temp, double Uw_temp)
{
	t_hat_peak = tp_temp;
	U_hat_peak = Up_temp;
	U_hat_width = Uw_temp;

	Array1D<double> zero_vec(n_nodes,0.0);
	zeta_hat = zero_vec.copy();
	zeta_intermediate = zero_vec.copy();
	zeta_last_timestep = zero_vec.copy();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this returns the U_star for a given t_hat based on a gaissian function
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double OneDImplicitHillslope::gaussian_uplift(double t_hat)
{
	double U_hat;
	U_hat = U_hat_peak*exp( - (t_hat-t_hat_peak)*(t_hat-t_hat_peak)/(U_hat_width*U_hat_width) );
	return U_hat;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs the model, using a gaussian uplift field, until it has reached the
// time indicated by t_ime
//
// this replaces the data in the E_star_modelled and R_star_modelled vectors
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::run_based_on_data_spacing
							(vector<double> t_hat_data, vector<double>& E_star_modelled,
							 vector<double>& R_star_modelled, double tolerance)
{

	int n_timesteps = t_hat_data.size();
	double target_time;
	double t_ime = 0;
	double dt_hat = 0.004;					// a default. The timestepping can change
											// so this will optimize after a few iterations

	for (int i = 0; i<n_timesteps; i++)
	{
		target_time = t_hat_data[i];
		advance_to_next_time_interval_gaussian_uplift(dt_hat, t_ime, target_time, tolerance);

		R_star_modelled[i] = calculate_R_star();
		E_star_modelled[i] = calculate_E_star();
		//cout << "LINE 311 time is: " << t_ime << " E*: " << E_star_modelled[i]
		//							 << " R*: " << R_star_modelled[i] << endl;
	}
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs the model, using a gaussian uplift field, until it has reached the
// time indicated by t_ime
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::advance_to_next_time_interval_gaussian_uplift
							(double& dt_hat, double& t_ime, double target_time, double tolerance)
{
	// this is just for getting the final step to advance to exactly the right time
	double last_dt_increment;

	// loop until you get to the timestep before the last timestep
	 while(t_ime < target_time-dt_hat)
	{
		hillslope_timestep_gaussian_uplift(dt_hat, t_ime, tolerance);
		//cout.precision(8);
		//cout << "LINE 332 time: " << t_ime << " dt_hat: " << dt_hat << " diff: " << target_time-dt_hat << endl;
	}

	//cout << endl << "NOW FOR THE LAST BIT, target time: " << target_time << " and time: " << t_ime << endl;
	// now increment up to the target time
	while(t_ime < target_time)
	{
		last_dt_increment = target_time-t_ime;
		//cout << "Looping, target_time: " << target_time << " time: " << t_ime;
		hillslope_timestep_gaussian_uplift(last_dt_increment, t_ime, tolerance);
		//cout << " and time after timestep" << t_ime << endl;
	}
	//cout << "Did the last bit, time: " << t_ime << endl;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function solves a timestep of the hillslope implicitly
// note uplift is considered to be spatially homogenous
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::hillslope_timestep(double& dt_hat, double& t_ime, double U_hat, double tolerance)
{
	int continue_switch = 0;
	int n_iter;
	double dt_hat_new = dt_hat;
	double dt_hat_old = dt_hat;
	//cout << "288 entering while loop" << endl;
	while (continue_switch == 0)
	{
		n_iter = hillslope_iterator(dt_hat_old, U_hat, tolerance);
		//cout << "292 iter: " << n_iter << " and dt hat new: " << dt_hat_new << endl;

		// if there is only one iteration, increase the timestep
		if (n_iter ==1)
		{

			dt_hat_new = dt_hat_old*1.5;
			continue_switch = 1;
			//cout << "SPEEDING UP TOTO!!!, dt is: " << dt_hat << endl;
		}
		else if(n_iter >=5)
		{
			dt_hat_new = dt_hat_new*0.5;
			zeta_hat = zeta_last_timestep.copy();
			zeta_intermediate = zeta_last_timestep.copy();
			dt_hat_old = dt_hat_new;
			//cout << "SLOWING DOWN TOTO!!!, dt is: " << dt_hat_new << endl;
		}
		else
		{
			continue_switch = 1;
			dt_hat_new = dt_hat_old;
		}
	}
	// the timestep advances by dt_hat_old, but the new timestep is dt_hat_new
	dt_hat = dt_hat_new;
	t_ime+= dt_hat_old;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function solves a timestep of the hillslope implicitly
// note uplift is considered to be spatially homogenous
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::hillslope_timestep_gaussian_uplift(double& dt_hat, double& t_ime, double tolerance)
{
	double U_hat;
	int continue_switch = 0;
	int n_iter;
	double dt_hat_new = dt_hat;
	double dt_hat_old = dt_hat;
	//cout << "288 entering while loop" << endl;
	//cout << "line 405, dt_hat: " << dt_hat << " and time: " << t_ime << endl;
	while (continue_switch == 0)
	{
		U_hat = gaussian_uplift(t_ime+dt_hat_old);
		n_iter = hillslope_iterator(dt_hat_old, U_hat, tolerance);
		//cout << "292 iter: " << n_iter << " and dt hat new: " << dt_hat_new << endl;

		// if there is only one iteration, increase the timestep
		if (n_iter ==1)
		{

			dt_hat_new = dt_hat_old*1.5;
			continue_switch = 1;
			//cout << "SPEEDING UP TOTO!!!, dt is: " << dt_hat << endl;
		}
		else if(n_iter >=5)
		{
			dt_hat_new = dt_hat_new*0.5;
			zeta_hat = zeta_last_timestep.copy();
			zeta_intermediate = zeta_last_timestep.copy();
			dt_hat_old = dt_hat_new;
			//cout << "SLOWING DOWN TOTO!!!, dt is: " << dt_hat_new << endl;
		}
		else
		{
			continue_switch = 1;
			dt_hat_new = dt_hat_old;
		}
	}
	// the timestep advances by dt_hat_old, but the new timestep is dt_hat_new
	dt_hat = dt_hat_new;
	t_ime+= dt_hat_old;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this iterates on the solution until a tolerance has been reached.
// there is a maximum number of iterations (5)
// If the model fails to converge within 5 iterations the function exits without updating
//  zeta.  A parent function can then change the timestep and try again
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
int OneDImplicitHillslope::hillslope_iterator(double dt_hat, double U_hat, double tolerance)
{
	// assemble the matrix using the last timestep
	double max_error = 10;
	int iterations = 0;

	// initialize arrays holding the last timestep and
	// an array for the intermediate timestep
	zeta_last_timestep = zeta_hat.copy();
	zeta_intermediate = zeta_hat.copy();

	// get the first intermediate zeta
	solve_for_zeta_intermediate(dt_hat, U_hat);

	//cout << endl << endl << "zeta is: " << zeta_intermediate << endl;
	zeta_hat = zeta_intermediate.copy();

	// now test intermediate zeta against further iterations
	while (max_error > tolerance && iterations <= 5)
	{
		iterations++;
		solve_for_zeta_intermediate(dt_hat, U_hat);
		max_error = get_zeta_rootsquare_error();
		zeta_hat = zeta_intermediate.copy();
		//cout << "iteration is: " << iterations << " and max error is: " << max_error << endl;
	}


	return iterations;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this assembles the coefficient matrix and solves for intermediate zeta
// it replaces the intermediate zeta provided to it
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::solve_for_zeta_intermediate(double dt_hat, double U_hat)
{
	double A_hat, B_hat;
	for(int i = 1; i<n_nodes-1; i++)
	{
		f[i] = dt_hat*0.5*U_hat+zeta_last_timestep[i];
		A_hat = dt_hat*A_hat_denom[i]/( 1 - A_slope_denom2[i]*
		             ( (zeta_intermediate[i+1]-zeta_intermediate[i])
		              *(zeta_intermediate[i+1]-zeta_intermediate[i]) ) );
		B_hat = dt_hat*B_hat_denom[i]/( 1 - B_slope_denom2[i]*
				     ( (zeta_intermediate[i]-zeta_intermediate[i-1])
		              *(zeta_intermediate[i]-zeta_intermediate[i-1]) ) );
		Coeff_matrix[i][i-1] = -B_hat;
		Coeff_matrix[i][i] = 1+A_hat+B_hat;
		Coeff_matrix[i][i+1] = -A_hat;
	}

	// Solve matrix equations using LU decomposition using the TNT JAMA package:
	// Coeff_matrix*zeta_hat = f, where coefs is the coefficients vector.

	//cout << endl << endl << endl << "f: " << f << endl;
	//cout << endl << endl << endl << "Coeff_matrix: " << Coeff_matrix << endl;

	LU<double> sol_CM(Coeff_matrix);  // Create LU object
	zeta_intermediate = sol_CM.solve(f);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets the r_star and E* values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double OneDImplicitHillslope::calculate_R_star()
{
	double R_star = zeta_hat[ridgetop_node];
	return R_star;
}

double OneDImplicitHillslope::analytical_R_star(double U_star)
{
	double x_with_zero_at_divide;
	double root_es_xs;
	double log_term;

	x_with_zero_at_divide = -1;
	root_es_xs = sqrt(1+ U_star*U_star*x_with_zero_at_divide*x_with_zero_at_divide);
	log_term = log( 0.5*(1+root_es_xs) );
	double z_zero =  (1/U_star)*(log_term+1-root_es_xs);


	return -z_zero;

}

double OneDImplicitHillslope::calculate_E_star()
{
	double E_star;
	int i = ridgetop_node;
	E_star = -2*( (zeta_hat[i+1] - 2*zeta_hat[i] + zeta_hat[i-1]) /
	              ( (x_hat[i]-x_hat[i-1])*(x_hat[i]-x_hat[i-1]) ) );
	return E_star;
}
// note there is no analytical e* since this is just U_star

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets the maximum difference between the intermediate zeta and the new zeta
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double OneDImplicitHillslope::get_zeta_rootsquare_error()
{
	double max_error = 0;
	for(int i = 1; i<n_nodes-1; i++)
	{
		if ( (zeta_hat[i]-zeta_intermediate[i])*(zeta_hat[i]-zeta_intermediate[i]) > max_error )
		{
			max_error = (zeta_hat[i]-zeta_intermediate[i])*(zeta_hat[i]-zeta_intermediate[i]);
		}
	}

	return sqrt(max_error);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// get the analaytical solution of C* as a function of x* from Roering et al
// equation 8a. Note that there must be a coordinate transformation since I have
// boundaries at x = 0 and x = 2 and Josh has boundaries at x = -1 and 1.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Array1D<double> OneDImplicitHillslope::get_analytical_SS_Cstar(double U_star)
{
	Array1D<double> analytical_Cstar(n_nodes,0.0);
	double x_with_zero_at_divide;
	double root_es_xs;
	double term1;
	double term2;
	double term3;

	for(int i = 0; i<n_nodes; i++)
	{
		x_with_zero_at_divide = x_hat[i]-1;
		if (x_with_zero_at_divide == 0)
		{
			analytical_Cstar[i] = U_star/2;
		}
		else
		{
			root_es_xs = sqrt(1+ U_star*U_star*x_with_zero_at_divide*x_with_zero_at_divide);
			term1 = (U_star/root_es_xs);
			term2 = (1/(U_star*x_with_zero_at_divide*x_with_zero_at_divide));
			term3 = term2*(root_es_xs-1);
			analytical_Cstar[i] = term1-term3;
		}
		//cout << "x["<<i<<"]: " << x_hat[i] << " \tC: " << analytical_Cstar[i]
		//     << "  \tterm1: " << term1 << " \tterm2: " << term2 << " \tterm3: " << term3 << endl;
	}

	return analytical_Cstar;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets the anayltical S* according to Roering et al
// equation 8b. Note that there must be a coordinate transformation since I have
// boundaries at x = 0 and x = 2 and Josh has boundaries at x = -1 and 1.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Array1D<double> OneDImplicitHillslope::get_analytical_SS_Sstar(double U_star)
{
	Array1D<double> analytical_Sstar(n_nodes,0.0);
	double x_with_zero_at_divide;
	double root_es_xs;

	for(int i = 0; i<n_nodes; i++)
	{
		x_with_zero_at_divide = x_hat[i]-1;
		root_es_xs = sqrt(1+ U_star*U_star*x_with_zero_at_divide*x_with_zero_at_divide);
		analytical_Sstar[i] = (1-root_es_xs)*(1/(U_star*x_with_zero_at_divide));
	}

	return analytical_Sstar;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets the anayltical z* according to Roering et al
// equation 8c. Note that there must be a coordinate transformation since I have
// boundaries at x = 0 and x = 2 and Josh has boundaries at x = -1 and 1.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Array1D<double> OneDImplicitHillslope::get_analytical_SS_zstar(double U_star)
{
	Array1D<double> analytical_zstar(n_nodes,0.0);
	double x_with_zero_at_divide;
	double root_es_xs;
	double log_term;

	x_with_zero_at_divide = -1;
	root_es_xs = sqrt(1+ U_star*U_star*x_with_zero_at_divide*x_with_zero_at_divide);
	log_term = log( 0.5*(1+root_es_xs) );
	double z_zero =  (1/U_star)*(log_term+1-root_es_xs);

	for(int i = 1; i<n_nodes-1; i++)
	{
		//cout << "zstar i: " << i << " x_hat: " << x_hat[i] << endl;
		x_with_zero_at_divide = x_hat[i]-1;
		root_es_xs = sqrt(1+ U_star*U_star*x_with_zero_at_divide*x_with_zero_at_divide);
		log_term = log( 0.5*(1+root_es_xs) );
		analytical_zstar[i] = (1/U_star)*(log_term+1-root_es_xs)-z_zero;
	}

	return analytical_zstar;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function prints the analytical solutions to screen
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void OneDImplicitHillslope::print_analytical_SS(double U_star)
{
	Array1D<double> analytical_zstar = get_analytical_SS_zstar(U_star);
	Array1D<double> analytical_Sstar = get_analytical_SS_Sstar(U_star);
	Array1D<double> analytical_Cstar = get_analytical_SS_Cstar(U_star);

	cout << "analytical solutions for U_star = " << U_star << endl;
	for(int i = 0; i<n_nodes; i++)
	{
		cout << "x+hat["<<i<<"]: " << x_hat[i] << " \t\t zeta_hat: " << analytical_zstar[i]
		     << " \t S_hat: " << analytical_Sstar[i] << " \t C_hat: " << analytical_Cstar[i] << endl;
	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function is simply used to test the steady solution
// it returns the root square of the maximum error between zeta_star modelled
// vs zeta_star analytical
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double OneDImplicitHillslope::test_steady(double dt_hat, double tolerance)
{

	double end_time = 4;
	double t_ime = 0;
	int iter = 0;
	double U_hat = U_hat_peak;

	// run the model to steady. Time is scaled by the relaxation time
	// so an end time of 10 should be plenty to get to steady state.
	//advance_to_next_time_interval(dt_hat, t_ime, end_time,tolerance);
	while (t_ime < end_time)
	{
		iter++;
		hillslope_timestep(dt_hat, t_ime, U_hat, tolerance);

		if(iter%1000 == 0)
		{
			cout << "LINE 577 t_ime is: " << t_ime << " and dt is: " << dt_hat << endl;
		}
	}
	cout << "Line 580 ended at time: " << t_ime << endl;


	// now compare the computational SS solution to the anayltical SS solution
	Array1D<double> analytical_zstar = get_analytical_SS_zstar(U_hat_peak);
	double max_error = 0;
	for(int i = 1; i<n_nodes-1; i++)
	{
		if ( (zeta_hat[i]-analytical_zstar[i])*(zeta_hat[i]-analytical_zstar[i]) > max_error )
		{
			//cout << "err["<<i<<"]: " << sqrt( (zeta_hat[i]-analytical_zstar[i])*
			//                                  (zeta_hat[i]-analytical_zstar[i]) ) << endl;
			max_error = (zeta_hat[i]-analytical_zstar[i])*(zeta_hat[i]-analytical_zstar[i]);
		}
	}

	for(int i = 1; i<n_nodes-1; i++)
	{
		cout << "x["<<i<<"]: " << x_hat[i] << " \t analytical zeta_hat: " << analytical_zstar[i]
		     << " \t zeta_hat: " << zeta_hat[i] << endl;
	}


	return sqrt(max_error);
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#endif

