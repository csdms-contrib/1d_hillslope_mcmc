//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// MCMCHillslopeDriver.cpp
//
// this drives the 1D hillslope model through an MCMC model
// Simon M. Mudd, University of Edinburgh, July 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include "normal.hpp"
#include "OneDImplicitHillslope.hpp"
#include "MCMCHillslopeDriver.hpp"
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// creates an MCMC object
// initiates the hillslope with default values
// dx = 1 (this yields 30 nodes and has errors of
// about 0.1% in E_star and R_star compared to
// analytical solutions: it is a good compromise
// between speed and accuracy
//
// t_hat_peak = 1
// U_hat_peak = 1
// U_hat_width = 1
//
// these choices shouldn't matter since we'll burn in through
// many 10k iterations along the chain
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void MCMCHillslopeDriver::create(char* data_fname)
{
	ifstream data_in;
	data_in.open(data_fname);

	double t_star_temp,E_star_temp,E_star_err_temp,R_star_temp,R_star_err_temp;
	vector<double> ts_temp;
	vector<double> Rs_temp;
	vector<double> Rse_temp;
	vector<double> Es_temp;
	vector<double> Ese_temp;

	while (data_in >> t_star_temp >> E_star_temp >> E_star_err_temp >> R_star_temp  >> R_star_err_temp)
	{
		ts_temp.push_back(t_star_temp);
		Rs_temp.push_back(R_star_temp);
		Rse_temp.push_back(R_star_err_temp);
		Es_temp.push_back(E_star_temp);
		Ese_temp.push_back(E_star_err_temp);
	}
	t_star_data = ts_temp;
	R_star_data = Rs_temp;
	R_star_error = Rse_temp;
	E_star_data = Es_temp;
	E_star_error = Ese_temp;

	E_star_modelled = E_star_data;
	R_star_modelled = R_star_data;

	n_data_elements = t_star_data.size();

	//double t_hat_peak = 1.0;
	//double U_hat_width = 1.0;
	//double U_hat_peak = 1.0;
	//double dx = 0.2;
	//OneDImplicitHillslope OneDHS(t_hat_peak,U_hat_peak,U_hat_width,dx);
	//Hillslope = OneDHS;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this calcualtes the combined likelihood using chi squared liklihoods
// of the E_star and R_star data
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double MCMCHillslopeDriver::calculate_likelihood()
{
	double sumEsSqExp = 1;
	double sumRsSqExp = 1;
	double Like;

	for (int i = 0; i<n_data_elements; i++)
	{
		//cout << "Es_data: " << E_star_data[i] << " Es mode: " << E_star_modelled[i]
		//     << " Rs data: " << R_star_data[i] << " Rs mod: " << R_star_modelled[i] << endl;
		sumEsSqExp = sumEsSqExp*exp(-0.5*(E_star_data[i]-E_star_modelled[i])*
							(E_star_data[i]-E_star_modelled[i])/E_star_error[i] );
		sumRsSqExp = sumRsSqExp*exp(-0.5*(R_star_data[i]-R_star_modelled[i])*
							(R_star_data[i]-R_star_modelled[i])/R_star_error[i] );
	}
	Like = sumEsSqExp*sumRsSqExp;
	return Like;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this runs a single iteration of the hillslope model for a given set of parameter values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double MCMCHillslopeDriver::run_one_HS_iteration(double t_hat_peak,double U_hat_peak,double U_hat_width)
{
	double Likelihood;			// the likelihood from the run
	double tolerance = 0.000001;

	// reset the hillslope
	Hillslope.reset_hillslope(t_hat_peak, U_hat_peak, U_hat_width);

	// now run the hillslope and retrieve the E_star_modelled and R_star_modelled vector
	Hillslope.run_based_on_data_spacing(t_star_data, E_star_modelled,R_star_modelled, tolerance);

	// calcualte the likelihood
	Likelihood = calculate_likelihood();
	//cout << "LINE 110 MCMC the likelihood is: " << Likelihood << endl;

	// return the likelihood
	return Likelihood;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this runs a single iteration of the hillslope model for a given set of parameter values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void MCMCHillslopeDriver::run_metropolis(int n_iterations, char* paramfilename, char* outfilename )
{
	// open the data file
	ofstream chain_out;
	chain_out.open(outfilename);
	chain_out << "parameter file is: " << paramfilename << endl;

	double test_likelihood;				// the last accepted likelihood
	double new_likelihood;				// the new likelihood
	double like_ratio;					// the ratio between these likelihoods
	int n_accepted = 0;					// the number of accepted parameters
	int n_rejected = 0;					// the number of rejected parameters
	double mean_change = 0.0;			// variables change based on a gaussian distribution
										// with a mean of 0 (that is, there is equal probability
										// of their value increasing or decreasing

	// parameters
	double thp_sd, Uhp_sd, Uhw_sd;		// standard deviations of the parameters
	ifstream param_in;
	param_in.open(paramfilename);
	string temp;
	param_in >> temp >> thp_sd >> temp >> Uhp_sd >> temp >> Uhw_sd;

	// set minimum and maximum values for parameters
	double thp_max,thp_min,Uhp_max,Uhp_min,Uhw_max ,Uhw_min;
	param_in >> temp >> thp_max>> temp >>thp_min>> temp >>Uhp_max>> temp
	         >>Uhp_min>> temp >>Uhw_max >> temp >>Uhw_min;

	double thp_init,Uh_init,Uhw_init;
	param_in >> temp >> thp_init >> temp >>Uh_init>> temp >> Uhw_init;

	param_in.close();

	temp = " ";
	cout << "parameters: " << endl;
	cout << temp << thp_sd << temp << Uhp_sd << temp << Uhw_sd << endl;
	cout << temp << thp_max<< temp <<thp_min<< temp <<Uhp_max<< temp
	         <<Uhp_min<< temp <<Uhw_max << temp <<Uhw_min << endl;

	int seed = time(NULL);       		// seed for random number generator
	double acceptance_probability;		// the new iteration is accepted if the ratio exceeds the
										// acceptance probability (which is determined by
										// a uniform distribution between zero and 1)
	// start the chain with a guess
	// this guess is a very coarse approximation of
	// what the 'real' values might be. The Metropolis algorithm
	// will sample around this
	double t_hat_peak_new = thp_init;
	double t_hat_peak_old = t_hat_peak_new;
	double U_hat_peak_new = Uh_init;
	double U_hat_peak_old = U_hat_peak_new;
	double U_hat_width_old = Uhw_init;
	double U_hat_width_new = U_hat_width_new;
	double dthp, dUhp, dUhw;				// change in parameter values
	double reflect;							// for enforcing minimum and maximum parameter values

	// run the model once to get the test_likelihood
	test_likelihood = run_one_HS_iteration(t_hat_peak_old, U_hat_peak_old,U_hat_width_old);

	// now do metropolis algorithm
	for (int i = 0; i<n_iterations; i++)
	{
		// vary variables
		// t_hat_peak
		dthp = r8_normal(mean_change,thp_sd,seed);
		t_hat_peak_new = t_hat_peak_old + dthp;
        if ( t_hat_peak_new < thp_min)
        {
        	reflect = thp_min -  t_hat_peak_new;
        	t_hat_peak_new = reflect+thp_min;
		}
        if ( t_hat_peak_new > thp_max)
        {
        	reflect = t_hat_peak_new - thp_max;
        	t_hat_peak_new = thp_max -  reflect;
		}
		// U_hat_peak
		dUhp = r8_normal(mean_change,Uhp_sd,seed);
		U_hat_peak_new = U_hat_peak_old + dUhp;
        if ( U_hat_peak_new < Uhp_min)
        {
        	reflect = Uhp_min -  U_hat_peak_new;
        	U_hat_peak_new = reflect+Uhp_min;
		}
        if ( U_hat_peak_new > Uhp_max)
        {
        	reflect = U_hat_peak_new - Uhp_max;
        	U_hat_peak_new = Uhp_max -  reflect;
		}
		// U_hat_width
		dUhw = r8_normal(mean_change,Uhw_sd,seed);
		U_hat_width_new = U_hat_width_old + dUhw;
        if ( U_hat_width_new < Uhw_min)
        {
        	reflect = Uhw_min -  U_hat_width_new;
        	U_hat_width_new = reflect+Uhw_min;
		}
        if ( U_hat_width_new > Uhw_max)
        {
        	reflect = U_hat_width_new - Uhw_max;
        	U_hat_width_new = Uhw_max -  reflect;
		}

		// run the model with the new parameters
		new_likelihood = run_one_HS_iteration(t_hat_peak_new, U_hat_peak_new,U_hat_width_new);

		// get the likelihood ratio
		like_ratio = new_likelihood/test_likelihood;

		// get the acceptance probability (this is set up so that occasional
		// guesses that are worse than the lst one get accepted so that
		// the chain can visit all of parameter space)
		acceptance_probability = r8_uniform_01(seed);

		// if accepted
    	if (like_ratio > acceptance_probability)
    	{
			test_likelihood = new_likelihood;
			n_accepted++;
			t_hat_peak_old = t_hat_peak_new;
			U_hat_peak_old = U_hat_peak_new;
			U_hat_width_old = U_hat_width_new;
		}
		else
		{
			n_rejected++;
		}

		chain_out << i << " " << t_hat_peak_new << " " << U_hat_peak_new << " " << U_hat_width_new << " "
		          << t_hat_peak_old << " " << U_hat_peak_old << " " << U_hat_width_old << " "
		          << new_likelihood << " " << test_likelihood << " " << n_accepted << " " << n_rejected << endl;
	}

	chain_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
