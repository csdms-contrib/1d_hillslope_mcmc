//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// MCMCHillslopeDriver.hpp
//
// this drives the 1D hillslope model through an MCMC model
// Simon M. Mudd, University of Edinburgh, July 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <fstream>
#include "OneDImplicitHillslope.hpp"
using namespace std;


#ifndef MCMCHillslopeDriver_H
#define MCMCHillslopeDriver_H

class MCMCHillslopeDriver
{
	public:

		MCMCHillslopeDriver(char* data_fname)		{create(data_fname); }

		double calculate_likelihood();			// calcualtes the likelihood using
												// measured and modelled data
		double run_one_HS_iteration(double t_hat_peak,double U_hat_peak,double U_hat_width);
												// runs a single iteration of the hillslope
												// model, then reports the likelihood of
												// the parameter suite
		void run_metropolis(int n_iterations, char* paramfilename, char* outfilename );
												// this runs the metropolis algorithm along
												// a chain with n_iterations
												// it prints to the file 'filename'

	protected:

		int n_data_elements;
		vector<double> t_star_data;
		vector<double> E_star_data;
		vector<double> E_star_error;
		vector<double> R_star_data;
		vector<double> R_star_error;
		vector<double> E_star_modelled;
		vector<double> R_star_modelled;

		OneDImplicitHillslope Hillslope;

	private:
		void create(char* data_fname);
};

#endif
