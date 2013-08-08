//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// DBRP_implicit_nonD_evolution.cpp
// this evolves a 1D hillslope in order to test if
// the evolution along the E* R* curve predicted by Roering et al 2007
// is consitent with field data
//
// The model just needs the name of the output filename.
// Simon M. Mudd, University of Edinburgh, July 2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "OneDImplicitHillslope.hpp"
#include "MCMCHillslopeDriver.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: not enough inputs. The program needs 1) number of iterations "
		     << " 2) the paramfilename filename and 3) the chainfile " << endl;
		exit(EXIT_SUCCESS);
	}

	// the name of the data
	char* data_name = "fixed_esRs.data";

	// the name of the chainfile
	char* chain_fname = argv[3];
	char* param_fname = argv[2];
	// the number of links in the chain
	int n_iterations = atoi(argv[1]);
	cout << "chain filename is: " << chain_fname << " and chain will have "
	     << n_iterations << " links, params are in: " << param_fname << endl;

	// load an MCMC driver object
	MCMCHillslopeDriver MCMC_driver(data_name);

	//now run the metropolis algorithm along a chain
	MCMC_driver.run_metropolis(n_iterations, param_fname, chain_fname);
}

