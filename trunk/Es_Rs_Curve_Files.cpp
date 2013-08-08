// this program takes a file that contains the mean and
// standard deviation of the three hillslope parameters
// t_star_peak
// U_star_peak
// U_star_width
// and prints out results from the mean, minimum
// and maximum versions of the model runs
// run from command line with something like
// EsRs_plotting.exe fitted_U_star_params.param DBPR_Estar_Rstar_data_SD fit_EsRs_trajectories
//
// Simon M. Mudd, University of Edinburgh, July 2012

#include <iostream>
#include <fstream>
#include <string>
#include "OneDImplicitHillslope.hpp"
using namespace std;

string itoa(int num)
{
    stringstream converter;
    converter << num;
    return converter.str();
}

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=4)
	{
		cout << "FATAL ERROR: not enough inputs. The program needs 1)fitted param filename" << endl;
		cout << "2) the fname of the data and 3) the prefix of the outfile" << endl;
		exit(EXIT_SUCCESS);
	}

	char* fitted_param_fname = argv[1];
	cout << "fitted parameter filename is: " << fitted_param_fname << endl;

	char* data_name = argv[2];
	cout << "data filename is: " << data_name << endl;

	char* outfile_prefix = argv[3];
	cout << "out file prefix is: " << outfile_prefix << endl;

	string of_prefix = outfile_prefix;
	string dot = ".";
	string num;
	string ext = "txt";
	string mean_str = "mean";

	// open the fitted_parameter fname
	ifstream fitted_param_in;
	fitted_param_in.open(fitted_param_fname);

	double tsp_mean;						// mean t_star_peak
	double tsp_sd_min;
	double tsp_sd_plus;			// standard deviation of t star peak
	double Usp_mean;						// mean U star peak
	double Usp_sd_min;
	double Usp_sd_plus;			// standard deviation of U star peak
	double Usw_mean;						// mean U star width
	double Usw_sd_min;
	double Usw_sd_plus;			// standard deviation of U star width
	//fitted_param_in >> tsp_mean >> tsp_sd_min >> tsp_sd_plus
	//                >> Usp_mean >> Usp_sd_min >> Usp_sd_plus
	//                >> Usw_mean >> Usw_sd_min >> Usw_sd_plus;
	fitted_param_in >> tsp_sd_min >> Usp_sd_min >> Usw_sd_min
	                >> tsp_mean >> Usp_mean >> Usw_mean
	                >> tsp_sd_plus >> Usp_sd_plus >> Usw_sd_plus;
	fitted_param_in.close();

	//cout << "fitted parameters are: " << endl;
	//cout << "t_s_mean: " << tsp_mean << " and sigma: " << tsp_sd << endl;
	//cout << "Usp_mean: " << Usp_mean << " and sigma: " << Usp_sd << endl;
	//cout << "Usw_mean: " << Usw_mean << " and sigma: " << Usw_sd << endl;

	// create a hillslope
	OneDImplicitHillslope Hillslope;

	// now load the data from the field
	ifstream data_in;
	data_in.open(data_name);

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
	cout << "loaded data " << endl;

	int n_data_elements = ts_temp.size();
	vector<double> mean_E_star(n_data_elements, 0.0);
	vector<double> max_E_star(n_data_elements, 0.0);
	vector<double> min_E_star(n_data_elements, 100.0);
	vector<double> mean_R_star(n_data_elements, 0.0);
	vector<double> max_R_star(n_data_elements, 0.0);
	vector<double> min_R_star(n_data_elements, 100.0);
	vector<double> mean_U_star(n_data_elements, 0.0);
	vector<double> max_U_star(n_data_elements, 0.0);
	vector<double> min_U_star(n_data_elements, 100.0);

	// run the mean values
	vector<double> E_star_modelled(n_data_elements, 0.0);
	vector<double> R_star_modelled(n_data_elements, 0.0);
	vector<double> U_star_modelled(n_data_elements, 0.0);
	double tolerance = 0.000001;

	Hillslope.reset_hillslope(tsp_mean, Usp_mean, Usw_mean);
	Hillslope.run_based_on_data_spacing(ts_temp, mean_E_star,mean_R_star, tolerance);

	ofstream data_out;
	string mean_fname = of_prefix+dot+mean_str+dot+ext;
	data_out.open(mean_fname.c_str());
	data_out << tsp_mean << " " << Usp_mean << " " << Usw_mean << " -9999" << endl;
	// now get the uplift
	for (int i = 0; i<n_data_elements; i++)
	{
		mean_U_star[i] = Hillslope.gaussian_uplift( ts_temp[i] );
		data_out << ts_temp[i] << " " << mean_E_star[i] << " " << mean_R_star[i] << " " << mean_U_star[i] << endl;
	}
	data_out.close();



	// now loop through the plus and minus standard deviations, collecting mean and min
	// values
	double tsp_this, Usp_this, Usw_this;
	int count = 0;
	for( int tsp_i = 0; tsp_i<2; tsp_i++)
	{
		for( int Usp_i = 0; Usp_i < 2; Usp_i++)
		{
			for( int Usw_i = 0; Usw_i < 2; Usw_i++)
			{
				count++;
				// see if this run is plus or minus one standard deviation about the mean parameter
				if(tsp_i == 0)
				{
					tsp_this = tsp_sd_min;
				}
				else
				{
					tsp_this = tsp_sd_plus;
				}

				if(Usp_i == 0)
				{
					Usp_this = Usp_sd_min;
				}
				else
				{
					Usp_this = Usp_sd_plus;
				}

				if(Usw_i == 0)
				{
					Usw_this = Usw_sd_min;
				}
				else
				{
					Usw_this = Usw_sd_plus;
				}

				cout << "model run, ts_peak: " << tsp_this << " Us_peak: " << Usp_this
				     << " and Us_width: " << Usw_this << endl;

				// run the hillslope model
				Hillslope.reset_hillslope(tsp_this, Usp_this, Usw_this);
				Hillslope.run_based_on_data_spacing(ts_temp, E_star_modelled,R_star_modelled, tolerance);

				// print results to file
				ofstream data_out2;
				num = itoa(count);
				string data_fname = of_prefix+dot+num+dot+ext;
				data_out2.open(data_fname.c_str());
				data_out2 << tsp_this << " " << Usp_this << " " << Usw_this << " -9999" << endl;
				for(int i = 0; i<n_data_elements; i++)
				{
					U_star_modelled[i] = Hillslope.gaussian_uplift( ts_temp[i] );
					data_out2 << ts_temp[i] << " " << E_star_modelled[i] << " " << R_star_modelled[i]
					          << " " << U_star_modelled[i] << endl;
				}
				data_out2.close();
			}
		}
	}
/*
	// now print the data to file
	ofstream data_out;
	data_out.open(outfile_fname);

	for(int i = 0; i< n_data_elements; i++)
	{
		data_out << ts_temp[i] << " " << Es_temp[i] << " " << Ese_temp[i] << " "
				 << max_E_star[i] << " " << mean_E_star[i] << " " << min_E_star[i] << " "
				 << Rs_temp[i] << " " << Rse_temp[i] << " "
				 << max_R_star[i] << " " << mean_R_star[i] << " " << min_R_star[i] << " "
				 << max_U_star[i] << " " << mean_U_star[i] << " " << min_U_star[i] << endl;
	}
	data_out.close();
	*/
}

