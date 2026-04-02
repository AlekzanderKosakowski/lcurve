#include "trm/subs.h"
#include "trm/constants.h"

/** Computes the Planck function Bnu = (2 h \nu^3/c^2)/(exp(h \nu/kT) - 1)
 *  as a function of wavelength and temperature. Output units are W/m**2/Hz/sr.
 *
 * \param wave wavelength in nanometres
 * \param temp temperature in K
 */

double Subs::planck(double wave, double temp){
  
    const double FAC1 = 2.e27*Constants::H*Constants::C;
    const double FAC2 = 1.e9*Constants::H*Constants::C/Constants::K;
    
    double efac = FAC2/(wave*temp);
    if(efac > 40.){
    	return FAC1*exp(-efac)/(wave*sqr(wave));
    }else{
    	return FAC1/(exp(efac)-1.)/(wave*sqr(wave));
    }
}

/** Computes the logarithmic derivative of the Planck function 
 * Bnu wrt wavelength (i.e. d ln(Bnu) / d ln(lambda)) as a function of wavelength and temperature
 * \param wave wavelength in nanometres
 * \param temp temperature in K
 */

double Subs::dplanck(double wave, double temp){
    
    const double FAC2 = 1.e9*Constants::H*Constants::C/Constants::K;
    
    double efac = FAC2/(wave*temp);
    return efac/(1.-exp(-efac)) - 3.;
}

/** Computes the logarithmic derivative of the Planck function 
 * Bnu wrt T (i.e. d ln(Bnu) / d ln(T)) as a function of wavelength and temperature
 * \param wave wavelength in nanometres
 * \param temp temperature in K
 */

double Subs::dlpdlt(double wave, double temp){
    
    const double FAC2 = 1.e9*Constants::H*Constants::C/Constants::K;
    
    double efac = FAC2/(wave*temp);
    return efac/(1.-exp(-efac));
}


// Given a filter name and temperature grid, return a vector<double> of planck fluxes integrated across the filter profile.
void Subs::integrate_filter(const std::vector<double>& temperature_array, std::vector<double>& planck_array, std::string filter){
    std::string filter_curve_path = "/trm_software/filter_curves/";
    std::string filename = filter_curve_path + filter + ".txt";

    std::vector<double> lambda;
    std::vector<double> R;

    double lam, r;
    std::string line;

    
    // Open the filter profile (two columns, space delimited, no header line).
    std::ifstream fin(filename);
    if(!fin.is_open()){
        throw std::runtime_error("Could not open filter file " + filename);
    }

    // Read in the filter profile line by line
    while(std::getline(fin, line)) {
        if(line.empty()){
            continue;
        }
        std::istringstream iss(line);
        if(!(iss >> lam >> r)) continue;
        lambda.push_back(lam);
        R.push_back(r);
    }
    fin.close();

    if(lambda.size() < 2){
        throw std::runtime_error("Too few properly-formatted lines in file " + filename);
    }
    
    for(size_t t=0; t < temperature_array.size(); ++t){
        double T = temperature_array[t];
        double numerator = 0.0;
        double denominator = 0.0;

        for(size_t k=0; k < lambda.size()-1; ++k){
            double dlambda = lambda[k+1] - lambda[k];
            double R_avg = 0.5 * (R[k] + R[k+1]);
            double B_avg = 0.5 * (Subs::planck(lambda[k], T) + Subs::planck(lambda[k+1], T));

            numerator += B_avg * R_avg * dlambda;
            denominator += R_avg * dlambda;
        }

        planck_array[t] = numerator / denominator;
        
    }
    
}

// 1D linear interpolator (within planck.cc to avoid changing the Makefile to link linterp.cc)
// Takes a 1D vector x_grid and the corresponding 1D vector y_grid
// Returns the interpolated y-value at the desired x position.
double Subs::interp1d(const std::vector<double>& x_grid, const std::vector<double>& y_grid, double x){

    int N = x_grid.size();

    if (x <= x_grid[0]){
        return y_grid[0];
    }
    
    if (x >= x_grid[N-1]){
        return y_grid[N-1];
    }

    // Binary search for lower/upper bound interval
    int i = 0, j = N-1;
    while(j - i > 1){
        
        int m = (i + j)/2;
        
        if(x_grid[m] > x){
            j = m;
        }else{
            i = m;
        }
    }
    // Linear interpolation
    double t = (x - x_grid[i]) / (x_grid[j] - x_grid[i]);
    return y_grid[i] + t * (y_grid[j] - y_grid[i]);
}