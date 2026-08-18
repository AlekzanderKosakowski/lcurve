import emcee
import numpy as np
import os; os.environ["OMP_NUM_THREADS"] = "1" # Avoid oversubscription of resources with OMP
import pandas as pd
import shutil # Used to copy chain.h5 as a backup

from datetime import datetime
from multiprocessing import Pool
from scipy.optimize import minimize

import sys; sys.path.append("/trm_software/wrapper/")
import lcurve_wrapper
from lcurve_model import *

# Define globals accessible to all of the walkers.
G = 6.673e-11    # mks
Msol = 1.9891e30 # mks
Rsol = 6.958e8   # mks
c = 299792458    # mks

worker_reference_model_parameters = None
worker_models = None # Each walker gets its own set of LCurve::Model models.

def init_worker(base_models_dict):
    """
    Initialize worker namespace
    Each worker is given its own set of N Lcurve_models to play with.
        where N is the number of filters.
        
    Clones the base LCurve::Model directly.
    This avoids having to create/destroy a new set of models for every MCMC iteration
    """
    global worker_models
    global worker_reference_model_parameters
    worker_models = {filter: model.clone() for filter, model in base_models_dict.items()}
    worker_reference_model_parameters = next(iter(base_models_dict.values())).summarize_model_parameters()


def test_parameters_match_between_models(base_models, fitted_parameter_names):
    '''
    Check non-fitted parameters in each model to confirm they are identical.
    This ensures that one model isn't using different values than another.

    Skips filter dependent parameters and fitted parameters.

    Inputs:
        base_models            : dictionary  {filter: Lcurve_model}
        fitted_parameter_names : 1D array of strings representing parameters to be fitted

    Raises:
        ValueError if not np.isclose() between two model parameter values.
    '''
    if len(base_models) == 1:
        return

    filter_dependent_parameters = ["ldc1_1", "ldc1_2", "ldc1_3", "ldc1_4",
                                   "ldc2_1", "ldc2_2", "ldc2_3", "ldc2_4",
                                   "ldc3_1", "ldc3_2", "ldc3_3", "ldc3_4",
                                   "beam_factor1", "beam_factor2",
                                   "gravity_dark1", "gravity_dark2",
                                   "slope", "quad", "cube",
                                   "wavelength", "filter"]

    mismatched_parameters = set()
    for i, model in enumerate(base_models.values()):
        
        if i==0:
            reference_parameter_dict = model.summarize_model_parameters()
            continue

        compare_parameter_dict = model.summarize_model_parameters()
        for parameter_name, reference_value in reference_parameter_dict.items():
            if parameter_name in list(fitted_parameter_names) + list(filter_dependent_parameters):
                continue

            compare_value = compare_parameter_dict[parameter_name]
            if not np.isclose(reference_value, compare_value):
                mismatched_parameters.add(parameter_name)

    if mismatched_parameters:
        raise ValueError(f"Mismatched parameters between input parameters files: {sorted(mismatched_parameters)}")
    
    return

def run_optimize(theta0, variable_parameters):
    '''
    Minimize chi-squared (-2*log_probability) of the combined multi-filter model

    Used with the optimize option.
    Replacement for simplex, with support for multi-filter and priors.

    Approximates convergence with Nelder-Mead algorithm using 
    approximately 1% change in optimal chi-squared between iterations as convergence.

    Input:
        theta0 : 2D np.ndarray of floats representing parameter values of each walker.
        variable_parameters : array of strings
                              representing parameter names in the same order as theta0
    '''    
    theta0 = np.average(theta0, axis=0)

    print("="*40)
    print("Starting parameters:")
    for i,k in enumerate(variable_parameters):
        print(f"  {k:<20s}  {theta0[i]:.5f}")

    N_data = sum(len(model.data) for model in worker_models.values())
    print(N_data)
    N_variable = len(variable_parameters)
    
    fxn = lambda x: -2*log_probability( dict(zip(variable_parameters, x)) )[0]
    results = minimize(fxn,
                       theta0,
                       method="Nelder-Mead",
                       options={"maxiter":2000,
                                "fatol":0.01*(N_data - N_variable),
                                "xatol":np.inf,
                                "disp":True,
                                }
                       )

    theta_best = results.x
    print("="*40)
    print("Final parameters:")
    for i,k in enumerate(variable_parameters):
        print(f"  {k:<20s}  {theta0[i]:.5f}  ->  {theta_best[i]:.5f}")
    print("="*40)
    
    return theta_best

def uniform_prior(parameter_value, lower, upper):
    '''
    Compares the input parameter value with lower and upper bounds.
    Returns:
        0.0    if within bounds
        np.inf if out of bounds
    '''
    if parameter_value < lower or parameter_value > upper:
        return np.inf
    return 0.0
    
def gaussian_prior_symmetric(walker_value, prior_value, prior_sigma):
    '''
    Gaussian prior with symmetric lower and upper sigma values

    Input:
        walker_value : float
                       Parameter value being tested by the walker.
        prior_value  : float
                       Parameter value from prior information.
        prior_sigma  : float
                       Parameter uncertainty from prior information.    
    '''
    return ((walker_value - prior_value)/prior_sigma)**2

def gaussian_prior_asymmetric(walker_value, prior_value, prior_sigma_lower, prior_sigma_upper):
    '''
    Piecewise Gaussian approximation to an asymmetric prior.

    Input:
        walker_value       : float
                             Parameter value being tested by the walker.
        prior_value        : float
                             Parameter value from prior information.
        prior_sigma_lower  : float
                             Parameter uncertainty from prior information.    
                             ~14-percentile from Gaussian distribution
        prior_sigma_upper  : float
                             Parameter uncertainty from prior information.    
                             ~86-percentile from Gaussian distribution
    '''
    if walker_value < prior_value:
        prior_sigma = prior_sigma_lower
    else:
        prior_sigma = prior_sigma_upper

    return ((walker_value - prior_value)/prior_sigma)**2

def log_probability(theta):
    '''
    Calculate log-likelihood for the set of parameters
        Update parameter values.
        Sums the chi-squared from lroche output across all filters
        Adds chi-squared from prior_chisq() call.
    
    Inputs:
        theta : Python dictionary of fitted parameters
                May contain filter-specific names, such as "ldc1_1_R" and "ldc1_1_V"

    Returns:
        -0.5 * chi_squared
    '''
    chi_squared = 0.0

    filters = worker_models.keys()
    
    for filter, model in worker_models.items(): 
        for parameter_name_F, suggested_value in theta.items():

            if (parameter_name_F.split("_")[-1] in filters) and (parameter_name_F.split("_")[-1] != model.filter):
                continue # Skip filter-dependent parameters when the loop iteration is not on that filter.
                
            elif (parameter_name_F.split("_")[-1] in filters) and (parameter_name_F_split("_")[-1] == model.filter):
                parameter_name = "_".join(parameter_name0.split("_")[:-1]) # Remove the "_{filter}" text from the parameter name

            else: # Not a filter-specific parameter, so use the parameter name as-is
                parameter_name = parameter_name_F

            getattr(model.parameters, parameter_name).value = suggested_value # Set each parameter
            
        try:
            lroche_output = model.lroche(scale=True)
            chi_squared += lroche_output["chisq"]
        except Exception as err:
            print(f"model.lroche() failed.\n{type(err).__name__}: {err}")
            chi_squared += np.inf

        if np.isinf(chi_squared):
            return -np.inf
    
    if not np.isinf(chi_squared):
        chi_squared += prior_chisq(theta, lroche_output)

    return (-0.5*chi_squared, 
            {
                "logg1":lroche_output["logg1"],
                "logg2":lroche_output["logg2"],
                "rvol1":lroche_output["rvol1"],
                "rvol2":lroche_output["rvol2"]
            }
           )

def prior_chisq(theta, lroche_output):
    '''
    Return chi-squared using prior information.

    Inputs:
        theta : Python dictionary of fitted parameters
                May contain filter-specific names, such as "ldc1_1_R" and "ldc1_1_V"
        lroche_output : dictionary
            keys = ['model_flux', 'chisq', 'wnok', 'wdwarf', 'logg1', 'logg2', 'rvol1', 'rvol2', 'ffac1', 'ffac2', 'sfac']
        
    Returns:
        chi_squared : float or np.inf
    '''
    chi_squared = 0.0

    if "iangle" in theta:
        chi_squared += uniform_prior(theta["iangle"], 0.0, 90.0)

    if "t1" in theta:
        chi_squared += gaussian_prior_symmetric(theta["t1"], 28_900, 400)

    if "t2" in theta:
        chi_squared += uniform_prior(theta["t2"], 500, 50_000)

    if np.isinf(chi_squared):
        return np.inf


    use_physical_priors = True
    if not use_physical_priors:
        return chi_squared
    
    # Priors based on physical constraints
    period_days = worker_reference_model_parameters["tperiod"]

    # Use fixed reference value from model parameters if not being fitted, else use the current iteration's fit value.
    q      = theta["q"] if "q" in theta else worker_reference_model_parameters["q"]
    vscale = theta["velocity_scale"] if "velocity_scale" in theta else worker_reference_model_parameters["velocity_scale"]
    iangle = theta["iangle"] if "iangle" in theta else worker_reference_model_parameters["iangle"]

    a  =  period_days*(86400)/(2*np.pi)*(vscale*1000)/Rsol               # Rsol; assumes circular orbit
    R1vol = lroche_output["rvol1"] * a
    R2vol = lroche_output["rvol2"] * a
    K1 =  q/(1+q)*vscale*np.sin(np.radians(iangle))                      # km/s
    K2 =  1./(1+q)*vscale*np.sin(np.radians(iangle))                     # km/s
    M1 =  1./(1+q)*(period_days*86400)/(2*np.pi*G)*(vscale*1000)**3/Msol # Msol
    M2 =  q/(1+q)*(period_days*86400)/(2*np.pi*G)*(vscale*1000)**3/Msol  # Msol

    chi_squared += gaussian_prior_symmetric(K1, 71.6, 1.7)
    chi_squared += gaussian_prior_symmetric(lroche_output["logg1"], 5.64, 0.02)
    
    return chi_squared
    

def main(profile=False):

    optimize = True
    
    ncores = 32
    nwalkers = 2*ncores # Use at least 2*ncores for efficiency with emcee.moves.RedBlueMove

    nsteps = 500
    output_filename = "chain.h5"

    fresh_mcmc = True

    init_param_filename = "init_parameters.txt"

    parameters_file = "parameters_files/parameters_simplex_{}.txt"
    data_file = "data_files/ec10246-2707_noise_model_{}.txt"
    filters = ["V"]


    # ================================================================================
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ================================================================================

    # Read the to-be fitted parameters and their lower and upper bounds for randomization.
    init_params = np.loadtxt(init_param_filename, dtype=[("parameter_name", "U20"),
                                                         ("lower_limit", np.float64),
                                                         ("upper_limit", np.float64)
                                                        ]
                            )
    ndim = len(init_params)
    print(f"Using {ndim} free parameters: {init_params['parameter_name']}")


    if not os.path.isfile(output_filename):
        fresh_mcmc = True

    # Generate random starting positions for all parameters for all walkers
    p0 = [
        [
            np.random.uniform(init_params["lower_limit"][i], init_params["upper_limit"][i]) for i in range(len(init_params))
        ] for n in range(nwalkers)
    ]
    state = p0 if fresh_mcmc else None    
    
    # Generate a set of models to use as a base for modification.
    # These depend on your provided parameters.txt files to allow different wavelength, LDC/GDC, and beaming values
    base_models = {filter: Lcurve_model(parameters_file.format(filter),
                                        data=data_file.format(filter),
                                        filter=filter)
                   for filter in filters
                  }
    
    # Confirm that the non-fitted parameters between filters match.
    test_parameters_match_between_models(base_models, init_params["parameter_name"])

    if optimize:
        print(f"Running simplex optimize with {len(base_models)} filters: {[k for k in base_models]}")
        init_worker(base_models)
        
        simplex_solution = run_optimize(p0, init_params["parameter_name"])
        sys.exit()

    # Generate workers, each with their own set of models to modify in memory.
    with Pool(ncores, initializer=init_worker, initargs=(base_models,)) as pool:

        # Define the backend for saving the outputs.
        backend = emcee.backends.HDFBackend(output_filename)
        if fresh_mcmc:
            if os.path.isfile(output_filename):
                print(f"Creating a back-up of {output_filename}")
                shutil.copy(output_filename, f"bk_{output_filename}")
            backend.reset(nwalkers, ndim)

        # Data type for the extra output stored with the walker positions.
        # Must be in the same order as the log_probability() output
        blobs_dtype = [
            ("logg1", float),
            ("logg2", float),
            ("rvol1", float),
            ("rvol2", float),
        ]
        
        # Initialize the sampler.
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        log_probability,
                                        parameter_names=init_params["parameter_name"],
                                        pool=pool,
                                        blobs_dtype=blobs_dtype,
                                        backend=backend,
                                       )
        
        # Run the sampler to completion.
        sampler.run_mcmc(state, nsteps, progress=False)


if __name__ == "__main__":

        main()
