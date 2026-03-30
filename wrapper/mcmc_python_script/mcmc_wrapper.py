import os; os.environ["OMP_NUM_THREADS"] = "1" # Avoid oversubscription of resources with OMP

import emcee
from multiprocessing import Pool
import numpy as np
import pandas as pd
import sys; sys.path.append("/trm_software/wrapper/")

# np.random.seed(75)

import lcurve_wrapper

# Define globals accessible to all of the walkers.
G = 6.673e-11
Msol = 1.9891e30
Rsol = 6.958e8
c = 299792458
defaults = {}        # Dictionary to contain values within taken from a provided parameters.txt file.
worker_models = None # Each walker gets its own set of LCurve::Model models. See the "init_worker()" function

# List of all adjustable model parameters in their required order.
all_params = ["q", "iangle", "r1", "r2", "cphi3", "cphi4", "spin1", "spin2", "t1", "t2",
              "ldc1_1", "ldc1_2", "ldc1_3", "ldc1_4", "ldc2_1", "ldc2_2", "ldc2_3", "ldc2_4",
              "velocity_scale", "beam_factor1", "beam_factor2",
              "t0", "period", "pdot", "deltat", "gravity_dark1", "gravity_dark2", "absorb",
              "slope", "quad", "cube", "third",
              "rdisc1", "rdisc2", "height_disc", "beta_disc", "temp_disc", "texp_disc", "lin_limb_disc", "quad_limb_disc", "temp_edge", "absorb_edge",
              "radius_spot", "length_spot", "height_spot", "expon_spot", "epow_spot", "angle_spot", "yaw_spot", "temp_spot", "tilt_spot", "cfrac_spot",
              "stsp11_long", "stsp11_lat", "stsp11_fwhm", "stsp11_tcen",
              "stsp12_long", "stsp12_lat", "stsp12_fwhm", "stsp12_tcen",
              "stsp13_long", "stsp13_lat", "stsp13_fwhm", "stsp13_tcen",
              "stsp21_long", "stsp21_lat", "stsp21_fwhm", "stsp21_tcen",
              "stsp22_long", "stsp22_lat", "stsp22_fwhm", "stsp22_tcen",
              "uesp_long1", "uesp_long2", "uesp_lathw", "uesp_taper", "uesp_temp",
             ]
param_map = {k:i for i,k in enumerate(all_params)}


def init_worker(base_models_dict):
    """
    When calling emcee with Pool:
        Initialize the worker, giving it its own set of filter-dependent models to play with

    Clones the base LCurve::Model directly.
    This avoids having to create/destroy a new set of models every iteration
    """
    global worker_models
    worker_models = {filter: model.clone() for filter, model in base_models_dict.items()}

class lcurve_model():
    def __init__(self, parameters_file, data_file, free_params, filter):
        '''
        Lcurve::Model  https://github.com/trmrsh/cpp-lcurve/blob/master/include/trm/lcurve.h#L249
        
        Build an LCurve::Model object with access to the legacy C++ methods:
            "compute_chisq", "get_model", "get_param_names", "get_params", "nvary", "set_params"

        See: "/trm_software/wrapper/wrapper.cpp" within the apptainer image for details.

        Parameters:
            parameters_file : string
                parameters.txt file
                
            data_file : string
                5-column space delimited LCURVE format light curve data file 
                
            free_params : numpy array of strings
                list of free parameter names   
                
            filter : string
                filter being used to build the model and represent the data

        Attributes:
            self.free_params : numpy array of strings
                Contains parameter names: ["q", "iangle", "t2", "height_disc", ...]
                
            self.filter : string
                Filter used with the data and model
                
            self.model : <class 'lcurve_wrapper.LCurveModel'>
                Legacy C++ LCurve::Model class object
                
            self.model_object : <class 'lcurve_wrapper.Model'>
                LCurve::Model wrapper. Used to communicate individual parameter value changes
                
            self.relevant_params_dict : dict{string: int}
                Relevant parameters taken from the full list of fitted parameters.
                We can fit filter-dependent parameters. Here we exlude parameters that use a filter not matching this model's self.filter attribute
                
            self.relevant_params_indices : Numpy array of integers
                A simple array of indices used to filter the theta vector given by emcee on each iteration
                self.relevant_params_dict.values()
        '''
        self.free_params = free_params
        self.filter = filter
        self.model = lcurve_wrapper.LCurveModel(parameters_file, data_file)
        self.model_object = self.model.get_model() # The model object used to modify parameters individually

        # Set all parameter.vary to False
        for param in all_params:
            getattr(self.model_object, param).vary = False

        # Extract list indices for relevant parameters. Used to set parameters later
        self.relevant_params_dict = {}
        for i, param in enumerate(free_params):
            param = str(param)
            if param.endswith(f"_{self.filter}"): # Filter-dependent parameter needs to have the "_filter" part removed now.
                param = param[:-(len(self.filter)+1)]
            
            if param in all_params:
                self.relevant_params_dict[param] = i
                getattr(self.model_object, param).vary = True

        self.relevant_indices = np.array(list(self.relevant_params_dict.values()))
        
    def get_chi_squared(self, theta): # __call__
        '''
        Given a set of model parameters:
            1) Update the model parameters using "set_params()"
            2) Call "compute_chisq()":  
                    *) Build the model light curve (C++)
                    *) Calculate chi-squared (C++)
                    *) Returns chi-squared, flux-weighted surface gravities, and flux2 contribution
            3) Return chi-squared, flux-weighted surface gravities, and flux2 contribution
        '''
        try:
            self.model.set_params(theta[self.relevant_indices])
        except RuntimeError as err:
            # print("Model set_params() failed: ",err)
            return np.inf, None, None, None

        # # Update parameters one by one, rather than all at once
        # for key, value in self.relevant_params_dict.items():
        #     # print(f"Updating  {key}  ->  {theta[value]}")
        #     getattr(self.model_object, key).value = theta[value]

        try:
            chi_squared, logg2, logg1, flux2_contribution = self.model.compute_chisq()
        except Exception as err:
            # print("Model compute_chisq() failed: ",err)
            return np.inf, None, None, None

        return chi_squared, logg2, logg1, flux2_contribution

    def clone(self):
        '''
        Create an independent clone of this class for model modifications.
        '''
        # Create a blank lcurve_model object with no attributes, skipping the __init__ call by using __new__
        cloned_model = lcurve_model.__new__(lcurve_model)

        # Manually assign attributes to the clone using attributes from the base model.
        cloned_model.model = self.model.clone() # Clones the C++ LCurve::Model, sharing a reference to the light curve data
        cloned_model.model_object = cloned_model.model.get_model() # Run get_model() on the cloned object, not cloning the original's get_model()
        cloned_model.free_params = self.free_params
        cloned_model.filter = self.filter
        cloned_model.relevant_params_dict = self.relevant_params_dict
        cloned_model.relevant_indices = self.relevant_indices
        
        return cloned_model

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
    
def gaussian_prior_symmetric(walker_value, observed_value, observed_sigma):
    '''
    Gaussian prior with symmetric lower and upper sigma values
    '''
    return ((walker_value - observed_value)/observed_sigma)**2

def gaussian_prior_asymmetric(walker_value, observed_value, observed_sigma_lower, observed_sigma_upper):
    '''
    Gaussian prior with different observed lower and upper sigma values
    '''
    if walker_value > observed_value:
        return ((walker_value - observed_value)/observed_sigma_upper)**2
    else:
        return ((walker_value - observed_value)/observed_sigma_lower)**2

def get_parameter(parameter_name, theta, variable_parameters):
    '''
    Extract the parameter value (float) from the theta vector provided by emcee
    Uses the variable_parameters array to determine the parameter's index

    Returns first satisfied of:
        1) The parameter value taken from theta (new parameter values suggested by emcee)
        2) The parameter value taken from the global "defaults" dictionary (static parameters.txt file)
        3) None
    '''
    global defaults

    if parameter_name in variable_parameters:
        index = np.where(variable_parameters==parameter_name)[0][0]
        return theta[index]
    elif parameter_name in defaults:
        return defaults[parameter_name]
    else:
        return None
    
def log_prior(theta, variable_parameters, logg2, logg1):
    '''
    Return chi-squared using prior information.
    
    TODO?: Make this a "log_prior" class with attributes for
            total chi-squared, methods for uniform/gaussian,
            __init__ takes theta and variable_parameters,
            use properties for masses, velocities, etc

    Inputs:
        theta : np.array of floats
            np.array of parameter values from emcee
        variable_parameters : np.array of strings
            variable parameter names
        logg2 : float
            flux-weighted surface gravity of star2 taken from the lroche output
        logg1 : float
            flux-weighted surface gravity of star1 taken from the lroche output

    Returns:
        chi_squared : float or +np.inf    
    '''
    chi_squared = 0.0

    q = get_parameter("q", theta, variable_parameters)
    chi_squared += uniform_prior(q, 0.0, 100.0)

    iangle = get_parameter("iangle", theta, variable_parameters)
    chi_squared += uniform_prior(iangle, 50.0, 130.0)

    vscale = get_parameter("velocity_scale", theta, variable_parameters)
    chi_squared += uniform_prior(vscale, 330.0, 1000.0)

    t2 = get_parameter("t2", theta, variable_parameters)
    chi_squared += uniform_prior(t2, 500.0, 60000.0)

    t1 = get_parameter("t1", theta, variable_parameters)
    chi_squared += uniform_prior(t1, 500.0, 100000.0)

    height_disc = get_parameter("height_disc", theta, variable_parameters)
    chi_squared += uniform_prior(height_disc, 0.0, 0.2)
    
    temp_disc = get_parameter("temp_disc", theta, variable_parameters)
    chi_squared += uniform_prior(temp_disc, 0, 100000.0)
    
    if np.isinf(chi_squared):
        return np.inf

    period = get_parameter("period", theta, variable_parameters)
        
    k1     =  q/(1+q)*vscale*np.sin(np.radians(iangle))                   # km/s
    k2     =  1./(1+q)*vscale*np.sin(np.radians(iangle))                  # km/s
    a      =  period*(86400)/(2*np.pi)*(vscale*1000)/Rsol                   # Rsol; assumes circular orbit
    M1     =  1./(1+q)*(period*86400)/(2*np.pi*G)*(vscale*1000)**3/Msol # Msol
    M2     =  q/(1+q)*(period*86400)/(2*np.pi*G)*(vscale*1000)**3/Msol  # Msol

    chi_squared += gaussian_prior_symmetric(k2, 328.68, 7.19)
    chi_squared += gaussian_prior_symmetric(t2, 35840.0, 640.0)
    chi_squared += gaussian_prior_symmetric(logg2, 5.5315, 0.09257)
    
    # Hard-coded limb and gravity darkening priors on the visible star by filter specific to BLG723:
    # TODO: Add Claret+2020 here (as part of the lcurve_model class?) to avoid hard-coding.
    try:
        darkening_constraints = {
               "i2": {"ldc2_1":[0.039690560, 0.040845147, 0.041944320], "ldc2_2":[0.164263072, 0.169791616, 0.175286048], "gravity_dark2":[0.354223616, 0.370227237, 0.385850432]},
               "g":  {"ldc2_1":[0.049968896, 0.051689723, 0.053434880], "ldc2_2":[0.240295232, 0.246197248, 0.252024992], "gravity_dark2":[0.375716864, 0.392724053, 0.409374528]},
               "R":  {"ldc2_1":[0.038720160, 0.040133488, 0.041517920], "ldc2_2":[0.174187328, 0.179642389, 0.185094048], "gravity_dark2":[0.369629216, 0.386514565, 0.402888992]},
               "B":  {"ldc2_1":[0.053280896, 0.054998688, 0.056750176], "ldc2_2":[0.255866976, 0.261733776, 0.267515008], "gravity_dark2":[0.376532192, 0.393527787, 0.410227328]},
               }
    
        for filter in darkening_constraints:
            for param in ["ldc2_1", "ldc2_2", "gravity_dark2"]:
                walker_position = get_parameter(f"{param}_{filter}", theta, variable_parameters)
                expected_mean = darkening_constraints[filter][param][1]     

                expected_upper = max(darkening_constraints[filter][param])
                sigma_upper = abs(expected_upper - expected_mean)

                expected_lower = min(darkening_constraints[filter][param])
                sigma_lower = abs(expected_mean - expected_lower)
    
                chi_squared += gaussian_prior_asymmetric(walker_position, expected_mean, sigma_lower, sigma_upper)
    except Exception as err:
        # print(f"Failed to use limb/gravity darkening priors: {err}")
        pass
        
    return chi_squared
    

def log_probability(theta, variable_parameters):
    '''
    Calculate log-probability for the set of parameters theta.
    Adds the chi-squared from the LCurve::Model evaluation with the chi-squared from log_prior() function.

    Inputs:
        theta : np.array of floats
            vector of suggested parameter values from emcee
        variable_parameters : np.array of strings
            variable parameter names

    Returns:
        -0.5 * chi_squared
    '''
    chi_squared = 0.0

    # Loop over filters
    # TODO: Add another Pool here to parallel the filters?
    for filter, model in worker_models.items():
        model_chi_squared, logg2, logg1, flux2_contribution = model.get_chi_squared(theta)
        chi_squared += model_chi_squared

    # Include priors
    if not np.isinf(chi_squared):
        chi_squared += log_prior(theta, variable_parameters, logg2, logg1)
    
    return (-0.5 * chi_squared,)

def get_parameter_order(parameter_name, filters):
    '''
    Using the legacy C++ set_param() method requires that the list of variable parameter be in a specific order.
    Here we find the index/position that each fitted parameter should be in to match the required order.
    This handles filter-dependent parameters by stripping them of the "{_filter}" before matching.
    '''
    for filter in filters:
        if parameter_name.endswith(f"_{filter}"):
            parameter_name = parameter_name[:-(len(filter)+1)]
            break
    return param_map[parameter_name]

def main(profile=False):
    
    ncores = 64
    nwalkers = 128 # Use at least 2*ncores for efficiency with emcee.moves.RedBlueMove

    nsteps = 2000
    output_filename = "chain.h5"

    fresh_mcmc = False # Reset the output file for a fresh chain

    init_param_filename = "init_parameters.txt"

    parameters_file = "../parameters_{}.txt"
    data_file = "../lcs/blg723_{}_lcurve2.txt"
    filters = ["B", "R", "i2", "g"]


    # ================================================================================
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ================================================================================
    
    
    # Populate the global defaults dictionary with values in the first of the provided parameters.txt files.
    with open(parameters_file.format(filters[0]), "r") as ifile:
        for line in ifile.readlines():
            if "=" in line:
                k = line.split()
                if len(k) < 4 and not line.startswith("tperiod"):
                    continue
                defaults[k[0]] = float(k[2])    

    # Read the to-be fitted parameters and their lower and upper bounds for randomization.
    init_params = np.loadtxt(init_param_filename, dtype=[("parameter_name", "U20"),
                                                         ("lower_limit", np.float64),
                                                         ("upper_limit", np.float64)]
                                                        )
    ndim = len(init_params)

    # Re-order init_params to match the required order for the legacy LCurve::Model::set_param() method
    init_params = init_params[
        np.argsort(
            [get_parameter_order(name, filters) for name in init_params["parameter_name"]]
        )
    ]
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
    
    # Generate a set of models to use as a basis for modification. These depend on your provided parameters.txt files for the "wavelength", LDC/GDC, and beaming values
    base_models = {filter: lcurve_model(parameters_file.format(filter), data_file.format(filter), init_params["parameter_name"], filter) for filter in filters}

    if profile:
        backend = emcee.backends.HDFBackend(output_filename)
        backend.reset(nwalkers, ndim)
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        log_probability,
                                        args=(init_params["parameter_name"],),
                                        backend=backend,
                                       )
        sampler.run_mcmc(p0, nsteps, progress=False)
        
        return
        
    
    # Generate workers, each with their own set of models to modify in memory.
    with Pool(ncores, initializer=init_worker, initargs=(base_models,)) as pool:

        # Define the backend for saving the outputs.
        backend = emcee.backends.HDFBackend(output_filename)
        if fresh_mcmc:
            backend.reset(nwalkers, ndim)
        
        # Initialize the sampler. TODO: include parameter_names -> update log_probability to accept a dict for convenience
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        log_probability,
                                        args=(init_params["parameter_name"],),
                                        pool=pool,
                                        backend=backend,
                                       )
        sampler.run_mcmc(state, nsteps, progress=False)


if __name__ == "__main__":

    profile = sys.argv[1].lower() in ["true", "1", "yes"]
    if profile:
        print("Running cProfile", flush=True)
        profile_filename = "profile_wrapper.prof"
        import cProfile
        cProfile.run("main(profile=True)", profile_filename)

        import pstats
        p = pstats.Stats(profile_filename)
        p.sort_stats("cumulative").print_stats(30)
        # Note: _thread.lock() is misleading since it includes all work done by the Pool, so we do not use Pool when profiling
        # backend writes appear to be at least 0.02s per step
        # Each compute_chisq() seems to be about 2.75s per call (5.5s for two filters, 11s for four filters). depends on grid resolution
    else:
        main()
