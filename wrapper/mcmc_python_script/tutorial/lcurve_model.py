ximport sys; sys.path.append("/trm_software/wrapper/")
import lcurve_wrapper
import tempfile
import os
import numpy as np
import pandas as pd
import numbers

fitted_params = ["q", "iangle",
              "r1", "r2", "r3", "cphi3", "cphi4",
              "spin1", "spin2",
              "t1", "t2", "t3",
              "ldc1_1", "ldc1_2", "ldc1_3", "ldc1_4", "ldc2_1", "ldc2_2", "ldc2_3", "ldc2_4", "ldc3_1", "ldc3_2", "ldc3_3", "ldc3_4",
              "velocity_scale",
              "beam_factor1", "beam_factor2",
              "t0", "period", "pdot", "deltat",
              "gravity_dark1", "gravity_dark2",
              "absorb",
              "slope", "quad", "cube",
              "rdisc1", "rdisc2", "height_disc", "beta_disc", "temp_disc", "texp_disc", "lin_limb_disc", "quad_limb_disc", "temp_edge", "absorb_edge",
              "radius_spot", "length_spot", "height_spot", "expon_spot", "epow_spot", "angle_spot", "yaw_spot", "temp_spot", "tilt_spot", "cfrac_spot",
              "stsp11_long", "stsp11_lat", "stsp11_fwhm", "stsp11_tcen",
              "stsp12_long", "stsp12_lat", "stsp12_fwhm", "stsp12_tcen",
              "stsp13_long", "stsp13_lat", "stsp13_fwhm", "stsp13_tcen",
              "stsp21_long", "stsp21_lat", "stsp21_fwhm", "stsp21_tcen",
              "stsp22_long", "stsp22_lat", "stsp22_fwhm", "stsp22_tcen",
              "stsp1i_long", "stsp1i_lat", "stsp1i_fwhm_long", "stsp1i_fwhm_lat", "stsp1i_escale_len", "stsp1i_tcen"
              ]
state_params = ["tperiod", "roche1", "roche2", "eclipse1", "eclipse2",
                "add_spot", "add_disc", "opaque", "wavelength", "filter",
                "finite_irr12", "finite_irr21", "third", "nlat1f", "nlat2f", "nlat1c", "nlat2c"
               ]

class Lcurve_model():
    def __init__(self, parameters_file):
        self.lcurve = lcurve_wrapper.LCurveModel(parameters_file) # Hosts the C++ model code.

        self.parameters = self.lcurve.get_model() # Reference to the model parameters
        self.data = None # Reference to the current cached data (NumPy array with shape (N, 6))

        
        for param in fitted_params:
            getattr(self.parameters, param).vary = False
    
    def create_parameters_file(self, output_filename=None):
        '''
        Create a new parameters file using the built-in Lcurve methods
        If no filename is specified, then saves to a temporary file in the tmp/ directory.
        '''
        if output_filename is None:
            if not os.path.isdir("tmp"):
                os.mkdir("./tmp")
            fd, tmp_filename = tempfile.mkstemp(dir="./tmp", suffix=".txt", prefix="parameters_")
            os.close(fd)
            self.lcurve.write_to_parameters_file(tmp_filename,)
            return tmp_filename
        else:
            self.lcurve.write_to_parameters_file(output_filename)
            return output_filename
        
    def update_parameters(self, new_params):
        '''
        Update the model parameter values.

        Inputs
            new_params : dictionary
        '''
        if not isinstance(new_params, dict):
            raise TypeError(f"Input parameters type must be dictionary with numerical values. Got {type(new_params).__name__}")
        
        print()
        print("Updating model parameters:")
        for param, value in new_params.items():
            
            if param not in fitted_params:
                if param in state_params:
                    raise ValueError(f"Model parameter \"{param}\" is a state parameter, not a fitted parameter.")
                raise ValueError(f"Unknown model parameter: {param}")

            if not isinstance(value, numbers.Real):
                raise TypeError(f"Parameter \"{param}\" must be numeric. Got {type(value).__name__}")
            
            old_value = getattr(self.parameters, param).value
            if float(old_value) != float(value):
                print(f"  {param}: {old_value} -> {value}")
            getattr(self.parameters, param).value = value

    def show_parameters(self):
        '''
        Print the values for all current binary parameters.
        '''
        for param in fitted_params:
            print(f"{param}: {getattr(self.parameters, param).value}")

    def summarize_model(self):
        '''
        Print the values for all parameters and flags used by the model.
        '''
        param_dict = {k:0 for k in fitted_params + state_params}
        for param_name in dir(self.parameters):
            if param_name.startswith("_"):
                continue

            value = getattr(self.parameters, param_name)

            if hasattr(value, "value"):
                param_dict[param_name] = value.value
            else:
                param_dict[param_name] = value

        for param_name, param_value in param_dict.items():
            print(f"{param_name:>20s}  {param_value}")

    @staticmethod
    def convert_data_to_array(data):
        '''
        Convert the given input data into the (N, 6) NumPy array required by C++ methods.

        If data is provided as a 1D numpy array, then assume those are "time" and all other columns are 1
        '''
        if not isinstance(data, (str, pd.DataFrame, np.ndarray)):
            raise TypeError("Input data must be NumPy array, Pandas dataframe, or existing filename string")

        # Provided filename string.
        if isinstance(data, str):
            assert os.path.isfile(data), f"File not found: {data}"
            data = np.loadtxt(data, unpack=False, dtype=np.float64)

        # Provided Pandas DataFrame
        elif isinstance(data, pd.DataFrame):

            if "time" in data.columns:
                time_col = "time"
            elif "phase" in data.columns:
                time_col = "phase" # Phase here is treated as "time" in C++. Be sure to set parameters.period.value=1.0 if using phase
            else:
                raise ValueError("DataFrame must contain either 'time' or 'phase' column. Uses to 'time' if both present.")
    
            if "exptime" in data.columns:
                exposure_col = "exptime"
            elif "exposure" in data.columns:
                exposure_col = "exposure"
            else:
                raise ValueError("DataFrame must contain either 'exptime' or 'exposure' column. Uses to 'exptime' if both present")

            if "flux_error" in data.columns:
                flux_error_col = "flux_error"
            elif "flux_err" in data.columns:
                flux_error_col = "flux_err"
            else:
                raise ValueError("DataFrame must contain either 'flux_error' or 'flux_err' column. Uses 'flux_error' if both present")

            required_columns = ["weight", "ndiv"]
            missing_columns = [column for column in required_columns if column not in data.columns]
            if missing_columns:
                raise ValueError(f"Input dataFrame is missing required column(s): {missing_columns}")
    
            # Convert to the six-column representation expected by C++
            data = data[[time_col, exposure_col, "flux", flux_error_col, "weight", "ndiv"]].to_numpy(dtype=np.float64)

        # Provided NumPy array
        elif isinstance(data, np.ndarray):
            data = np.asarray(data, dtype=np.float64)

            if data.ndim == 1:
                time = data
                other_col = np.ones_like(time)
                data = np.asarray([time,
                                   other_col,
                                   other_col,
                                   other_col,
                                   other_col,
                                   other_col
                                  ]).T

        
        if data.ndim != 2:
            raise ValueError(f"Data must be 2-dimensional. Data shape = {data.shape}")
        if data.shape[1] != 6:
            raise ValueError(
                "Data must have 6 columns: "
                "[time, exptime, flux, flux_error, weight, ndiv]. "
                f"Data shape = {data.shape}"
            )
        if data.shape[0] <= 1:
            raise ValueError(
                f"Data must contain more than 1 observation. "
                f"Data shape = {data.shape}"
            )
    
        return data
        
    def set_data(self, data):
        '''
        Assign given data to the model object.
        Stores an LCurve::Data object on the C++ side
        Stores a NumPy array of shape (N, 6) on the Python side.

        Inputs:
            data : string or NumPy array or Pandas DataFrame
        '''
        data = self.convert_data_to_array(data)

        self.lcurve.set_data(data) # Stores an LCurve::Data object
        self.data = data           # Stores a numpy array
    
    def lroche(self, data=None, scale=True):
        '''
        Runs compute_light_curve() from the C++ code.
        Uses parameters from the stored model.parameters
        If data is provided, then uses timestamps and flux from that data
        If data is not provided, then uses timestamps and flux from the cached data.
        
        This is equivalent to running lroche

        Inputs:
            data : str, NumPy array, or Pandas DataFrame
                   None (Uses the cached data from set_data())
            scale : boolean
                    True:  Use SVD to scale the model flux to match data flux (for chi-squared calculations)
                    False: Use scale_factor = 1.0 for all fits (will output model_flux in B_lambda units)

        Outputs:
            lroche_output : dictionary containing typical lroche output keys:
                 model_flux (array), chisq, wnok, wdwarf, logg1, logg2, rvol1, rvol2, ffac1, ffac2, sfac (array)

        '''
        if data is None: # No data provided; use cached data
            return self.lcurve.compute_light_curve(scale=scale)
            
        else: # New data provided; convert it to a numpy array then use it temporarily
            data = self.convert_data_to_array(data)
            return self.lcurve.compute_light_curve(data, scale=scale)
    
    
    def clone(self):
        '''
        Create an independent clone of this model object.
        '''
        # Create a blank lcurve_model object with no attributes, skipping the __init__ call by using __new__
        cloned_model = Lcurve_model.__new__(Lcurve_model)

        # Manually assign attributes to the clone using attributes from the base model.
        cloned_model.lcurve = self.lcurve.clone() # Clones the C++ LCurve::Model, sharing a reference to the light curve data
        cloned_model.parameters = cloned_model.lcurve.get_model() # Run get_model() on the cloned object, not cloning the original's get_model()
        cloned_model.data = self.data
        
        return cloned_model
