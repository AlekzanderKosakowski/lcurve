import sys; sys.path.append("/trm_software/wrapper/")
import lcurve_wrapper
import tempfile
import os

all_params = ["q", "iangle",
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
              "stsp1i_long", "stsp1i_lat", "stsp1i_fwhm_long", "stsp1i_fwhm_lat", "stsp1i_escale_len", "stsp1i_tcen",
             ]

class Lcurve_model():
    def __init__(self, parameters_file, data_file,):
        self.model = lcurve_wrapper.LCurveModel(parameters_file, data_file)
        self.parameters = self.model.get_model()

        for param in all_params:
            getattr(self.parameters, param).vary = False

    def get_chi_squared(self,):
        try:
            chi_squared, logg2, logg1, flux2_contribution = self.model.compute_chisq()
        except Exception as err:
            print("Model compute_chisq() failed: ",err)
            return np.inf, None, None, None

        return chi_squared, logg2, logg1, flux2_contribution

    def save_parameters_file(self, output_filename=None):
        if output_filename is None:
            fd, tmp_filename = tempfile.mkstemp(dir="./tmp", suffix=".txt", prefix="parameters_")
            os.close(fd)
            self.model.write_to_parameters_file(tmp_filename,)
            return tmp_filename
        else:
            self.model.write_to_parameters_file(output_filename)
            return output_filename
            
    def set_parameters(self, new_params):
        print()
        print("Updating model parameters:")
        for param, value in new_params.items():
            old_value = getattr(self.parameters, param).value
            if float(old_value) != float(value):
                print(f"  {param}: {old_value} -> {value}")
            getattr(self.parameters, param).value = value

    def show_parameters(self,):
        for param in all_params:
            print(f"{param}: {getattr(self.parameters, param).value}")

    def summarize(self,):
        param_dict = {k:0 for k in all_params + ["tperiod", "roche1", "roche2", "eclipse1", "eclipse2", "add_spot", "add_disc", "opaque", "wavelength", "filter", "finite_irr12", "finite_irr21", "third", "nlat1f", "nlat2f", "nlat1c", "nlat2c"]}
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
    
    def clone(self):
        '''
        Create an independent clone of this class for model modifications.
        '''
        # Create a blank lcurve_model object with no attributes, skipping the __init__ call by using __new__
        cloned_model = Lcurve_model.__new__(Lcurve_model)

        # Manually assign attributes to the clone using attributes from the base model.
        cloned_model.model = self.model.clone() # Clones the C++ LCurve::Model, sharing a reference to the light curve data
        cloned_model.parameters = cloned_model.model.get_model() # Run get_model() on the cloned object, not cloning the original's get_model()

        return cloned_model
