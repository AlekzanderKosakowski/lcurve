// wrapper.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "trm/lcurve.h"
#include "trm/subs.h"
#include <iostream>
#include <string>

namespace py = pybind11;


// This is a wrapper for very specific parts of lcurve to be accessed from python
// It includes the Lcurve::Model class, but not all of the methods inside
//   It only include the methods we define here, such as get_params, set_params, compute_chisq()
class LCurveWrapper {
public:
    std::unique_ptr<Lcurve::Model> model_;
    std::shared_ptr<Lcurve::Data> data_;
    Subs::Buffer1D<double> sfac_;

    Subs::Array1D<double> fit_;
    Subs::Array1D<double> param_buffer_;

    // Constructor: initialize model from a file
    LCurveWrapper(const std::string& model_file,
                  const std::string& data_file) {
        try {

            // const char* dir = std::getenv("LCURVE_DIR");
            // const char* env = std::getenv("LCURVE_ENV");
            // std::cout << "DEBUG: LCURVE_DIR = " << (dir ? dir : "not set") << std::endl;
            // std::cout << "DEBUG: LCURVE_ENV = " << (env ? env : "not set") << std::endl;

            // char cwd[1024];
            // if (getcwd(cwd, sizeof(cwd)) != nullptr)
            //     std::cout << "DEBUG: Current working directory = " << cwd << std::endl;

            // std::ifstream fin_model(model_file.c_str());
            // if (!fin_model) {
            //     std::cerr << "DEBUG: Tried to open file: " << model_file << std::endl;
            //     throw std::runtime_error("Cannot open model file for Lcurve::Model");
            // }
            // std::ifstream fin_data(data_file.c_str());
            // if (!fin_data) {
            //     std::cerr << "DEBUG: Tried to open file: " << data_file << std::endl;
            //     throw std::runtime_error("Cannot open data file for Lcurve::Data");
            // }

            model_ = std::make_unique<Lcurve::Model>(model_file);
            data_ = std::make_shared<Lcurve::Data>(data_file);


            // Scale factors buffer
            sfac_.resize(4);
            for (int i = 0; i < 4; i++)
                sfac_[i] = 1.0;
        }
        catch (const Lcurve::Lcurve_Error& e) {
            throw std::runtime_error(std::string("TRM Lcurve_Error during constructor: ") + e);
        }
        catch (const std::exception& e) {
            throw std::runtime_error(std::string("Standard exception during constructor: ") + e.what());
        }
        catch (...) {
            throw std::runtime_error("Unknown exception during constructor");
        }
    }

    // -------------------------
    // Custom copy constructor
    // -------------------------
    LCurveWrapper(const LCurveWrapper& other)
        : model_(std::make_unique<Lcurve::Model>(*other.model_)), // deep copy
          data_(other.data_),                                     // shared pointer
          sfac_(other.sfac_),                                     // copy buffer
          fit_(other.fit_),
          param_buffer_(other.param_buffer_)
    {
        // Nothing else needed
    }


    // -------------------------
    // Clone method for emcee
    // -------------------------
    std::unique_ptr<LCurveWrapper> clone() const {
        // Use copy constructor — NO file I/O
        return std::make_unique<LCurveWrapper>(*this);
    }


    Lcurve::Model& get_model() {
        return *model_;
    }


    // Update the values for the variable parameters in bulk
    // Provide a list of floats for each parameter (order matters)
    void set_params(const std::vector<double>& params) {
        try {
            if ((int)params.size() != model_->nvary())
                throw std::runtime_error("Parameter vector size does not match model nvary()");

            if (static_cast<size_t>(param_buffer_.size()) != params.size())
                param_buffer_.resize(params.size());
        
            for (size_t i = 0; i < params.size(); i++)
                param_buffer_[i] = params[i];
        
            if (model_->is_not_legal(param_buffer_))
                throw std::runtime_error("Attempted to set illegal parameter values");
        
            model_->set_param(param_buffer_);
        }
        catch (const Lcurve::Lcurve_Error& e) {
            throw std::runtime_error(std::string("TRM Lcurve_Error in set_params: ") + e);
        }
        catch (const std::exception& e) {
            throw std::runtime_error(std::string("Standard exception in set_params: ") + e.what());
        }
        catch (...) {
            throw std::runtime_error("Unknown exception in set_params");
        }
    }

    // Compute chi-squared through Lcurve::light_curve_comp()
    std::vector<double> compute_chisq() {
        py::gil_scoped_release release;
        try {
            double wdwarf, chisq, wnok, logg1, logg2, rv1, rv2;

            constexpr bool scale = true;        // autoscale
            constexpr bool do_copy = false;     // no copy data
            constexpr bool add_noise = false;   // no noise

            // Call TRM light curve computation
            Lcurve::light_curve_comp(*model_, *data_, scale, do_copy, add_noise,
                                     sfac_, fit_, wdwarf, chisq, wnok,
                                     logg1, logg2, rv1, rv2);

            return {chisq, logg2, logg1, wdwarf};
        }
        catch (const Lcurve::Lcurve_Error& e) {
            throw std::runtime_error(std::string("TRM Lcurve_Error in compute_chisq: ") + e);
        }
        catch (const std::exception& e) {
            throw std::runtime_error(std::string("Standard exception in compute_chisq: ") + e.what());
        }
        catch (...) {
            throw std::runtime_error("Unknown exception in compute_chisq");
        }
    }

    // Get the curent values of the variable parameters
    std::vector<double> get_params() {
        try {
            Subs::Array1D<double> v = model_->get_param();
            std::vector<double> vec(v.size());
            for (int i = 0; i < v.size(); i++)
                vec[i] = v[i];
            return vec;
        }
        catch (const Lcurve::Lcurve_Error& e) {
            throw std::runtime_error(std::string("TRM Lcurve_Error in get_params: ") + e);
        }
        catch (const std::exception& e) {
            throw std::runtime_error(std::string("Standard exception in get_params: ") + e.what());
        }
        catch (...) {
            throw std::runtime_error("Unknown exception in get_params");
        }
    }

    // Get the names of the variable parameters
    std::vector<std::string> get_param_names() {
        try {
            int n = model_->nvary();
            std::vector<std::string> names;
            for (int i = 0; i < n; i++)
                names.push_back(model_->get_name(i));
            return names;
        }
        catch (const Lcurve::Lcurve_Error& e) {
            throw std::runtime_error(std::string("TRM Lcurve_Error in get_param_names: ") + e);
        }
        catch (const std::exception& e) {
            throw std::runtime_error(std::string("Standard exception in get_param_names: ") + e.what());
        }
        catch (...) {
            throw std::runtime_error("Unknown exception in get_param_names");
        }
    }

    // C++ nvary() function to update how many parameters are set to vary.
    // Required to use the bulk set_params() method
    int nvary() const {
        return model_->nvary();
    }

};

// Define the stuff you want Python to access here.
PYBIND11_MODULE(lcurve_wrapper, m) {
    py::class_<Lcurve::Pparam>(m, "Pparam") // Lcurve::Model stores variables in a PParam. We want this and the ability to modify param.value, param.vary, param.defined
        .def_readwrite("value", &Lcurve::Pparam::value)
        .def_readwrite("vary", &Lcurve::Pparam::vary)
        .def_readwrite("defined", &Lcurve::Pparam::defined);

    py::class_<Lcurve::Model>(m, "Model") // We want the ability to modify all the parameters, do we include them all here explicitly.
        .def_readwrite("q", &Lcurve::Model::q)
        .def_readwrite("iangle", &Lcurve::Model::iangle)
        .def_readwrite("r1", &Lcurve::Model::r1)
        .def_readwrite("r2", &Lcurve::Model::r2)
        .def_readwrite("cphi3", &Lcurve::Model::cphi3)
        .def_readwrite("cphi4", &Lcurve::Model::cphi4)
        .def_readwrite("spin1", &Lcurve::Model::spin1)
        .def_readwrite("spin2", &Lcurve::Model::spin2)
        .def_readwrite("t1", &Lcurve::Model::t1)
        .def_readwrite("t2", &Lcurve::Model::t2)
        .def_readwrite("ldc1_1", &Lcurve::Model::ldc1_1)
        .def_readwrite("ldc1_2", &Lcurve::Model::ldc1_2)
        .def_readwrite("ldc1_3", &Lcurve::Model::ldc1_3)
        .def_readwrite("ldc1_4", &Lcurve::Model::ldc1_4)
        .def_readwrite("ldc2_1", &Lcurve::Model::ldc2_1)
        .def_readwrite("ldc2_2", &Lcurve::Model::ldc2_2)
        .def_readwrite("ldc2_3", &Lcurve::Model::ldc2_3)
        .def_readwrite("ldc2_4", &Lcurve::Model::ldc2_4)
        .def_readwrite("velocity_scale", &Lcurve::Model::velocity_scale)
        .def_readwrite("beam_factor1", &Lcurve::Model::beam_factor1)
        .def_readwrite("beam_factor2", &Lcurve::Model::beam_factor2)
        .def_readwrite("t0", &Lcurve::Model::t0)
        .def_readwrite("period", &Lcurve::Model::period)
        .def_readwrite("pdot", &Lcurve::Model::pdot)
        .def_readwrite("deltat", &Lcurve::Model::deltat)
        .def_readwrite("gravity_dark1", &Lcurve::Model::gravity_dark1)
        .def_readwrite("gravity_dark2", &Lcurve::Model::gravity_dark2)
        .def_readwrite("absorb", &Lcurve::Model::absorb)
        .def_readwrite("slope", &Lcurve::Model::slope)
        .def_readwrite("quad", &Lcurve::Model::quad)
        .def_readwrite("cube", &Lcurve::Model::cube)
        .def_readwrite("third", &Lcurve::Model::third)
        .def_readwrite("rdisc1", &Lcurve::Model::rdisc1)
        .def_readwrite("rdisc2", &Lcurve::Model::rdisc2)
        .def_readwrite("height_disc", &Lcurve::Model::height_disc)
        .def_readwrite("beta_disc", &Lcurve::Model::beta_disc)
        .def_readwrite("temp_disc", &Lcurve::Model::temp_disc)
        .def_readwrite("texp_disc", &Lcurve::Model::texp_disc)
        .def_readwrite("lin_limb_disc", &Lcurve::Model::lin_limb_disc)
        .def_readwrite("quad_limb_disc", &Lcurve::Model::quad_limb_disc)
        .def_readwrite("temp_edge", &Lcurve::Model::temp_edge)
        .def_readwrite("absorb_edge", &Lcurve::Model::absorb_edge)
        .def_readwrite("radius_spot", &Lcurve::Model::radius_spot)
        .def_readwrite("length_spot", &Lcurve::Model::length_spot)
        .def_readwrite("height_spot", &Lcurve::Model::height_spot)
        .def_readwrite("expon_spot", &Lcurve::Model::expon_spot)
        .def_readwrite("epow_spot", &Lcurve::Model::epow_spot)
        .def_readwrite("angle_spot", &Lcurve::Model::angle_spot)
        .def_readwrite("yaw_spot", &Lcurve::Model::yaw_spot)
        .def_readwrite("temp_spot", &Lcurve::Model::temp_spot)
        .def_readwrite("tilt_spot", &Lcurve::Model::tilt_spot)
        .def_readwrite("cfrac_spot", &Lcurve::Model::cfrac_spot)
        .def_readwrite("stsp11_long", &Lcurve::Model::stsp11_long)
        .def_readwrite("stsp11_lat", &Lcurve::Model::stsp11_lat)
        .def_readwrite("stsp11_fwhm", &Lcurve::Model::stsp11_fwhm)
        .def_readwrite("stsp11_tcen", &Lcurve::Model::stsp11_tcen)
        .def_readwrite("stsp12_long", &Lcurve::Model::stsp12_long)
        .def_readwrite("stsp12_lat", &Lcurve::Model::stsp12_lat)
        .def_readwrite("stsp12_fwhm", &Lcurve::Model::stsp12_fwhm)
        .def_readwrite("stsp12_tcen", &Lcurve::Model::stsp12_tcen)
        .def_readwrite("stsp13_long", &Lcurve::Model::stsp13_long)
        .def_readwrite("stsp13_lat", &Lcurve::Model::stsp13_lat)
        .def_readwrite("stsp13_fwhm", &Lcurve::Model::stsp13_fwhm)
        .def_readwrite("stsp13_tcen", &Lcurve::Model::stsp13_tcen)
        .def_readwrite("stsp21_long", &Lcurve::Model::stsp21_long)
        .def_readwrite("stsp21_lat", &Lcurve::Model::stsp21_lat)
        .def_readwrite("stsp21_fwhm", &Lcurve::Model::stsp21_fwhm)
        .def_readwrite("stsp21_tcen", &Lcurve::Model::stsp21_tcen)
        .def_readwrite("stsp22_long", &Lcurve::Model::stsp22_long)
        .def_readwrite("stsp22_lat", &Lcurve::Model::stsp22_lat)
        .def_readwrite("stsp22_fwhm", &Lcurve::Model::stsp22_fwhm)
        .def_readwrite("stsp22_tcen", &Lcurve::Model::stsp22_tcen)
        .def_readwrite("uesp_long1", &Lcurve::Model::uesp_long1)
        .def_readwrite("uesp_long2", &Lcurve::Model::uesp_long2)
        .def_readwrite("uesp_lathw", &Lcurve::Model::uesp_lathw)
        .def_readwrite("uesp_taper", &Lcurve::Model::uesp_taper)
        .def_readwrite("uesp_temp", &Lcurve::Model::uesp_temp)
        .def_readonly("add_spot", &Lcurve::Model::add_spot)
        .def_readonly("add_disc", &Lcurve::Model::add_disc);
    
    py::class_<LCurveWrapper>(m, "LCurveModel")
        .def(py::init<const std::string&, const std::string&>(),
             py::arg("model_file"),
             py::arg("data_file"))
        .def("set_params", &LCurveWrapper::set_params, py::arg("params"))
        .def("compute_chisq", &LCurveWrapper::compute_chisq)
        .def("get_params", &LCurveWrapper::get_params)
        .def("get_param_names", &LCurveWrapper::get_param_names)
        .def("get_model", &LCurveWrapper::get_model, py::return_value_policy::reference)
        .def("nvary", &LCurveWrapper::nvary)
        .def("clone", &LCurveWrapper::clone);
}

