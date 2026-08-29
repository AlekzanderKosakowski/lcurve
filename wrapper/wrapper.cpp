// wrapper.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "trm/lcurve.h"
#include "trm/subs.h"
#include "trm/input.h"


#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;


struct lroche_output {
    // Expected outputs from lroche (compute_light_curve() call)
    std::vector<double> model_flux;

    double chisq;
    double wnok;
    double wdwarf;

    double logg1;
    double logg2;

    double rvol1;
    double rvol2;

    double ffac1;
    double ffac2;

    std::vector<double> sfac;
};


std::vector<std::array<double, 3>> project_visible(
    const Subs::Buffer1D<Lcurve::Point>& object,
    const Subs::Vec3& earth,
    const Subs::Vec3& cofm,
    const Subs::Vec3& xsky,
    const Subs::Vec3& ysky,
    double phase)
{
    std::vector<std::array<double, 3>> points;

    for (int i = 0; i < object.size(); ++i) {

        if (Subs::dot(earth, object[i].dirn) > 0. &&
            object[i].visible(phase)) {

            Subs::Vec3 r = object[i].posn - cofm;

            double xpos = Subs::dot(r, xsky);
            double ypos = Subs::dot(r, ysky);

            points.push_back({
                xpos,
                ypos,
                object[i].temp
            });
        }
    }

    return points;
}


// This is a wrapper for very specific parts of lcurve to be accessed from python
// It includes the Lcurve::Model class, but not all of the methods inside
class LCurveWrapper {
public:
    std::unique_ptr<Lcurve::Model> model_;
    std::shared_ptr<Lcurve::Data> data_;
    Subs::Buffer1D<double> sfac_;

    Subs::Array1D<double> model_flux_;
    Subs::Array1D<double> param_buffer_;

    // Constructor: initialize model from a file
    LCurveWrapper(const std::string& model_file) {
        try {

            model_ = std::make_unique<Lcurve::Model>(model_file);

            sfac_.resize(5);
            for (int i = 0; i < 5; i++)
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
        : model_(std::make_unique<Lcurve::Model>(*other.model_)), // deep copy (gets unique model)
          data_(other.data_),                                     // shared data between cloned models
          sfac_(other.sfac_),
          model_flux_(other.model_flux_),
          param_buffer_(other.param_buffer_)
    {
        // Nothing else needed
    }


    // -------------------------
    // Clone method for emcee
    // -------------------------
    std::unique_ptr<LCurveWrapper> clone() const {
        // Use copy constructor (allows us to skip file I/O to build a new model object)
        return std::make_unique<LCurveWrapper>(*this);
    }


    // Update the model object's stored LCurve::Data object using NumPy inputs
    void set_data(const py::array_t<double, py::array::c_style | py::array::forcecast>& data){
        auto buf = data.request();
    
        if (buf.ndim != 2) // Input must be of shape (N, M)
            throw std::runtime_error(
                "set_data: data must be a 2-dimensional NumPy array"
            );
    
        if (buf.shape[1] != 6) // Each row must have 6 columns
            throw std::runtime_error(
                "set_data: data must have exactly 6 columns"
            );

        const ssize_t n = buf.shape[0];
    
        // Pointer to the first element of the contiguous NumPy array.
        const double* ptr = static_cast<const double*>(buf.ptr);
        
        auto new_data = std::make_shared<Lcurve::Data>(
            static_cast<int>(n)
        );
    
        for (ssize_t i = 0; i < n; ++i) {

            // Each row is:
            // time, exposure, flux, flux_error, weight, ndiv
            const double* row = ptr + i * 6;
    
            (*new_data)[i].time   = row[0];
            (*new_data)[i].expose = row[1];
            (*new_data)[i].flux   = row[2];
            (*new_data)[i].ferr   = row[3];
            (*new_data)[i].weight = row[4];

            if (row[5] < 1 || row[5] != std::floor(row[5]))
                throw std::runtime_error("Light curve data column 6 (ndiv) must contain positive integers");
            
            (*new_data)[i].ndiv = static_cast<int>(row[5]); // NumPy made this a float. Now we have C++ fix it back to integer
        }
    
        data_ = std::move(new_data); // Update the model's data_ object
    }
    
    
    Lcurve::Model& get_model() {
        return *model_;
    }

    // Write the model parameters to a new parameters.txt file.
    void write_to_parameters_file(const std::string& ofilename) {
        if (!model_) {
            throw std::runtime_error("model_ is null");
        }
        model_->wrasc(ofilename);    
    }

    // Validate that the model's parameters are legal.
    // See the original trm_lcurve.cc::is_not_legal()
    bool is_legal() const {

        if(model_->q.vary){
            if(model_->q.value < 0. || model_->q.value > 100.) return false;
        }

        if(model_->iangle.vary){
            if(model_->iangle.value < 0. || model_->iangle.value > 90.) return false;
        }
    
        if(model_->use_radii){
            

            if(model_->r1.vary){
                double xl11 = Roche::xl11(model_->q.value, model_->spin1.value);
                if(model_->r1.value <= 0. || model_->r1.value >= xl11) return false;
            }
            if(model_->r2.vary){
                double xl12 = 1.0 - Roche::xl12(model_->q.value, model_->spin2.value);
                if(model_->r2.value <= 0. || model_->r2.value > xl12) return false;
            }
        }else{
            if(model_->cphi3.vary){
                if(model_->cphi3.value <= 0. || model_->cphi3.value > 0.25) return false;
            }
            if(model_->cphi4.vary){
                if(model_->cphi4.value <= 0. || model_->cphi4.value > 0.25) return false;
            }
        }

        if(model_->r3.vary) {
            if(model_->r3.value < 0.) return false;
        }


        if(model_->spin1.vary) {
            if(model_->spin1.value <= 0. || model_->spin1.value > 1000.) return false;
        }
    
        if(model_->spin2.vary) {
            if(model_->spin2.value <= 0. || model_->spin2.value > 1000.) return false;
        }
    
        if(model_->t1.vary) {
            if(model_->t1.value <= 500. || model_->t1.value > 1.e6) return false;
        }
    
        if(model_->t2.vary) {
            if(model_->t2.value <= 500. || model_->t2.value > 1.e6) return false;
        }
    
        if(model_->t3.vary) {
            if(model_->t3.value <= 500. || model_->t3.value > 1.e6) return false;
        }



        if(model_->ldc1_1.vary){
            if(model_->ldc1_1.value < -2. || model_->ldc1_1.value > 2.) return false;
        }
    
        if(model_->ldc1_2.vary){
            if(model_->ldc1_2.value < -2. || model_->ldc1_2.value > 2.) return false;
        }
    
        if(model_->ldc1_3.vary){
            if(model_->ldc1_3.value < -2. || model_->ldc1_3.value > 2.) return false;
        }
    
        if(model_->ldc1_4.vary){
            if(model_->ldc1_4.value < -2. || model_->ldc1_4.value > 2.) return false;
        }
    
        if(model_->ldc2_1.vary){
            if(model_->ldc2_1.value < -2. || model_->ldc2_1.value > 2.) return false;
        }
    
        if(model_->ldc2_2.vary){
            if(model_->ldc2_2.value < -2. || model_->ldc2_2.value > 2.) return false;
        }
    
        if(model_->ldc2_3.vary){
            if(model_->ldc2_3.value < -2. || model_->ldc2_3.value > 2.) return false;
        }
    
        if(model_->ldc2_4.vary){
            if(model_->ldc2_4.value < -2. || model_->ldc2_4.value > 2.) return false;
        }
    
        if(model_->ldc3_1.vary){
            if(model_->ldc3_1.value < -2. || model_->ldc3_1.value > 2.) return false;
        }
    
        if(model_->ldc3_2.vary){
            if(model_->ldc3_2.value < -2. || model_->ldc3_2.value > 2.) return false;
        }
    
        if(model_->ldc3_3.vary){
            if(model_->ldc3_3.value < -2. || model_->ldc3_3.value > 2.) return false;
        }
    
        if(model_->ldc3_4.vary){
            if(model_->ldc3_4.value < -2. || model_->ldc3_4.value > 2.) return false;
        }

        if(model_->period.vary) {
            if(model_->period.value <= 1.e-3 || model_->period.value > 100.) return false;
        }

        if(model_->deltat.vary) {
            if(model_->deltat.value <= -1. || model_->deltat.value > 1.) return false;
        }

        if(model_->gravity_dark1.vary){
            if(model_->gravity_dark1.value < -1. || model_->gravity_dark1.value > 1.) return false;
        }
    
        if(model_->gravity_dark2.vary){
            if(model_->gravity_dark2.value < -1. || model_->gravity_dark2.value > 1.) return false;
        }
    
        if(model_->absorb.vary){
            if(model_->absorb.value < 0. || model_->absorb.value > 10.) return false;
        }
    
        if(model_->slope.vary){
            if(model_->slope.value < -2. || model_->slope.value > 2.) return false;
        }
    
        if(model_->quad.vary){
            if(model_->quad.value < -2. || model_->quad.value > 2.) return false;
        }
    
        if(model_->cube.vary){
            if(model_->cube.value < -2. || model_->cube.value > 2.) return false;
        }

        if(model_->add_disc){

            if(model_->rdisc1.vary){
                if(model_->rdisc1.value < 0. || model_->rdisc1.value > 1.) return false;
            }

            if(model_->rdisc2.vary){
                if(model_->rdisc2.value < 0. || model_->rdisc2.value > 1.) return false;
            }

            if(model_->height_disc.vary){
                if(model_->height_disc.value < 0.) return false;
            }

            if(model_->beta_disc.vary){
                if(model_->beta_disc.value < 1. || model_->beta_disc.value > 100.) return false;
            }

            if(model_->temp_disc.vary){
                if(model_->temp_disc.value < 500. || model_->temp_disc.value > 1.e6) return false;
            }

            if(model_->texp_disc.vary){
                if(model_->texp_disc.value < -100. || model_->texp_disc.value > 100.) return false;
            }

            if(model_->lin_limb_disc.vary){
                if(model_->lin_limb_disc.value < 0. || model_->lin_limb_disc.value > 1.) return false;
            }

            if(model_->quad_limb_disc.vary){
                if(model_->quad_limb_disc.value < 0. || model_->quad_limb_disc.value > 1.) return false;
            }

            if(model_->temp_edge.vary){
                if(model_->temp_edge.value <= 0. || model_->temp_edge.value > 1.e6) return false;
            }

            if(model_->absorb_edge.vary){
                if(model_->absorb_edge.value < 0. || model_->absorb_edge.value > 10.) return false;
            }
        }

        if(model_->add_spot){
    
            if(model_->radius_spot.vary){
                if(model_->radius_spot.value < 0. || model_->radius_spot.value > 1.) return false;
            }
    
            if(model_->length_spot.vary){
                if(model_->length_spot.value < 0. || model_->length_spot.value > 1.) return false;
            }
    
            if(model_->height_spot.vary){
                if(model_->height_spot.value < 0. || model_->height_spot.value > 1.) return false;
            }
    
            if(model_->expon_spot.vary){
                if(model_->expon_spot.value <= 0. || model_->expon_spot.value > 100.) return false;
            }
    
            if(model_->epow_spot.vary){
                if(model_->epow_spot.value <= 0. || model_->epow_spot.value > 10.) return false;
            }
    
            if(model_->angle_spot.vary){
                if(model_->angle_spot.value < -1000. || model_->angle_spot.value > 1000.) return false;
            }
    
            if(model_->yaw_spot.vary){
                if(model_->yaw_spot.value < -1000. || model_->yaw_spot.value > 1000.) return false;
            }
    
            if(model_->temp_spot.vary){
                if(model_->temp_spot.value < 0. || model_->temp_spot.value > 1.e10) return false;
            }
    
            if(model_->tilt_spot.vary){
                if(model_->tilt_spot.value < -1000. || model_->tilt_spot.value > 1000.) return false;
            }
    
            if(model_->cfrac_spot.vary){
                if(model_->cfrac_spot.value < 0. || model_->cfrac_spot.value > 1.) return false;
            }
        }
    
        if(model_->stsp11_long.defined && model_->stsp11_long.vary){
            if(model_->stsp11_long.value < -400. || model_->stsp11_long.value > 400.) return false;
        }
        if(model_->stsp11_lat.defined && model_->stsp11_lat.vary){
            if(model_->stsp11_lat.value < -90. || model_->stsp11_lat.value > 90.) return false;
        }
        if(model_->stsp11_fwhm.defined  && model_->stsp11_fwhm.vary){
            if(model_->stsp11_fwhm.value <= 0. || model_->stsp11_fwhm.value > 180.) return false;
        }
        if(model_->stsp11_tcen.defined  && model_->stsp11_tcen.vary){
            if(model_->stsp11_tcen.value <= 0.) return false;
        }

        if(model_->stsp12_long.defined && model_->stsp12_long.vary){
            if(model_->stsp12_long.value < -400. || model_->stsp12_long.value > 400.) return false;
        }
        if(model_->stsp12_lat.defined && model_->stsp12_lat.vary){
            if(model_->stsp12_lat.value < -90. || model_->stsp12_lat.value > 90.) return false;
        }
        if(model_->stsp12_fwhm.defined  && model_->stsp12_fwhm.vary){
            if(model_->stsp12_fwhm.value <= 0. || model_->stsp12_fwhm.value > 180.) return false;
        }
        if(model_->stsp12_tcen.defined  && model_->stsp12_tcen.vary){
            if(model_->stsp12_tcen.value <= 0.) return false;
        }

        if(model_->stsp13_long.defined && model_->stsp13_long.vary){
            if(model_->stsp13_long.value < -400. || model_->stsp13_long.value > 400.) return false;
        }
        if(model_->stsp13_lat.defined && model_->stsp13_lat.vary){
            if(model_->stsp13_lat.value < -90. || model_->stsp13_lat.value > 90.) return false;
        }
        if(model_->stsp13_fwhm.defined  && model_->stsp13_fwhm.vary){
            if(model_->stsp13_fwhm.value <= 0. || model_->stsp13_fwhm.value > 180.) return false;
        }
        if(model_->stsp13_tcen.defined  && model_->stsp13_tcen.vary){
            if(model_->stsp13_tcen.value <= 0.) return false;
        }
        
        if(model_->stsp21_long.defined && model_->stsp21_long.vary){
            if(model_->stsp21_long.value < -400. || model_->stsp21_long.value > 400.) return false;
        }
        if(model_->stsp21_lat.defined && model_->stsp21_lat.vary){
            if(model_->stsp21_lat.value < -90. || model_->stsp21_lat.value > 90.) return false;
        }
        if(model_->stsp21_fwhm.defined  && model_->stsp21_fwhm.vary){
            if(model_->stsp21_fwhm.value <= 0. || model_->stsp21_fwhm.value > 180.) return false;
        }
        if(model_->stsp21_tcen.defined  && model_->stsp21_tcen.vary){
            if(model_->stsp21_tcen.value <= 0.) return false;
        }

        if(model_->stsp22_long.defined && model_->stsp22_long.vary){
            if(model_->stsp22_long.value < -400. || model_->stsp22_long.value > 400.) return false;
        }
        if(model_->stsp22_lat.defined && model_->stsp22_lat.vary){
            if(model_->stsp22_lat.value < -90. || model_->stsp22_lat.value > 90.) return false;
        }
        if(model_->stsp22_fwhm.defined  && model_->stsp22_fwhm.vary){
            if(model_->stsp22_fwhm.value <= 0. || model_->stsp22_fwhm.value > 180.) return false;
        }
        if(model_->stsp22_tcen.defined  && model_->stsp22_tcen.vary){
            if(model_->stsp22_tcen.value <= 0.) return false;
        }

        if(model_->stsp1i_long.defined && model_->stsp1i_long.vary) {
            if(model_->stsp1i_long.value < -400. || model_->stsp1i_long.value > 400.) return false;
        }

        if(model_->stsp1i_lat.defined && model_->stsp1i_lat.vary) {
            if(model_->stsp1i_lat.value < -90. || model_->stsp1i_lat.value > 90.) return false;
        }

        if(model_->stsp1i_fwhm_long.defined && model_->stsp1i_fwhm_long.vary) {
            if(model_->stsp1i_fwhm_long.value < 0. || model_->stsp1i_fwhm_long.value > 180.) return false;
        }

        if(model_->stsp1i_fwhm_lat.defined && model_->stsp1i_fwhm_lat.vary) {
            if(model_->stsp1i_fwhm_lat.value <= 0. || model_->stsp1i_fwhm_lat.value > 180.) return false;
        }

        if(model_->stsp1i_escale_len.defined && model_->stsp1i_escale_len.vary) {
            if(model_->stsp1i_escale_len.value < 0. || model_->stsp1i_escale_len.value > 400.) return false;
        }

        if(model_->stsp1i_tcen.defined && model_->stsp1i_tcen.vary) {
            if(model_->stsp1i_tcen.value <= 0.) return false;
        }

        return true;
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

    // Get the current values of the variable parameters
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

    // Give LCurveWrapper object a visualise command at a specific phase. test
    py::dict visualise_data(double phase)
    {
        // Prepare radii
        double r1, r2, rdisc1=0., rdisc2=0.;
        model_->get_r1r2(r1, r2);
        double rl1 = Roche::xl11(model_->q,model_->spin1);
        if(r1 <= 0){
            r1 = 0.99999999999*rl1;
        }
        double rl2 = 1.-Roche::xl12(model_->q,model_->spin2);
        if(r2 <= 0){
            r2 = 0.99999999999*rl2;
        }
        
        // One object per plotted component.
        Subs::Buffer1D<Lcurve::Point> star1, star2, disc, outer_edge, inner_edge, bspot, stream;

        // Generate grids for each component
        Lcurve::set_star_grid(*model_, Roche::PRIMARY, true, star1);
        Lcurve::set_star_grid(*model_, Roche::SECONDARY, true, star2);

        // --------------------------------------------------
        // Temperature grid for filter integration
        // --------------------------------------------------
        double temperature_grid_min = 100.0;
        double temperature_grid_max = 100000.0;
        double temperature_grid_step = 200.0;
    
        int N_temperatures = static_cast<int>(
            std::ceil(
                (temperature_grid_max - temperature_grid_min) / temperature_grid_step)
            ) + 1;
    
        std::vector<double> temperature_array(N_temperatures);
        std::vector<double> planck_array(N_temperatures);
    
        for (int i = 0; i < N_temperatures; ++i){
            temperature_array[i] = temperature_grid_min + temperature_grid_step * i;
    
            if (temperature_array[i] > temperature_grid_max){
                temperature_array[i] = temperature_grid_max;
            }
            planck_array[i] = 0.0;
        }
    
        bool integrate_filter = !(
            Subs::tolower(model_->filter) == "none" ||
            Subs::tolower(model_->filter) == "false" ||
            Subs::tolower(model_->filter) == "n" ||
            Subs::tolower(model_->filter) == "no" ||
            Subs::tolower(model_->filter) == "0"
        );
    
        if (integrate_filter){
            Subs::integrate_filter(temperature_array, planck_array, model_->filter);
        }
    
        Lcurve::LDC ldc1 = model_->get_ldc1();
        Lcurve::LDC ldc2 = model_->get_ldc2();
    
        set_star_continuum(*model_, star1, star2, integrate_filter, temperature_array, planck_array, ldc1, ldc2);

        if(model_->add_disc){

            rdisc1 = model_->rdisc1 > 0. ? model_->rdisc1 : r1;
            rdisc2 = model_->rdisc2 > 0. ? model_->rdisc2 : model_->radius_spot;

            Lcurve::set_disc_grid(*model_, disc);
            Lcurve::set_disc_edge(*model_, true, outer_edge);
            Lcurve::set_disc_edge(*model_, false, inner_edge);


            std::vector<std::pair<double,double> > eclipses;
            if(model_->opaque){
	    
                // Apply eclipse by disc to star 1
                for(int i=0; i<star1.size(); i++){
                    eclipses =  Roche::disc_eclipse(model_->iangle, rdisc1, rdisc2, model_->beta_disc, model_->height_disc, star1[i].posn);
                    for(size_t j=0; j<eclipses.size(); j++)
                        star1[i].eclipse.push_back(eclipses[j]);
                }
	    
                // Apply eclipse by disc to star 2
                for(int i=0; i<star2.size(); i++){
                    eclipses =  Roche::disc_eclipse(model_->iangle, rdisc1, rdisc2, model_->beta_disc, model_->height_disc, star2[i].posn);
                    for(size_t j=0; j<eclipses.size(); j++)
                        star2[i].eclipse.push_back(eclipses[j]);
                }
            }

            // Set the surface brightness of the disc
            set_disc_continuum(rdisc2, model_->temp_disc, model_->texp_disc, model_->wavelength, disc, integrate_filter, temperature_array, planck_array);

            // Set the surface brightness of outer edge, accounting for irradiation by star 2
            set_edge_continuum(model_->temp_edge, r2, std::abs(model_->t2), model_->absorb_edge, model_->wavelength, outer_edge, integrate_filter, temperature_array, planck_array);
        }

        if(model_->add_spot){
            Subs::Vec3 dir(1,0,0), posn, v;
            Lcurve::set_bright_spot_grid(*model_, bspot, integrate_filter, temperature_array, planck_array);

            double rl1_spot = Roche::xl1(model_->q);
      
            // Calculate a reference radius and potential for the two stars
            double rref1, pref1, ffac1 = r1/rl1_spot;
            Roche::ref_sphere(model_->q, Roche::PRIMARY, model_->spin1, ffac1, rref1, pref1);
      
            double rref2, pref2, ffac2 = r2/rl2;
            Roche::ref_sphere(model_->q, Roche::SECONDARY, model_->spin2, ffac2, rref2, pref2);
      
            dir.set(0,0,1);
            Roche::strinit(model_->q, posn, v);
      
            Lcurve::Point::etype eclipses, edisc;
            Lcurve::star_eclipse(model_->q, r1, model_->spin1, ffac1, model_->iangle, posn, model_->delta_phase, model_->roche1, Roche::PRIMARY,   eclipses);
            Lcurve::star_eclipse(model_->q, r2, model_->spin2, ffac2, model_->iangle, posn, model_->delta_phase, model_->roche2, Roche::SECONDARY, eclipses);
            stream.push_back(Lcurve::Point(posn, dir, 0., 1., eclipses));
		       
            const int NSTREAM = int((rl1_spot-model_->radius_spot)/0.001);
            double radius;
            for(int i=0; i<NSTREAM; i++){
                radius = rl1_spot + (model_->radius_spot-rl1_spot)*(i+1)/NSTREAM;
                Roche::stradv(model_->q, posn, v, radius, 1.e-10, 1.e-3);
                eclipses.clear();
                Lcurve::star_eclipse(model_->q, r1, model_->spin1, ffac1, model_->iangle, posn, model_->delta_phase, model_->roche1, Roche::PRIMARY,   eclipses);
                Lcurve::star_eclipse(model_->q, r2, model_->spin2, ffac2, model_->iangle, posn, model_->delta_phase, model_->roche2, Roche::SECONDARY, eclipses);
                if(model_->add_disc){
                    edisc = Roche::disc_eclipse(model_->iangle, rdisc1, rdisc2, model_->beta_disc, model_->height_disc, posn);
                    for(size_t j=0; j<edisc.size(); j++)
                        eclipses.push_back(edisc[j]);
                }
                stream.push_back(Lcurve::Point(posn, dir, 0., 1., eclipses));
            }	
        }

        // --------------------------------------------------
        // Geometry -- copied from legacy visualise
        // --------------------------------------------------
    
        Subs::Vec3 cofm(
            model_->q / (1. + model_->q),
            0.,
            0.
        );
    
        Subs::Vec3 earth = Roche::set_earth(model_->iangle, phase);

        double cosp =std::cos(Constants::TWOPI * phase);
        double sinp = std::sin(Constants::TWOPI * phase);
    
        Subs::Vec3 xsky(sinp, cosp, 0.);
        Subs::Vec3 ysky = Subs::cross(earth, xsky);
    
        // --------------------------------------------------
        // Project visible points
        // --------------------------------------------------
    
        auto points1 = project_visible(star1, earth, cofm, xsky, ysky, phase);    
        auto points2 = project_visible(star2, earth, cofm, xsky, ysky, phase);
        auto points_bspot = project_visible(bspot, earth, cofm, xsky, ysky, phase);
        auto points_disc = project_visible(disc, earth, cofm, xsky, ysky, phase);
        auto points_outer_edge = project_visible(outer_edge, earth, cofm, xsky, ysky, phase);
        auto points_stream = project_visible(stream, earth, cofm, xsky, ysky, phase);


        // --------------------------------------------------
        // Return to Python
        // --------------------------------------------------
    
        py::dict result;

        // Output order used for generic Python dict.items() loops.
        // Alignment will place the later output points on top of the first where overlapping.
        result["stream"] = points_stream;
        result["outer_edge"] = points_outer_edge;
        result["bspot"] = points_bspot;
        result["disc"] = points_disc;
        result["star2"] = points2;
        result["star1"] = points1;
    
        return result;
    }

    lroche_output calculate_light_curve(const Lcurve::Data& data, bool scale){
        if (!model_)
            throw std::runtime_error("Model is null");
    
        if (data.empty())
            throw std::runtime_error("Light curve data is empty");
    
        lroche_output result; // Define the data type of "result"
        
        try {
            double wdwarf;
            double chisq;
            double wnok;
            double logg1;
            double logg2;
            double rvol1;
            double rvol2;
            double ffac1;
            double ffac2;
    
            constexpr bool rdata = true;
            constexpr bool info = false;

            // Reset scale factors.
            for (int i = 0; i < 5; ++i){
                sfac_[i] = 1.0;
            }
            
            {
                // Release GIL since we aren't doing Python work for a bit.
                py::gil_scoped_release release;

                Lcurve::light_curve_comp(*model_, data, scale, rdata, info, sfac_, model_flux_,
                    wdwarf, chisq, wnok, logg1, logg2, rvol1, rvol2, ffac1, ffac2);
            }
            
            // Copy the calculated light curve out of model_flux_
            // into an ordinary C++ vector.
            result.model_flux.resize(model_flux_.size());
    
            for (int i = 0; i < model_flux_.size(); ++i)
                result.model_flux[i] = model_flux_[i];

            // Copy scale factors.
            result.sfac.resize(5);
            for (int i = 0; i < 5; ++i)
                result.sfac[i] = sfac_[i];
            
            result.chisq  = chisq; // TODO: Add compute_chisq() method to skip unnecessary outputs for emcee
            result.wnok   = wnok;
            result.wdwarf = wdwarf;
    
            result.logg1 = logg1;
            result.logg2 = logg2;
    
            result.rvol1 = rvol1;
            result.rvol2 = rvol2;
    
            result.ffac1 = ffac1;
            result.ffac2 = ffac2;
    
        }
        catch (const Lcurve::Lcurve_Error& e) {
            throw std::runtime_error(std::string("TRM Lcurve_Error in calculate_light_curve: ") + e);
        }
        catch (const std::exception& e) {
            throw std::runtime_error(std::string("Standard exception in calculate_light_curve: ") + e.what());
        }
        catch (...) {
            throw std::runtime_error("Unknown exception in calculate_light_curve");
        }
    
        return result;
    }

    py::dict compute_light_curve(bool scale)
    {
        if (!data_)
            throw std::runtime_error("No cached data set has been loaded. Use model.set_data(data_array)");
    
        lroche_output result = calculate_light_curve(*data_, scale);
        
        return make_lroche_output_dict(result);
    }

    py::dict compute_light_curve(const py::array_t<double, py::array::c_style | py::array::forcecast>& input_data, bool scale){
        auto buf = input_data.request();
    
        if (buf.ndim != 2)
            throw std::runtime_error("Light curve data must be a 2D NumPy array");
    
        if (buf.shape[1] != 6)
            throw std::runtime_error("Light curve data must have 6 columns");
    
        const size_t ndata = buf.shape[0];
    
        const double* ptr = static_cast<const double*>(buf.ptr);
    
        Lcurve::Data data(ndata);
    
        for (size_t i = 0; i < ndata; ++i) {
    
            const double* row = ptr + i * 6;
    
            data[i].time   = row[0];
            data[i].expose = row[1];
            data[i].flux   = row[2];
            data[i].ferr   = row[3];
            data[i].weight = row[4];
            data[i].ndiv   = static_cast<int>(row[5]);
        }

        lroche_output result = calculate_light_curve(data, scale);
    
        return make_lroche_output_dict(result);
    }


    py::dict make_lroche_output_dict(const lroche_output& result){
        py::dict output;
    
        output["model_flux"] = result.model_flux;
        output["chisq"]      = result.chisq;
        output["wnok"]       = result.wnok;
        output["wdwarf"]     = result.wdwarf;
    
        output["logg1"]      = result.logg1;
        output["logg2"]      = result.logg2;
    
        output["rvol1"]        = result.rvol1;
        output["rvol2"]        = result.rvol2;
    
        output["ffac1"]      = result.ffac1;
        output["ffac2"]      = result.ffac2;

        output["sfac"]       = result.sfac;
        
        return output;
    }

};

// Define the stuff you want Python to access here.
PYBIND11_MODULE(lcurve_wrapper, m) {
    py::class_<Lcurve::Pparam>(m, "Pparam") // Lcurve::Model stores variables in a PParam. We want this and the ability to modify param.value, param.vary, param.defined
        .def_readwrite("value", &Lcurve::Pparam::value)
        .def_readwrite("range", &Lcurve::Pparam::range)
        .def_readwrite("dstep", &Lcurve::Pparam::dstep)
        .def_readwrite("vary", &Lcurve::Pparam::vary)
        .def_readwrite("defined", &Lcurve::Pparam::defined);

    py::class_<Lcurve::Model>(m, "Model") // We want the ability to modify all the parameters, do we include them all here explicitly.
        .def_readwrite("q", &Lcurve::Model::q)
        .def_readwrite("iangle", &Lcurve::Model::iangle)
        .def_readwrite("r1", &Lcurve::Model::r1)
        .def_readwrite("r2", &Lcurve::Model::r2)
        .def_readwrite("r3", &Lcurve::Model::r3)
        .def_readwrite("cphi3", &Lcurve::Model::cphi3)
        .def_readwrite("cphi4", &Lcurve::Model::cphi4)
        .def_readwrite("spin1", &Lcurve::Model::spin1)
        .def_readwrite("spin2", &Lcurve::Model::spin2)
        .def_readwrite("t1", &Lcurve::Model::t1)
        .def_readwrite("t2", &Lcurve::Model::t2)
        .def_readwrite("t3", &Lcurve::Model::t3)
        .def_readwrite("ldc1_1", &Lcurve::Model::ldc1_1)
        .def_readwrite("ldc1_2", &Lcurve::Model::ldc1_2)
        .def_readwrite("ldc1_3", &Lcurve::Model::ldc1_3)
        .def_readwrite("ldc1_4", &Lcurve::Model::ldc1_4)
        .def_readwrite("ldc2_1", &Lcurve::Model::ldc2_1)
        .def_readwrite("ldc2_2", &Lcurve::Model::ldc2_2)
        .def_readwrite("ldc2_3", &Lcurve::Model::ldc2_3)
        .def_readwrite("ldc2_4", &Lcurve::Model::ldc2_4)
        .def_readwrite("ldc3_1", &Lcurve::Model::ldc3_1)
        .def_readwrite("ldc3_2", &Lcurve::Model::ldc3_2)
        .def_readwrite("ldc3_3", &Lcurve::Model::ldc3_3)
        .def_readwrite("ldc3_4", &Lcurve::Model::ldc3_4)
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
        .def_readwrite("stsp1i_long", &Lcurve::Model::stsp1i_long)
        .def_readwrite("stsp1i_lat", &Lcurve::Model::stsp1i_lat)
        .def_readwrite("stsp1i_fwhm_long", &Lcurve::Model::stsp1i_fwhm_long)
        .def_readwrite("stsp1i_fwhm_lat", &Lcurve::Model::stsp1i_fwhm_lat)
        .def_readwrite("stsp1i_escale_len", &Lcurve::Model::stsp1i_escale_len)
        .def_readwrite("stsp1i_tcen", &Lcurve::Model::stsp1i_tcen)
        .def_readwrite("tperiod", &Lcurve::Model::tperiod)
        .def_readwrite("eclipse1", &Lcurve::Model::eclipse1)
        .def_readwrite("eclipse2", &Lcurve::Model::eclipse2)
        .def_readwrite("roche1", &Lcurve::Model::roche1)
        .def_readwrite("roche2", &Lcurve::Model::roche2)
        .def_readwrite("nlat1f", &Lcurve::Model::nlat1f)
        .def_readwrite("nlat2f", &Lcurve::Model::nlat2f)
        .def_readwrite("nlat1c", &Lcurve::Model::nlat1c)
        .def_readwrite("nlat2c", &Lcurve::Model::nlat2c)
        .def_readwrite("nrad", &Lcurve::Model::nrad)
        .def_readwrite("nspot", &Lcurve::Model::nspot)
        .def_readwrite("add_spot", &Lcurve::Model::add_spot)
        .def_readwrite("add_disc", &Lcurve::Model::add_disc)
        .def_readwrite("opaque", &Lcurve::Model::opaque)
        .def_readwrite("filter", &Lcurve::Model::filter)
        .def_readwrite("wavelength", &Lcurve::Model::wavelength)
        .def_readwrite("finite_irr12", &Lcurve::Model::finite_irr12)
        .def_readwrite("finite_irr21", &Lcurve::Model::finite_irr21)
        .def_readwrite("third", &Lcurve::Model::third);

    py::class_<LCurveWrapper>(m, "LCurveModel")
        .def(
            py::init<const std::string&>(),
            py::arg("model_file")        
        )
        .def("set_data",
            &LCurveWrapper::set_data,
            py::arg("data")
        )
        .def("set_params",
            &LCurveWrapper::set_params,
            py::arg("params")
        )
        .def("is_legal",
            &LCurveWrapper::is_legal
        )
        .def("get_params",
            &LCurveWrapper::get_params
        )
        .def("get_param_names",
            &LCurveWrapper::get_param_names
        )
        .def("get_model",
            &LCurveWrapper::get_model, 
            py::return_value_policy::reference
        )
        .def("nvary",
            &LCurveWrapper::nvary
        )
        .def("clone",
            &LCurveWrapper::clone
        )
        .def("write_to_parameters_file",
            &LCurveWrapper::write_to_parameters_file
        )
        .def("visualize_data",
            &LCurveWrapper::visualise_data,
            py::arg("phase")
        )
        .def(
            "compute_light_curve",
            py::overload_cast<bool>(
                &LCurveWrapper::compute_light_curve
            ),
            py::arg("scale")=true
        )                
        .def(
            "compute_light_curve",
            py::overload_cast<const py::array_t<double, py::array::c_style | py::array::forcecast>&,bool>(
                &LCurveWrapper::compute_light_curve
            ),
            py::arg("data"),
            py::arg("scale")=true
        );
            
}