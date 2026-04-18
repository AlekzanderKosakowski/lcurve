/*

!!begin
!!title    Plots the appearance of Roche-distorted stars in a binary
!!author   T.R.Marsh 
!!created  17 Sep 2003
!!revised  29 May 2011
!!descr    Plots the appearance of Roche-distorted stars in a binary
!!index    visualise
!!root     visualise
!!css      style.css
!!class    Model
!!head1    visualise -- plots the appearance of Roche-distorted stars in a binary

!!emph{visualise} computes grids of points representing two stars and optionally a disc
and plots those visible at a particular phase projected onto the plane of the sky.
It can only handle stars less than or equal to their Roche lobes. This allows visualisation
and checking of the eclipse computations.

!!head2 Command invocation

visualise model nphase nphase (phase)/(phase1 phase2) device x1 x2 y1 y2

!!head2 Arguments

!!table
!!arg{model}{Parameter file defining the model, as used e.g. by !!ref{lroche.html}{lrcohe}. }
!!arg{nphase}{Number of orbital phases to display.}
!!arg{phase}{Orbital phase if nphase=1}
!!arg{phase1}{First orbital phase if nphase>1}
!!arg{phase2}{Last orbital phase if nphase>1}
!!arg{device}{Plot device}
!!arg{x1}{left-hand X limit}
!!arg{x2}{right-hand X limit}
!!arg{y1}{lower Y limit}
!!arg{y2}{upper Y limit}
!!arg{width}{width of plot in inches}
!!arg{colormap}{select from predefined sequential colormaps}
!!arg{colorscale}{select from log or linear scaling of flux -> colormap}
!!arg{reverse}{reverse the chosen colormap}
!!arg{ncolors}{number of colors used to represent data}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include "cpgplot.h"
#include "trm/subs.h"
#include "trm/plot.h"
#include "trm/vec3.h"
#include "trm/input.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

int    Lcurve::Fobj::neval = 0;
double Lcurve::Fobj::chisq_min;
Subs::Buffer1D<double> Lcurve::Fobj::scale_min;

// Polynomial approximation of various colormaps
void set_colormap(std::string colormap, bool reverse, int ncolors) {

    std::string cmap = Subs::tolower(colormap);
    
    for(int i = 0; i < ncolors; i++){

        float t = reverse ? 1.0f - i/float(ncolors-1) : i/float(ncolors - 1);

        float r, g, b;

        // Sequential colormaps
        if(cmap=="viridis"){
            r = 0.2f + 0.8f * std::pow(t, 1.5f);
            g = 0.3f + 0.7f * std::pow(t, 1.0f);
            b = 0.5f - 0.5f * std::pow(1.0f - t, 1.2f);
        }else if(cmap=="inferno"){
            r = std::pow(t, 0.5f);
            g = std::pow(t, 1.5f);
            b = std::pow(t, 3.0f);
        }else if(cmap=="magma"){
            r = 1.0f * std::pow(t, 0.8f);
            g = 0.6f * std::pow(t, 1.4f);
            b = 0.3f * std::pow(t, 2.0f);
        }else if(cmap=="plasma"){
            r = 0.1f + 0.9f * std::pow(t, 1.2f);
            g = 0.0f + 0.8f * std::pow(t, 1.5f);
            b = 0.4f + 0.6f * std::pow(1.0f - t, 0.8f);
        }else if(cmap=="cividis"){
            r = 0.0f + 1.0f * std::pow(t, 0.9f);
            g = 0.1f + 0.9f * std::pow(t, 1.2f);
            b = 0.3f + 0.7f * std::pow(t, 0.8f);
            
        // Diverging colormaps
        }else if(cmap=="vanimo"){
            if(t < 0.5f){ // purple to white
                float u = t*2.0f;
                r = 0.4f + 0.6f*u;
                g = 0.0f + 1.0f*u;
                b = 0.6f + 0.4f*u;
            }else{ // white to green
                float u = (t-0.5f)*2.0f;
                r = 1.0f - 1.0f*u;
                g = 1.0f;
                b = 1.0f - 1.0f*u;
            }
        }else if(cmap=="seismic"){
            if(t < 0.5f){ // blue to white
                float u = t*2.0f;
                r = std::pow(u, 1.0f);
                g = std::pow(u, 1.0f);
                b = 1.0f;
            }else{ // white to red
                float u = (t-0.5f)*2.0f;
                r = 1.0f;
                g = 1.0f - std::pow(u, 1.0f);
                b = 1.0f - std::pow(u, 1.0f);
            }
            
        // Other colormaps
        }else if(cmap=="redblue"){
            r = t;
            g = 0.0;
            b = 1.0 - t;
            
        // Single color colormaps
        }else if(cmap=="black"){
            r = 0.0f;
            g = 0.0f;
            b = 0.0f;
        }
        
        r = std::min(1.0f, std::max(0.0f, r));
        g = std::min(1.0f, std::max(0.0f, g));
        b = std::min(1.0f, std::max(0.0f, b));
        
        cpgscr(16 + i, r, g, b); 
    }
}



// Main program
int main(int argc, char* argv[]){

    // Defined at the end
    void plot_visible(const Subs::Buffer1D<Lcurve::Point>& object, const Subs::Vec3& earth, const Subs::Vec3& cofm, const Subs::Vec3& xsky, const Subs::Vec3& ysky, double phase, double fmin, double fmax, std::string colorscale, int ncolors, int plt_marker);  
  
    try{
    
        // Construct Input object
        Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);
    
        // sign-in input variables
        input.sign_in("model",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
        input.sign_in("nphase",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("phase",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("phase1",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("phase2",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("device",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
        input.sign_in("x1",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("x2",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("y1",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("y2",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("colormap", Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("colorscale", Subs::Input::LOCAL,  Subs::Input::PROMPT); // log,linear
        input.sign_in("reverse", Subs::Input::LOCAL,  Subs::Input::PROMPT);    // false,true
        input.sign_in("ncolors", Subs::Input::LOCAL,  Subs::Input::PROMPT);    // < 240 (higher is not always better)
        input.sign_in("width",    Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
        input.sign_in("marker",    Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    
        std::string smodel;
        input.get_value("model", smodel, "model", "model file of parameter values");
        Lcurve::Model model(smodel);
    
        int nphase;
        input.get_value("nphase", nphase, 1, 1, 1000, "number of orbital phases");
        double phase, phase1, phase2;
        if(nphase == 1){
            input.get_value("phase", phase, 0., -10., 10., "orbital phase of interest");
        }else{
            input.get_value("phase1", phase1, 0., -10., 10., "first orbital phase of interest");
            input.get_value("phase2", phase2, 0., -10., 10., "last orbital phase of interest");
        }

        std::string device;
        input.get_value("device", device, "/xs", "plot device");
        float x1;
        input.get_value("x1", x1, -2.f, -100.f, 100.f, "left-hand X limit");
        float x2;
        input.get_value("x2", x2,  2.f, -100.f, 100.f, "right-hand X limit");
        float y1;
        input.get_value("y1", y1, -2.f, -100.f, 100.f, "lower Y limit");
        float y2;
        input.get_value("y2", y2,  2.f, -100.f, 100.f, "upper Y limit");
        std::string colormap;
        input.get_value("colormap", colormap, "inferno", "[viridis, inferno, magma, plasma, cividis, seismic, vanimo, redblue, black]");
        std::string colorscale;
        input.get_value("colorscale", colorscale, "log", "[linear, log]");
        bool reverse;
        input.get_value("reverse", reverse, false, "[yes, no]");
        float ncols;
        input.get_value("ncolors", ncols, 32.f, 16.f, 239.f, "color grid resolution (< 240)");
        float marker;
        input.get_value("marker", marker,  1.f, -1.f, 20.f, "surface element marker style");
        float width;
        input.get_value("width", width,  8.f, 0.f, 100.f, "width of the plot (inches)");
        
        input.save();

        int ncolors = std::floor(ncols); // Force integer ncolors
        int plt_marker = static_cast<int>(marker);
        
        // Complain if the user provides an invalid colormap
        if (colormap != "viridis" &&
            colormap != "inferno" &&
            colormap != "magma" &&
            colormap != "plasma" &&
            colormap != "cividis" &&
            colormap != "vanimo" &&
            colormap != "seismic" &&
            colormap != "redblue" &&
            colormap != "black") {
            std::cerr << "Invalid colormap. Try one of [viridis, inferno, magma, plasma, cividis, seismic, vanimo, redblue, black]" << std::endl;
            exit(1);
        }
        
        // Complain if the user provides an invalid colorscale
        if (colorscale != "linear" &&
            colorscale != "log") {
            std::cerr << "Invalid colorscale. Try one of [linear, log]" << std::endl;
            exit(1);
        }

        // Complain if the user provides an invalid number of colors
        if (ncolors >= 240) {
            std::cerr << "Color grid resolution is too high. Use ncolors < 240" << std::endl;
            exit(1);
        }

        double r1, r2, rdisc1=0., rdisc2=0.;
        model.get_r1r2(r1, r2);
        double rl1 = Roche::xl11(model.q,model.spin1);
        
        if(r1 <= 0){
            r1 = 0.99999999999*rl1;
        }
        
        double rl2 = 1.-Roche::xl12(model.q,model.spin2);
        if(r2 <= 0){
            r2 = 0.99999999999*rl2;
        }
    
        // Generate arrays over each star's face. 
        Subs::Buffer1D<Lcurve::Point> star1, star2, disc, outer_edge, inner_edge, bspot, stream;
        Lcurve::set_star_grid(model, Roche::PRIMARY, true, star1);
        Lcurve::set_star_grid(model, Roche::SECONDARY, true, star2);

        // Apply star continuum for coloring.
        ////////////////////////////////////////////////////////////////////////////
        double temperature_grid_min = 100.0;
        double temperature_grid_max = 100000.0;
        double temperature_grid_step = 200.0;
        int N_temperatures = static_cast<int>(std::ceil((temperature_grid_max - temperature_grid_min)/temperature_grid_step)) + 1;
    
        std::vector<double> temperature_array(N_temperatures);
        std::vector<double> planck_array(N_temperatures);
        for(int i=0; i < N_temperatures; ++i){
            temperature_array[i] = temperature_grid_min + temperature_grid_step*i;
            if (temperature_array[i] > temperature_grid_max) {
                temperature_array[i] = temperature_grid_max;
            }
            planck_array[i] = 0.0;
        }
    
        bool integrate_filter = !(
                                Subs::tolower(model.filter) == "none" ||
                                Subs::tolower(model.filter) == "false" ||
                                Subs::tolower(model.filter) == "n" ||
                                Subs::tolower(model.filter) == "no" ||
                                Subs::tolower(model.filter) == "0"
                                );
        
        if(integrate_filter){
            Subs::integrate_filter(temperature_array, planck_array, model.filter);
        }

        Lcurve::LDC ldc1 = model.get_ldc1();
        Lcurve::LDC ldc2 = model.get_ldc2();
        
        set_star_continuum(model, star1, star2, integrate_filter, temperature_array, planck_array, ldc1, ldc2);
        ////////////////////////////////////////////////////////////////////////////
        // Find the min and max temperature across the entire system.
        // Using "flux" adds an artificial gradient from equator to pole due to grid element area
        // TODO: Add more loops for disc and accretion spot.
        double min_stat=1e20, max_stat=1e-20;
        for(int i=0; i < star1.size(); i++){
            if(star1[i].temp < min_stat){min_stat = star1[i].temp;}
            if(star1[i].temp > max_stat){max_stat = star1[i].temp;}
        }
        for(int i=0; i < star2.size(); i++){
            if(star2[i].temp < min_stat){min_stat = star2[i].temp;}
            if(star2[i].temp > max_stat){max_stat = star2[i].temp;}
        }
        
        
        if(model.add_disc){

            rdisc1 = model.rdisc1 > 0. ? model.rdisc1 : r1;
            rdisc2 = model.rdisc2 > 0. ? model.rdisc2 : model.radius_spot;

            // note that the inner radius of the disc is set equal to that of the white dwarf if rdisc1 <= 0
            // while the outer disc is set equal to the spot radius
            Lcurve::set_disc_grid(model, disc);
            Lcurve::set_disc_edge(model, true, outer_edge);
            Lcurve::set_disc_edge(model, false, inner_edge);
	
            std::vector<std::pair<double,double> > eclipses;

            if(model.opaque){
	    
                // Apply eclipse by disc to star 1
                for(int i=0; i<star1.size(); i++){
                    eclipses =  Roche::disc_eclipse(model.iangle, rdisc1, rdisc2, model.beta_disc, model.height_disc, star1[i].posn);
                    for(size_t j=0; j<eclipses.size(); j++)
                        star1[i].eclipse.push_back(eclipses[j]);
                }
	    
                // Apply eclipse by disc to star 2
                for(int i=0; i<star2.size(); i++){
                    eclipses =  Roche::disc_eclipse(model.iangle, rdisc1, rdisc2, model.beta_disc, model.height_disc, star2[i].posn);
                    for(size_t j=0; j<eclipses.size(); j++)
                        star2[i].eclipse.push_back(eclipses[j]);
                }
            }
        }

        if(model.add_spot){
            //      double wave = 5000, temp_spot = 10000, cfrac_spot = 0.5, height_spot = 0.01;
            //      Lcurve::set_bright_spot_grid(q, iangle, r1, r2, roche1, roche2, eclipse1, eclipse2, delta_phase, radius_spot, 
            //			   length_spot, height_spot, expon_spot, epow_spot, angle_spot, yaw_spot, temp_spot, tilt_spot, 
            //			   cfrac_spot, nspot, wave, bspot);
            Subs::Vec3 dir(1,0,0), posn, v;
	  
            double rl1 = Roche::xl1(model.q);
      
            // Calculate a reference radius and potential for the two stars
            double rref1, pref1, ffac1 = r1/rl1;
            Roche::ref_sphere(model.q, Roche::PRIMARY, model.spin1, ffac1, rref1, pref1);
      
            double rref2, pref2, ffac2 = r2/rl2;
            Roche::ref_sphere(model.q, Roche::SECONDARY, model.spin2, ffac2, rref2, pref2);
      
            dir.set(0,0,1);
            Roche::strinit(model.q, posn, v);
      
            Lcurve::Point::etype eclipses, edisc;
            Lcurve::star_eclipse(model.q, r1, model.spin1, ffac1, model.iangle, posn, model.delta_phase, model.roche1, Roche::PRIMARY,   eclipses);
            Lcurve::star_eclipse(model.q, r2, model.spin2, ffac2, model.iangle, posn, model.delta_phase, model.roche2, Roche::SECONDARY, eclipses);
            stream.push_back(Lcurve::Point(posn, dir, 0., 1., eclipses));
		       
            const int NSTREAM = int((rl1-model.radius_spot)/0.001);
            double radius;
            for(int i=0; i<NSTREAM; i++){
                radius = rl1 + (model.radius_spot-rl1)*(i+1)/NSTREAM;
                Roche::stradv(model.q, posn, v, radius, 1.e-10, 1.e-3);
                eclipses.clear();
                Lcurve::star_eclipse(model.q, r1, model.spin1, ffac1, model.iangle, posn, model.delta_phase, model.roche1, Roche::PRIMARY,   eclipses);
                Lcurve::star_eclipse(model.q, r2, model.spin2, ffac2, model.iangle, posn, model.delta_phase, model.roche2, Roche::SECONDARY, eclipses);
                if(model.add_disc){
                    edisc = Roche::disc_eclipse(model.iangle, rdisc1, rdisc2, model.beta_disc, model.height_disc, posn);
                    for(size_t j=0; j<edisc.size(); j++)
                        eclipses.push_back(edisc[j]);
                }
                stream.push_back(Lcurve::Point(posn, dir, 0., 1., eclipses));
            }	
        }

        // Plot
        Subs::Plot plot(device);

        cpgpap(width, (y2-y1)/(x2-x1));
        
        // User-selected colormap
        set_colormap(colormap, reverse, ncolors);
        
        // Make stars orbit around centre of mass of system
        const Subs::Vec3 cofm(model.q/(1.+model.q),0.,0.);
        Subs::Vec3 r, earth;

        for(int np=0; np<nphase; np++){

            if(nphase > 1) phase = phase1 + (phase2-phase1)*np/double(nphase-1);

            // cpgenv(x1, x2, y1, y2, 1, -2);
            cpgpage(); // Manually start a new Postscript page, otherwise all pages get stacked into one.
            cpgsvp(0.01, 0.99, 0.01, 0.99); // X and Y limits for the plotting area (think Python's plt.subplots_adjust(left, right, bottom, top)
            cpgswin(x1, x2, y1, y2);    // Axis limits similar to ax.set_xlim() and ax.set_ylim()
            // cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0); // Similar to ax.tick_params(), but more detailed and with extra features.
                                                      // BCNST is split into X and Y sets, which is why there are two in this line.
                                                      //     B = show bottom(left) axis
                                                      //     C = show top(right) axis
                                                      //     N = Display numeric tick labels
                                                      //     S = Display minor tick markers
                                                      //     T = Display major tick markers
                                                      // 2nd and 3rd arguments defined major and minor tick spacings
                                                      //     using "0" means "auto" to let the plot determine tick spacings itself
        
            
            earth = Roche::set_earth(model.iangle, phase);

            // Compute sky basis vectors
            double cosp = cos(Constants::TWOPI*phase);
            double sinp = sin(Constants::TWOPI*phase);

            Subs::Vec3 xsky(sinp,cosp,0.);
            Subs::Vec3 ysky = Subs::cross(earth, xsky);

            cpgbbuf();

            cpgsch(2);

            // star 1
            cpgsci(4);
            plot_visible(star1, earth, cofm, xsky, ysky, phase, min_stat, max_stat, colorscale, ncolors, plt_marker);

            // star 2
            cpgsci(2);
            plot_visible(star2, earth, cofm, xsky, ysky, phase, min_stat, max_stat, colorscale, ncolors, plt_marker);

            if(model.add_disc){

                // disc surface
                cpgsci(3);
                plot_visible(disc, earth, cofm, xsky, ysky, phase, min_stat, max_stat, colorscale, ncolors, plt_marker);

                // edges
                cpgsci(1);
                plot_visible(outer_edge, earth, cofm, xsky, ysky, phase, min_stat, max_stat, colorscale, ncolors, plt_marker);
                plot_visible(inner_edge, earth, cofm, xsky, ysky, phase, min_stat, max_stat, colorscale, ncolors, plt_marker);
            }

            if(model.add_spot){
                cpgsci(2);
                plot_visible(stream, earth, cofm, xsky, ysky, phase, min_stat, max_stat, colorscale, ncolors, -1);
                cpgsci(2);
                double cosbs = Subs::dot(earth, stream[stream.size()-1].posn);
                if(cosbs > 0. && stream[stream.size()-1].visible(phase)){
                    r = stream[stream.size()-1].posn - cofm;
                    cpgslw(6);
                    cpgsch(0.5+3.5*cosbs);
                    cpgpt1(Subs::dot(r, xsky), Subs::dot(r, ysky), 18);
                }
            }
            cpgebuf();
        }
    }
    catch(const std::string& err){
        std::cerr << err << std::endl;
        exit(EXIT_FAILURE);
    }
}


// plots visible points
void plot_visible(const Subs::Buffer1D<Lcurve::Point>& object, const Subs::Vec3& earth, const Subs::Vec3& cofm, const Subs::Vec3& xsky, const Subs::Vec3& ysky, double phase, double fmin, double fmax, std::string colorscale, int ncolors, int plt_marker){
    Subs::Vec3 r;

    bool log_scaling = false;
    if(Subs::tolower(colorscale)=="log"){
        log_scaling = true;
    }

    bool object_scaling=false; //  TODO, user input for "object" or "system" scaling.
    if(object_scaling){
        fmin = 1e30;
        fmax = -1e30;
        for(int i=0; i<object.size(); i++){
            fmin = std::min<double>(fmin, object[i].temp);
            fmax = std::max<double>(fmax, object[i].temp);
        }
    }

    
    double log_fmin = log10(fmin + 1e-20);
    double log_fmax = log10(fmax + 1e-20);

    // Allowed up to 255 color indices with /cps (color postscript)
    // The first 16 are saved for PGPLOT defaults, such as background color
    for(int i=0; i<object.size(); i++){
        if(Subs::dot(earth, object[i].dirn) > 0. && object[i].visible(phase)){

            r = object[i].posn - cofm;

            double point_flux = object[i].temp;
            double log_point_flux = log10(point_flux + 1e-20);

            double norm; // Normalize the point's flux.
            if(log_scaling){
                norm = (log_point_flux - log_fmin) / (log_fmax - log_fmin);
            }else{
                norm = (point_flux - fmin) / (fmax - fmin);
            }
            norm = std::min(1.0, std::max(0.0, norm)); // Ensure normalization between 0.0 and 1.0    
            
            int color = 16 + int(norm * (ncolors - 1) + 0.5); // Set the color index based on the normalization.
            cpgsci(color);
            
            // Plot one data point as a small dot.
            cpgpt1(Subs::dot(r, xsky), Subs::dot(r, ysky), plt_marker); 
        }
    }
}