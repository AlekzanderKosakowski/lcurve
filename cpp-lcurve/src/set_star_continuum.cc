#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/lcurve.h"

/** set_star_continuum computes the continuum face-on brightness*area of each
 * element of the two stars assuming a black-body relation . The actual
 * contribution to the light-curve is the brightness times the area times (1 -
 * lin_limb*(1-mu) - quad_limb*(1-mu)**2), following Wade & Rucinski (1985)
 * where lin_limb and quad_limb are the linear and quadratic limb darkening
 * coefficients and mu is the cosine of the face angle.  This leads to a
 * correction factor of 1/(1-a/3-b/6) to the brightness to get the right
 * effective temperature.
 *
 * Irradiation is computed by adding the contribution of the other star to the
 * one being irradiated; the process is not iterated to compute
 * 'back-heating'. The finite size of the source is included only
 * approximately.
 *
 * \param mdl          Model defining parameters
 * \param star1        grid of elements over star 1, modified on exit
 * \param star2        grid of elements over star 2, modified on exit.
 */

void Lcurve::set_star_continuum(const Model& mdl,
                                Subs::Buffer1D<Lcurve::Point>& star1,
                                Subs::Buffer1D<Lcurve::Point>& star2){

    double r1, r2;
    mdl.get_r1r2(r1, r2);

    double rl1 = Roche::xl1(mdl.q);
    if(r1 < 0) r1 = rl1;

    double rl2 = 1.-rl1;
    if(r2 < 0) r2 = rl2;

    double mu, r, temp;
    Subs::Vec3 vec;
    const Subs::Vec3 cofm2(1.,0.,0.);
    
    // First compute irradiation of star 1 by star 2. 'geom' is the
    // geometrical factor by which the radiation from star 2 is reduced
    // by the time it hits the element in terms of flux per unit area
    // compared to its level as it leaves star 2
    double geom;

    // modify the gravity darkening coefficient to allow for two possibilities
    // the gravity darkening is implemented by modifying the temperature and
    // then calculating the flux in a BB approx. The 'bolometric' method does
    // this directly; the 'filter integrated' method modifies the exponent
    // used to give the desired behaviour of flux with gravity.
    const double GDCBOL1 = mdl.gdark_bolom1 ? mdl.gravity_dark1 :
        mdl.gravity_dark1 / Subs::dlpdlt(mdl.wavelength, mdl.t1);

    int nelem1 = star1.size();

    // compute direction of star spot 11, 12, 13, and the direct impact starspot (stsp1i)
    Subs::Vec3 spot11, spot12, spot13, spot1i;
    bool is_spot11 = mdl.stsp11_long.defined && mdl.stsp11_lat.defined &&
        mdl.stsp11_fwhm.defined && mdl.stsp11_tcen.defined;
    bool is_spot12 = mdl.stsp12_long.defined && mdl.stsp12_lat.defined &&
      mdl.stsp12_fwhm.defined && mdl.stsp12_tcen.defined;
    bool is_spot13 = mdl.stsp13_long.defined && mdl.stsp13_lat.defined &&
        mdl.stsp13_fwhm.defined && mdl.stsp13_tcen.defined;
    
    bool is_spot1i = mdl.stsp1i_long.defined && mdl.stsp1i_lat.defined &&
                    mdl.stsp1i_fwhm_lat.defined && mdl.stsp1i_fwhm_long1.defined && mdl.stsp1i_fwhm_long2.defined &&
                    mdl.stsp1i_tcen.defined;
    
    double clong11=0., slong11=0., clat11=0., slat11=0.;
    double clong12=0., slong12=0., clat12=0., slat12=0.;
    double clong13=0., slong13=0., clat13=0., slat13=0.;
    double clong1i=0., slong1i=0., clat1i=0., slat1i=0.;
    if(is_spot11){
        clong11 = std::cos(Subs::deg2rad(mdl.stsp11_long.value));
        slong11 = std::sin(Subs::deg2rad(mdl.stsp11_long.value));
        clat11 = std::cos(Subs::deg2rad(mdl.stsp11_lat.value));
        slat11 = std::sin(Subs::deg2rad(mdl.stsp11_lat.value));
    	spot11 = Subs::Vec3(clat11*clong11, clat11*slong11, slat11);
    }

    if(is_spot12){
        clong12 = std::cos(Subs::deg2rad(mdl.stsp12_long.value));
        slong12 = std::sin(Subs::deg2rad(mdl.stsp12_long.value));
        clat12  = std::cos(Subs::deg2rad(mdl.stsp12_lat.value));
        slat12  = std::sin(Subs::deg2rad(mdl.stsp12_lat.value));
    	spot12 = Subs::Vec3(clat12*clong12, clat12*slong12, slat12);
    }

    if(is_spot13){
        clong13 = std::cos(Subs::deg2rad(mdl.stsp13_long.value));
        slong13 = std::sin(Subs::deg2rad(mdl.stsp13_long.value));
        clat13  = std::cos(Subs::deg2rad(mdl.stsp13_lat.value));
        slat13  = std::sin(Subs::deg2rad(mdl.stsp13_lat.value));
    	spot13 = Subs::Vec3(clat13*clong13, clat13*slong13, slat13);
    }

    if(is_spot1i){ // Star1 impact spot
        clong1i = std::cos(Subs::deg2rad(mdl.stsp1i_long.value));
        slong1i = std::sin(Subs::deg2rad(mdl.stsp1i_long.value));
        clat1i = std::cos(Subs::deg2rad(mdl.stsp11_lat.value));
        slat1i = std::sin(Subs::deg2rad(mdl.stsp11_lat.value));
    	spot1i = Subs::Vec3(clat1i*clong1i, clat1i*slong1i, slat1i); // (x,y,z) cartesian position of impact starspot
    }

    // Build star1 grid details
    for(int i=0; i<nelem1; i++){
        vec = cofm2 - star1[i].posn;
        r = vec.length();
        mu = Subs::dot(star1[i].dirn, vec)/r;

        // All grid point start with the same temperature t1
        double t1 = mdl.t1;

        // Adjust temperature from Gaussian spots based on angular distance from center.
        if(is_spot11){
        	double dist = Subs::rad2deg(std::acos(Subs::dot(star1[i].posn, spot11)/star1[i].posn.length()));
        	t1 += (mdl.stsp11_tcen-mdl.t1)*std::exp(-Subs::sqr(dist/(mdl.stsp11_fwhm/Constants::EFAC))/2.);
        }
    	if(is_spot12){
        	double dist = Subs::rad2deg(std::acos(Subs::dot(star1[i].posn, spot12)/star1[i].posn.length()));
            t1 += (mdl.stsp12_tcen-mdl.t1)*std::exp(-Subs::sqr(dist/(mdl.stsp12_fwhm/Constants::EFAC))/2.);
        }
    	if(is_spot13){
        	double dist = Subs::rad2deg(std::acos(Subs::dot(star1[i].posn, spot13)/star1[i].posn.length()));
    	    t1 += (mdl.stsp13_tcen-mdl.t1)*std::exp(-Subs::sqr(dist/(mdl.stsp13_fwhm/Constants::EFAC))/2.);
        }
        if(is_spot1i && mdl.stsp1i_fwhm_long1 > 0 && mdl.stsp1i_fwhm_long2 > 0 && mdl.stsp1i_fwhm_lat > 0){
            // Direct-impact starspot with approximate advective tail in the spin direction.
            
            // Latitude and longitudes all in radians here for trig calculations.
            double surface_longitude = std::atan2(star1[i].posn.y(), star1[i].posn.x());
            double surface_latitude = std::asin(star1[i].posn.z() / star1[i].posn.length());
            
            double impact_longitude = Subs::deg2rad(mdl.stsp1i_long);
            double impact_latitude = Subs::deg2rad(mdl.stsp1i_lat);

            double fwhm_upstream   = Subs::deg2rad(mdl.stsp1i_fwhm_long1);
            double fwhm_downstream = Subs::deg2rad(mdl.stsp1i_fwhm_long2);
            double fwhm_latitude   = Subs::deg2rad(mdl.stsp1i_fwhm_lat);
            
            
            // Latitude calculations
            double latitude_offset = Subs::abs(impact_latitude - surface_latitude);
            double latitude_decay_term = std::exp(-Subs::sqr(latitude_offset/(fwhm_latitude/Constants::EFAC))/2.);


            // Longitude calculations
            double delta_phi = surface_longitude - impact_longitude; // Raw difference in longitude (radians). Positive means surface element is downstream from impact center
            delta_phi = std::fmod(delta_phi + Constants::PI, 2.0*Constants::PI) - Constants::PI;
            double delta_phi_phys = delta_phi * (1+0*std::cos(surface_latitude)); // Apply spherical correction to fix decay at high latitudes (energy transport stuff)

            double upstream_limit = -5.0 * fwhm_upstream / Constants::EFAC; // Hard cutoff at 5 sigma for "upstream" direction limits.
            bool upstream = (delta_phi_phys >= upstream_limit) && (delta_phi_phys <= 0.0);

            double longitude_decay_term;
            if(upstream){
                longitude_decay_term = std::exp(-Subs::sqr(delta_phi_phys/(fwhm_upstream/Constants::EFAC))/2.);
            }else{
                double delta_phi_downstream = std::fmod(surface_longitude - impact_longitude, 2.0*Constants::PI);
                if(delta_phi_downstream < 0.0){
                    delta_phi_downstream += 2.0*Constants::PI;
                }
                double delta_phi_downstream_phys = delta_phi_downstream * (1+0*std::cos(surface_latitude));
                
                double fraction_core = std::exp(-Subs::sqr(delta_phi_phys/(fwhm_upstream/Constants::EFAC))/2.);
                double fraction_tail = std::exp(-delta_phi_downstream_phys/fwhm_downstream) * (1.0 - fraction_core);
                longitude_decay_term = fraction_core + fraction_tail;
            }
            t1 += (mdl.stsp1i_tcen - mdl.t1) * longitude_decay_term * latitude_decay_term;                
        }



        
        // Irradiation handled by treating star2 as a point source with no starspots.
        if(mu >= r2){
            // Full tilt irradiation
            geom = Subs::sqr(r2/r)*mu;
            temp = pow(pow(t1*pow(double(star1[i].gravity),GDCBOL1),4) + mdl.absorb*pow(mdl.t2,4)*geom, 0.25);
        }else if(mu > -r2){
            // 'sunset' case
            double x0 = -mu/r2;
            // The following factor is a weighted version of 'mu' as the
            // secondary sets as far as this element is concerned.  When x0 =
            // -1 it equals r2 = mu. when x0 = 0 it equals 2*r2/(3*Pi) as
            // opposed to zero which it would be in the point source case.
            geom = Subs::sqr(r2/r)*r2*(sqrt(1.-x0*x0)*(2+x0*x0)/3 - x0*(Constants::PI/2-asin(x0)))/Constants::PI;
            temp = pow(pow(t1*pow(double(star1[i].gravity),GDCBOL1),4) + mdl.absorb*pow(mdl.t2,4)*geom, 0.25);
        }else{
            // No irradiation
            geom = 0.;
            temp = t1*pow(double(star1[i].gravity),GDCBOL1);
        }

        star1[i].temp = temp;
            
        // At this stage also add in a directly reflected part too
        star1[i].flux  = star1[i].area*Subs::planck(mdl.wavelength, temp);

        if(mdl.mirror){
            star1[i].flux  += star1[i].area*geom*Subs::planck(mdl.wavelength, Subs::abs(double(mdl.t2)));
        }
    }

    const Subs::Vec3 cofm1(0.,0.,0.);

    // See comments on GDCBOL1
    const double GDCBOL2 = mdl.gdark_bolom2 ? mdl.gravity_dark2 :
        mdl.gravity_dark2 / Subs::dlpdlt(mdl.wavelength, Subs::abs(double(mdl.t2)));

    int nelem2 = star2.size();

    // compute direction of star spot 21
    Subs::Vec3 spot21, spot22;
    
    bool is_spot21 = mdl.stsp21_long.defined && mdl.stsp21_lat.defined &&
        mdl.stsp21_fwhm.defined && mdl.stsp21_tcen.defined;
    bool is_spot22 = mdl.stsp22_long.defined && mdl.stsp22_lat.defined &&
      mdl.stsp22_fwhm.defined && mdl.stsp22_tcen.defined;
    
    double clong21=0., slong21=0., clat21=0., slat21=0.;
    double clong22=0., slong22=0., clat22=0., slat22=0.;
    if(is_spot21){
        clong21 = std::cos(Subs::deg2rad(mdl.stsp21_long.value));
        slong21 = std::sin(Subs::deg2rad(mdl.stsp21_long.value));
        clat21  = std::cos(Subs::deg2rad(mdl.stsp21_lat.value));
        slat21  = std::sin(Subs::deg2rad(mdl.stsp21_lat.value));
    	spot21 = Subs::Vec3(clat21*clong21, clat21*slong21, slat21);
    }

    if(is_spot22){
        clong22 = std::cos(Subs::deg2rad(mdl.stsp22_long.value));
        slong22 = std::sin(Subs::deg2rad(mdl.stsp22_long.value));
        clat22  = std::cos(Subs::deg2rad(mdl.stsp22_lat.value));
        slat22  = std::sin(Subs::deg2rad(mdl.stsp22_lat.value));
    	spot22 = Subs::Vec3(clat22*clong22, clat22*slong22, slat22);
    }


    
    // Star1 surface elements have been completely built with starspots included by now.
    // Now use star1's surface elements in a nested loop with star2's surface elements to find star2's flux, including irradiation from starspots.
    for(int i=0; i<nelem2; i++){

        // Begin with a single temperature for all surface elements
        double t2 = Subs::abs(double(mdl.t2));

        // Increase temperature based on starspot parameters.
        if(is_spot21){
            Subs::Vec3 off = star2[i].posn-cofm2;
            double dist = Subs::rad2deg(std::acos(Subs::dot(off, spot21)/off.length()));
            t2 += (mdl.stsp21_tcen-t2)*std::exp(-Subs::sqr(dist/(mdl.stsp21_fwhm/Constants::EFAC))/2.);
        }
    	if(is_spot22){
            Subs::Vec3 off = star2[i].posn-cofm2;
            double dist = Subs::rad2deg(std::acos(Subs::dot(off, spot22)/off.length()));
            t2 += (mdl.stsp22_tcen-t2)*std::exp(-Subs::sqr(dist/(mdl.stsp22_fwhm/Constants::EFAC))/2.);
        }
    
        if(mdl.finite_irr12 && mdl.absorb > 0){ // Use finite surface elements from star1 to irradiate star2.
                                                // Allows the donor to be irradiated by starspots on the surface of the accretor (think direct-impact accretion).
        
            // Adjust temperature based on gravity darkening
            double T_intrinsic = t2 * pow(double(star2[i].gravity), GDCBOL2);
    
            // Begin loop over star1's elements for precise irradiations
            double F_irradiation = 0.0;  // Accumulate irradiation flux across all star1 surface elements
            for(int j=0; j<nelem1; j++){
                // Vector from star2 element to this star1 element
                Subs::Vec3 vec12 = star1[j].posn - star2[i].posn;
                double r12 = vec12.length(); // For geometric dilution
        
                // Cosine of angle between star1 element normal and direction to star2
                double mu1 = Subs::dot(star1[j].dirn, -vec12) / r12;
        
                // Cosine of angle between star2 element normal and direction to star1 element
                double mu2 = Subs::dot(star2[i].dirn, vec12) / r12;
        
                // Only include contribution if both elements face each other (ignores edge cases like rotation/beaming)
                if(mu1 > 0.0 && mu2 > 0.0){
                    
                    // Geometric factor: projected area / distance^2
                    double geom_factor = mu1 * mu2 * star1[j].area / (r12*r12);
        
                    // Add to irradiation flux (absorbed fraction)
                    F_irradiation += mdl.absorb * pow(star1[j].temp, 4) * geom_factor;
                }
            }

            // Total effective temperature including irradiation
            double T_eff = pow(pow(T_intrinsic, 4) + F_irradiation, 0.25);
            // std::cout << i << "/" << nelem2 << " " << T_intrinsic << " " << F_irradiation << " " << T_eff << std::endl;
            star2[i].temp = T_eff;
                
            // Compute monochromatic flux for this element
            star2[i].flux = star2[i].area * Subs::planck(mdl.wavelength, T_eff);

        }else{ // Original behavior: treat star1 as a point source without starspots

            vec = cofm1 - star2[i].posn;
            r   = vec.length();
            mu  = Subs::dot(star2[i].dirn,vec)/r;

            // compute unirradiated temperature allowing for
            // offset from spot centre
            double t2 = Subs::abs(double(mdl.t2));
            if(is_spot21){
                Subs::Vec3 off = star2[i].posn-cofm2;
                double dist =
                    Subs::rad2deg(std::acos(Subs::dot(off, spot21)/
                                            off.length()));
                t2 += (mdl.stsp21_tcen-t2)*
                    std::exp(-Subs::sqr(dist/(mdl.stsp21_fwhm/Constants::EFAC))/2.);
            }
    
            if(mu >= r1){
                geom = Subs::sqr(r1/r)*mu;
                temp = pow(pow(t2*pow(double(star2[i].gravity),GDCBOL2),4) + mdl.absorb*pow(mdl.t1,4.)*geom, 0.25);
    
            }else if(mu > -r1){
                double x0 = -mu/r1;
                geom = Subs::sqr(r1/r)*r1*(sqrt(1.-x0*x0)*(2+x0*x0)/3 - x0*(Constants::PI/2-asin(x0)))/Constants::PI;
                temp = pow(pow(t2*pow(double(star2[i].gravity),GDCBOL2),4) + mdl.absorb*pow(mdl.t1,4.)*geom, 0.25);
    
            }else{
                geom = 0.;
                temp = t2*pow(double(star2[i].gravity), GDCBOL2);
            }
    
            star2[i].flux = star2[i].area*Subs::planck(mdl.wavelength, temp);
    
            if(mdl.mirror){
                star2[i].flux += star2[i].area*geom*Subs::planck(mdl.wavelength, mdl.t1);
            }
        }        
    }

}
// Do not iterate to calculate star1's flux based on star2's surface element flux.
// TODO: Consider systems with starspot and disc. The disc may block irradiation.
//     Check that the vector between star1[j].posn and star2[i].posn doesn't intersect disc[k].posn

