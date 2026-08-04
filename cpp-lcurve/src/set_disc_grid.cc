#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

/**
 * set_disc_grid sets up the elements needed to define an accretion
 * disc around the primary star in a co-rotating binary system.  For
 * each element set_disc_grid computes the position within the binary,
 * the direction perpendicular to the element, the area of the element
 * and whether it is eclipsed by the stars. The disc is modelled as
 * having a power law variation of height with radius. Elements are
 * only computed for the top surface of the disc. The area calculation
 * ignores the angle of the elements.
 *
 * All numbers assume a separation = 1. The primary star is centred at
 * (0,0,0), the secondary at (1,0,0). The y-axis points in the
 * direction of motion of the secondary star.
 *
 * \param mdl    the Model
 * \param disc      array representing the disc
 * \exception Exceptions are thrown if the specified radii over-fill the
 * Roche lobes.
 */

void Lcurve::set_disc_grid(const Model& mdl,
			   Subs::Buffer1D<Lcurve::Point>& disc){

    const double EFAC = 1.0000001;

    double r1, r2;
    mdl.get_r1r2(r1, r2);

    double rl1 = Roche::xl11(mdl.q, mdl.spin1);
    if(r1 < 0)
        r1 = rl1;
    else if(r1 > rl1)
        throw Lcurve_Error("set_disc_grid: the primary star is larger than its Roche lobe!");

    double rl2 = 1.- Roche::xl12(mdl.q, mdl.spin2);
    if(r2 < 0)
        r2 = rl2;
    else if(r2 > rl2)
        throw Lcurve_Error("set_disc_grid: the secondary star is larger than its Roche lobe!");

    // note that the inner radius of the disc is set equal to that of
    // the white dwarf if rdisc1 <= 0 while the outer disc is set
    // equal to the spot radius
    double rdisc1 = mdl.rdisc1 > 0. ? mdl.rdisc1 : r1;
    double rdisc2 = mdl.rdisc2 > 0. ? mdl.rdisc2 : mdl.radius_spot;

    // Calculate a reference radius and potential for the two stars
    double rref1, pref1, ffac1 = r1/rl1;
    Roche::ref_sphere(mdl.q, Roche::PRIMARY, mdl.spin1, ffac1, rref1, pref1);

    double rref2, pref2, ffac2 = r2/rl2;
    Roche::ref_sphere(mdl.q, Roche::SECONDARY, mdl.spin2, ffac2, rref2, pref2);

    // Centre of masses of the stars
    const Subs::Vec3 cofm1(0.,0.,0.), cofm2(1.,0.,0.);

    int nface = 0, ntheta = 0;

    // drrad here introduced to avoid ntheta getting silly if rdisc1 and rdisc2 
    // are close to each other.
    double drad = (rdisc2-rdisc1)/mdl.nrad, drrad = rdisc2/mdl.nrad;
    for(int i=0; i<mdl.nrad; i++){
        ntheta = int(ceil(Constants::TWOPI*(rdisc1 + (rdisc2-rdisc1)*(i+0.5)/mdl.nrad)/drrad));
        ntheta = ntheta > 8 ? ntheta : 8;
        nface += ntheta;
    }

    disc.resize(nface);

    Subs::Vec3 dirn, posn;
    double theta, rad, h, area, sint, cost, tanp, sinp, cosp;

    int i;
    Lcurve::Point::etype eclipses;
    for(i=0, nface=0; i<mdl.nrad; i++){
        rad     = rdisc1 + (rdisc2-rdisc1)*(i+0.5)/mdl.nrad;
        ntheta  = int(ceil(Constants::TWOPI*rad/drrad));
        ntheta  = ntheta > 8 ? ntheta : 8;
        h       = EFAC*mdl.height_disc*pow(rad, mdl.beta_disc);
        tanp    = mdl.beta_disc*h/rad; // tan of tilt angle of element
        cosp    = 1./sqrt(1+Subs::sqr(tanp));
        sinp    = tanp*cosp;

        // NB: no accounting for the angle of the face
        area    = Constants::TWOPI/ntheta*rad*drad;

        for(int j=0; j<ntheta; j++, nface++){

            theta = Constants::TWOPI*j/ntheta;
            sint  = sin(theta);
            cost  = cos(theta);
            posn.set(rad*cost, rad*sint, h);
            dirn.set(-cost*sinp, -sint*sinp, cosp);

            if(mdl.opaque)
                eclipses =  Roche::disc_eclipse(mdl.iangle, rdisc1, rdisc2,
						mdl.beta_disc, mdl.height_disc,
						posn);
            else
                eclipses.clear();
            if(mdl.eclipse1) Lcurve::star_eclipse(mdl.q, r1, mdl.spin1, ffac1,
						  mdl.iangle, posn,
						  mdl.delta_phase,
                                                  mdl.roche1, Roche::PRIMARY,
						  eclipses);
            if(mdl.eclipse2) Lcurve::star_eclipse(mdl.q, r2, mdl.spin2, ffac2,
						  mdl.iangle, posn,
						  mdl.delta_phase,
                                                  mdl.roche2,
						  Roche::SECONDARY, eclipses);

            disc[nface] = Lcurve::Point(posn, dirn, area, 1., eclipses);

        }
    }
}

/**
 * set_disc_edge sets up the elements needed to define a rim of a
 * cyclindrical accretion disc. The elements are assumed to face
 * radially outwards or inwards depending upon whether it is an outer
 * or inner edge, except for elements on the upper edges of the rim
 * which are set to face upwards to guarantee that they will be
 * displayed unless they are eclipsed.
 *
 * All numbers assume a separation = 1. The primary star is centred at
 * (0,0,0), the secondary at (1,0,0). The y-axis points in the
 * direction of motion of the secondary star.
 *
 *
 * \param mdl       input model
 * \param outer     true for the outer edge, false for the inner edge
 * \param edge      array representing the outer edge
 * \param visial    array representing the outer edge
 * \exception Exceptions are thrown if the specified radii over-fill the
 * Roche lobes.
 */

void Lcurve::set_disc_edge(const Model& mdl, bool outer,
			   Subs::Buffer1D<Lcurve::Point>& edge,
			   bool visual){

    const double EFAC = 1.0000001;

    double r1, r2;
    mdl.get_r1r2(r1, r2);

    double rl1 = Roche::xl11(mdl.q, mdl.spin1);
    if(r1 < 0)
        r1 = rl1;
    else if(r1 > rl1)
        throw Lcurve_Error("set_disc_grid: the primary star is larger than its Roche lobe!");

    double rl2 = 1.- Roche::xl12(mdl.q, mdl.spin2);
    if(r2 < 0)
        r2 = rl2;
    else if(r2 > rl2)
        throw Lcurve_Error("set_disc_grid: the secondary star is larger than its Roche lobe!");

    // note that the inner radius of the disc is set equal to that of the white dwarf if rdisc1 <= 0
    // while the outer disc is set equal to the spot radius
    double rdisc1 = mdl.rdisc1 > 0. ? mdl.rdisc1 : r1;
    double rdisc2 = mdl.rdisc2 > 0. ? mdl.rdisc2 : mdl.radius_spot;
    double size = rdisc2/mdl.nrad;

    // Calculate a reference radius and potential for the two stars
    double rref1, pref1, ffac1 = r1/rl1;
    Roche::ref_sphere(mdl.q, Roche::PRIMARY, mdl.spin1, ffac1, rref1, pref1);

    double rref2, pref2, ffac2 = r2/rl2;
    Roche::ref_sphere(mdl.q, Roche::SECONDARY, mdl.spin2, ffac2, rref2, pref2);

    // Centre of masses of the stars
    const Subs::Vec3 cofm1(0.,0.,0.), cofm2(1.,0.,0.);

    const double RAD = outer ? rdisc2 : rdisc1;
    const double H = mdl.height_disc*pow(RAD, mdl.beta_disc);
    const int NOUT = int(2*H/size) + 2;
    const int NTHETA = int(ceil(Constants::TWOPI*RAD/size)) > 8 ? int(ceil(Constants::TWOPI*RAD/size)) : 8;

    edge.resize(NTHETA*NOUT);

    Subs::Vec3 dirn, posn;
    double theta, area, sint, cost;

    // Assign half-areas for top and bottom elements
    area = Constants::TWOPI/NTHETA*RAD*2*H/(NOUT-1);

    // Note nface carries on here as we add more elements to the end
    Lcurve::Point::etype eclipses;
    for(int i=0, nface=0; i<NTHETA; i++){

        theta = Constants::TWOPI*i/NTHETA;
        sint = sin(theta);
        cost = cos(theta);

        // Upper rim element
        if(outer)
            posn.set(EFAC*RAD*cost, EFAC*RAD*sint, EFAC*H);
        else
            posn.set(RAD*cost/EFAC, RAD*sint/EFAC, EFAC*H);

        // Direction set so that it will always be visible if is not
        // eclipsed if for visualisation
	if(visual)
	  dirn.set(0, 0, 1);
	else
	  dirn.set(cost, sint, 0);


        if(mdl.opaque)
            eclipses = Roche::disc_eclipse(mdl.iangle, rdisc1, rdisc2, mdl.beta_disc, mdl.height_disc, posn);
        else
            eclipses.clear();

        // Primary star can expand to wipe out inner disc; trap but report such errors
        if(mdl.eclipse1){
            try{
                Lcurve::star_eclipse(mdl.q, r1, mdl.spin1, ffac1, mdl.iangle, posn, mdl.delta_phase,
                                     mdl.roche1, Roche::PRIMARY, eclipses);
            }
            catch(const Roche::Roche_Error& err){
                std::cerr << "Trapped following error:\n" << err.what() << std::endl;
                eclipses.push_back(std::make_pair(0.,1.));
            }
        }

        if(mdl.eclipse2) Lcurve::star_eclipse(mdl.q, r2, mdl.spin2, ffac2, mdl.iangle, posn, mdl.delta_phase,
                                              mdl.roche2, Roche::SECONDARY, eclipses);
        edge[nface++] = Lcurve::Point(posn, dirn, area/2., 1., eclipses);

        // All lower elements
        for(int j=0; j<NOUT-1; j++){

            if(outer){
                posn.set(EFAC*RAD*cost, EFAC*RAD*sint, -H+2*H*j/(NOUT-1));
                dirn.set(cost, sint, 0);
            }else{
                posn.set(RAD*cost/EFAC, RAD*sint/EFAC, -H+2*H*j/(NOUT-1));
                dirn.set(-cost, -sint, 0);
            }

            if(mdl.opaque)
                eclipses =  Roche::disc_eclipse(mdl.iangle, rdisc1, rdisc2, mdl.beta_disc, mdl.height_disc, posn);
            else
                eclipses.clear();
            if(mdl.eclipse1) Lcurve::star_eclipse(mdl.q, r1, mdl.spin1, ffac1, mdl.iangle, posn, mdl.delta_phase,
                                                  mdl.roche1, Roche::PRIMARY, eclipses);
            if(mdl.eclipse2) Lcurve::star_eclipse(mdl.q, r2, mdl.spin2, ffac2, mdl.iangle, posn, mdl.delta_phase,
                                                  mdl.roche2, Roche::SECONDARY, eclipses);

            edge[nface++] = Lcurve::Point(posn, dirn, area, 1., eclipses);

        }
    }
}




/**
 * Finite-radius irradiation between stars introduces the problem of
 * an accretion disc intercepting rays between stars, preventing
 * irradiation of certain parts of the receiving star. This function 
 * tries to handle that by determining if a ray is intercepted on its path
 * between two surface elements.
 *
 * This function is called once per surface element pair during finite irradiation
 * calculations. This means it is only called if finite_irr12 is enabled and
 * only if mu1>0 and mu2>0.
 *
 * This approximates the disc interception by checking the ray height (z-position)
 * where it crosses the outer disc radius. This assumes the outer rim of the disc
 * is the disc's maximum vertical extent (a simplified concave flared disc model).
 *    
 * This is done by solving the quadratic for the eqn of a circle with r = rdisc2
 * for the intersection points of the irradiation with the disc edges.
 * If the irradiation passed through the disc radial extent, then compare its z-posn
 * to the height of the disc at that (x, y) position to confirm occultation.
 *
 * Outputs :
 *         True  : Boolean
 *                 The elements can not see each other
 *                 The disc intercepts the ray
 *         False : Boolean
 *                 The elements can see each other
 *                 The disc does not intercept the ray
 *
 *
 * Inputs :
 *         ray_origin : Lcurve::Point object
 *                      Starting location of the irradiation (surface element on star)
 *         ray_vec    : Subs::Vec3 object
 *                      Full displacement vector connecting the ray's start and end locations
                          starting on one surface element and ending on another
 *         mdl        : Lcurve::Model object
 *                      Contains necessary disc geometry parameters
 *
 */
bool Lcurve::disc_blocks_ray(const Lcurve::Point& ray_origin,
                             const Subs::Vec3& ray_vec,
                             const Model& mdl
                            ){
    if(!mdl.add_disc){
        return false;
    }
    
    double a = ray_vec.x()*ray_vec.x() + ray_vec.y()*ray_vec.y();
    if(a < 1e-12){
        return false;
    }

    double b = 2.0*(ray_origin.posn.x()*ray_vec.x() +
                    ray_origin.posn.y()*ray_vec.y());

    double c = ray_origin.posn.x()*ray_origin.posn.x() +
               ray_origin.posn.y()*ray_origin.posn.y() -
               mdl.rdisc2*mdl.rdisc2;


    double discriminant = b*b - 4.0*a*c;
    if (discriminant < 0.0){
        if (discriminant > -1e-12){
            discriminant = 0.0;
        }
        else{
            throw Lcurve_Error("Lcurve::disc_blocks_ray(): Unexpected ray geometry: Negative discriminant");
        }
    }
    
    double root = sqrt(discriminant);

    // Two roots (the points the line intersects the circle)
    // "t" is the fractional distance between the surface elements in direction of propagation
    double t1 = (-b-root)/(2.0*a);
    double t2 = (-b+root)/(2.0*a);

    // Since we solve for these roots as fractional distances between the start/end of the ray's path,
    //   the only physically correct root is the one with value between 0.0 and 1.0
    // If the root is t < 0.0, then the disc-crossing point is behind the starting point
    //   This will be one of the roots if the irradiation comes from the accretor
    // If the root is t > 0.0, then the disc-crossing point is behind the ending point
    //   This will be one of the roots if the irradiation comes from the donor

    // Check the first root.
    double t_hit = t1;
    if(t_hit < 0.0 || t_hit > 1.0)
        t_hit = t2;

    if(t_hit < 0.0 || t_hit > 1.0){ // One of the roots must be between 0.0 and 1.0 since the disc must surround one of the stars.
        throw Lcurve_Error("Lcurve::disc_blocks_ray(): Unexpected ray geometry: No roots found between 0.0 and 1.0");
    }

    // We now know the (x,y,z) position of the ray when it reaches a distance
    // rdisc2 from the accretor, which is the tallest part of the disc in this model.
    // Now we check if the ray's z-position is lower than the disc's z-position at this (x,y) position
    
    Subs::Vec3 p = ray_origin.posn + t_hit*ray_vec;

    double R = sqrt(p.x()*p.x() + p.y()*p.y());
    double H = mdl.height_disc * pow(R, mdl.beta_disc);

    if(fabs(p.z()) < H)
        return true;

    return false;
}
