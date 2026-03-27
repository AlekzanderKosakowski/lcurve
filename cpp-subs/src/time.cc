#include <iomanip>
#include <math.h>
#include <sstream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "slalib.h"
#include "trm/time.h"
#include "trm/vec3.h"
#include "trm/telescope.h"

Subs::Time::Time(int day_, Date::Month month_, int year_, double hour) : Date(day_,month_,year_) {

    if(hour < 0. || hour >= 24.)
        throw Time_Error("Subs::Time::Time(int, Date::Month, int, double): hour out of range");

    hour_  = hour;
}

Subs::Time::Time(int day_, int month_, int year_, double hour) : Date(day_,month_,year_) {

    if(hour < 0. || hour >= 24.)
        throw Time_Error("Subs::Time::Time(int, int, int, double): hour out of range");

    hour_  = hour;
}

Subs::Time::Time(int day_, Date::Month month_, int year_, int hour, int min, double sec) : Date(day_,month_,year_) {

    try{
        valid_time(hour,min,sec);
    }
    catch(const Time_Error& msg){
        throw Time_Error(std::string("Subs::Time::Time(int, Date::Month, int, int, int, double):\n") + msg);
    }

    hour_  = hour + min/60. + sec/3600.;
}

Subs::Time::Time(int day_, int month_, int year_, int hour, int min, double sec) : Date(day_,month_,year_) {

    try{
        valid_time(hour,min,sec);
    }
    catch(const Time_Error& msg){
        throw Time_Error(std::string("Subs::Time::Time(int, int, int, int, int, double):\n") + msg);
    }

    hour_  = hour + min/60. + sec/3600.;
}

Subs::Time::Time(const std::string& time){
    try{
        set(time);
    }
    catch(const Time_Error& msg){
        throw Time_Error(std::string("Subs::Time::Time(const std::string&):\n") + msg);
    }
}

Subs::Time::Time(const Date& date, double hour) : Date(date) {

    if(hour < 0. || hour >= 24.)
        throw Time_Error("Subs::Time::Time(Date, double, float): hour out of range");
    hour_  = hour;
}

Subs::Time::Time(double mjd){
    set(mjd);
}

Subs::Time::HMS Subs::Time::hms() const{
    HMS t;
    t.hour  = (short int)(floor(hour_));
    t.min   = (short int)(floor(60.*(hour_-t.hour)));
    t.sec   = (short int)(floor(60.*(60.*(hour_-t.hour)-t.min)));
    t.fsec  = 60.*(60.*(hour_-t.hour)-t.min)-t.sec;

    // avoid irritating round off problem without compromising the time too much
    if(t.fsec > 0.9999999){
        t.sec++;
        t.fsec = 0.f;
        if(t.sec == 60){
            t.min++;
            t.sec = 0;
            if(t.min == 60){
                t.hour++;
                t.min = 0;
                if(t.hour == 24) t.hour = 0;
            }
        }
    }
    return t;
}

/** Returns time as a modified JD, i.e. JD - 2400000.5. That is they
 * start at midnight UT. They are also more accurate than JDs when round-off
 * is considered.
 * \return Modified Julian dates are Julian dates - 2400000.5
 */

double Subs::Time::mjd() const{
    return Date::mjd() + hour()/24.;
}

/** Returns the Julian epoch corresponding to the Time
 */

double Subs::Time::jepoch() const{
    // Back to sla 25/11/07, TRM
    return slaEpj(mjd());
}


/** Returns Terrestrial Time (TT) minus the UTC
 * \return Returns TT-UTC in seconds
 */
double Subs::Time::dtt() const {
    // Changed to pure SLA 23/11/2007 TRM
    return slaDtt(mjd());
}

/** Returns Terestial Time (TT) in the form of MJD
 * \return Returns TT as modified Julian date.
 */
double Subs::Time::tt() const {
    // Changed to pure SLA 23/11/2007 TRM
    double mj = this->mjd();
    return mj + slaDtt(mj)/86400.;
}

/** Returns Barycentric Dynamical Time (TDB) minus TT in seconds
 * \param  tel the telescope where the observations were taken.
 * \return Returns TDB-TT in seconds
 */
double Subs::Time::dtdb(const Subs::Telescope& tel) const {

    // Switched back to pure SLA. Note that old
    // wl had reverse sign, I assume because of the SOFA
    // routine used in place of slaRCC. TRM, 23/11/2007
    double wl, u, v;
    slaGeoc(tel.latituder(),tel.height(),&u, &v);
    wl  = -tel.longituder(); // Longitude, radians west
    u  *= Constants::AU/1000.0; // km from spin axis
    v  *= Constants::AU/1000.0; // km north of equator
    return slaRcc(tt(),hour()/24.,wl,u,v);

}

/** Returns Barycentric Dynamical Time (TDB) in the form of MJD (TDB). Note
 * this does not apply any light travel time corrections.
 * \param  tel the telescope where the observations were taken.
 * \return Returns TDB as an MJD
 */
double Subs::Time::tdb(const Subs::Telescope& tel) const {

    // Switched back to pure SLA. Note that old
    // wl had reverse sign, I assume because of the SOFA
    // routine used in place of slaRCC. TRM, 23/11/2007
    double wl, u, v;
    slaGeoc( tel.latituder(), tel.height(), &u, &v);
    wl  = -tel.longituder(); // Longitude, radians west
    u  *= Constants::AU/1000.0; // km from spin axis
    v  *= Constants::AU/1000.0; // km north of equator
    double tti = tt();
    return tti + slaRcc( tti, hour()/24., wl, u, v)/86400.0;

}

/** Computes the position of the observatory in barycentric
 * coordinates with respect to the barycentric reference frame. This is the
 * precursor to finding light-travel time corrections and is based upon the
 * SLA routine 'slaEpv'. The quoted error between this routine and the JPL
 * DE405 ephemeris is 4.6km RMS, 13.4 km max over the period 1900 to 2100.
 * \param tel the telescope where the observations were taken.  
 * \return Barycentric position, units of metres, relative to the BCRS
 */
Subs::Vec3 Subs::Time::earth_pos_bar(const Telescope& tel) const {

    double td = tdb(tel);
    
    // Compute position and velocity of Earth's centre in the BCRS
    // Switched to this routine, 23/11/2007, TRM
    double ph[3], pb[3], vh[3], vb[3];
    slaEpv(td, ph, vh, pb, vb);

    Vec3 position(pb);
    
    // Add on extra part from Earth's centre to observatory. Note
    // that slaGmst should take UT1 not UT but the error is only
    // about 1 second at most and no part of the Earth travels 
    // far in 1 second. 
    double last  = slaGmst(mjd()) + tel.longituder() + slaEqeqx(td);
    double pv[6];
    
    // this routine from slalib gives the telescope position and velocity relative 
    // to the centre of the Earth in a coordinate frame defined by the true equator
    // and equinox corresponding to 'last'.
    slaPvobs(tel.latituder(), tel.height(), last, pv);

    // Convert this position to the BCRS.
    // First compute the precession-nutation matrix to go from GCRS to true equator
    // equinox equivalent to 'td', then applies in iverse form with slaDimxv
    double rnpb[3][3];
    slaPneqx(td, rnpb);
    slaDimxv(rnpb, pv, pv);

    Vec3 pextra(pv);
    position += pextra;
    position *= Constants::AU;

    return position;

}

/** Computes the position and velocity of the Earth in heliocentric coordinates. 
 * This is the precursor to finding light-travel time corrections and is based 
 * is based upon the SLA routine 'slaEpv'. The quoted error between this
 * routine and the JPL DE405 ephemeris is 4.6km RMS, 13.4 km max over the period 1900 to 2100.
 * \param tel the telescope where the observations were taken.  
 * \return Heliocentric position, units of metres, relative to the BCRS
 */
Subs::Vec3 Subs::Time::earth_pos_hel(const Telescope& tel) const {

    double td = tdb(tel);
    
    // Compute position and velocity of Earth's centre in the BCRS
    // Switched to this routine, 23/11/2007, TRM
    double ph[3], pb[3], vh[3], vb[3];
    slaEpv(td, ph, vh, pb, vb);

    Vec3 position(ph);
    
    // Add on extra part from Earth's centre to observatory. Note
    // that slaGmst should take UT1 not UT but the error is only
    // about 1 second at most and no part of the Earth travels 
    // far in 1 second. 
    double last  = slaGmst(mjd()) + tel.longituder() + slaEqeqx(td);
    double pv[6];
    
    // this routine from slalib gives the telescope position and velocity relative 
    // to the centre of the Earth in a coordinate frame defined by the true equator
    // and equinox corresponding to 'last'.
    slaPvobs(tel.latituder(), tel.height(), last, pv);

    // Convert this position to the BCRS.
    // First compute the precession-nutation matrix to go from GCRS to true equator
    // equinox equivalent to 'td', then applies in iverse form with slaDimxv
    double rnpb[3][3];
    slaPneqx(td, rnpb);
    slaDimxv(rnpb, pv, pv);

    Vec3 pextra(pv);
    position += pextra;
    position *= Constants::AU;

    return position;

}

/** Computes the position and velocity of the observatory in barycentric and heliocentric
 * coordinates with respect to the barycentric reference frame. This is based on the
 * same principles as \c earth_pos_bar and \c earth_pso_hel, but is faster than 
 * separate calls to each of those and gives velocities too.
 * \param tel the telescope where the observations were taken.  
 * \param ph position of Earth wrt heliocentre, units of metres, returned.
 * \param vh velocity of Earth wrt heliocentre, units of metres/sec, returned.
 * \param pb position of Earth wrt barycentre, units of metres, returned.
 * \param vb velocity of Earth wrt barycentre, units of metres/sec, returned.
 */
void Subs::Time::earth(const Telescope& tel, Vec3& ph, Vec3& vh, Vec3& pb, Vec3& vb) const {

    double td = tdb(tel);
    
    // Compute position of observatory relative to centre of Earth. Note
    // that slaGmst should take UT1 not UT but the error is only
    // about 1 second at most and no part of the Earth travels 
    // far in 1 second. 
    double last  = slaGmst(mjd()) + tel.longituder() + slaEqeqx(td);
    double pv[6];
    
    // this routine from slalib gives the telescope position and velocity relative 
    // to the centre of the Earth in a coordinate frame defined by the true equator
    // and equinox corresponding to 'last'.
    slaPvobs(tel.latituder(), tel.height(), last, pv);

    // Convert the position and velocity to the BCRS.
    // First compute the precession-nutation matrix to go from GCRS to true equator
    // equinox equivalent to 'td', then applies in iverse form with slaDimxv
    double rnpb[3][3];
    slaPneqx(td, rnpb);
    slaDimxv(rnpb, pv,   pv);
    slaDimxv(rnpb, pv+3, pv+3);

    Vec3 padd(pv), vadd(pv+3);
    
    // Compute position and velocity of Earth's centre in the BCRS
    double tph[3], tpb[3], tvh[3], tvb[3];
    slaEpv(td, tph, tvh, tpb, tvb);

    ph.set(tph);
    vh.set(tvh);
    pb.set(tpb);
    vb.set(tvb);

    ph += padd;
    vh += vadd;
    pb += padd;
    vb += vadd;

    ph *= Constants::AU;
    vh *= Constants::AU;
    pb *= Constants::AU;
    vb *= Constants::AU;
    vh *= Constants::AU;

}

std::string Subs::Time::str() const {
    std::ostringstream ost;
    HMS t = hms();
    ost << Date::str() << ", " 
        << std::setw(2) << std::setfill('0') << t.hour << ":" 
        << std::setw(2) << std::setfill('0') << t.min  << ":" 
        << std::setw(2) << std::setfill('0') << t.sec  << "." 
        << std::setw(3) << std::setfill('0') << int(floor(1000.*t.fsec+0.5))
        << std::setfill(' ');
    return ost.str();
}

void Subs::Time::set_hour(double hour){
  
    if(hour < 0. || hour >= 24.)
        throw Time_Error("Subs::Time::set_hour(double): hours out of range");

    hour_ = hour;
}

void Subs::Time::set(int day_, Date::Month month_, int year_, int hour, int min, double sec){
  
    Date::set(day_,month_,year_);
    try{
        valid_time(hour,min,sec);
    }
    catch(const Time_Error& msg){
        throw Time_Error(std::string("void Subs::Time::set(int, Date::Month month_, int, int, int, double): ") + msg);
    }
    hour_  = hour + min/60. + sec/3600.;
}

/** Sets a time from individual components.
 * \param day_   the day of the month
 * \param month_ the month of the year
 * \param year_  the year
 * \param hour   the hour (0-23)
 * \param min    the number of minutes (0-59)
 * \param sec    the number of seconds (0-59.9999..)
 */
void Subs::Time::set(int day_, int month_, int year_, int hour, int min, double sec){
  
    Date::set(day_,month_,year_);
    try{
        valid_time(hour,min,sec);
    }
    catch(const Time_Error& msg){
        throw Time_Error(std::string("void Subs::Time::set(int, int, int, int, int, double): ") + msg);
    }
    hour_  = hour + min/60. + sec/3600.;
}

/** Sets a time from a string of the form "Date, 13:04:56.345" where
 * "Date" is a supported string of the Subs::Date class. 
 * \param time the string containing the time.
 */
void Subs::Time::set(const std::string& time){
    Date::set(time);
    std::istringstream ist(time.substr(11));
    char c;
    int hour_s, mn;
    double sc;
    ist >> c >> hour_s >> c >> mn >> c >> sc;
    if(!ist) throw Time_Error("Subs::Time::set(const std::string&, float): error reading in time = " + time);
    valid_time(hour_s,mn,sc); 

    hour_  = hour_s + mn/60. + sc/3600.;
}

void Subs::Time::set(double mjd){
    int imjd = int(floor(mjd));
    Date::set(imjd);

    hour_ = 24.*(mjd-imjd);
}

void Subs::Time::set(const Date& date, double hour){
    *this  = date;
    hour_  = hour;
}

// set to UTC

void Subs::Time::set(){
    time_t t;
    time(&t);
    tm* gt = gmtime(&t);
    Date::set(short(gt->tm_mday),Month(gt->tm_mon+1),1900+gt->tm_year);
    hour_  = gt->tm_hour+gt->tm_min/60.+gt->tm_sec/3600.;
}

void Subs::Time::add_hour(double hour){
    set(mjd() + hour/24.);
}

void Subs::Time::add_second(double second){
    set(mjd() + second/24./3600.);
}

/** Returns Greenwich Mean Sidereal Time (GMST) in radians.
 * Not quite right here in that I am using UTC rather than UT1 so this could
 * be up to 0.9 seconds out.
 */

double Subs::Time::GMST() const {
    return slaGmst(mjd());
}

/** Converts from GMT to GST */

void Subs::Time::GMTtoGST(){
    Date tdate(1,Jan,year());
    int mjs  = tdate.mjd();
    int days = Date::mjd() - mjs; 

    double t = (mjs - 15019.5)/36525.;
    double r = 6.6460656 + (0.00002581*t+2400.051262)*t;
    double b = 24.-(r-24.*(year()-1900));

    double t0 = 0.0657098*days-b;
    t0 += 1.002738*hour();
    t0  = t0 > 24. ? t0 - 24. : t0;
    t0  = t0 < 0.  ? t0 + 24. : t0;
    set_hour(t0);
}

void Subs::Time::read(std::ifstream& s, bool swap_bytes){
    Date::read(s, swap_bytes);

    double h;
    s.read((char*)&h,sizeof(double));
    if(!s) throw Time_Error("void Subs::Time::read(std::ifstream&): error reading hour of day");
    if(swap_bytes) h = Subs::byte_swap(h);
    if(h < 0. || h >= 24.) 
        throw Time_Error("hour out of range in Subs::Time::read(std::ifstream&)");

    hour_ = h;
}

void Subs::Time::write(std::ofstream& s) const {
    Date::write(s);
    s.write((char*)&hour_,sizeof(double));
    if(!s) throw Time_Error("void Subs::Time::write(std::ofstream&) const: error writing hour of day");
}

void Subs::Time::skip(std::ifstream& s){
    Date::skip(s);
    s.ignore(sizeof(double));
    if(!s) throw Time_Error("void Subs::Time::skip(std::ifstream&): error skipping hour of day");
}

/** Inputs a zero padded time e.g. "17 Nov 1961, 01:03:04.02"
 * \relates Subs::Time
 */
std::istream& Subs::operator>>(std::istream& ist, Subs::Time& time){

    if(!ist) return ist;

    Date date;
    ist >> date;
    if(!ist) return ist;

    char c;
    std::string hms;
    int hour_s, mn;
    double sc;

    ist >> c >> hour_s >> c >> mn >> c >> sc;  
    if(!ist) return ist;

    try{
        Time::valid_time(hour_s,mn,sc);
    }
    catch(...){
        ist.setstate(std::ios::failbit);
        return ist;
    }

    double hour  = hour_s + mn/60. + sc/3600.;
    time.set(date, hour);
    return ist;
}

/** Outputs a time with zero padding e.g. 01:03:04.02
 * \relates Subs::Time
 */
std::ostream& Subs::operator<<(std::ostream& ost, const Time& time){

    if(!ost) return ost;

    Subs::Time::HMS t = time.hms();
    ost.setf(std::ios::right);
    ost << (Date&)time << ", "
        << std::setw(2) << std::setfill('0') << t.hour << ":" 
        << std::setw(2) << std::setfill('0') << t.min  << ":" 
        << std::setw(2) << std::setfill('0') << t.sec  << "." 
        << std::setw(5) << std::setfill('0') << int(floor(100000.*t.fsec))
        << std::setfill(' ');
    if(!ost)
        throw Time::Time_Error("ostream& Subs::operator<<(std::ostream&, Subs::Time&): error output");

    return ost;
}

// defines allowable times

void Subs::Time::valid_time(int hour, int minute, double second){
    if(hour < 0 || hour > 23) 
        throw Time_Error("void Subs::Time::valid_time(int, int, double): hour = " + Subs::str(hour) + " is out of allowable range 0 to 23");

    if(minute < 0 || minute > 59)
        throw Time_Error("void Subs::Time::valid_time(int, int, double): minute = " + Subs::str(minute) + " is out of allowable range 0 to 59");
  
    if(second < 0. || second >= 60.)
        throw Time_Error("void Subs::Time::valid_time(int, int, double): second = " + Subs::str(second) + " is out of allowable range 0. to <60.");
}

bool Subs::operator>(const Time& t1, const Time& t2){
    const Date &d1 = t1, &d2 = t2;
    if(d1 > d2){
        return true;
    }else if(d1 == d2){
        return (t1.hour_ > t2.hour_);
    }else{
        return false;
    }
}

bool Subs::operator<(const Time& t1, const Time& t2){
    const Date &d1 = t1, &d2 = t2;
    if(d1 < d2){
        return true;
    }else if(d1 == d2){
        return (t1.hour_ < t2.hour_);
    }else{
        return false;
    }
}

bool Subs::operator>=(const Time& t1, const Time& t2){
    const Date &d1 = t1, &d2 = t2;
    if(d1 > d2){
        return true;
    }else if(d1 == d2){
        return (t1.hour_ >= t2.hour_);
    }else{
        return false;
    }
}

bool Subs::operator<=(const Time& t1, const Time& t2){
    const Date &d1 = t1, &d2 = t2;
    if(d1 < d2){
        return true;
    }else if(d1 == d2){
        return (t1.hour_ <= t2.hour_);
    }else{
        return false;
    }
}

bool Subs::operator==(const Time& t1, const Time& t2){
    const Date &d1 = t1, &d2 = t2;
    return (d1 == d2 && t1.hour_ == t2.hour_);
}

double Subs::operator-(const Time& t1, const Time& t2){
    const Date &d1 = t1, &d2 = t2;
    return (86400.*(d1-d2) + 3600.*(t1.hour_ - t2.hour_));
}








