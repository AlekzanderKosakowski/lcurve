#include "trm/plot.h"

/** Constructor to open a plot device
 * \param device name of plot device to open
 */
Subs::Plot::Plot(const std::string& device) {
	idev = cpgopen(device.c_str());
	if(idev < 1){
		idev = 0;
		throw Plot_Error("Failed to open plot device = " + device + " in Plot::Plot(const std::string&)");
	}else{
		devname = device;
	}
}

Subs::Plot::~Plot(){
	if(idev){
		cpgslct(idev);
		cpgclos();
	}
}


/** Closes a plot
 */
void Subs::Plot::close(){
	if(idev){
		cpgslct(idev);
		cpgclos();
		idev = 0;
	}
}


/** Opens a plot, closing it first if it is already
 * open. The focus of next plotting commands switches to
 * this one
 * \param device plot device to be opened
 */
void Subs::Plot::open(const std::string& device) {
	close();
	idev    = cpgopen(device.c_str());
	if(idev < 1){
		idev = 0;
		throw Plot_Error("Failed to open plot device = " + device + " in Plot::open(const std::string&)");
	}
	devname = device;
}

/** Moves the focus of the next plotting commands to the
 * plot 
 */ 
void Subs::Plot::focus() const {
	if(idev) cpgslct(idev);
}

/** Tests whether a plot is in focus or not
 */
bool Subs::Plot::is_in_focus() const {
	if(idev){
		int current;
		cpgqid(&current);
		return (current == idev);
	}else{
		return false;
	}
}

/** Sets the ranges of the plot panel
 * \param x_1 left-hand X limit
 * \param x_2 right-hand X limit
 * \param y_1 lower Y limit
 * \param y_2 upper Y limit
 * \param scale true for equal physical scales
 */    
void Subs::Plot::Panel::set(float x_1, float x_2, float y_1, float y_2, bool scale) {
	x1   = x_1;
	x2   = x_2;
	y1   = y_1;
	y2   = y_2;
	just = scale;
	stan = true;
}

/** Sets the ranges of the plot panel
 * \param x_1 left-hand X limit
 * \param x_2 right-hand X limit
 * \param y_1 lower Y limit
 * \param y_2 upper Y limit
 * \param scale true for equal physical scales
 * \param xv_1 left-hand X limit of viewport
 * \param xv_2 right-hand X limit of viewport
 * \param yv_1 lower Y limit of viewport
 * \param yv_2 upper Y limit of viewport
 */    
void Subs::Plot::Panel::set(float x_1, float x_2, float y_1, float y_2, bool scale, float xv_1, 
							float xv_2, float yv_1, float yv_2){
	x1   = x_1;
	x2   = x_2;
	y1   = y_1;
	y2   = y_2;
	just = scale;
	stan = false;
	xv1  = xv_1;
	xv2  = xv_2;
	yv1  = yv_1;
	yv2  = yv_2;
}

void Subs::Plot::Panel::focus() const {
	if(stan){
		cpgvstd();
	}else{
		cpgsvp(xv1, xv2, yv1, yv2);
	}
	if(just){
		cpgwnad(x1, x2, y1, y2);
	}else{
		cpgswin(x1, x2, y1, y2);
	}
}

/** Plots a graph with the UT axis running 21, 22, 23, 0, 1 etc if
 * it has to cross 24. Plot must be open, viewport, colours set etc beforehand.
 * \param t1 start time
 * \param t2 end time
 * \param y1 lower Y limit
 * \param y2 upper Y limit
 * \param sides true to plot axes along the sides
 * \param top true to plot axis along the top
 */

void Subs::ut_plot(float t1, float t2, float y1, float y2, bool sides, bool top){
	
	if(t2 > 24.){
		float xv1, xv2, yv1, yv2;
		cpgqvp(0,&xv1, &xv2, &yv1, &yv2);
		float xvt1, xvt2;

		// left hand box
		xvt2 = xv1 + (xv2-xv1)*(24.-t1)/(t2-t1);
		cpgsvp(xv1,xvt2,yv1,yv2);
		cpgswin(t1,23.9999,y1,y2);
		if(sides && top){
			cpgbox("BCNST",2.,2,"BNST",0.,0);
		}else if(sides){
			cpgbox("BNST",2.,2,"BNST",0.,0);
		}else if(top){
			cpgbox("BCNST",2.,2," ",0.,0);
		}else{
			cpgbox("BNST",2.,2," ",0.,0);
		}

		// right-hand box
		xvt1 = xvt2;
		cpgsvp(xvt1,xv2,yv1,yv2);
		cpgswin(0.,t2-24.,y1,y2);
		if(sides && top){
			cpgbox("BCNST",2.,2,"CST",0.,0);
		}else if(sides){
			cpgbox("BNST",2.,2,"CST",0.,0);
		}else if(top){
			cpgbox("BCNST",2.,2," ",0.,0);
		}else{
			cpgbox("BNST",2.,2," ",0.,0);
		}
		cpgsvp(xv1,xv2,yv1,yv2);
		cpgswin(t1,t2,y1,y2);
	}else{
		cpgswin(t1,t2,y1,y2);
		if(sides && top){
			cpgbox("BCNST",2.,2,"BNST",0.,0);
		}else if(sides){
			cpgbox("BNST",2.,2,"BNST",0.,0);
		}else if(top){
			cpgbox("BCNST",2.,2," ",0.,0);
		}else{
			cpgbox("BNST",2.,2," ",0.,0);
		}
	}
}

Subs::Plots::Plots(int nplot) : nplot_(nplot), nfocus_(0) {
	if(nplot < 1)
		throw Subs::Plot::Plot_Error("Subs::Plots(): need at least 1 plot");
	plot_ = new Plot[nplot];
}

/** Opens the plot now in focus, if it is not open already. If
 * the focus is set < 0, then it opens the first one that is not already
 * open.
 * \param device the plot device to open.
 */
void Subs::Plots::open(const std::string& device){
	if(nfocus_ >= 0){
		if(plot_[nfocus_].is_open())
			plot_[nfocus_].focus();
		else
			plot_[nfocus_].open(device);

	}else{
		bool ok = false;
		for(int i=0; i<nplot_; i++){
			if(!plot_[i].is_open()){
				plot_[i].open(device);
				nfocus_ = i;
				ok = true;
				break;
			}
		}
		if(!ok)
			throw Plot::Plot_Error("Subs::Plots::open(const std::string&): no more plots available. You should close some");
	}
}

/** Explicitly close a particular plot.
 * \param nclose the plot to close. -1 to close the current plot.
 */
void Subs::Plots::close(int nclose){
	if(nclose >= 0 && nclose < nplot_){
		plot_[nclose].close();
		nfocus_ = -1;
	}else if(nclose == -1){
		if(nfocus_ >= 0)
			plot_[nfocus_].close();
	}else{
		throw Plot::Plot_Error("Subs::Plots::close(int): plot number out of range -1 to " + Subs::str(nplot_-1));
	}
}

void Subs::Plots::close(){
	for(int i=0; i<nplot_; i++)
		plot_[i].close();
	nfocus_ = -1;
}

bool Subs::Plots::is_open() const {
	if(nfocus_ >= 0 && nfocus_ < nplot_){
		return plot_[nfocus_].is_open();
	}else{
		return false;
	}
}

void Subs::Plots::focus(){
	if(nfocus_ >= 0)
		plot_[nfocus_].focus();
}

void Subs::Plots::focus(int i){
	if(i >= 0 && i < nplot_){
		if(plot_[i].is_open())
			plot_[i].focus();
		nfocus_ = i;
	}else if(i < 0){
		nfocus_ = -1;
	}else{
		throw Plot::Plot_Error("void Subs::Plots::focus(int): attempt to focus on non existent plot number " + Subs::str(i+1));
	}
}
   
/** Lists the plots and those which are open
 * \param ostr the output stream to send the list to
 */
void Subs::Plots::list(std::ostream& ostr){
	for(int i=0; i<nplot_; i++){
		if(plot_[i].is_open()){
			ostr << "Plot " << i+1 << " was opened with device name = " << plot_[i].get_name() << std::endl;
		}else{
			ostr << "Plot " << i+1 << " is not open."  << std::endl;
		}
	}
} 

Subs::Plots::~Plots(){
	for(int i=0; i<nplot_; i++)
		plot_[i].close();
	delete[] plot_;
}


