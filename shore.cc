/******************************************************************************
 $Id: flux_stationary.cc,v 1.27 2005/04/13 18:33:12 duran Exp $
 ******************************************************************************/

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <limits>

#include "physics_const.h"
#include "globals.h"
#include "shore.h"



//*****************************************************************************
//  class flux3d_stationary

/*!  Constructor: reads relevant parameters and creates auxiliary arrays for
 the velocity and the density of entrained sand rho and its gradient.  */

shore3d::shore3d( const dunepar& parameters ): dunedata(parameters)
{
    
    m_x_periodic= duneglobals::periodic_x();
    m_y_periodic= duneglobals::periodic_y();
    
    //!! BEACH
    m_shore_HMWL = duneglobals::HMWL(); //parameters.getdefault("shore.HMWL", 0.0);
    m_watertable = duneglobals::MSL(); //parameters.getdefault("shore.sealevel", 0.0);
    m_slope = duneglobals::slope();//parameters.getdefault("beach.angle", 0.0);
    m_swash_normal = 0.5; //duneglobals::swash_norm();

        
    m_sealevelrise = parameters.getdefault("shore.sealevelrise", 0.0); // m/yr
    m_sealevelrise /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec
    
//    m_shoreface_lenght = parameters.getdefault("shore.facelenght", 1000.0);
    /* m_grad_alongshore = parameters.getdefault("shore.alongshore_grad", 0.0); // m/yr
    m_grad_alongshore *= m_slope;
    m_grad_alongshore /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec
    */

    // shoreline change parameters
    m_grad_alongshore_type = parameters.getrequired<string>("Init-shore.alongshore_grad");
    if(m_grad_alongshore_type.compare("constant") == 0) {
		m_grad_alongshore = parameters.getdefault("shore.alongshore_grad", 0.0); // m/yr
		m_grad_alongshore *= m_slope;
		m_grad_alongshore /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec
	} else if (m_grad_alongshore_type.compare("init_h") == 0) {
		string m_filename = parameters.getrequired<string>("shore.alongshore_grad.init_h.file" );
		double m_grad_alongshore_stepsize = parameters.getdefault("shore.alongshore_grad_timestep", 1.0); // units: years
		double m_grad_alongshore_nsteps = parameters.getdefault("shore.alongshore_grad_nsteps", 50); // units: years
		m_grad_alongshore_sim_stepsize = m_grad_alongshore_stepsize*duneglobals::secyear()*duneglobals::timefrac(); //units: seconds
		// load array
		ifstream file(m_filename);
			if(file.is_open())
			{
				int n = m_grad_alongshore_nsteps;

				for(int i = 0; i < n; ++i)
				{
					file >> m_grad_alongshore_array[i]; //read in value from file
					m_grad_alongshore_array[i] *= m_slope;
					m_grad_alongshore_array[i] /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec
				}
			} else {
				cout << "gradalongshore::create: FATAL ERROR: init_h file `" << m_filename << "\' does not exist for parameter `Init-shore.alongshore_grad\' !" << endl;
				exit(1);
			}

	} else {
		cout << "gradalongshore::create: FATAL ERROR: Unknown value `" << m_grad_alongshore_type << "\' for parameter `Init-shore.alongshore_grad\' !" << endl;
		cout << "  Valid values are: constant, init_h" << endl;
		exit(1);
	}

    m_shoreline = 0;

    // 
    
}


/*!  The destructor currently does nothing.  */
shore3d::~shore3d()
{
}

void shore3d::shorelinecal(const TFktScal& h)
{
    // CALCULATION OF THE MIN SHORELINE POSITION
    double shorelinepos = duneglobals::nx();
    int xaux = 0;
    for(int y = 0; y< duneglobals::ny(); ++y ){
        for(int x = 0; x< duneglobals::nx(); ++x ){
            if (h(x, y) >= m_shore_HMWL){
                xaux = x - 1;
                break;
            }
        }
        // take smaller position
        shorelinepos = (shorelinepos > xaux ? xaux : shorelinepos);
    }   
    m_shoreline = shorelinepos; // store shoreline position
}

// correct shoreline geometry
void shore3d::restoreshoreface(TFktScal& h)
{
    int x0 = 0;//m_shoreline-5;
    for( int y= 0; y< duneglobals::ny(); ++y ){
        for( int x= x0; x < m_shoreline; ++x ){
            h(x, y) = h(x0, y) + m_slope*duneglobals::dx()*(x-x0);
        }
    }
}

void shore3d::restoreberm(TFktScal& h)
{
	int x0 = m_shoreline;
	int x1 = (int) ceil(m_swash_normal/m_slope);
    for( int y= 0; y< duneglobals::ny(); ++y ){
        for( int x= x0; x < x0+x1+1; ++x ){
            h(x, y) = (h(x, y) < (h(x0, y) + m_slope*duneglobals::dx()*(x-x0)) ? h(x0, y) + m_slope*duneglobals::dx()*(x-x0) : h(x,y));
        }
        bool cont = true;
        for( int x= (x1+1); (x < (x1+(x1-x0))) && cont; ++x ){
        	// if next point is at an elevation lower than berm peak - slope*dist
        	if(h(x, y) < h(x1, y) - m_slope*duneglobals::dx()*(x-x1)) {
        		h(x, y) = h(x1, y) - m_slope*duneglobals::dx()*(x-x1);
        	} else {
        		cont = false;
        	}
		}
    }
}

/*
void shore3d::gradalongshore( const dunepar& parameters)
{
	string m_grad_alongshore_type = parameters.getrequired<string>("Init-shore.alongshore_grad");
	if(m_grad_alongshore_type == "constant") {
	    double m_grad_alongshore = parameters.getdefault("shore.alongshore_grad", 0.0); // m/yr
	    m_grad_alongshore *= m_slope;
	    m_grad_alongshore /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec
	} else if (m_grad_alongshore_type == "init_h") {
		string m_filename = parameters.getrequired<string>("shore.alongshore_grad.init_h.file" );
		double m_grad_alongshore_stepsize = parameters.getdefault("shore.alongshore_grad_timestep", 1.0); // units: years
		double m_grad_alongshore_nsteps = parameters.getdefault("shore.alongshore_grad_timestep", 50); // units: years
		double m_grad_alongshore_sim_stepsize = m_grad_alongshore_stepsize*duneglobals::secyear()*duneglobals::timefrac(); //units: seconds
		// load array
		ifstream file(m_filename);
			if(file.is_open())
			{
				int n = m_grad_alongshore_nsteps;
				double m_grad_alongshore_array[n];

				for(int i = 0; i < n; ++i)
				{
					file >> m_grad_alongshore_array[i]; //read in value from file
					m_grad_alongshore_array[i] *= m_slope;
					m_grad_alongshore_array[i] /= duneglobals::secyear()*duneglobals::timefrac(); // convert to m/sec
				}
			}

	} else {
		cout << "gradalongshore::create: FATAL ERROR: Unknown value `" << m_grad_alongshore_type << "\' for parameter `Init-shore.alongshore_grad\' !" << endl;
		cout << "  Valid values are: constant, init_h" << endl;
		exit(1);
	}
}
*/

int shore3d::shorefacemotion(TFktScal& h, double timestep, double time)
{
	double m_grad_alongshore_current;
	if(m_grad_alongshore_type.compare("constant") == 0) {
		m_grad_alongshore_current = m_grad_alongshore;
	} else if (m_grad_alongshore_type.compare("init_h") == 0) {
		int i = (int) floor(time / m_grad_alongshore_sim_stepsize);
		m_grad_alongshore_current = m_grad_alongshore_array[i];
	} else {
		m_grad_alongshore_current=0;
	}
    int mean_shoreshift = 0;
    double shoreshift_rate = m_sealevelrise + m_grad_alongshore_current;
    
    // CONSTANT SLOPE!
    int xend = (m_sealevelrise > 0 ? duneglobals::nx() : m_shoreline);
    for( int y= 0; y< duneglobals::ny(); ++y ){
        for( int x= 0; x < xend; ++x ){
            
            h(x, y) += - shoreshift_rate * timestep; // * (h(x, y) < m_shore_HMWL ? 1 : 1);
            
        }
        if (shoreshift_rate > std::numeric_limits<double>::epsilon()) {
            // EROSION
            mean_shoreshift = (h(1,0) < 0 ? 1 : 0);
            
            // correct shoreline geometry
            restoreshoreface(h);
            
        } else {
            // ACCRETION
            mean_shoreshift = (h(0,0) >= m_slope * duneglobals::dx() ? -1 : 0);

            // correct shoreline geometry
            restoreshoreface(h);
        }
    }
    // CORRECT SHORELINE
    
    double shift = h(0,0)/m_slope;
        
    cout << h(0,0) << " " << m_shoreline << " " << m_slope << ' ' << mean_shoreshift << ' ' << m_grad_alongshore_current << " # SHORE 1" << endl;
    
    return mean_shoreshift;
}


/*!  Saves the arrays  */

void shore3d::save_arrays()
{
}

