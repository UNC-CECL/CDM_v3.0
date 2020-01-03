/******************************************************************************
 $Id: flux_stationary.h,v 1.14 2005/04/13 17:47:53 duran Exp $
 ******************************************************************************/

#ifndef SHORE_H
#define SHORE_H

#include "vec.h"
#include "func.h"
#include "globals.h"

class dunepar;

/*!  \brief Class defining the main shore processes.  */

class shore3d : public dunedata
{
public:
    shore3d( const dunepar& parameters );
    virtual ~shore3d();

    /*!  Returns m_shoreline.  */
    virtual double shorelinepos() { return m_shoreline; }

    virtual void shorelinecal(const TFktScal& h);
    virtual void restoreshoreface(TFktScal& h);
    virtual void restoreberm(TFktScal& h);
    virtual int shorefacemotion(TFktScal& h, double timestep, double time);
    //virtual void gradalongshore( const dunepar& parameters );
    
    virtual void save_arrays();
    
protected:
    
    bool m_x_periodic, m_y_periodic;
    
    double m_watertable, m_shore_HMWL, m_slope, m_shoreline, m_swash_normal;

    double m_sealevelrise;
    //double m_shoreface_lenght;
    double m_shoreface_aspectratio;
    double m_grad_alongshore;
    double m_grad_alongshore_sim_stepsize;
    string m_grad_alongshore_type, m_filename;
    double m_grad_alongshore_array [500];

};


#endif // SHORE_H
