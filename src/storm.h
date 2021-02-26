/******************************************************************************
 $Id: flowbeach.h,v 1.5 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#ifndef STORM_H
#define STORM_H

#include "globals.h"
#include "func.h"
#include "PTG_FileNames.h"
#include "vec.h"

class dunepar;
class avalanche;

/*! Continuum aproach of sand relaxation by coasts (including coasts over solid substrates) */

class storm : public dunedata
{

public:
    // construction
    storm(const dunepar& p);
    virtual ~storm() {}
        
    void CalcGradUp(TFktScal &h);
    void Step(TFktScal &h, TFktScal &overwash, double &m_Tsurge_eff);

    virtual int dcrest_ident(TFktScal &h, int y);
    virtual double impact(double shoreline, double time, double timestep, TFktScal &h, TFktScal &h_nonerod, TFktScal &overwash );
    virtual void stop( double time, double timestep, bool &calc_storm);

    virtual void calc( TFktScal &h, TFktScal &overwash );
  
    virtual void save_arrays();
        
private:
 
    /*!  Avalanche relaxation object.  */
    avalanche *m_avalanche;

    TFktScal m_sflux;
    TFktScal m_hst;
    
    int m_storm_iter;
    double m_Smax;
    double m_Sdt;
    double m_Q, m_scalefactor, m_shape, m_scalefactor0, m_shore_HMWL, m_watertable, m_surge, m_Tsurge;
    double m_slope_storm, m_slope_storm_angle, m_shoreface_slope, m_repose_dyn, Kc, Kb, m_g, R_eff, Cls;
    double m_season_scalefactor, m_storm_season_start, m_storm_season_end;

    int m_shoreline;
    int m_storm_morphodynamics;
    int m_storm_seed;
    string m_storm_type;
    int m_tmax, m_dtmax;

    double surge[5000];
    int stormindex;

};


#endif

