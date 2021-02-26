/******************************************************************************
 $Id: dune_evolution.h,v 1.18 2005/04/26 14:26:54 duran Exp $
 ******************************************************************************/

#ifndef DUNE_EVOLUTION_H
#define DUNE_EVOLUTION_H

#include "func.h"
#include "evolution.h"

class dunepar;
class wind;
class arrayinit;
class shear;
class flux3d_stationary;
class avalanche;
class analyze;
class influx;
class shore3d;
class vegetation;
class storm;

///*!  Base class of all classes handling the top level of the simulation and
//  containing the height array.  */
//
//class dune_evolution : public evolution
//{
//public:
//  /*!  This constructor takes only pointers to an initialisation object and a
//    wind object as arguments which are written to the member variables m_init
//    and m_wind.  The type of the classes for shear stress and sand transport
//    calculations may depend on the derived class and is therefore left to be
//    defined in the subclass constructor argument list.  */
//  dune_evolution(const dunepar &par, arrayinit *i, wind *w) :
//    evolution(par), m_init(i), m_wind(w) {}
//  virtual ~dune_evolution() {}
//
//protected:
//  /*!  Height profile initialisation object.  */
//  arrayinit *m_init;
//  /*!  Object giving wind speed and direction.  */
//  wind *m_wind;
//};


/*!  Dune evolution in 3-dimensional simulation.  */

class dune_evol_3d : public evolution
{
public:
    dune_evol_3d(const dunepar& par);
    virtual ~dune_evol_3d();

    virtual void init(const dunepar& par);

    virtual void save_arrays();
    
    virtual double step_implementation();
    
protected:
    virtual double update_height(double halfmeanLength);
    virtual void volume_correction(double timestep);
    virtual int shiftback(int shift);

    virtual void update_dhdt();

    /*!  Object giving wind speed and direction.  */
    wind *m_wind;
    /*!  Shear stress calculation object.  */
    shear *m_calcshear;
    /*!  Object which sets the sand influx.  */
    influx *m_influx;
    /*!  Sand flux calculation object.  */
    flux3d_stationary *m_calcflux;
    /*!  Avalanche relaxation object.  */
    avalanche *m_avalanche;
    /*!  Shore calculation object.  */
    shore3d *m_shore;
    /*!  vegetation calculation object.  */
    vegetation *m_grass;
    /*!  storm calculation object.  */
    storm *m_storm;
    
    /*!  Analyze object.  */
    analyze *m_analyze;
    
    bool m_x_periodic, m_y_periodic;
    /*!  Maximal time step for step_implementation().  */
    double m_dtmax;
    /*!  If !=0, the number of iterations before the height profile is copied to
     the non-erodable height profile.  Parameter update.fix_every in parameter
     file.  */
    int m_fix_every;
    /*!  If true, the height profile is moved back if necessary so that its
     centre of mass stays in approximately the same place.  This is the
     parameter calc.shift_back in the parameter file.  */
    bool m_shiftback;
    /*! If true, the reference position of the dune is shifted back to the initial position of its center of mass, otherwise it is shifted back to the center of the field*/
    bool m_shiftback_target_cm0;
    /*! If true, the reference position of the dune to be shifted back is the actual center of mass, otherwise it is a new center calculated by the analyze.center function*/
    bool m_shiftback_centrex_cmT;
    /*!  If true, the volume conservation is fixed by hand, otherwise an small error appears.  This is the
     parameter calc.volume.correction in the parameter file.  */
    bool m_vol_correct;
    /*! If true, the timestep, heigth, width, large, volume, influx and outflux are computed and copyied in the
     'time.dat' file*/
    bool m_calc_analyze;
    /*! If true, the vegetation evolution is computed*/
    bool m_calc_veget;
    //bool m_calc_season0, m_calc_season;
    
    /*! Shore parameters*/
    double m_shoreline, m_shorelinechange;
    bool m_shoremotion, m_calc_shore, m_calc_berm;
    
    bool m_calc_storm0, m_calc_storm;
    double m_surge;
    
    /*!  Dune height profile.  */
    TFktScal m_h;
    TFktScal m_h0;  // store initial profile
    TFktScal m_hprev;  // store previous profile
    /*!  Non-erodable height profile under the sand.  */
    TFktScal m_h_nonerod;
    /*!  Height change per time.  */
    TFktScal m_dh_dt;
    /*!  Shear stress calculated by m_calcshear.  */
    TFktVec m_tau;
    /* to define where there is no erosion nor deposition */
    TFktScal m_gamma;
    /*!  Sand flux computed by m_calcflux.  */
    TFktVec m_flux;
    /* sand influx along wind direction x- */
    TFktScal m_flux_in;
    /* Vegetation cover */
    TFktScal m_rho_veget;
    /* Overwash */
    TFktScal m_overwash;

    /*!  Cumulative amount by which the height profile has been shifted.  This is
     the x/y component of a vector, where the x direction is the origin of
     the wind direction angle.  */
    double m_shift_dist_x, m_shift_dist_y;
    double m_wind_factor;
    double m_Satflux_upwind, m_maxchange, m_qin, m_qout;
};

#endif // DUNE_EVOLUTION_H

