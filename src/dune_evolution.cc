/******************************************************************************
 $Id: dune_evolution.cc,v 1.32 2005/04/26 14:26:53 duran Exp $
 ******************************************************************************/

#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "globals.h"
#include "dune_evolution.h"
#include "shear.h"
#include "wind.h"
#include "flux_stationary.h"
#include "initsurf.h"
#include "avalanche.h"
#include "influx.h"
#include "shear_hlr.h"
#include "vegetation.h"
#include "storm.h"
#include "analyze_new.h"

#include "func.h"
#include "shore.h"


//*****************************************************************************
//  class dune_evol_3d


/*!  The pointers to calculation objects are written to member variables.
 Space for all two-dimensional arrays is allocated, and m_h is initialised
 with \a i and m_h_nonerod with i_nonerod.  It is ensured that h is nowhere
 below h_nonerod.  */


dune_evol_3d::dune_evol_3d(const dunepar& par) : evolution(par)
{
    
    if( !duneglobals::sim3d() ) {
        cerr << "dune_evol_3d Constructor: FATAL ERROR: This class mustn't be used in 2D simulations!" << endl;
        exit(1);
    }
    m_calcshear= new shearHLR(par);
    m_wind= wind::create(par);
    m_influx= influx::create(par);
    m_calcflux= new flux3d_stationary(par);
    m_avalanche= avalanche::create(par);
    
    m_shore= new shore3d(par);
    m_grass= new vegetation(par);
    m_storm = new storm(par);
    
    m_analyze = new analyze(par);
    
    /* Selecting boundary conditions*/
    m_x_periodic= duneglobals::periodic_x();
    m_y_periodic= duneglobals::periodic_y();
    
    m_dtmax = par.getdefault("dt_max", 1000.0);
    m_shiftback = par.getdefault("calc.shift_back", true);
    m_shiftback_target_cm0= par.getdefault("shift_back.cm0", false);
    m_shiftback_centrex_cmT= par.getdefault("shift_back.cmT", true);

    m_shift_dist_x= m_shift_dist_y= 0.0;
    m_vol_correct = par.getdefault("calc.volume.correction", true);
    m_calc_analyze = par.getdefault("calc.analyze", true);
    m_calc_veget = par.getdefault("veget.calc", false);
    
    m_wind_factor = par.getdefault("wind.fraction", 1.0);
    
    /* Shore parameters */
    m_calc_shore = par.getdefault("calc.shore", true);
    m_calc_berm = par.getdefault("calc.berm", false);
    m_shoremotion = par.getdefault("shore.motion", 1.0);
    
    /* Storm parameters */
    m_calc_storm0 = par.getdefault("calc.storm", false);
    m_calc_storm = m_calc_storm0;   // temporal storm variable

    /*Functions creation*/
    m_h.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    m_h0.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    m_hprev.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    m_h_nonerod.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    m_dh_dt.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    m_tau.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), vec2(0.0, 0.0));
    m_flux.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), vec2(0.0, 0.0));
    m_flux_in.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    m_gamma.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    
    m_rho_veget.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);

    m_overwash.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);

    /* Veg parameters */
    /*m_calc_season0 = par.getdefault("calc.season", false);
    m_calc_season = m_calc_season0;   // temporal storm variable
    */

    /* Surface initialization*/
    init(par);
    m_grass->init(par);
    m_grass->getcover(m_rho_veget, evolution::time(), m_dtmax);

    /* Auxiliar */
    m_Satflux_upwind = (m_wind->u_star()==0? 0 : m_calcflux->Satflux_upwind(m_wind->u_star()));
    cout << "!!!Flux0 = " << m_Satflux_upwind/duneglobals::rho_sand()*duneglobals::secyear() << ", " << m_Satflux_upwind << endl;
    
    double SatLength_upwind = (m_wind->u_star()==0? 0 : m_calcflux->SatLength_upwind(m_wind->u_star()));
    cout << "!!!ls = " << SatLength_upwind << endl;
    
    if (m_calc_shore) {
        /* Shoreline calculation*/
        m_shore->shorelinecal( m_h );
        /* return shoreline */
        m_shoreline = m_shore->shorelinepos();
    } else {
        m_shoreline = 0;
    }
    /* init shoreline change in m */
    m_shorelinechange = 0;
    cout << "!!! shoreline = " << m_shoreline*duneglobals::dx() << endl;
}


/*!  Deletes all objects performing parts of the calculation, which were
 created in the constructor.  */

dune_evol_3d::~dune_evol_3d()
{
    delete m_avalanche;
    delete m_calcflux;
    delete m_influx;
    delete m_wind;
    delete m_calcshear;
    delete m_shore;
    delete m_analyze;
}

void dune_evol_3d::init(const dunepar& par){
    
    arrayinit *init_h;
    
    init_h= arrayinit::create(par);
    init_h->init_2d_scal( m_h );
    delete init_h;
    
    if( par.exists("nonerod.Init-Surf") )
        init_h= arrayinit::create(par, "nonerod.");
    else
        init_h= new CInitSurfPlain(0.0);
    
    init_h->init_2d_scal( m_h_nonerod );
    delete init_h;
    for( int x= 0; x< duneglobals::nx(); ++x )
        for( int y= 0; y< duneglobals::ny(); ++y )
            if( m_h(x, y) < m_h_nonerod(x, y) )
                m_h(x, y)= m_h_nonerod(x, y);

    m_h0 = m_h; // store initial profile
    
}

/*!  Saves the two-dimensional arrays m_h, m_tau and m_flux.  The difference
 between h and the non-erodable surface is computed and saved as h_deposit.  */

void dune_evol_3d::save_arrays()
{
    save_2d_scalarray( "dhdt", m_dh_dt );
    save_2d_scalarray( "h", m_h );
    save_2d_vecarray( "shear", m_tau );
    save_2d_vecarray( "flux", m_flux );
}


/*!  This is the function which computes the actual change in the surface
 profile.  The method for time evolution is forward Euler.  First, the shear
 stress is computed with the help of m_calcshear.  Then the sand flux is
 calculated with m_calcflux.  The rate of change of the height profile is
 calculated as the divergence of the flux divided by the sand density.  m_h is
 changed according to this rate with a time step which is at most m_dtmax
 (parameter salt.dt_h_max in the .par file) and which is small enough that the
 change in height is at most m_dhmax (parameter salt.dh_max).  */

double dune_evol_3d::step_implementation()
{
    double newwind, angle, timestep, halfmeanLength = 0, m_ustar0, factor=0;
    
    m_wind->advance( evolution::time() );
    newwind= m_wind->direction();
    m_ustar0 = m_wind->u_star();
    cout << "dune_evol_3d::step_implementation: step " << steps()+1 <<
	": wind direction " << newwind*360 << " deg., u* " << m_ustar0 << " m/s\n";
	
    if(m_ustar0 > 0.6*duneglobals::u_star_ft()){
        m_Satflux_upwind = (m_ustar0 > duneglobals::u_star_ft() ? m_calcflux->Satflux_upwind(m_ustar0) : 0.0);
        
        m_calcshear->set_ustar(m_ustar0);
        if(m_calc_veget) halfmeanLength = m_calcshear->Calc( m_h, m_tau, m_rho_veget );
        else halfmeanLength = m_calcshear->Calc( m_h, m_tau);
        
        m_influx->set(m_flux_in, m_flux, m_Satflux_upwind, angle );
        m_gamma.SetAll(0.0);
//       m_calcflux->waterlevel_factor(evolution::time()); // get temporal variation of wet-dry level
        m_calcflux->calc( m_flux_in, m_flux, m_h, m_h_nonerod, m_tau, m_gamma );
        
        m_hprev = m_h; // copy previous profile;
        timestep= update_height(halfmeanLength);
        m_avalanche->calc(m_h, m_h_nonerod);
        
        update_dhdt();
        
    }else{
        timestep = m_dtmax;
    }
    
    /*STORMS*/
    if (m_calc_storm0) {
        m_surge = 0;
        if (m_calc_storm) {
            m_surge = m_storm->impact( m_shoreline, evolution::time(), timestep, m_h, m_h_nonerod, m_overwash );
            /* no innundation condition */
            if (m_surge < 0)
            {
                /* code */
                m_h = m_h0; // restore initial profile
                m_surge *= -1.0;
                m_surge += 1000; // identify inundation
            }
        }
        m_storm->stop( evolution::time(), timestep, m_calc_storm);
    }
    
    // SHIFT AND SHORE
    if (m_calc_shore) {
        // Shoreline calculation
        // Shore motion
        int shoreshift = (m_shoremotion ? m_shore->shorefacemotion(m_h, timestep, evolution::time()) : 0);
        m_shorelinechange += shoreshift * duneglobals::dx();
        if(abs(shoreshift) > 0)
        	m_shift_dist_x = shiftback(shoreshift);
        // Restore shoreface
        m_shore->restoreshoreface( m_h );
        if(m_calc_berm)
        	m_shore->restoreberm( m_h );
    } else {
        m_shift_dist_x = shiftback(0);
    }
    
    // VEGETATION
    int m_veget_X0;
    if (m_calc_veget) {
    	/*if(m_calc_season0) {
			m_grass->season(evolution::time(), timestep, m_calc_season);
		}*/
        m_veget_X0 = m_grass->evol(m_rho_veget, evolution::time(), timestep, m_shoreline, m_h, m_dh_dt, m_overwash);
    }
    
    // PROCESS DATA
    if(m_calc_analyze){

        // Rescaling
        m_flux.rescale(m_Satflux_upwind);

        int steps = evolution::time() / timestep;
        int process = (steps % 1000 == 0 ? 1 : 0);
        if (process || m_surge > 0) {
            // // Rescaling
            // m_flux.rescale(m_Satflux_upwind);
            m_analyze->Calc(steps, evolution::time(), m_shift_dist_x*duneglobals::dx(), m_shoreline, m_shorelinechange, m_veget_X0, m_qin/m_Satflux_upwind/duneglobals::ny(), m_qout/m_Satflux_upwind/duneglobals::ny(), m_ustar0, m_surge, m_h, m_rho_veget);
            // prevent saving storms many times
        }
    }
    

    return timestep;
}



/*!  Computes the height profile change from the divergence of the sand flux. The return value is the time step.  */

double dune_evol_3d::update_height(double halfmeanLength)
{
    double timestep;
    int x, y;
    
    timestep = m_dtmax;

    m_maxchange= 0.0;
    for( y= 0; y< duneglobals::ny(); ++y ){
        for( x= 0; x< duneglobals::nx(); ++x ) {

            m_dh_dt(x, y) = -m_gamma(x, y) / duneglobals::rho_sand();

            if( fabs(m_dh_dt(x, y)) > m_maxchange )	m_maxchange= fabs(m_dh_dt(x, y));
        }
    }
    
    // Flux calculation
    m_qin=0.0;
    m_qout= 0.0;
    for( y=0; y< duneglobals::ny(); y++) {
        m_qin+= (m_shoreline>0? m_flux(m_shoreline,y)[0] : m_flux_in(0,y));
        m_qout+= m_flux(duneglobals::nx()-1,y)[0];
    }
    
    // Volume correction
    if(m_vol_correct && m_maxchange> 1e-10){
        volume_correction( timestep);
    }
    
    // height update
    // ---> periodic boundary conditions
    for( y= 0; y< duneglobals::ny(); ++y ){
        for( x= m_shoreline; x< duneglobals::nx(); ++x )
        {
            m_h(x, y) += timestep * (x < m_shoreline && m_dh_dt(x, y) < 0 ? 0 : m_dh_dt(x, y));
            // HIT ROCK
            if( m_h(x, y) < m_h_nonerod(x, y) ){
                // NO EROSION
                m_h(x, y) = m_h_nonerod(x, y);
            }

        }
    }
    
    return timestep;
}

// Update dhdt
void dune_evol_3d::update_dhdt()
{
    double timestep = m_dtmax;
    
    for( int y= 0; y< duneglobals::ny(); ++y ){
        for( int x= 0; x< duneglobals::nx(); ++x ) {
            m_dh_dt(x, y) = (m_h(x, y) - m_hprev(x, y)) / timestep;
        }
    }
}

void dune_evol_3d::volume_correction(double timestep)
{
    double dVReal, dV= 0, dVh= 0, dVh2= 0, dV_ex_correction= 0, Dex, dV_correction, dV_correction2, dx2= duneglobals::dx()*duneglobals::dx();
    int x,y;
    
  	dVReal = (m_qin - m_qout) * duneglobals::dx()/duneglobals::rho_sand();
	
    dV = m_dh_dt.Integrate(m_shoreline);
    
	for( y= 0; y< duneglobals::ny(); ++y ){
        for( x= m_shoreline; x< duneglobals::nx(); ++x ){
            Dex= (m_h_nonerod(x,y)-m_h(x,y))/timestep - m_dh_dt(x,y);
            if(Dex > 0) dV_ex_correction-= Dex;
            if(m_h(x,y) - m_h_nonerod(x,y) > timestep*m_maxchange && fabs(m_dh_dt(x,y))>1e-2*m_maxchange){
                if((dVReal-dV)*m_dh_dt(x,y) < 0)
                    dVh+= m_h(x,y) - m_h_nonerod(x,y);
                if(m_h(x,y) - m_h_nonerod(x,y) < 10.0*timestep*m_maxchange)
                    dVh2++;
            }
        }
    }
    dV_ex_correction*= dx2;
    dVh2*= dx2;
    dVh*= dx2;
    dV_correction = (dVh > 0 ? (dVReal-dV)/dVh : 0);
    dV_correction2 = (dVh2 > 0 ? dV_ex_correction/dVh2 : 0);
    for ( y=0; y < duneglobals::ny(); y++) {
        for ( x=m_shoreline; x < duneglobals::nx(); x++) {
            if(m_h(x,y) - m_h_nonerod(x,y) > timestep*m_maxchange && fabs(m_dh_dt(x,y))>1e-2*m_maxchange){
                if((dVReal-dV)*m_dh_dt(x,y) < 0)
                    m_dh_dt(x,y) += dV_correction * (m_h(x,y) - m_h_nonerod(x,y));
                
                if(m_h(x,y) - m_h_nonerod(x,y) < 10.0*timestep*m_maxchange)
                    m_dh_dt(x,y) += dV_correction2;
            }
        }
    }
}

/*!  Computes the shiftback.  */

int dune_evol_3d::shiftback(int shift)
{
    /*if (m_shiftback) {
        double m_centrex= m_h.CenterOfMassX()/duneglobals::dx();
        double m_targetcentrex= duneglobals::nx()*0.5;
        if (m_centrex - m_targetcentrex > 1) {
            //           m_shift_dist_x++;
            shift = 1;
        }
    }*/
    // Shift
    m_h.ShiftOne(shift);
    m_h_nonerod.ShiftOne(shift);
    if (m_calc_veget) {
        m_grass->shiftback(shift);
    }
    m_shift_dist_x += shift;
    
    return m_shift_dist_x;
}
