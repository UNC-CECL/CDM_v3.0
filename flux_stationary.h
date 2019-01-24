/******************************************************************************
 $Id: flux_stationary.h,v 1.14 2005/04/13 17:47:53 duran Exp $
 ******************************************************************************/

#ifndef FLUX_STAT_H
#define FLUX_STAT_H

#include "vec.h"
#include "func.h"
#include "globals.h"

class dunepar;

/*!  \brief Class which calculates equilibrium sand flux from shear stress.
 For 3-dimensional simulations.
 
 The sand flux is assumed to be stationary, that is no time transients are
 taken into account.  No aerodynamic entrainment is assumed to occur.  Spatial
 transient effects, however, are accounted for.
 
 This class just implements the calculation of the effective grain velocity and
 the calculation of the flux from the sand density rho.  The calculation of rho
 is left to be implemented in a subclass to allow different algoriths, for
 instance for taking into account vegetation.  */

class flux3d_stationary : public dunedata
{
public:
    flux3d_stationary( const dunepar& parameters );
    virtual ~flux3d_stationary();
    
    virtual void calc( TFktScal& influx, TFktVec& flux, const TFktScal& h, const TFktScal& h_nonerod, const TFktVec& tau, TFktScal& gamma );
    
//    virtual void waterlevel_factor( const double time );

    virtual double Satflux_upwind(double u_star);
    virtual double SatLength_upwind(double u_star);
    
    virtual void save_arrays();
    
protected:
    bool m_x_periodic, m_y_periodic;
    
    void calc_u( double *u_x, double *u_y, const double& tau_x, const double& tau_y,
  				double gradh_x, double gradh_y );
    
    double u_star_at( double u_star);
    
    void calc_rho( const TFktScal& h, const TFktScal& h_nonerod, const TFktVec& tau, TFktScal& gamma );

    /*! Shear stress threshold : beach */
    double tau_t_factor( double h );
    /*! tau_t recovery lenght */
    bool m_calc_shore;
    double m_tau_t_L;
    /*! watertable distance to the MSL and watertable height*/
    double m_watertable, m_watertable0;
    double m_wind_factor;
    
    /*!  Density of saltating grains.  */
    TFktScal m_rho;
    TFktScal m_rhoaux;
    /*!  input density of saltating grains.  */
    TFktScal m_inrho;
    /*!  Effective grain velocity. (= velocity at height m_z1)  */
    TFktVec m_u;
    TFktVec m_flux_sat;
    TFktVec m_TempVec;
    
    /*!  Typical grain diametre (in metres).  */
    double m_d_grain;
    /*!  Fluid viscosity, needed for computing the drag coefficient if not given.  */
    double m_fluid_viscosity;
    /*!  Gravity acceleration.  Retained as a parameter for the time when we'll
     simulate dunes on Mars ;).  */
    double m_g;
    /*!  Diffusion constant.  Describes short-range smearing of the height
     profile by turbulence and other processes we don't model.  */
    double m_C_diff;
    /*!  Fluid entrainment coefficient (parameter beta in parameter file).  */
    double m_beta;
    /*!  Sustained entrainment coefficient (parameter gamma in parameter file).
     Proportional to linear coefficient of expansion of the number of grains
     dislodged by each impacting grain in terms of the air shear stress at
     ground level.  Determines spatial saturation transients.  */
    double m_gamma;
    /*!  True roughness length of the sand surface (in metres).  */
    double m_z0;
    /*!  Mean saltation height (in metres).  */
    double m_zm;
    /*!  Intermediate height between m_z0 and m_zm.  This is the height at which
     the wind speed is the effective wind velocity (see v_eff).  It has to be
     determined by fitting experimental data.  */
    double m_z1;
    /*!  Drag coefficient of a sand grain.  */
    double m_C_drag;
    /*!  Proportionality factor between the vertical component of the rebounding
     velocity and the horizontal velocity loss of a grain impacting on the sand
     bed.  To be fitted to experiment or computed from a splash function.  */
    double m_alpha;
    /*!  Shear velocity threshold for sustaining grain entrainment.  (Parameter
     u*t) */
    double m_u_star_t;
    /*!  Shear velocity threshold for initial grain entrainment by the fluid
     alone.  (Parameter u*_ft)  */
    double m_u_star_ft;
    
    /*! Settling velocity / 2*alpha acts as a characteristic relative wind velocity*/
    double m_DeltaU;
    
    
    /*!  Diffusion constant divided by grid spacing.  Needed frequently by
     calc_rho.  */
    double m_Cd_dx;
    /*! M = ln(Z/zm)/b(Z), ln(z1/z0) and b(z1), Used in calc_u*/
    double m_M, m_log_z1z0, m_b_z1;
    /*!  (2*alpha/g) and (gamma/(2*alpha/g)^2) and (beta/(2*alpha/g * tau_ft))  Used in calc_rho.  */
    double m_2alpha_g, m_gamma_2alpha_g2, m_beta_2alpha_g_tau_ft;
    /*!  Threshold shear stress derived from m_u_star_t.  */
    double m_tau_t;
    /*!  Threshold shear stress for direct entrainment derived from m_u_star_ft.  */
    double m_tau_ft;
    
private:  
    double solve_for_rho(const double u,
                         const double tau, const double splash, const double A, const double B, const double tau_t_fact);
};


#endif // FLUX_STAT_H
