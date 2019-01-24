/******************************************************************************
 $Id: flux_stationary.cc,v 1.27 2005/04/13 18:33:12 duran Exp $
 ******************************************************************************/

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <iostream>

#include "physics_const.h"
#include "globals.h"

#include "flux_stationary.h"


//*****************************************************************************
//  class flux3d_stationary

/*!  Constructor: reads relevant parameters and creates auxiliary arrays for
 the velocity and the density of entrained sand rho and its gradient.  */

flux3d_stationary::flux3d_stationary( const dunepar& parameters ): dunedata(parameters)
{
    
    m_rho.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );
    m_rho.SetAll(0.0);
	m_inrho.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );
	m_flux_sat.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );
    
    m_rhoaux.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );
    m_rhoaux.SetAll(0.0);
    
    m_u.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );
    
    m_x_periodic= duneglobals::periodic_x();
    m_y_periodic= duneglobals::periodic_y();
    
    // Get the parameters:
    m_d_grain = duneglobals::d_grain();
    m_u_star_ft = duneglobals::u_star_ft();
    m_u_star_t = duneglobals::u_star_t();
    m_tau_ft = duneglobals::tau_ft();
    m_tau_t = duneglobals::tau_t();
    
    m_DeltaU = duneglobals::DeltaU();
    m_M = duneglobals::M();
    m_b_z1 = duneglobals::b_z1();
    m_Cd_dx= duneglobals::Cd_dx();
    m_alpha = duneglobals::alpha();
    m_gamma = duneglobals::gamma();
    m_2alpha_g = duneglobals::alpha2_g();
    m_log_z1z0 = duneglobals::log_z1z0();
    m_gamma_2alpha_g2 = duneglobals::gamma_2alpha_g2();
    m_beta_2alpha_g_tau_ft = duneglobals::beta_2alpha_g_tau_ft();

    // time factor
    m_wind_factor =  duneglobals::timefrac();

    //!! BEACH
    m_calc_shore = parameters.getdefault("calc.shore", false);
    m_tau_t_L = parameters.getdefault( "beach.tau_t_L", 0.05);  // beach param!!!
    m_watertable0 = duneglobals::MSL();//parameters.getdefault("shore.sealevel", 0.0);
    m_watertable = m_watertable0;
}

/*!  The destructor currently does nothing.  */

flux3d_stationary::~flux3d_stationary()
{
}

double flux3d_stationary::Satflux_upwind(double u_star)
{
    double u_s = 2.5 *u_star *(m_log_z1z0 - (1.- u_star_at(u_star)/u_star)*m_b_z1) - m_DeltaU;
    double flux = m_2alpha_g * m_tau_t*((u_star/m_u_star_t)*(u_star/m_u_star_t) - 1) * u_s;
    return flux;
}

double flux3d_stationary::SatLength_upwind(double u_star)
{
    double u_s = 2.5 *u_star *(m_log_z1z0 - (1.- u_star_at(u_star)/u_star)*m_b_z1) - m_DeltaU;
    double L_sat = m_2alpha_g * u_s*u_s /(m_gamma * ((u_star/m_u_star_t)*(u_star/m_u_star_t) - 1));
    return L_sat;
}

double flux3d_stationary::u_star_at( double u_star )
{
    double u_t = m_M * m_u_star_t - (m_M - 1.) * u_star;
    return u_t;//(u_t > 0 ? u_t : 0);
}

// void flux3d_stationary::waterlevel_factor( const double time )
// {
//     double period = 60*60*24*365. * m_wind_factor;  // convert to years
    
//     // growth step function
//     double factor = (sin(2*M_PI * time / period) > 0 ? 1 : 0);
     
//     m_watertable = m_watertable0 + 0.3*factor;
// }

double flux3d_stationary::tau_t_factor( double h )
{
	double tau_t = (m_calc_shore ? 10.+(1.-10.)*(1-exp(-(h > m_watertable ? h - m_watertable : 0) / m_tau_t_L)) : 1);
	return tau_t;
}

/*!  Calculates the effective wind velocity.  The drag force under this
 velocity is the mean force affecting the sand grains in the air.  The
 effective velocity is the wind speed at the height m_z1.
 
 In addition to the simplifications used for all this class, this calculation
 neglects the influence of even spatial transient effects on v_eff since they
 were seen to be small (Veit Schwaemmle, Diploma Thesis).
 
 V_eff is calculated according to the formula:
 
 v_eff = tau/|tau| * ( sqrt( [4/kappa^2 * z1/zm * 1/rho_fluid] * tau
 + [(1 - z1/zm) * u_star_t^2] ) + [(log(z1/z0) - 2) * u_star_t/kappa )
 
 This formula is equivalent to equation 3.27 in Schwaemmle's thesis.  The
 three expressions in square brackets are computed in the constructor
 (m_veff_sqrt_slope, m_veff_sqrt_offset and m_veff_offset, respectively).  */

void flux3d_stationary::calc_u( double *u_x, double *u_y,
                               const double& tau_x, const double& tau_y, double gradh_x, double gradh_y )
{
    double abs_tau = sqrt( tau_x*tau_x + tau_y*tau_y );
    double u_star = sqrt(abs_tau/duneglobals::rho_fluid());
    double u_star_a_t = u_star_at(u_star);
    if(abs_tau <= m_tau_t) {
        *u_x= *u_y= 0.0;
        return;
    }
    
    //------efective wind velocity-----------
    
    double abs_u_wind = 2.5 *u_star *(m_log_z1z0 - (1.- u_star_a_t/u_star)*m_b_z1);
    
    if(abs_u_wind < m_DeltaU){
        *u_x= *u_y= 0.0;
        return;
    }
    
    double windx= (abs_u_wind/m_DeltaU)*tau_x/abs_tau;
    double windy= (abs_u_wind/m_DeltaU)*tau_y/abs_tau;
    double Dhx= 2.*m_alpha*gradh_x;
    double Dhy= 2.*m_alpha*gradh_y;
    double Gx, Gy, G, Fx, Fy, F, u, du, delta_du= 1.0;
    
    double dux= tau_x/abs_tau+Dhx;
    double duy= tau_y/abs_tau+Dhy;
    du= sqrt(dux*dux+duy*duy);
    du= sqrt(du);
    int j=0;
    while(fabs(delta_du) > 1e-3){
        delta_du= du;
        Fx= du*windx-Dhx;
        Fy= du*windy-Dhy;
        F= sqrt(Fx*Fx+Fy*Fy);
        if(du>0) u = (F-1.0)/du;
        else u=0;
        Gx= u*Dhx+windx;
        Gy= u*Dhy+windy;
        G= sqrt(Gx*Gx+Gy*Gy);
        if(u>0) du= 0.5*(sqrt(1.0+4.0*u*G)-1.0)/u;
        else break;
        
        delta_du-=du;
        delta_du= delta_du/du;
    }
    double ux, uy;
    ux= u*Fx/F;
    uy= u*Fy/F;
    //ux= windx - du;
    //uy= windy - du;
    if( (ux*tau_x+uy*tau_y) > 0){
        *u_x= ux*m_DeltaU;
        *u_y= uy*m_DeltaU;
    }else{
        *u_x= *u_y= 0;
    }
    
}

/*!  Computes the sand flux.  To give the initial condition, the column x=0 of
 flux has to be filled with the influx (which will be projected onto the
 direction of the air at those points).  */

void flux3d_stationary::calc(  TFktScal& influx, TFktVec& flux, const TFktScal& h, const TFktScal&
                             h_nonerod, const TFktVec& tau, TFktScal& gamma )
{
    double gradh_x, gradh_y, qin, qout;
    int x, y, prevx, nextx, prevy, nexty;
    
    // calculate surface gradients under periodic bc...
    for( x= 0; x< duneglobals::nx(); ++x ) {
        prevx = (x==0)?(duneglobals::nx()-1):(x-1);
        nextx = (x==duneglobals::nx()-1)?(0):(x+1);
        for( y= 0; y< duneglobals::ny(); ++y ) {
            gradh_x= 0.5*(h(nextx, y) - h(prevx, y))/(duneglobals::dx());
            
            prevy = (y==0)?(duneglobals::ny()-1):(y-1);
            nexty = (y==duneglobals::ny()-1)?(0):(y+1);
            
            gradh_y= (nexty==prevy? 0: (0.5*(h(x,nexty)-h(x,prevy))/(duneglobals::dx())));
            calc_u( &m_u(x, y)[0], &m_u(x, y)[1], tau(x, y)[0], tau(x, y)[1], gradh_x, gradh_y );
        }
    }
    
    // now compute the influx in the direction of the shear velocity and the
    // initial sand density.
	for( y= 0; y< duneglobals::ny(); ++y ){
		if( fabs(m_u(0, y)[0]) > 0 ) {
			m_inrho(0, y)= influx(0,y)/m_u(0,y)[0];
			if( !isfinite(m_rho(0, y)) )
				m_inrho(0, y)= 0.0;   // very small u_abs might cause the result infinity
		}
		else
			m_inrho(0, y)= 0.0;
    }
    
    // Saturated flux= satuated density * saturated velocity
	for( x= 0; x< duneglobals::nx(); ++x )
		for( y= 0; y< duneglobals::ny(); ++y ) {
            // for beach:
            double tau_t_fact = tau_t_factor(h(x, y));
			double abs_tau = sqrt( tau(x,y)[0]*tau(x,y)[0] + tau(x,y)[1]*tau(x,y)[1] );
			double rho_sat = (abs_tau > m_tau_t * tau_t_fact ? (abs_tau - m_tau_t * tau_t_fact) * m_2alpha_g : 0);
			m_flux_sat(x, y)[0] = rho_sat * m_u(x, y)[0];
			m_flux_sat(x, y)[1] = rho_sat * m_u(x, y)[1];
		}
    
    calc_rho( h, h_nonerod, tau, gamma);
    
    // Flux= density * saturated velocity
    for( x= 0; x< duneglobals::nx(); ++x )
        for( y= 0; y< duneglobals::ny(); ++y ) {
            flux(x, y)[0] = m_rho(x, y) * m_u(x, y)[0];
            flux(x, y)[1] = m_rho(x, y) * m_u(x, y)[1];
        }
    
    for( y= 0; y< duneglobals::ny(); ++y ){
        for( x= 0; x< duneglobals::nx(); ++x ){
            if(gamma(x, y)== -1.0) {
                
                double qx = flux(x,y)[0];
                double qy = flux(x,y)[1];
                double qxprevx = (x==0?influx(0,y):flux(x-1,y)[0]);
                //                prevx = (x==0)?((m_x_periodic)?(duneglobals::nx()-1):0):(x-1);
                nextx = (x==duneglobals::nx()-1)?((m_x_periodic)?0:x):(x+1);
                prevy = (y==0)?((m_y_periodic)?(duneglobals::ny()-1):0):(y-1);
                nexty = (y==duneglobals::ny()-1)?((m_y_periodic)?0:y):(y+1);
                
                double qxnextx = flux(nextx,y)[0];
                double qyprevy = flux(x,prevy)[1];
                double qynexty = flux(x,nexty)[1];
                
                qout= fabs(qx + qy);
                qin =   (qxprevx > 0 ? qxprevx : 0) +
                (qxnextx < 0 ? -qxnextx : 0) +
                (qyprevy > 0 ? qyprevy : 0) +
                (qynexty < 0 ? -qynexty : 0);
                gamma(x, y)= (qout - qin) / duneglobals::dx();
            }
        }
    }
    
}
/*!  Computes the saltating sand density.  The erosion/deposition rate gamma is
 not computed from the density rho but determined separately while computing
 rho.  */

void flux3d_stationary::calc_rho( const TFktScal& h, const TFktScal& h_nonerod,
                                 const TFktVec& tau, TFktScal& gamma )
{
    double A, B;
    
    double tau_abs, u_abs, div_u, splash;
    int x, y, miny, maxy, prevx, nextx, prevy, nexty;
    bool yascend, firsty;

    double tau_t_fact = 1; //beach

    for(y=0;y<duneglobals::ny();++y) {
        m_rhoaux(0,y) = m_rho(0,y);
    }
    for( x= 0; x< duneglobals::nx(); ++x )
    {
        
        maxy= -1;
        
        while( true)
        {
            y= maxy+1;
            if( y >= duneglobals::ny() )	break;
            
            for( ; y< duneglobals::ny() && m_u(x, y)[0]==0.0 && m_u(x, y)[1]==0.0; ++y )
            {  m_rho(x, y)= 0.0; gamma(x,y)= -1.0; m_rhoaux(x,y) = 0.0;}
            if( y >= duneglobals::ny() )	break;
            miny= y;
            for( ; y< duneglobals::ny() && m_u(x, y)[1]==0 && m_u(x, y)[0]!=0; ++y );
            if( m_u(x, y)[1] < 0 ) {
                yascend= false;
                while( y< duneglobals::ny() &&
                      (m_u(x, y)[1] < 0 || (m_u(x, y)[1]==0 && m_u(x, y)[0]!=0)) )
                    ++y;
                --y;
                maxy= y;
            }
            else {
                yascend= true;
                while( y< duneglobals::ny() &&
                      (m_u(x, y)[1] > 0 || (m_u(x, y)[1]==0 && m_u(x, y)[0]!=0)) )
                    ++y;
                maxy= y-1;
                y= miny;
            }
            for(/* firsty= true*/; yascend? y<=maxy : y>=miny;/* firsty= false,*/
                yascend? ++y: --y )
            {
                if (m_x_periodic) {
                    
                    prevx = (x==0)?(duneglobals::nx()-1):(x-1);
                    nextx = (x==duneglobals::nx()-1)?(0):(x+1);
                    
                    if( m_u(x+1, y)[0] == 0.0)
                        div_u= m_u(x, y)[0] - m_u(x-1, y)[0];
                    else if(m_u(x-1, y)[0] == 0.0)
                        div_u= m_u(x+1, y)[0] - m_u(x, y)[0];
                    else
                        div_u = 0.5 * (m_u(nextx,y)[0] - m_u(prevx,y)[0]);
                    
                } else {
                    if( x==duneglobals::nx()-1 || m_u(x+1, y)[0] == 0.0)
                        div_u= m_u(x, y)[0] - m_u(x-1, y)[0];
                    else if( x==0 || m_u(x-1, y)[0] == 0.0)
                        div_u= m_u(x+1, y)[0] - m_u(x, y)[0];
                    else
                        div_u= 0.5 * (m_u(x+1, y)[0] - m_u(x-1, y)[0]);
                }
                
                prevy = (y==0)?(duneglobals::ny()-1):(y-1);
                nexty = (y==duneglobals::ny()-1)?(0):(y+1);
                
                div_u += 0.5 * (m_u(x,nexty)[1] - m_u(x,prevy)[1]);
                
                A= m_u(x, y)[0] + div_u;
                
                B= m_u(x, y)[0] * (x==0 ? m_inrho(0, y) : m_rhoaux(x-1, y));
                if( yascend ) {
                    A += m_u(x, y)[1];
                    B += m_u(x, y)[1] * m_rhoaux(x, prevy);
                }
                else {
                    A += - m_u(x, y)[1];
                    B += -m_u(x, y)[1] * m_rhoaux(x,nexty);
                }
                A /= duneglobals::dx();
                B /= duneglobals::dx();
                //  Erosion can happen only if the sand layer is thick enough.
                splash= (h(x, y) - h_nonerod(x, y)) / m_d_grain;
                if( splash > 1.0 )
                    splash= 1.0;
                else if( splash < /*0.001*/100*m_d_grain )
                    splash = 0;
                
                //  Diffusion term:
                if(splash == 0 && fabs(m_u(x, y)[1]/m_u(x, y)[0])<1e-2){
                    
                    if (m_x_periodic) {
                        A -= -m_Cd_dx;
                        
                        prevx = (x==0?duneglobals::nx()-1:x-1);
                        double prevxx = (x==0?duneglobals::nx()-2:(x==1?duneglobals::nx()-1:x-2));
                        
                        B += 0*m_Cd_dx * (m_rhoaux(prevxx, y) - 2.0 * m_rhoaux(prevx, y));
                        
                        prevy = (y==0)?((m_y_periodic)?(duneglobals::ny()-1):0):(y-1);
                        nexty = (y==duneglobals::ny()-1)?((m_y_periodic)?0:y):(y+1);
                        
                        B += m_Cd_dx * (0.5 *(m_rhoaux(x, nexty) + m_rhoaux(x, prevy))-0*m_rhoaux(prevx,y));
                        //                        B += m_Cd_dx * (0.5 *(m_rhoaux(prevx, nexty) + m_rhoaux(prevx, prevy))-m_rhoaux(prevx,y));
                        
                    } else {
                        if( x>1 ) {
                            A -= -m_Cd_dx;
                            
                            B += 0*m_Cd_dx * (m_rho(x-2, y) - 2.0 * m_rho(x-1, y));
                            
                            prevy = (y==0)?((m_y_periodic)?(duneglobals::ny()-1):0):(y-1);
                            nexty = (y==duneglobals::ny()-1)?((m_y_periodic)?0:y):(y+1);
                            
                            B += m_Cd_dx * (0.5 *(m_rhoaux(x, nexty) + m_rhoaux(x, prevy))-0*m_rhoaux(x,y));
                        }
                        
                    }
                    
                }
                
                u_abs= sqrt( m_u(x, y)[0]*m_u(x, y)[0] + m_u(x, y)[1]*m_u(x, y)[1] );
                tau_abs= sqrt( tau(x, y)[0]*tau(x, y)[0] + tau(x, y)[1]*tau(x, y)[1] );
                // for beach:
                tau_t_fact = tau_t_factor(h(x, y));

                if(B < 0 || A <= 0 || !u_abs){
                    m_rho(x,y)= 0.0;
                    gamma(x,y)= -1.0;
                }else{
                    m_rho(x,y) = solve_for_rho(u_abs, tau_abs, splash, A, B, tau_t_fact);
                    
                    gamma(x,y)= A*m_rho(x,y)-B;
                }
                
            }
        }
        for(y=0;y<duneglobals::ny();++y) {
            m_rhoaux(x,y) = m_rho(x,y);
        }
    }
    
    
    
}

/*!  Determines rho from a quadratic equation which depends on the
 circumstances.  */

double flux3d_stationary::solve_for_rho(const double u,
                                        const double tau, const double splash, const double A, const double B, const double tau_t_fact)
{
    double rho, rho_s, rho_fs, A0_1, A0, A2, A1, rho_dep, rho_direct, rho_eros;
    double u_star = sqrt(tau/duneglobals::rho_fluid());
    double u_star_a_t = u_star_at(u_star);
    double tau_t = duneglobals::rho_fluid()*u_star_a_t*u_star_a_t;
    // Saturated rho
    rho_s= (tau > m_tau_t*tau_t_fact ? (tau - tau_t*tau_t_fact) * m_2alpha_g : 0);
    // Saturated rho for direct entrainment
    rho_fs= (tau > m_tau_ft*tau_t_fact ? (tau - m_tau_ft*tau_t_fact) * m_2alpha_g : 0);
    // ---Auxiliar const
    // for aeolian entrainment
    A0_1= m_beta_2alpha_g_tau_ft * splash;
    A0= A0_1 * rho_fs;
    // for aeolian deposition
    A2= m_gamma_2alpha_g2 / (u * tau_t);
    A1= A2 * rho_s;
    
    // deposition:
    // calculation of rho with the assumtion that rho > rho_sat
    rho_dep=(B==0 ? 0 : 0.5 * (sqrt((A - A1)*(A - A1) + 4* A2 * B) - (A - A1)) / A2);
    // verifying the assumtion
    if( rho_dep < rho_s ){
        // erosion:
        if( splash > 0 ){
            // calculation of rho without direct entrainment ( rho > rho_fs)
            A2*= splash;
            A1= A2 * rho_s;
            rho_eros=(B==0 ? 0 : 0.5 * (sqrt((A - A1)*(A - A1) + 4* A2 * B) - (A - A1)) / A2);
            
            // verifying the assumtion
            if( rho_eros < rho_fs ){
                // calculation of rho with direct entrainment ( rho < rho_fs)
                A1-= A0_1;
                rho_direct= 0.5 * (sqrt((A - A1)*(A - A1) + 4* A2 * (B + A0)) - (A - A1)) / A2;
                
                // verifying the assumtion
                if( rho_eros > rho_fs )    cout << "no solution for rho!!!!!" << endl;
                rho= rho_direct;
            }else	rho= rho_eros;
        }else	rho= B / A;
    }else	rho= rho_dep;
    
    if(rho > 0) return rho;
    else return 0;
}


/*!  Saves the arrays m_u and m_rho.  */

void flux3d_stationary::save_arrays()
{
    save_2d_vecarray( "u", m_u );
    save_2d_vecarray( "flux_s", m_flux_sat );
    save_2d_scalarray( "rho", m_rho );
}

