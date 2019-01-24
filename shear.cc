#include <time.h>

#include "globals.h"
//#include "save.h"
#include "shear.h"
#include "shear_hlr.h"
#include "sepbubble.h"

//*****************************************************************************
//  class shear
//

shear::shear(const dunepar & par) :
dunedata(par),
m_rho_veget_smooth(duneglobals::nx(), duneglobals::ny(), duneglobals::dx()),
m_hSepBub(duneglobals::nx(), duneglobals::ny(), duneglobals::dx()),
m_hSepBub_aux(duneglobals::nx(), duneglobals::ny(), duneglobals::dx()),
m_TauP(duneglobals::nx(), duneglobals::ny(), duneglobals::dx())
{
    m_pSepBub= sepbubble::create(par);
    
    m_tau_sepbub = par.getdefault("sepbub.tau", 0.05);
    
    m_bTauY = par.getdefault( "hlr.tau_y", 1.);
    m_stall.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    m_rho_veget.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    
    m_addsealevel = par.getdefault("shore.correct.profile", true);
    m_sealevel = par.getdefault("shore.sealevel", 0.0);
    m_sealevel = 0;
//    m_sealevel = par.getdefault("shore.HMWL", 0.0);
    
    m_calc_veget = par.getdefault("veget.calc", false);
    m_veget_m= par.getdefault("veget.m", 0.16);
    m_veget_beta_sigma= par.getdefault("veget.beta", 150.0)/par.getdefault("vegetation.sigma", 1.5);
    h_limit= 1e-3;
    
    double m_rho_max = par.getdefault("veget.rho.max", 1.0);
    double max_corr = (1 - m_veget_m * m_rho_max)*(1 + m_veget_m*m_veget_beta_sigma * m_rho_max);

//    m_calc_veget = P.getdefault("calc.vegetation", false);
//    m_rhofactor= P.getdefault("vegetation.concent", 1.0)*P.getdefault("vegetation.sigma", 1.0)*P.getdefault("vegetation.sigma", 1.0);//(duneglobals::dx()*duneglobals::dx());
//    m_veget_m= P.getdefault("vegetation.m", 1.0);
//    m_veget_beta_sigma= P.getdefault("vegetation.beta", 100.0)/P.getdefault("vegetation.sigma", 1.0);
    h_limit= 1e-3;
    
    cout << "!!! RHO_MAX = " << m_rho_max << endl;
    cout << "!!! CORR_MAX = " << max_corr << endl;

}


shear::~shear()
{
}

/*!  Calls calc(... TFktVec* veget = NULL).  */
double shear::Calc(const TFktScal& h, TFktVec& tau)
{
    return Calc( h, tau, NULL);
}

/*!  Calls calc(... TFktVec* veget).  */
double shear::Calc(const TFktScal& h, TFktVec& tau, const TFktScal& rho_veget)
{
    return Calc( h, tau, &rho_veget);
}

void shear::set_ustar( double u_star )
{
    m_u_star= u_star;
    m_dTau0= u_star * u_star * duneglobals::rho_fluid();
}

double shear::Calc(const TFktScal& h, TFktVec& tau, const TFktScal *rho_veget)
{
    clock_t clocktime;
    
    if( duneglobals::timing() )
        clocktime= clock();

    if(m_calc_veget){
        //m_rho_veget_smooth.Smooth(*rho_veget);
        m_rho_veget_smooth = *rho_veget;
        //        m_stall = m_rho_veget_smooth;
        for (int y=0; y<duneglobals::ny(); y++)
      		for (int x=0; x<duneglobals::nx(); x++){
                m_stall(x,y) -= h_limit;
            }
    }else	m_stall.SetAll(-1);
    m_pSepBub->Calc(m_hSepBub, m_stall, h);

    m_hSepBub_aux = m_hSepBub;
    if (m_addsealevel) {
        for (int y=0; y<duneglobals::ny(); y++) {
            for (int x=0; x<duneglobals::nx(); x++) {
                if (m_hSepBub(x,y) < m_sealevel)
                {
                    m_hSepBub_aux(x,y) = m_sealevel;
                }
            }
        }
    }

    // calc shear stress pertubation
    double L = CalcPertTau(m_hSepBub_aux, m_TauP);
    
    // PARTELI TEST
    //  if(duneglobals::check_error())
    //     cout << "shear.cc; double L = CalcPertTau(m_hSepBub, m_TauP); done!" << endl;
    
    // calc shear stress
    // The shear stress will be consider 'null' only INSIDE the separation region, at the border
    // the actual values are conserved.
    
    // tau
    for (int y=0; y<duneglobals::ny(); y++) {
        for (int x=0; x<duneglobals::nx(); x++) {
            if(m_calc_veget){
                double factor;

                factor= (1 - m_veget_m * m_rho_veget_smooth(x,y))*(1 + m_veget_m*m_veget_beta_sigma * m_rho_veget_smooth(x,y));
                factor= (factor > 0? 1.0/factor : 1e20);
                
                tau(x,y)[0]= tau(x,y)[1]= factor;
            }else{
                tau(x,y)[0]= tau(x,y)[1]= 1.0;
            }
        }
    }
    
    const double slope= tan(duneglobals::repose_dyn()*M_PI/180.0)*duneglobals::dx();
    const double delta=1.0/(m_tau_sepbub*slope);
    for (int y=0; y<duneglobals::ny(); y++)
        for (int x=0; x<duneglobals::nx(); x++) {
            // Reducing the shear stress below separation surface (to mimic the
            //turbulence effects)
            double h_delta= 1 - delta*(m_hSepBub(x,y)-h(x,y));
            if(h_delta < 0) h_delta= 0.0;
            else if(h_delta > 1.0) h_delta= 1.0;
            
            // TEST
            // m_TauP(x,y)[0] *= 0.5;
            // m_TauP(x,y)[1] *= 0.5;

            // Acotation
            if(m_TauP(x,y)[0] > 10) m_TauP(x,y)[0] = 10;
            else if(m_TauP(x,y)[0] < -10) m_TauP(x,y)[0] = -10;
            
            if(m_TauP(x,y)[1] > 10) m_TauP(x,y)[1] = 10;
            else if(m_TauP(x,y)[1] < -10) m_TauP(x,y)[1] = -10;

            // Shear stress calculation
            tau(x,y)[0] *= m_dTau0 * fabs( 1 + 0.5*m_TauP(x,y)[0])*( 1 + 0.5*m_TauP(x,y)[0]) * h_delta;
            tau(x,y)[1] *= m_bTauY * m_dTau0 * fabs( 1 + 0.5*m_TauP(x,y)[0])*m_TauP(x,y)[1] * h_delta;
            
            if(tau(x,y)[0] < 0) tau(x,y)[0] = 0.;
            
            //redefinition of tau_p in order to print the wind speed-up
//            m_TauP(x,y)[0]= 1 + 0.5*m_TauP(x,y)[0];
//            m_TauP(x,y)[0]= (m_TauP(x,y)[0] > 0 ? m_TauP(x,y)[0] : 0);
//            m_TauP(x,y)[1]= (m_TauP(x,y)[1] < 0 ? -1 : 1)*sqrt(m_TauP(x,y)[0]*fabs(m_TauP(x,y)[1]));
        }
    
    if( duneglobals::messages() )
        if( duneglobals::timing() ) {
            clocktime= clock() - clocktime;
            cout << ", tau: " << clocktime*1000/CLOCKS_PER_SEC << "ms) ";
        }
        else
            cout << ", tau) ";
    
    return L;
}
/*!  Saves the arrays m_hSepBub, m_stall and m_TauP.  */

void shear::save_arrays()
{
    save_2d_scalarray("h_sep", m_hSepBub);
    save_2d_scalarray("stall", m_stall);
    save_2d_vecarray("shear_pert", m_TauP );
    //if(m_calc_veget)	save_2d_scalarray("rho_veget", m_rho_veget);
}

