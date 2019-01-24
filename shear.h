#ifndef __SHEAR_H__
#define __SHEAR_H__


#include "vec.h"
#include "func.h"
#include "globals.h"

class sepbubble;


////////////////////////////////////////////////////////////
// shear
//

class shear : public dunedata
{
public:
    shear(const dunepar& P);
    virtual ~shear();
    
    virtual double Calc(const TFktScal& h, TFktVec& tau);
    virtual double Calc(const TFktScal& h, TFktVec& tau, const TFktScal& rho_veget);
    
    //virtual bool Save() { return false; }
    
    virtual void save_arrays();
    
    void set_ustar( double u_star );
    
protected:
    virtual double Calc(const TFktScal& h, TFktVec& tau, const TFktScal *rho_veget);
    virtual double CalcPertTau(TFktScal& h, TFktVec& tau) = 0;
    //virtual void CalcPertTau(const TFktScal& h, const TFktVec *veget, TFktVec& tau) = 0;
    /*!  The shear velocity on a flat plain corresponding to the current wind
     speed.  */
    double m_u_star;
    /*!  The unperturbed shear stress corresponding to m_u_star.  */
    double m_dTau0;
    
private:
    
    TFktScal m_stall;
    TFktScal m_rho_veget;
    TFktScal m_rho_veget_smooth;
    
    TFktScal m_hSepBub;
    TFktScal m_hSepBub_aux;
    sepbubble* m_pSepBub;
    
    TFktVec m_TauP;
    
    bool m_bTauYZero;
    double m_bTauY;
    bool m_bSepbub;
    /*!  This seems to give the decrease in shear stress after the brink of the
     slip face.  That would determine the length of the region after the brink
     where the sand blown over the brink is deposited.  This variable is given
     by the parameter sepbub.tau.  */
    double m_tau_sepbub;
    
    bool m_calc_veget, m_addsealevel;
    double m_sealevel;
    double m_rhofactor, m_veget_m, m_veget_beta_sigma, h_limit;
    
    //  bool m_calc_veget;
    //  double m_rhofactor, m_veget_m, m_veget_beta_sigma, h_limit;
};

#endif
