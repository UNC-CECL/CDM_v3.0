/******************************************************************************
 $Id: vegetation.h,v 1.5 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#ifndef VEGETATION_H
#define VEGETATION_H

#include "globals.h"
#include "func.h"
//#include "PTG_FileNames.h"
#include "vec.h"

class arrayinit;

class vegetation : public dunedata
{

public:
    // construction
    vegetation(const dunepar& par);
    virtual ~vegetation() {}
        
    //  Functions reimplemented from coast:
    virtual void init(const dunepar& par);
    //virtual void season( double time, double timestep, bool &calc_veg_seasonality);

    virtual int evol(TFktScal& rho_veget, const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dhdt, TFktScal& overwash);
    void getcover(TFktScal& rho_veget, const double time, const double timestep);

    virtual void shiftback(const int plusminus);

    virtual void save_arrays();

private:
    
    virtual int evolspec(const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dhdt, TFktScal& overwash, int species);

    /* auxiliary function for vegetation */
    TFktVec m_veget;    // vegetation cover fraction for two species

    /*! Vegetation parameters (see default.par)*/
    int m_xmin, m_xmin0; // vegetation minimum distance to shoreline
    int m_veg_type;
    double m_Lveg; // vegetation limit = m_xmin * dx
    double m_zmin; // veget minimum height above sealevel
    double m_veget_init0, m_veget_init1, m_veget_init2_1, m_veget_init2_2;
    double m_Tveg, m_rho_max, m_rho_min, m_Hveg, m_sens, m_veg_season, m_wind_factor, m_Vlateral_factor, m_angle_ref;
    double m_veg1_r, m_veg1_C, m_veg1_seed_deterministic, m_veg1_seed_probabilistic, m_veg1_seed_zmax, m_veg1_seed_wrack_pct, m_veg1_K, m_veg1_r_sand, m_veg1_r_sand2, m_overwash_sens1, m_veg1_season_start, m_veg1_season_end, m_veg1_season_senescence_prop, m_veg1_Hveg, m_veg1_alpha21;
    double m_veg2_r, m_veg2_C, m_veg2_seed_deterministic, m_veg2_seed_probabilistic, m_veg2_seed_zmax, m_veg2_seed_wrack_pct, m_veg2_K, m_veg2_r_sand, m_veg2_r_sand2, m_overwash_sens2, m_veg2_season_start, m_veg2_season_end, m_veg2_season_senescence_prop, m_veg2_Hveg, m_veg2_alpha12;
    double m_wrack_max, m_wrack_ht, m_wrack_seed_probabilistic;
    double m_Sdt;
    bool m_spec1, m_spec2, m_lateral;
    bool m_survive;
    bool m_wrack;
    bool m_veg2, m_veg1_annual, m_veg2_annual, m_veg12_interactions;
    int seasonindex;
    string m_veget_init_surf;
    string m_veg1_seed_owregion;
    string m_veg2_seed_owregion;

};


#endif //  VEGETATION_H

