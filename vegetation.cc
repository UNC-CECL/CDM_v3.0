/******************************************************************************
 $Id: vegetation.cc,v 1.9 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#include <math.h>
#include <random>
#include <chrono>

#include "globals.h"
#include "func.h"
#include "initsurf.h"
#include "vegetation.h"

//*****************************************************************************
//  class vegetation

vegetation::vegetation(const dunepar& p) : dunedata(p)
{
    m_veget.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );
    
    /*Vegetation parameters*/
    m_veg_type = p.getdefault("veget.type", 1);

    m_Lveg = p.getdefault("veget.xmin", 0.0);
    m_xmin0 = m_Lveg/duneglobals::dx();
    m_xmin = m_xmin0;
    m_zmin = p.getdefault("veget.zmin", 0.0);

    m_Tveg = p.getdefault("veget.Tveg", 1.0) * duneglobals::secday();       // typical generic plant growth time (sec)
    cout << "grass constructor: Tveg = " << m_Tveg << endl;
    
    // initial value
    m_veget_init0 = p.getdefault("veget.0.init", 1e-2); // fix to 10%
    m_veget_init1 = p.getdefault("veget.1.init", 1e-5); // fix to 10
    m_veget_init2_1 = p.getdefault("veget1.2.init", 0.0); // fix to 10%
    m_veget_init2_2 = p.getdefault("veget2.2.init", 0.0); // fix to 10%
    
    // erosion/acc
    m_sens = p.getdefault("veget.erosion.sensitivity", 1.0);       // species sensitivity for erosion/acc. rate
    m_overwash_sens1 = p.getdefault("veget1.overwash.sensitivity", 1.0);       // proportional death with overwash
    m_overwash_sens2 = p.getdefault("veget2.overwash.sensitivity", 1.0);       // proportional death with overwash

    m_Hveg = p.getdefault("veget.Hveg", 0.3);       // typical dune-forming plant height  (~ 0.3 m)
    
    // Veg Logistic Growth Parameterization
    bool m_calc_season = p.getdefault("calc.season", false);
	if(m_calc_season) {
		m_veg1_season_start =  p.getdefault("veget1.seasonstart", 0.0); // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg1_season_end =  p.getdefault("veget1.seasonend", 1.0); // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg1_season_senescence_prop = p.getdefault("veget1.senescence", 0.0); // proportional reduction in percent cover due to leaf senescence
		m_veg2_season_start =  p.getdefault("veget2.seasonstart", 0.0); // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg2_season_end =  p.getdefault("veget2.seasonend", 1.0); // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg2_season_senescence_prop = p.getdefault("veget2.senescence", 0.0); // proportional reduction in percent cover due to leaf senescence

	} else {
		m_veg1_season_start =  0.0; // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg1_season_end =  1.0; // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg1_season_senescence_prop = 0.0; // proportional reduction in percent cover due to leaf senescence
		m_veg2_season_start =  0.0; // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg2_season_end =  1.0; // beginning of the growing season (as a proportion of entire year). 0 = 1/1. 1=12/31
		m_veg2_season_senescence_prop = 0.0; // proportional reduction in percent cover due to leaf senescence
	}

    m_veg1_annual = p.getdefault("veget1.annual", false);
    m_veg1_r = p.getdefault("veget1.r", 1.)/(m_veg1_season_end-m_veg1_season_start);       // typical dune-forming plant intrinsic growth rate  (~ 1)
    m_veg1_C = p.getdefault("veget1.C", 1.)/(m_veg1_season_end-m_veg1_season_start);       // typical dune-forming plant lateral growth rate  (~ 1 m/yr)
    m_veg1_K = p.getdefault("veget1.K", 1.);       // carrying capacity (percent cover where dRho/dt = 0)
    m_veg1_seed_deterministic = p.getdefault("veget1.seed.deterministic", 0.);
    m_veg1_seed_zmax = p.getdefault("veget1.seed.zmax", 10000.);
    m_veg1_seed_probabilistic = p.getdefault("veget1.seed.probabilistic", 0.1); // Poisson rate for random occurrence of % cover (mean % cover from seed/yr)
    m_veg1_seed_owregion = p.getrequired<string>("veget1.seed.owregion"); // high_water_line or all
    m_veg1_r_sand = p.getdefault("veget1.sand", 0.0); // modification of r from sand burial (dhdt)
	m_veg1_r_sand2 = p.getdefault("veget1.sand2", 0.0); // modification of r from (sand burial)^2 (dhdt)
	m_veg1_seed_wrack_pct = p.getdefault("veget1.seedfactor", 0.0); // relative proportion of propagules that are derived from wrack (e.g., rhizomes) vs. other sources (e.g., random seeding)
	m_veg1_Hveg = p.getdefault("veget1.Hveg", 0.3);       // typical dune-forming plant height  (~ 0.3 m)

	m_veg2 = p.getdefault("veget2", false);
	m_wrack = p.getdefault("calc.wrack", false);
	if(m_veg2 && m_wrack==false) {
		// veg characteristics
		m_veg2_annual = p.getdefault("veget2.annual", false);
	    m_veg2_r = p.getdefault("veget2.r", 0.);       // typical dune-forming plant intrinsic growth rate  (~ 1)
	    m_veg2_C = p.getdefault("veget2.C", 0.);       // typical dune-forming plant lateral growth rate  (~ 1 m/yr)
	    m_veg2_K = p.getdefault("veget2.K", 1.);       // carrying capacity (percent cover where dRho/dt = 0)
	    m_veg2_seed_deterministic = p.getdefault("veget2.seed.deterministic", 0.);
	    m_veg2_seed_zmax = p.getdefault("veget2.seed.zmax", 10000.);
		m_veg2_seed_probabilistic = p.getdefault("veget2.seed.probabilistic", 0.1); // Poisson rate for random occurrence of % cover (mean % cover from seed/yr)
		m_veg2_seed_owregion = p.getrequired<string>("veget2.seed.owregion"); // high_water_line or all
		m_veg2_r_sand = p.getdefault("veget2.sand", 0.0); // modification of r from sand burial (dhdt)
		m_veg2_r_sand2 = p.getdefault("veget2.sand2", 0.0); // modification of r from (sand burial)^2 (dhdt)
		m_veg2_seed_wrack_pct = p.getdefault("veget2.seedfactor", 0.0); // relative proportion of propagules that are derived from wrack (e.g., rhizomes) vs. other sources (e.g., random seeding)
		m_veg2_Hveg = p.getdefault("veget2.Hveg", 0.3);       // typical dune-forming plant height  (~ 0.3 m)
	}
	if(m_veg2 && m_wrack) {
		// wrack characteristics
	    m_wrack_max = p.getdefault("wrack.max", 0.);       // carrying capacity (percent cover where dRho/dt = 0)
	    m_wrack_ht = p.getdefault("wrack.ht", 0.0); // Poisson rate for random occurrence of % cover (mean % cover from seed/yr)
	    m_wrack_seed_probabilistic  = p.getdefault("wrack.seed.probabilistic", 0.1); // Bernoulli rate for random occurrence of % cover (mean % cover from seed/yr)
	}


	// Species Interactions
	m_veg12_interactions = p.getdefault("veget12.interactions", false);       // Do species 1 and 2 interact?
	if(m_veg12_interactions) {
		m_veg1_alpha21 = p.getdefault("veget1.alpha21", 0.0);       // effect of species 2 on species 1
		m_veg2_alpha12 = p.getdefault("veget2.alpha12", 0.0);       // effect of species 1 on species 2
	} else {
		m_veg1_alpha21 = 0.0;
		m_veg2_alpha12 = 0.0;
	}



    m_angle_ref = p.getdefault("veget.max.slope", 15);  // highest slope for lateral vegetation growth (degrees)

    // time conversion
    m_wind_factor =  duneglobals::timefrac();
    
    // different species
    m_spec1 = p.getdefault("veget.spec.1", 1);
    m_spec2 = p.getdefault("veget.spec.2", 0);
    
    // lateral growth
    m_lateral = p.getdefault("veget.lateralgrowth", true);

    if (m_veg_type == 0)
    {
        m_Vlateral_factor = p.getdefault("veget.Llateral", 10);
    } else {
        m_Vlateral_factor = p.getdefault("veget.Vlateral.factor", 100);
    }
    
    // extra
    m_rho_max = p.getdefault("veget.rho.max", 1.0);
    m_rho_min = p.getdefault("veget.rho.min", 0.0);
    
    m_survive = p.getdefault("veget.survive", false);
    
    // seasonality
    // Frequency
    const int m_freq = 1; // freq. per year
    m_Sdt = duneglobals::secyear()*duneglobals::timefrac()/(double)m_freq; // adjusted freq. per year (accounting for m_freq)

    seasonindex = 0;

    // error check for seeding function parameterization
    if(m_veg1_seed_probabilistic > 0 && !(m_veg1_seed_owregion == "all" || m_veg1_seed_owregion == "hwl" || m_veg1_seed_owregion == "none")) {
        cout << "gradalongshore::create: FATAL ERROR: Unknown value `" << m_veg1_seed_owregion << "\' for parameter `veget1.seed.owregion\' !" << endl;
		cout << "  Valid values are: all, hwl, none" << endl;
		exit(1);

    }
}

/*! Initialize vegetation */
void vegetation::init(const dunepar& par)
{
    arrayinit *init_veget;
    
    if( par.exists("veget.Init-Surf") )
        init_veget= arrayinit::create(par, "veget.");
    else
        init_veget= new CInitSurfPlain(0.0);
    
    init_veget->init_2d_vec( m_veget );
    
    delete init_veget;
    
    // INIT VEGET
    if (m_veg_type == 0 && 0)
    {
        for( int x= 0; x< duneglobals::nx(); ++x ){
        // int x = 0.5 * duneglobals::nx();
            for( int y= 0; y< duneglobals::ny(); ++y ){
                m_veget(x,y)[0]= 0.0;
            }
        }
    }

    string strType = par.getrequired<string>("veget.Init-Surf");
    if (m_veg_type == 2 && strType == "plain")
    {
        for( int x= 0; x< duneglobals::nx(); ++x ){
        // int x = 0.5 * duneglobals::nx();
            for( int y= 0; y< duneglobals::ny(); ++y ){
                m_veget(x,y)[0]= m_veget_init2_1;
                m_veget(x,y)[1]= m_veget_init2_2;
            }
        }
    }
}

// kills off annual vegetation seasonally/annually
/*void vegetation::season( double time, double timestep, bool &calc_season)
{
    if( calc_season > 0 ){
    	if(m_veg1_annual) {
    		for( int x= 0; x< duneglobals::nx(); ++x ){
				for( int y= 0; y< duneglobals::ny(); ++y ){
					m_veget(x,y)[0]= 0.0;
				}
			}
    	}
    	if(m_veg2_annual) {
    		for( int x= 0; x< duneglobals::nx(); ++x ){
				for( int y= 0; y< duneglobals::ny(); ++y ){
					m_veget(x,y)[1]= 0.0;
				}
			}
    	}
    	calc_season = 0;
    } else {
        double tstep = time / timestep;
        double Sstep = m_Sdt / timestep;
        int tnextseason = (int) tstep % (int) Sstep; //750; //1000; //500;
        calc_season = (tnextseason == 0 && tstep >= 0 ? 1 : 0);

        if (calc_season > 0) {
            seasonindex++;
        }
        cout << "!! Season = " << calc_season << ' ' << tstep << ' ' << seasonindex << ' ' << tnextseason << endl;
	}
}*/

/*! return cover fraction*/
void vegetation::getcover(TFktScal& rho_veget, const double time, const double timestep)
{
	// --- TIME ----
	double realtimestep = timestep/duneglobals::secyear()/duneglobals::timefrac(); // years per timestep (1000 sec/timestep)*(timefrac timesteps / timestep)*(secyear year / sec) ]
	double time_yr = time/timestep*realtimestep; // time/timestep = number of timesteps; time/timestep*realtimestep = number of years elapsed;
	double m_season = fmod(time_yr, 1.0); // month of year, represented as a decimal (i.e., 0.00 = 1/1, 0.998 = 12/31)

    // NORMALIZATION
    for( int x= 0; x< duneglobals::nx(); ++x )
        for( int y= 0; y< duneglobals::ny(); ++y ){
            // added density:
            rho_veget(x,y) = m_veget(x,y)[0]*((m_season >= m_veg1_season_start && m_season <= m_veg1_season_end) ? 1 : (1 - m_veg1_season_senescence_prop)) + m_veget(x,y)[1];
            if (rho_veget(x,y) > m_rho_max) {
                rho_veget(x,y) = m_rho_max;
            }
        }
    
}

/*!  Computes the evolution of the vegetation.  */
int vegetation::evol(TFktScal& rho_veget, const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dh_dt, TFktScal& overwash)
{
    int veget_X0 = 0;
    // SPEC 1
    veget_X0 = evolspec(time, timestep, shoreline, h, dh_dt, overwash, 0);
    // SPEC 2
    if(m_veg2 == true) {
    	int veget_X1 = 0;
        veget_X1 = evolspec(time, timestep, shoreline, h, dh_dt, overwash, 1);
    }
    
    getcover(rho_veget, time, timestep);
    
    // set overwash to zero
    overwash.SetAll(0.0);
    return veget_X0;
}

/*!  Computes the evolution of each species.  */
int vegetation::evolspec(const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dh_dt, TFktScal& overwash, int species)
{
	// Stationary shoreline position
	/*if(m_Lveg_change = stationary) {
		m_xmin = m_xmin0 + m_shorelinechange;
	} else {
		m_xmin = m_xmin0;
	}*/


    // VEGET GRAD (for lateral propagation)
    TFktVec grad_veget;
    grad_veget.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx(),vec2(0.0,0.0));
    
    TFktScal veget_aux;
    veget_aux.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    // calculates an effective gradient based on interactions
    if(m_veg_type == 2 && m_veg12_interactions) {
    	if(species == 0) {
        	for( int y= 0; y< duneglobals::ny(); ++y ){
        	        for( int x= 0; x< duneglobals::nx(); ++x ){
        	            veget_aux(x,y) = m_veget(x,y)[0]; // + m_veg1_alpha21 * m_veget(x,y)[1];
        	        }
        	    }
    	}
    	if(species == 1) {
        	for( int y= 0; y< duneglobals::ny(); ++y ){
        	        for( int x= 0; x< duneglobals::nx(); ++x ){
        	            veget_aux(x,y) = m_veget(x,y)[1]; // + m_veg2_alpha12 * m_veget(x,y)[0];
        	        }
        	    }
    	}
    } else {
    	for( int y= 0; y< duneglobals::ny(); ++y ){
    	        for( int x= 0; x< duneglobals::nx(); ++x ){
    	            veget_aux(x,y) = m_veget(x,y)[species];
    	        }
    	    }
    }
    grad_veget.GradMin(veget_aux);
    
    // VEGET LAPLACIAN
    TFktVec laplacian_veget;
    laplacian_veget.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx(),vec2(0.0));
	laplacian_veget.DiscreteLaplacian(veget_aux);



    // CALCULATION OF VEG LIMIT RELAXATION

    // repose angle for normalization
    //    double m_angle_ref= 15.;
    double m_tan_ref= tan(m_angle_ref * M_PI / 180.0);


    // calculate gradient
    TFktVec grad_h;
    grad_h.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h.GradMid(h);

    // General variables
    double growthrate, proprate, reprod_rate, mortal_rate;
    // int vegetpoints = 0;
    
    //double m_veg1_D = (m_veg1_C * m_veg1_C / (4*m_veg1_r)); // vegetation diffusion rate (m2/yr) for 1m grid cells
    //double m_veg2_D = (m_veg2_C * m_veg2_C / (4*m_veg2_r)); // vegetation diffusion rate (m2/yr) for 1m grid cells

    // double m_shore_HMWL = duneglobals::HMWL();
    
    // --- TIME ----
	double realtimestep = timestep/duneglobals::secyear()/duneglobals::timefrac(); // years (1000 timesteps)*(timefrac sec / timestep)*(secyear sec / year) ]
	double time_yr = time/timestep*realtimestep;
	double m_season = fmod(time_yr, 1.0); // month of year, represented as a decimal (i.e., 0.00 = 1/1, 0.998 = 12/31)

        
    for(int y = 0; y< duneglobals::ny(); ++y ){
        for(int x = 0; x< duneglobals::nx()-1; ++x ){
            
            // AUXILIAR
            double dhdt = dh_dt(x,y); // * m_wind_factor; // average erosion rate
            double dhdtimestep = dhdt*timestep;
            double erosion = (dhdt < 0. ? 1 : 0);
            
            // limiting factors for vegetation growth
            double shorefactor = (x < shoreline + m_xmin ? 0 : 1);  // 0: maximum effect; 1: no effect
            double watertable = (h(x,y) > m_zmin? 1 : 0);

            // lateral propagation only if slope is lower than m_tan_ref
			double abs_grad_h = sqrt(grad_h(x,y)[0]*grad_h(x,y)[0]+grad_h(x,y)[1]*grad_h(x,y)[1]);
			double dhdxfactor = (1 - abs_grad_h/m_tan_ref);

            if (m_veg_type == 0 || m_veg_type == 1)
			{
				// absolute value gradient
				double abs_grad_veget = sqrt(grad_veget(x,y)[0]*grad_veget(x,y)[0]+grad_veget(x,y)[1]*grad_veget(x,y)[1]);

				// GROWTH RATE FOR GENERIC VEGETATION (no feedback with accretion rate & no lateral propagation)
				double V_gen = m_veget(x,y)[species] * (1 - m_veget(x,y)[species]) / (m_Tveg * m_wind_factor);

				// GROWTH RATE FOR AMMOPHILIA-LIKE VEGETATION (feedback with accretion rate & lateral propagation)
				// using continuous logistic growth curve (deltaRho = 1/(1+((1-rho0)/rho0)*exp(-1*r*deltaT))-rho0, assuming rhoMax=1
				double V_ammoph = (1. - erosion) * m_veget(x,y)[species] * (1 - m_veget(x,y)[species]) * dhdt / m_Hveg;

				// LATERAL PROPAGATION
				double proprate_ammoph =  (1. - erosion) * dhdt * (dhdxfactor > 0 ? 1 : 0) * m_Vlateral_factor * abs_grad_veget * (m_veget(x,y)[species] < 1 && dhdt > 0 ? 1 : 0);

				// TOTAL GROWTH RATES
				proprate    = watertable * shorefactor * (m_veg_type == 0 ? 0 : proprate_ammoph);
				reprod_rate = watertable * shorefactor * (m_veg_type == 0 ? V_gen : V_ammoph);
				mortal_rate = (m_veg_type == 0 ? 1 : erosion) * m_sens * m_veget(x,y)[species] * fabs(dhdt);

				growthrate = proprate + reprod_rate - mortal_rate;

				// evolution of cover fraction
				m_veget(x,y)[species] += timestep * growthrate;
			}
            if (m_veg_type == 2)
            {
            	if (species == 0)
            	{
            		// vertical and lateral growth
            		double Vert_ammoph;
            		double Seed_ammoph;

            		double m_veg_r_eff1 = m_veg1_r + m_veg1_r_sand * dhdt / realtimestep + m_veg1_r_sand2  * pow(dhdt / realtimestep,2) ;
            		double m_veg1_D = (m_veg1_C * m_veg1_C / (4*m_veg_r_eff1)); // vegetation diffusion rate (m2/yr) for 1m grid cells
            		//double m_veget_eff1 = m_veg1_D * realtimestep * laplacian_veget(x,y)[species] + m_veget(x,y)[species];

					if(m_veg12_interactions) {
	            		Vert_ammoph = m_veget(x,y)[species] + m_veg_r_eff1  * m_veget(x,y)[species] * (m_veg1_K - m_veg1_alpha21 * m_veget(x,y)[1] - m_veget(x,y)[species]) / m_veg1_K * realtimestep;
					} else {
						Vert_ammoph =  (m_veg1_K / (1 + ( (m_veg1_K - m_veget(x,y)[species]) / m_veget(x,y)[species]) * exp(-1 * m_veg_r_eff1 * realtimestep)));
					}
					double Lat_ammoph = m_veg1_D * realtimestep * laplacian_veget(x,y)[0]; // * m_veget(x,y)[species];

					// colonization
					if (m_veg1_seed_probabilistic > 0.0) {
						// random seedling/rhizome fragment growth
						std::default_random_engine generator;
						generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

						double m_veg_seed_mod;

						// specifies region of overwash-facilitated seeding. Either at high water line or throughout overwash area
						bool m_ow_seed = 0;
						if(m_veg1_seed_owregion == "all") {
							 m_ow_seed = (overwash(x,y)>0);
						} else if(m_veg1_seed_owregion == "hwl") {
							bool m_ow_seed = (overwash(x,y)>0 && ((overwash(x+1,y)==0) || x==(duneglobals::nx()-1)));
						} else if (m_veg1_seed_owregion == "none") {
							// m_ow_seed=0
						}

						m_veg_seed_mod = m_veg1_seed_probabilistic * (1.-m_veg1_seed_wrack_pct) * realtimestep; // scaled by number of timesteps per year
						if(m_ow_seed) {
							m_veg_seed_mod += m_veg1_seed_probabilistic * m_veg1_seed_wrack_pct / 183.; // scaled by number of storm events per year
						}



						std::bernoulli_distribution distribution(m_veg_seed_mod);

						int number = distribution(generator);
						Seed_ammoph = (double)number/100;

					} else {
						Seed_ammoph = 0.0;
					}

					// evolution of cover fraction
					if(m_season >= m_veg1_season_start && m_season <= m_veg1_season_end) {
						// m_veget(x,y)[species] = (m_veget(x,y)[species] + Vert_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) * ((dhdt / realtimestep) < -1.0 ? 0 : 1) + Seed_ammoph;
						//m_veget(x,y)[species] = (m_veget(x,y)[species] + Vert_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) + Seed_ammoph;
						 m_veget(x,y)[species] = (Vert_ammoph + Lat_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) + Seed_ammoph;
						// m_veget(x,y)[species] = (m_veget(x,y)[species] + Vert_ammoph + Lat_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) + Seed_ammoph;
					} else {
						m_veget(x,y)[species] = m_veget(x,y)[species] * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) * ((dhdt / realtimestep) < -1.0 ? 0 : 1) ; //(dhdt >= 0 ? (1-dhdt/m_veg1_Hveg) : 1) *
					}

					if(m_veg1_seed_deterministic>0 && m_veget(x,y)[species]<m_veg1_seed_deterministic && h(x,y) <= m_veg1_seed_zmax) {
						m_veget(x,y)[species] = m_veg1_seed_deterministic;
					}
            	}

            	if (species == 1)
            	{
            		if(m_veg2 && m_wrack==false) {
            			// vertical and lateral growth
						double Vert_ammoph;
						double Seed_ammoph;

						double m_veg_r_eff2 = m_veg2_r + m_veg2_r_sand * dhdt / realtimestep + m_veg2_r_sand2  * pow(dhdt / realtimestep,2) ;
						double m_veg2_D = (m_veg2_C * m_veg2_C / (4*m_veg_r_eff2)); // vegetation diffusion rate (m2/yr) for 1m grid cells
						//double m_veget_eff1 = m_veg1_D * realtimestep * laplacian_veget(x,y)[species] + m_veget(x,y)[species];

						if(m_veg12_interactions) {
							Vert_ammoph = m_veget(x,y)[species] + m_veg_r_eff2  * m_veget(x,y)[species] * (m_veg2_K - m_veg2_alpha12 * m_veget(x,y)[0] - m_veget(x,y)[species]) / m_veg2_K * realtimestep;
						} else {
							Vert_ammoph =  (m_veg2_K / (1 + ( (m_veg2_K - m_veget(x,y)[species]) / m_veget(x,y)[species]) * exp(-1 * m_veg_r_eff2 * realtimestep)));
						}
						double Lat_ammoph = m_veg2_D * realtimestep * laplacian_veget(x,y)[0]; // * m_veget(x,y)[species]; laplacian_veget index 0 is overwritten for each species, so index 0 is correct for species 2

						// colonization
						if (m_veg2_seed_probabilistic > 0.0) {
							// random seedling/rhizome fragment growth
							std::default_random_engine generator;
							generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

							double m_veg_seed_mod;

							// specifies region of overwash-facilitated seeding. Either at high water line or throughout overwash area
							bool m_ow_seed = 0;
							if(m_veg2_seed_owregion == "all") {
								 m_ow_seed = (overwash(x,y)>0);
							} else if(m_veg2_seed_owregion == "hwl") {
								bool m_ow_seed = (overwash(x,y)>0 && ((overwash(x+1,y)==0) || x==(duneglobals::nx()-1)));
							} else if (m_veg2_seed_owregion == "none") {
								// m_ow_seed=0
							}

							m_veg_seed_mod = m_veg2_seed_probabilistic * (1.-m_veg2_seed_wrack_pct) * realtimestep; // scaled by number of timesteps per year
							if(m_ow_seed) {
								m_veg_seed_mod += m_veg2_seed_probabilistic * m_veg2_seed_wrack_pct / 183.; // scaled by number of storm events per year
							}


							std::bernoulli_distribution distribution(m_veg_seed_mod);

							int number = distribution(generator);
							Seed_ammoph = (double)number/100;

						} else {
							Seed_ammoph = 0.0;
						}

    					// evolution of cover fraction
						if(m_season >= m_veg2_season_start && m_season <= m_veg2_season_end) {
							// m_veget(x,y)[species] = (m_veget(x,y)[species] + Vert_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) * ((dhdt / realtimestep) < -1.0 ? 0 : 1) + Seed_ammoph;
							//m_veget(x,y)[species] = (m_veget(x,y)[species] + Vert_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) + Seed_ammoph;
							 m_veget(x,y)[species] = (Vert_ammoph + Lat_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens2) : 1.0) + Seed_ammoph;
							// m_veget(x,y)[species] = (m_veget(x,y)[species] + Vert_ammoph + Lat_ammoph) * (overwash(x,y) > 0 ? (1.0-m_overwash_sens1) : 1.0) + Seed_ammoph;
						} else {
							m_veget(x,y)[species] = m_veget(x,y)[species] * (overwash(x,y) > 0 ? (1.0-m_overwash_sens2) : 1.0) * ((dhdt / realtimestep) < -1.0 ? 0 : 1) ; //(dhdt >= 0 ? (1-dhdt/m_veg1_Hveg) : 1) *
						}

						if(m_veg2_seed_deterministic>0 && m_veget(x,y)[species]<m_veg2_seed_deterministic && h(x,y) <= m_veg2_seed_zmax) {
							m_veget(x,y)[species] = m_veg2_seed_deterministic;
						}
            		}

            		if(m_veg2 && m_wrack) {
            			if(overwash(x,y)>0) {
            				double p = m_wrack_seed_probabilistic / 183.;
            				double p1 = p*10.;
            				std::default_random_engine generator;
							generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

							if((overwash(x+1,y)==0) || x==(duneglobals::nx()-1)) { // high water line
								std::bernoulli_distribution d(p1);
								bool m_wrack_presence = d(generator);
	    						m_veget(x,y)[1] = (double) m_wrack_presence * m_wrack_max;
							} else { // other overwashed areas
								std::bernoulli_distribution d(p);
								bool m_wrack_presence = d(generator);
	    						m_veget(x,y)[1] = (double) m_wrack_presence * m_wrack_max;
							}
						} else {
							if( m_veget(x,y)[1] > 0 ) {
								if(dhdtimestep > 0) {
									double m_wrack_eff_ht = (m_veget(x,y)[1] / m_wrack_max) * m_wrack_ht;
									m_veget(x,y)[1] = m_wrack_max * (m_wrack_eff_ht - dhdtimestep)/m_wrack_ht;
								}
								// implicitly, if dhdt < 0, then m_veget(x,y)[1] remains the same
							}
						}
            		}
				}
            }
            
            
            // limiting conditions
            if(m_veget(x,y)[species] > 1 ){
                m_veget(x,y)[species] = 1;
            }

            if(m_veg_type == 2 && m_wrack) {
				if(m_veget(x,y)[0] < 0 || shorefactor == 0 || watertable == 0){
					m_veget(x,y)[0] = 0;
				}
				if(m_veget(x,y)[1] < 0 || watertable == 0){ // remove Lveg limiting condition for wrack
					m_veget(x,y)[1] = 0;
				}
			} else {
				if(m_veget(x,y)[species] < 0 || shorefactor == 0 || watertable == 0){
					m_veget(x,y)[species] = 0;
				}
            }

            // initial condition
            if(m_veg_type == 0 || m_veg_type == 1) {
                if (m_veget(x,y)[species] == 0 && shorefactor * dhdxfactor * watertable > 0 && dhdt >= 0){
                    m_veget(x,y)[species] = (m_veg_type == 0 ? m_veget_init0 : m_veget_init1);
                }
            }
        }
    }
    
    return m_xmin;
}

/* Shift back*/
void vegetation::shiftback(const int plusminus)
{
    m_veget.ShiftOne(plusminus);
}
/*!  Saves the arrays m_u and m_rho.  */
void vegetation::save_arrays(){
    save_2d_vecarray( "veget", m_veget );
}
