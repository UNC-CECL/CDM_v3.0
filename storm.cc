/******************************************************************************
 $Id: storm.cc,v 1.9 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#include <math.h>
#include <random>
#include <iostream>
//#include <Rmath.h>
//#include <R.h>
#include <algorithm>
#include <chrono>
#include <array>

using namespace std;

#include "globals.h"
#include "func.h"
#include "storm.h"
#include "avalanche.h"


//*****************************************************************************
//  class storm

storm::storm(const dunepar& p) : dunedata(p)
{
	bool calc_storm = p.getdefault("calc.storm", false);
	if(calc_storm) {
		m_avalanche = avalanche::create(p);

		m_sflux.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
		m_hst.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);

		m_storm_iter = p.getdefault("storm.totaliter", 100);
		m_storm_seed = p.getdefault("storm.seed", 0);
		m_Q = p.getdefault("storm.step", 0.1);
		m_scalefactor = p.getdefault("storm.scalefactor", 1.0);
		m_season_scalefactor = p.getdefault("storm.seasonscalefactor", 1.0);
		m_storm_season_start = p.getdefault("storm.seasonsstart", 0.);
		m_storm_season_end = p.getdefault("storm.seasonend", 1.);
		m_shape = p.getdefault("storm.gamma.shapeparameter", 3.0);
		m_scalefactor0 = p.getdefault("storm.gamma.scaleparameter", 0.4);
		Kc = 3.0e-3;
		Kb = 0.005;
		m_g = 9.8;
		Cls = 0.12;

		m_shore_HMWL = duneglobals::HMWL();
		m_watertable = duneglobals::MSL();
		m_slope_storm_angle = p.getdefault("storm.beach.slope", 1.0); // duneglobals::slope();
		m_slope_storm = tan(m_slope_storm_angle * M_PI / 180.);
		m_shoreface_slope = duneglobals::slope();
		m_repose_dyn = duneglobals::repose_dyn();

		m_Smax = p.getdefault("storm.Smax", 10000.0); // max surge?

		// storm parameters
		m_storm_morphodynamics = p.getdefault("storm.morphodynamics", 1);
		m_storm_type = p.getrequired<string>("Init-storm.timeseries");

		// storm Frequency
		const int m_freq = 183; //p.getdefault("storm.freq", 5); // freq. per year // 183 = HWE every 2 days
		m_Sdt = duneglobals::secyear()*duneglobals::timefrac()/(double)m_freq; // adjusted freq. per year (accounting for m_freq)
		m_tmax = p.getrequired<int>("Nt");
		m_dtmax = p.getdefault("dt_max", 1000.0);
		double m_simlength_yrs = m_tmax * m_dtmax /(duneglobals::secyear()*duneglobals::timefrac());

		// storm sequence array
		if(m_storm_type.compare("annual_distribution") == 0) {
			// Storm initialization from distribution
			/// ERLANG (GAMMA) DISTRIBUTION!!!
			double shape = m_shape;
			double scalefactor0 = m_scalefactor0;

			default_random_engine generator;
			if(m_storm_seed == 0) {
				generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
			} else {
				generator.seed(m_storm_seed);
			}

			gamma_distribution<double> distribution(shape,scalefactor0);

			const int days_year = 365;
			array<double,days_year> daily_TWL;
			array<double,m_freq> largest_annual_TWLs;
			int counter=0;
			int k_max = round(m_simlength_yrs) + 1;
			for(int k = 0; k < k_max; k++) {
				for(int j = 0; j < m_freq; j++) {
					if(j==0) {
						for(int i = 0; i < days_year; i++){
							daily_TWL[i] = distribution(generator); // randomly draw daily maximum TWL from gamma dist.
						}
						sort(daily_TWL.begin(),daily_TWL.end(), greater<int>()); // sort daily max TWL by decreasing
					}
					largest_annual_TWLs[j] = daily_TWL[j];
				}
				for(int j = 0; j < m_freq; j++) {
					if(j==0) {
						random_shuffle ( largest_annual_TWLs.begin(), largest_annual_TWLs.end() );
					}
					surge[counter] = largest_annual_TWLs[j]; // extract highest m_freq storms from distribution
					counter+=1;
				}
			}

		} else if (m_storm_type.compare("init_h") == 0) {
			string m_filename = p.getrequired<string>("storm.morphodynamics.init_h.file" );
			double m_storm_nsteps = m_freq * m_simlength_yrs; // (number of storms)*(total_sim_length(sec)/(sec/yr*timefrac)

			// load array
			ifstream file(m_filename);
				if(file.is_open())
				{
					int n = round(m_storm_nsteps);

					for(int i = 0; i < n; ++i)
					{
						file >> surge[i]; //read in value from file
					}
				} else {
					cout << "surge::create: FATAL ERROR: init_h file `" << m_filename << "\' does not exist for parameter `Init-storm.timeseries\' !" << endl;
					exit(1);
				}

		} else {
			cout << "gradalongshore::create: FATAL ERROR: Unknown value `" << m_storm_type << "\' for parameter `Init-storm.timeseries\' !" << endl;
			cout << "  Valid values are: annual_distribution, init_h" << endl;
			exit(1);
		}
    }
    
    //for(int i = 0; i < 3000; i++) {
    //    surge[i] = distribution(generator); // generate a sample of 3000 surge elevations from the gamma distribution for subsequent sampling
        //        cout << surge[i] << "# StormT" << endl;
    //}
    stormindex = 0;
    
}

double storm::impact( double shoreline, double time, double timestep, TFktScal &h, TFktScal &h_nonerod, TFktScal &overwash )
{
    // read the shoreline position
    m_shoreline = shoreline;
    
    overwash.SetAll(0.0);
    // --- TIME ----
	double realtimestep = timestep/duneglobals::secyear()/duneglobals::timefrac(); // years (1000 timesteps)*(timefrac sec / timestep)*(secyear sec / year) ]
	double time_yr = time/timestep*realtimestep;
	double m_season = fmod(time_yr, 1.0); // month of year, represented as a decimal (i.e., 0.00 = 1/1, 0.998 = 12/31)

 //   int iterstorm0 = 1; // 10;
	double m_season_scalefactor1 = 1.;
	if(m_season >= m_storm_season_start && m_season <= m_storm_season_end) {
		double m_season_scalefactor1 = m_season_scalefactor;
	}
    double isurge = m_scalefactor * m_season_scalefactor1 * surge[stormindex];//iterstorm0 - 100 + rand() % 300; // sampled surge * scalar

    // avoid inundation
    m_Tsurge = m_shore_HMWL + (isurge < m_Smax ? isurge : m_Smax); // actual surge. m_Smax is the maximum possible surge
    
//    for (int i=0; i<iterstorm0; i++) {
        calc(h, overwash);
        m_avalanche->calc(h, h_nonerod);

//    }

    double HMax = h.GetMax();
    if (HMax < m_shore_HMWL + 0.3 && false)
    {
        m_Tsurge *= -1.0;
    }

    return m_Tsurge;
}

void storm::stop( double time, double timestep, bool &calc_storm)
{
    
    if( calc_storm > 0 ){
        calc_storm = 0;
    } else {
        double tstep = time / timestep;
        double Sstep = m_Sdt / timestep;
        int tnextstorm = (int) tstep % (int) Sstep; //750; //1000; //500;
        calc_storm = (tnextstorm == 0 && tstep >= 0 ? 1 : 0);

        if (calc_storm > 0) {
            stormindex++;
        }
        
        cout << "!! HWE = " << calc_storm << ' ' << tstep << ' ' << stormindex << ' ' << tnextstorm << endl;
    }
}

void storm::calc( TFktScal &h, TFktScal &overwash )
{
	double itermax = (double)m_storm_iter * m_Tsurge * m_Tsurge;  //m_storm_iter;
	int iterstorm0 = 10; // 1;
    for (int iter=0; iter < itermax; iter++) { // why is the storm iterated at m_storm_iter*surge_elev^2???
    	// double m_Tsurge_scalar = dnorm((double)iter, itermax/2., itermax/4., 0) / dnorm(itermax/2., itermax/2., itermax/4, 0);
    	double m_Tsurge_scalar = 1;
    	double m_Tsurge_eff = (m_Tsurge - m_shore_HMWL) * m_Tsurge_scalar + m_shore_HMWL;
    	//double m_Tsurge_eff = m_Tsurge * m_Tsurge_scalar + m_shore_HMWL;
    	for (int i=0; i<iterstorm0; i++) {
    		Step(h, overwash, m_Tsurge_eff);
    	}
    }
}

int storm::dcrest_ident(TFktScal &h, int y)
{
	// initialize variables
	int dcrest_loc = duneglobals::nx();
	int loc_max = 1;
	bool loc = false;
	int n = 5;

	// loop to find first local maximum (region where h exceeds h of the 4 adjacent grids on either side)
	for (int x = n-1; x < duneglobals::nx()-n-1; x++)
	{
		if(loc == false) {
			int loc_ct = 0;
			for (int k = 1; k < n; k++) {
				if((h(x,y) > h(x-k,y)) && (h(x,y) > h(x+k,y)))
					{loc_ct += 1;}
			}
			if (loc_ct == (n-1)) {
				loc = true;
				dcrest_loc = x;
			}
		}
		if(h(x,y) > h(loc_max, y)) {
			loc_max = x;
		}
	}
	if(dcrest_loc == duneglobals::nx()) {
		dcrest_loc = loc_max;
	}
	return dcrest_loc;
}


void storm::Step(TFktScal &h, TFktScal &overwash, double &m_Tsurge_eff)
{
    double hnext, hi, hprev, Sfactor, hx, hxx, divq;
    double dx = duneglobals::dx();
    int m_swl = (int) floor((m_shore_HMWL-m_watertable)/m_shoreface_slope); // still water level location

    for (int y=0; y < duneglobals::ny(); y++) {
        bool cont = true;
        //int dcrest_loc = dcrest_ident(h, y);

        // identify backshore slope
        int Bf_range = 20;
        double Bf_obs = (h(m_swl+Bf_range,y)-h(m_swl,y))/(Bf_range*dx);
        double Bf_ref = m_slope_storm;
        double H0L0 = pow(m_Tsurge_eff/(1.1*(0.35*Bf_ref+sqrt(0.563*pow(Bf_ref,2)+0.004)/2)),2);
        // double R_obs = m_Tsurge_eff;
        double R_obs = 1.1*(0.35 * Bf_obs * sqrt(H0L0) + sqrt(H0L0 * (0.563 * pow(Bf_obs,2) + 0.004))/2);


        for (int x=m_shoreline; x < duneglobals::nx() && cont; x++) {

            hi = h(x,y);
            // definition of B.C
            // left: h = MHWL
            hprev = (x==m_shoreline? m_shore_HMWL : h(x-1,y));
            // right: depend on storm surge
            // surge > h -> h(x+1) = h(x) (hx = 0)
            hnext = (x==duneglobals::nx()-1? m_shore_HMWL : h(x+1,y));
            if (h(x+1,y) > R_obs) {
                // surge < h -> h(x+1) = surge and STOP
                hnext = R_obs;
                cont = false;
            }

            //if(x == dcrest_loc)
            //	cont = false;


			// Auxiliar
			hx = 0.5 * (hnext - hprev)/dx;
			hxx = (hnext - 2*hi + hprev)/dx/dx;
			Sfactor = (R_obs - hi);

			// Flux
			// in
		 //   m_sflux(x,y) = (m_slope - hx) * Sfactor * Sfactor;
			// div Q
			divq = Sfactor * ( hxx * Sfactor + 2 * (m_slope_storm - hx) * hx );

			// Evol
			if(m_storm_morphodynamics)
				h(x,y) += m_Q * divq / m_Tsurge / m_Tsurge;

            // Submerge index
			if(hi < R_obs) {
				overwash(x,y) = 1;
			}
        }
    }

//    cout << "!! SURGE = " << m_Tsurge << ' ' << x << endl;

}


/*void storm::Step(TFktScal &h, TFktScal &overwash, double &m_Tsurge_eff)
{
    double hnext, hi, hprev, hx, divq;
    double qS, qD, hx_xS, scalar;
    double dx = duneglobals::dx();
    int m_swl = (int) floor((m_shore_HMWL-m_watertable)/m_shoreface_slope)+1;

    for (int y=0; y < duneglobals::ny(); y++) {
        bool cont = true;
        //bool surfzone = true;
        bool swashzone = true;
        bool backdunezone = false;
        bool dcrest = false;
        int dcrest_loc = dcrest_ident(h, y);
        double dc = h(dcrest_loc,y);
        double qtot = 0;
        for (int x=0; x<m_swl; x++) {
        	overwash(x,y) = 1;
        }
        for (int x=m_swl; x < (duneglobals::nx()-1) && cont; x++) {

        	hi = h(x,y);
        	// definition of B.C
            // left: h = MHWL
            hprev = (x==m_swl ? (m_shore_HMWL-m_slope_storm) : h(x-1,y));
            // right: depend on storm surge
            // surge > h -> h(x+1) = h(x) (hx = 0)
            hnext = (x==duneglobals::nx()-1? m_shore_HMWL : h(x+1,y));            
            if (hnext > m_Tsurge_eff) {
                // surge < h -> h(x+1) = surge and STOP
                hnext = m_Tsurge_eff;
                cont = false;
            }
            

			// Auxiliar
			//hx = 0.5 * (hnext - hprev)/dx;
			hx = (hi - hprev)/dx;
            //hxx = (hnext - 2*hi + hprev)/dx/dx;
			//Sfactor = (m_Tsurge_eff - hi);
			R_eff = m_Tsurge_eff - m_shore_HMWL;

			// Flux
			// in
		 //   m_sflux(x,y) = (m_slope_storm - hx) * Sfactor * Sfactor;
			// div Q
			if (m_Tsurge_eff <= dc) { // collision regime
				//if(surfzone) {
				//	divq = Kc*2.*sqrt(2.*m_g)*pow(R_eff,3./2.)*(hx - m_slope_storm);
				//}

				if(swashzone) {
					divq = Kc*2.*sqrt(2.*m_g)*pow(R_eff,3./2.)*pow(1-((hi-m_shore_HMWL)/R_eff),2)*(hx - m_slope_storm);
				}

				//if(surfzone && hnext > m_shore_HMWL) {
				//	surfzone = false;
				//	swashzone = true;
				//}
				//divq = Sfactor * ( hxx * Sfactor + 2 * (m_slope_storm - hx) * hx );
			}

			if (m_Tsurge_eff > dc) { //overwash regime
				//if(surfzone) {
				//	divq = Kc*2.*sqrt(2.*m_g)*pow(R_eff,3./2.)*(hx - m_slope_storm);
				//}

				if(swashzone) {
					divq = qD + (qS-qD)*pow(1.-(hi - m_shore_HMWL )/(dc - m_shore_HMWL),3./2.); //
				}

				if(dcrest) {
					divq = qD;
					dcrest = false;
					backdunezone = true;
				}

		        if(x<=dcrest_loc) {
		          qtot += divq;
		        }

				if(backdunezone && x!=dcrest_loc) {
					scalar = qtot/(qD/Cls*(log(Cls*(duneglobals::nx()-dcrest_loc)+1)));
					//divq = -1*scalar*qD/(1.+(x-dcrest_loc)*Cls);
					divq = 0;
				}

				//if(surfzone && hnext > m_shore_HMWL) {
				//	surfzone = false;
				//	swashzone = true;
				if(x == m_swl) {
					hx_xS = (hnext - hi)/dx;
				    qS = Kc*2.*sqrt(2.*m_g)*pow(R_eff,3./2.)*(hx_xS - m_slope_storm);
					qD = Kb*2.*sqrt(2.*m_g/R_eff)*pow(R_eff-(dc-m_shore_HMWL),2);
				}

				if(swashzone && (x+1)==dcrest_loc) {
					swashzone = false;
					dcrest = true;
				}
			}

			// Evol
			h(x,y) -= m_Q * divq;
			if(h(x,y) < m_shore_HMWL) {
				h(x,y) = m_shore_HMWL;
			}
			//h(x,y) += m_Q * divq / m_Tsurge / m_Tsurge;

            // Submerge index
			if(hi < m_Tsurge_eff) {
				overwash(x,y) = 1;
			}
        }
    }

//    cout << "!! SURGE = " << m_Tsurge << ' ' << x << endl;

} */

/*!  Saves the arrays m_u and m_rho.  */

void storm::save_arrays()
{
    //	save_2d_vecarray( "grad_h", m_grad_h_up );
//    save_2d_scalarray( "h_st", m_hst);
//    save_2d_scalarray( "", m_div_q );
//	save_2d_scalarray( "flux_beach", m_sflux );
}
