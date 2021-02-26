/******************************************************************************
 $Id: wind.cc,v 1.5 2004/12/22 10:21:57 schatz Exp $
 ******************************************************************************/

#include <random>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <fstream>
#include <ctype.h>

#include "globals.h"
#include "wind.h"
#include "initsurf.h"


//*****************************************************************************
//  class wind

/*! Named constructor which creates the wind object requested by the parameter file.  */

wind *wind::create(const dunepar &p)
{
    string windtype;
    
    windtype= p.getdefault<string>( "wind", "const" );
    if( windtype=="const" )
        return new wind_const(p);
    else if( windtype=="real" )
        return new wind_real(p);
    else if( windtype=="flatrand" )
        return new wind_flatrand(p);
    else if( windtype=="sine" )
        return new wind_sine(p);
    else if( windtype=="bi" )
        return new wind_bi(p);
    else {
        cerr << "wind::create: ERROR: illegal value `" << windtype << "' for "
        << "parameter wind.  Valid values are `const', `sine' or `flatrand'.\n";
        exit(1);
        return NULL;
    }
}


//*****************************************************************************
//  class wind_const

/*!  Reads the member variables m_ustar and m_dir from the parameter file.  If
 the shear velocity parameter constwind.u* does not exist, the old parameter
 name u* is tried.  */

wind_const::wind_const(const dunepar& par)
{
    if( par.exists("constwind.u") )
        m_ustar= par.getrequired<double>( "constwind.u" );
    else
        m_ustar= par.getrequired<double>( "u" );
    m_dir= par.getdefault( "constwind.direction", 0.0 ) / 360.0;
}

//*****************************************************************************
//  class wind_sine

/*!  Reads the wind data from file.  */

wind_real::wind_real(const dunepar& par)
{
    arrayinit *init_wind;
    m_wind_data.Create( 2, 182, 1 );
    init_wind= arrayinit::create(par, "wind.");
    init_wind->init_2d_scal( m_wind_data );
    delete init_wind;
    
    t0 = 0;
    n = 0;
}

void wind_real::advance( double time )
{
    int t;
    if(t0==m_wind_data.SizeY()){
        t0= 0;
        n++;
    }
    for(t= t0; t < m_wind_data.SizeY(); t++){
        //cout << m_wind(0,t) << "	" << m_wind(1,t) << endl;
        if(time < m_wind_data(0,t) + 31536000*n){
            m_ustar= m_wind_data(1,t);
            //m_ustar= 0.4;
            break;
        }
    }
    t0= t;
    m_dir= 0;
}


//*****************************************************************************
//  class wind_flatrand

/*!  Reads the member variables m_ustar and m_dir from the parameter file.  If
 the shear velocity parameter flatwind.u* does not exist, the old parameter
 name u* is tried.  */

wind_flatrand::wind_flatrand(const dunepar& par)
{
    if( par.exists("flatwind.u") )
        m_ustar0= par.getrequired<double>( "flatwind.u" );
    else
        m_ustar0= par.getrequired<double>( "u" );
    m_dustar= par.getrequired<double>("flatwind.du");
    if( m_ustar0 - m_dustar < 0.0 )
    {
        m_ustar0= (m_ustar0 + m_dustar)/2.0;
        m_dustar= m_ustar0;
    }
    m_dir0= par.getdefault( "flatwind.avgdir", 0.0 ) / 360.0;
    m_ddir= par.getrequired<double>( "flatwind.ddir" ) / 360.0;
    
    // initstate( time(NULL), m_statearray, 256 );
    srand(1973);
}


/*!  Writes new random values to m_ustar and m_dir in the range
 [m_ustar0-m_dustar, m_ustar0+m_dustar] and [m_dir0-m_ddir, m_dir0+m_ddir],
 respectively.  */

void wind_flatrand::advance( double )
{
    char *prevrngstate;
    
    // prevrngstate = setstate(m_statearray);
    srand(1973);
    if( m_ddir==0.0 )
        m_dir= m_dir0;
    else
        m_dir= m_dir0 + m_ddir * 2.0 * ((double)rand()/(double)RAND_MAX - 0.5);
    if( m_dustar==0.0 )
        m_ustar= m_ustar0;
    else
        m_ustar= m_ustar0 +
	    m_dustar * 2.0 * ((double)rand()/(double)RAND_MAX - 0.5);
    // setstate(prevrngstate);
    srand(1973);
}


//*****************************************************************************
//  class wind_sine

/*!  Reads the member variables from the parameter file.  If the shear velocity
 parameter sinewind.u* does not exist, the old parameter name u* is tried.  */

wind_sine::wind_sine(const dunepar& par)
{
    m_dtpeak= par.getrequired<double>( "sinewind.dtpeak" );
    m_period= par.getrequired<double>( "sinewind.period" );
    if( par.exists("sinewind.u") )
        m_ustar0= par.getrequired<double>( "sinewind.u" );
    else
        m_ustar0= par.getrequired<double>( "u" );
    m_dustar= par.getrequired<double>("sinewind.du");
    /*if( m_ustar0 < m_dustar )
     {
     cerr << "wind_sine constructor: Error in parameter file: sinewind.du* must not be larger than sinewind.u*\n";
     exit(1);
     }*/
    m_dir0= par.getdefault( "sinewind.avgdir", 0.0 ) / 360.0;
    m_ddir= par.getrequired<double>( "sinewind.ddir" ) / 360.0;
}


void wind_sine::advance( double time )
{
    double value = ((1-m_dtpeak)-cos(M_PI*2.0*time/m_period))/(2-m_dtpeak);
//    m_ustar= m_ustar0 + m_dustar * (value > 0 ? 1. : 0);
    m_ustar= m_ustar0 + m_dustar * (value > 0 ? value : 0);
//      m_ustar= m_ustar0 + m_dustar * sin(M_PI*2.0*time/m_period)*sin(M_PI*2.0*time/m_period);
    m_dir= m_dir0 + m_ddir * sin(M_PI*2.0*time/m_period);
}


//*****************************************************************************
//  class wind_bi

wind_bi::wind_bi(const dunepar &p)
{
    m_interval1= p.getrequired<double>("biwind.interval1");
    m_interval2= p.getrequired<double>("biwind.interval2");
    m_dir1= p.getrequired<double>("biwind.dir1") / 360.0;
    m_dir2= p.getrequired<double>("biwind.dir2") / 360.0;
    if( p.exists("biwind.u*1") || p.exists("biwind.u*2") || !p.exists("u*") )
    {
        m_ustar1= p.getrequired<double>("biwind.u*1");
        m_ustar2= p.getrequired<double>("biwind.u*2");
    }
    else
        m_ustar1= m_ustar2= p.getrequired<double>("u*");
    
    m_dir= m_dir1;
    m_ustar= m_ustar1;
    m_currentregime= 1;
    m_nextchange= duneglobals::starttime() + m_interval1;
}


/*!  Changes the current wind direction and speed every time \a time exceeds
 the end of the current regime, m_nextchange.  */

void wind_bi::advance( double time )
{
    if( time < m_nextchange )
        return;
    if( m_currentregime==1 ) {
        m_currentregime= 2;
        m_dir= m_dir2;
        m_ustar= m_ustar2;
        m_nextchange += m_interval2;
    }
    else {
        m_currentregime= 1;
        m_dir= m_dir1;
        m_ustar= m_ustar1;
        m_nextchange += m_interval1;
    }
}

