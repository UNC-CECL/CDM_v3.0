/******************************************************************************
  $Id: influx.cc,v 1.3 2005/04/26 14:28:07 duran Exp $
******************************************************************************/

#include "influx.h"

//*****************************************************************************
//  class influx


/*!  Named constructor which chooses the right influx class according to the
  parameter file, the data from which are obtained from \a p.  The type of
  influx is determined by the parameter "influx" and can currently be "const"
  or "outflux".  */

influx *influx::create(const dunepar& p)
{
  string influxtype;

  influxtype= p.getdefault<string>("influx", "const");
  if( influxtype=="const" )
    return new influx_const(p);
  else if( influxtype=="outflux" )
    return new influx_outflux(p);
  else {
    cerr << "influx::create: ERROR: illegal value `" << influxtype << 
     "' for parameter influx.  Valid values are `const' or `outflux'." << endl;
    exit(1);
    return NULL;  // to keep the compiler happy...
  }
}




//*****************************************************************************
//  class influx_const

/*!  Reads the constant influx from the parameter "q_in" and stores it in
  m_influx.  */

influx_const::influx_const(const dunepar &p )
{
  m_influx= p.getdefault("q_in", 0.0);
  m_y_range= p.getdefault("y_range", 10000000.0);
}


/*!  Sets the x component of the column x=0 of \a flux to the constant value
  m_influx.  */

void influx_const::set(TFktScal& influx, TFktVec &flux, double &fluxin, double &angle)
{
  double y1, y2;
  y1= 0.5*duneglobals::ny()*(1.-m_y_range*cos(angle))-0.5*duneglobals::nx()*sin(angle);
  y2= 0.5*duneglobals::ny()*(1.+m_y_range*cos(angle))-0.5*duneglobals::nx()*sin(angle);
    
  for( int y= 0; y< duneglobals::ny(); ++y ) {
    influx(0, y)= (y > y1 && y< y2 ? m_influx*fluxin : 0);
  }
}


//*****************************************************************************
//  class influx_outflux

/*!  Gets values of the member variables m_average, m_factor and m_yoff from
  the parameter file.  */

influx_outflux::influx_outflux(const dunepar &p )
{
  m_x_periodic = duneglobals::periodic_x();
  m_average= p.getdefault("qinout.avg", true);
  m_factor= p.getdefault("qinout.fact", 1.0);
  //m_yoff= p.getdefault("qinout.yoff", 0);
}


void influx_outflux::set( TFktScal& influx, TFktVec &flux, double &fluxin, double &angle )
{
  if( m_average && !m_x_periodic )
  {
    double avgout;
    int y;
    
    avgout= 0;
    //  We play it safe, guard against infinities and NAN in flux
    for( y= 0; y< duneglobals::ny(); ++y ) {
      if( isfinite(flux(duneglobals::nx()-1, y)[0]) )
        avgout += flux(duneglobals::nx()-1, y)[0];
    }
    avgout *= m_factor/duneglobals::ny();
    if( !isfinite(avgout) )
      avgout= 0.0;
    for( y= 0; y< duneglobals::ny(); ++y ) {
      influx(0, y)= (avgout > 0 ? avgout : 0);
    }
  }
  else
    for( int y= 0; y< duneglobals::ny(); ++y )
	  influx(0, y)= flux(duneglobals::nx()-1, y)[0];
}

