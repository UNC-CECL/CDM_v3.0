/******************************************************************************
  $Id: influx.h,v 1.1 2004/09/22 14:44:21 schatz Exp $
******************************************************************************/

#ifndef INFLUX_H
#define INFLUX_H

#include "globals.h"
#include "func.h"

/*!  Parent class of all classes for setting the sand influx at the upwind
  boundary of the simulation area.  */

class influx
{
public:
  influx() {}
  virtual ~influx() {}
  
  /*!  Function which should set the first column (x=0) of the flux variable
    field \a flux.  (To be implemented by subclasses.)  */
  virtual void set( TFktScal& influx, TFktVec& flux , double &fluxin, double &angle ) = 0;
  
  static influx *create(const dunepar& p);
};


/*!  Constant influx according to the parameter "q_in".  */

class influx_const : public influx
{
public:
  influx_const(const dunepar &p);
  virtual ~influx_const() {}

  virtual void set( TFktScal& influx, TFktVec& flux , double &fluxin, double &angle );

private:
  /*!  Constant influx to use.  Parameter q_in in parameter file.  */
  double m_influx;
  /*!  Y-range for no-null influx.  Parameter y_range in parameter file.  */
  double m_y_range;
};


/*!  Influx proportional to outflux.  Bear in mind that the outflux from the
  previous simulation step is used even if the wind direction has changed in
  the mean time.  */

class influx_outflux : public influx
{
public:
  influx_outflux(const dunepar &p);
  virtual ~influx_outflux() {}

  virtual void set( TFktScal& influx, TFktVec& flux, double &fluxin, double &angle );

private:
  bool m_x_periodic;
  /*!  If true, the outflux is averaged over the width of the simulation region
    before being used as influx.  Parameter "qinout.avg" in the parameter file
    (defaults to true).  */
  bool m_average;
  /*!  Constant proportionality factor influx / outflux.  Parameter qinout.fact
    in parameter file (defaults to 1).  */
  double m_factor;
  /*!  Offset of influx relative to outflux in y direction.  Parameter
    "qinout.yoff" in parameter file (defaults to 0).  The influx values which
    move out of range when shifted by this offset are set to 0.  This offset is
    ignored if "qinout.avg" is true.  */
  int m_yoff;
};


#endif

