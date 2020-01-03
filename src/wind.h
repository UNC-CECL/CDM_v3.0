/******************************************************************************
  $Id: wind.h,v 1.5 2004/10/18 11:40:29 schatz Exp $
******************************************************************************/

#ifndef WIND_H
#define WIND_H

#include "globals.h"

class arrayinit;

/*! Base class of classes giving the wind speed and direction depending on
  time.  Defines interface for obtaining wind speed and direction and advancing
  to the next time step.  */

class wind 
{
public:
  wind() {}
  virtual ~wind() {}
  
  /*!  Should return a number between 0 and 1.  0 is the initial x direction of
    the arrays, and ascending angles mean a counterclockwise change in wind
    direction (as is the convention in mathematics).  */
  virtual double direction()= 0;
  /*!  Should return the shear velocity corresponding to the current wind
    speed.  */
  virtual double u_star()= 0;

  /*!  This function is called by classes using this class to notify the wind
    object that a time step has been completed and the wind speed and direction
    for the time \a newtime should now be returned by u_star() and direction().
    This allows these functions to be called several times per time step.  */
  virtual void advance( double newtime )= 0;

  static wind *create(const dunepar &p);
};


/*!  Class describing constant wind speed and direction.  */

class wind_const : public wind
{
public:
  wind_const(const dunepar& par);
  ~wind_const() {}
  
  /*!  Returns the constant.  */
  virtual double direction() { return m_dir; }
  /*!  Returns the constant m_ustar, "u*" from the parameter file.  */
  virtual double u_star() { return m_ustar; }
  
  /*!  Does nothing since the wind doesn't change with time.  */
  virtual void advance( double ) {}

private:
  /*!  Constant wind shear stress.  Parameter constwind.u* in parameter file.  */
  double m_ustar;
  /*!  Constant wind direction in units of 2 pi, that is a right angle
    corresponds to 0.25.  0 is the direction into which the height profile is
    rotated before being saved to disk.  This variable is read from the
    parameter constwind.dir in the parameter file and defaults to 0.  */
  double m_dir;
};

/*!  Real wind field.  Wind strength is given by an external data file which reproduces a real wind data.  */

class wind_real : public wind
{
public:
  wind_real(const dunepar &p);
  virtual ~wind_real() {}

  /*!  Returns m_dir.  */
  virtual double direction() { return m_dir; }
  /*!  Returns m_ustar.  */
  virtual double u_star() { return m_ustar; }
  
  virtual void advance( double );

private:
  /*!  Wind data.  */
  TFktScal m_wind_data;
  
  int t0, n;
  /*!  Current wind shear velocity.  */
  double m_ustar;
  /*!  Current wind direction.  */
  double m_dir;
  
  string m_filename;
};


/*!  Class producing uncorrelated statistical variation in wind speed and
  direction.  */

class wind_flatrand : public wind
{
public:
  wind_flatrand(const dunepar &p);
  virtual ~wind_flatrand() {}
  
  /*!  Returns m_dir.  */
  virtual double direction() { return m_dir; }
  /*!  Returns m_ustar.  */
  virtual double u_star() { return m_ustar; }
  
  virtual void advance( double );

private:
  /*!  Average shear velocity.  Parameter flatwind.u* in parameter file.  */
  double m_ustar0;
  /*!  Amount of variation of wind shear velocity.  Parameter flatwind.du*.
    The actual shear velocity is in the range [m_ustar0-m_dustar,
    m_ustar0+m_dustar].  */
  double m_dustar;
  /*!  Average wind direction dividecd by 2 pi.  0 is the direction into which
    the height profile is rotated before being saved to disk.  Parameter
    flatwind.avgdir.  Defaults to 0.  */
  double m_dir0;
  /*!  Variation of wind direction to both sides.  Parameter flatwind.ddir.  */
  double m_ddir;
  /*!  Current wind shear stress.  */
  double m_ustar;
  /*!  Current wind direction.  */
  double m_dir;
  /*!  State array for the random number generator.  */
  char m_statearray[256];
};


/*!  Oscillating wind.  Direction and strength change with the same sine oscillation.  */

class wind_sine : public wind
{
public:
  wind_sine(const dunepar &p);
  virtual ~wind_sine() {}

  /*!  Returns m_dir.  */
  virtual double direction() { return m_dir; }
  /*!  Returns m_ustar.  */
  virtual double u_star() { return m_ustar; }
  
  virtual void advance( double );

private:
  /*!  Period of the sine in seconds.  */
  double m_period, m_dtpeak;
  /*!  Average shear velocity.  Parameter sinewind.u* in parameter file.  */
  double m_ustar0;
  /*!  Amplitude of the sine variation of wind shear velocity.  Parameter
    sinewind.du*.  */
  double m_dustar;
  /*!  Average wind direction divided by 2 pi.  0 is the direction into which
    the height profile is rotated before being saved to disk.  Parameter
    sinewind.avgdir.  Defaults to 0.  */
  double m_dir0;
  /*!  Amplitude of sine variation of wind direction to both sides.  In units
    of 2 pi / 360 degrees.  Parameter sinewind.ddir.  */
  double m_ddir;
  /*!  Current wind shear velocity.  */
  double m_ustar;
  /*!  Current wind direction.  */
  double m_dir;
};


/*!  Bidirectional wind.  */

class wind_bi : public wind
{
public:
  wind_bi(const dunepar &p);
  ~wind_bi() {}
  
  /*!  Returns m_dir.  */
  virtual double direction() { return m_dir; }
  /*!  Returns m_ustar.  */
  virtual double u_star() { return m_ustar; }
  
  virtual void advance( double );

private:
  /*!  Time intervals for the two wind regimes.  */
  double m_interval1, m_interval2;
  /*!  Wind speeds.  */
  double m_dir1, m_dir2;
  /*!  Wind directions.  */
  double m_ustar1, m_ustar2;
  /*!  Current wind shear velocity.  */
  double m_ustar;
  /*!  Current wind direction.  */
  double m_dir;
  /*!  Time of next change of wind.  */
  double m_nextchange;
  /*!  Current wind regime, 1 or 2.  */
  int m_currentregime;
};


#endif //  WIND_H


