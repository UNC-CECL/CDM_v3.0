/******************************************************************************
  $Id: evolution.h,v 1.4 2004/10/07 14:10:55 schatz Exp $

New time evolution classes which will replace CPde.
******************************************************************************/

#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <fstream>
#include <string>

#include "globals.h"

/*!  Parent class of all time evolution classes.  This class just defines code
  for time bookkeeping and step counting.  */

class evolution : public dunedata
{
public:
  evolution(const dunepar &par);
  virtual ~evolution();
  
  /*!  Advances m_steps and m_time and calls step_implementation().  */
  void step() { m_started= true; m_time += step_implementation(); ++m_steps; }
  
  /*!  Number of evolution steps done so far.  */
  static int steps() { if( instance ) return instance->m_steps; 
						else return 0; }
  /*!  Time elapsed since start of evolution.  */
  static double time() { if( instance ) return instance->m_time;
  						else return 0.0; }
						 
protected:
  /*!  Implementation of the time evolution step.  To be implemented by
    subclass.  Must return time interval to be added to m_time.  */
  virtual double step_implementation()=0;
  
private:
  /*!  Flag indicating whether a simulation has started already.  Set to false
    in the constructor and to true in step().  */
  bool m_started;
  /*!  Number of evolution steps done so far.  */
  int m_steps;
  /*!  Time elapsed since start of evolution.  */
  double m_time;

  /*!  Pointer to what ist meant to be the only instance of this class.  See
    the constructor.  This exists to give static functions in the dunedata
    class a chance to obtain the evolution step (which is part of the file
    names under which the arrays are saved).  */
  static evolution *instance;
};




#endif //  EVOLUTION_
