/******************************************************************************
  $Id: avalanche.h,v 1.5 2004/09/23 12:41:36 schatz Exp $
******************************************************************************/

#ifndef AVALANCHE_H
#define AVALANCHE_H

#include "globals.h"
#include "func.h"
#include "PTG_FileNames.h"

#include "vec.h"

/*!  Base class of all classes computing sand relaxation by way of avalanches.  */

class avalanche
{
public:
  avalanche();
  virtual ~avalanche() {}

  virtual void calc( TFktScal &h, TFktScal &h_nonerod )= 0;
  
protected:
  /*!  Static angle of repose.  */
  double m_angle_stat;
  /*!  Tangent of the static angle of repose.  */
  double m_tan_stat;
  /*!  Grid spacing * tangent of static angle of repose.  The height difference
    between neighbouring sites the slope between which forms the static angle
    of repose.  */
  double m_dh_stat;
  /*!  Dynamic angle of repose.  */
  double m_angle_dyn;
  /*!  Tangent of the dynamic angle of repose.  */
  double m_tan_dyn;
  /*!  Grid spacing * tangent of dynamic angle of repose.  The height
    difference between neighbouring sites the slope between which forms the
    dynamic angle of repose.  */
  double m_dh_dyn;

public:
  static avalanche *create(const dunepar &p);
};


/*!  Avalanche class which does nothing.  For testing.  */

class dummyaval : public avalanche
{
public:
  dummyaval() {}
  virtual ~dummyaval() {}
  virtual void calc( TFktScal& h, TFktScal &h_nonerod ) {}
};

/*! Continuum aproach of sand relaxation by avalanches (including avalanches over solid substrates) */

class flowaval : public avalanche
{
  TFktScal *m_h;
  TFktScal *m_h_nonerod;
  
  /*!  Number of iterations.  */
  int m_n_iter;
  bool m_x_periodic, m_y_periodic;

  TFktVec  m_grad_h_down;
  TFktVec  m_flux_down;
  TFktScal m_div_q;

  double m_E;

public:
  // construction

  flowaval(const dunepar& p);

  virtual ~flowaval() {}

  // helper functions

  double CalcGradDown();
  double Step(double max_slope);

  //  Functions reimplemented from avalanche:
  virtual void calc( TFktScal &h , TFktScal &h_nonerod );
};

/*!  Cellular automaton implementation of sand relaxation by avalanches.  */

class cellaval : public avalanche
{
public:
  cellaval(const dunepar &p);
  virtual ~cellaval();
  
  virtual void calc( TFktScal &h, TFktScal &h_nonerod );

private:
  void set_sliporder( TFktScal &h );

  /*!  Convenience function for accessing the two-dimensional array m_done.  */
  bool &done(int x, int y) { return m_done[duneglobals::nx()*y + x]; }
  /*!  Auxiliary variable field for keeping track of which values have already
    been computed.  */
  bool *m_done;
  /*!  Maximal number of passes over the grid.  Parameter aval.cell.maxiter in parameter file.  */
  int m_maxiter;
  /*!  Pointer to an array in which the coordinate pairs of all sites are
    stored in the order in which they will be processed.  */
  int *m_sliporder;

  static const TFktScal *hpointer;

  static int cmp_coordpair(const void *a, const void *b);
};


#endif //  AVALANCHE_H

