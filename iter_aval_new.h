/******************************************************************************
  $Id: iter_aval_new.h,v 1.6 2005/01/20 17:26:26 schatz Exp $
******************************************************************************/

#ifndef __ITER_AVAL_NEW_H__
#define __ITER_AVAL_NEW_H__

#include "func.h"
#include "globals.h"
#include "avalanche.h"

class CSaveFields;

class CIterAvalNew : public dunedata, public avalanche
{
  TFktScal *m_h;
  TFktScal *m_h_nonerod;

  /*!  Number of iterations.  */
  int m_n_iter;

  CBoundary* m_pactBoundary;

  TFktVec  m_grad_h_down;
  TFktVec  m_flux_down;

  TFktScal m_div_q;

  double m_tan_angle_repose_stat;
  double m_tan_angle_repose_dyn;
  double m_E;

  int m_largefluxwarned;

public:
  // construction

  CIterAvalNew(const dunepar& P, CBoundary* boundary);

  virtual ~CIterAvalNew() {}

  // helper functions

  double CalcGradDown();
  double Step(double max_slope);

  //  To allow using this class with the new version of the software:
//  void set_and_calc( TFktScal *h, TFktScal *stall );

  //  Reimplemented from dunedata:
  void save_arrays();

  //  Functions reimplemented from avalanche:
  virtual void calc( TFktScal &h , TFktScal &h_nonerod );
};

#endif
