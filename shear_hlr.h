/******************************************************************************
  $Id: shear_hlr.h,v 1.20 2005/07/01 16:26:15 duran Exp $
******************************************************************************/

#ifndef __SHEAR_HLR__
#define __SHEAR_HLR__

#include "globals.h"
#include "vec.h"
#include "func.h"
#include "shear.h"
#include "rfftw12d.h"

class sepbubble;


////////////////////////////////////////////////////////////
// CShearHLR
//


class shearHLR : public shear
{
public:
  shearHLR(const dunepar& P);
  virtual ~shearHLR();

protected:
  virtual double CalcPertTau(TFktScal& h, TFktVec& tau);
private:
  double innerlayer_height(double factor);
  double middlelayer_height(double factor);
  double J(double p, double p0, double dp);
  
  double m_z0;

  /*! 2 k_cutoff^2 for smoothing of h (zero means no smoothing)  */
  double m_h_cut;
  
  int iNx, iNy;
  /*! x size of the FFT matrices (larger than original matrices)  */
  int m_fftxsize;
  /*! y size of the FFT matrices (larger than original matrices)  */
  int m_fftysize;

  /*! grid spacing in position space  */
  double m_dx;
  /*! wave number discretisation  */
  double m_dkx, m_dky;

  /*! 2k^2  */
  double m_inner_const;

};

#endif
