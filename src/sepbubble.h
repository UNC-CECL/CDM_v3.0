/******************************************************************************
  $Id: sepbubble.h,v 1.16 2005/07/05 14:11:25 duran Exp $
******************************************************************************/

#ifndef __SEPBUBBLE_H__
#define __SEPBUBBLE_H__

#include "rfftw12d.h"
#include "vec.h"
#include "func.h"
#include "globals.h"


/*!  Base class for classes computing the surface seen by the flow, which
  includes the recirculation bubble.  */

class sepbubble
{
public:
  sepbubble() {}
  virtual ~sepbubble() {}

  /*!  Should determine the height profile including the separation bubble and
    write it to \a h_sepbub.  Besides, for backward compatibility stall should
    be set to a negative value outside the slip face (not the bubble) and a
    positive one inside.  */
  virtual void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)= 0;

  static sepbubble *create(const dunepar& p);
};


/*!  Dummy separation bubble class which simulates the absence of flow
  separation.  For testing.  */

class nosepbub : public sepbubble
{
public:
  nosepbub() {}
  virtual ~nosepbub() {}
  
  virtual void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h);
};


/*!  Separation bubble smaller than the one traditionally chosen, more in line
  with the separation bubble obtained by Sauermann's FLUENT simulation.  */

class sepbubsmall : public sepbubble
{
public:
  sepbubsmall();
  virtual ~sepbubsmall();
  
  virtual void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h);

private:
  /*!  Slope in x direction which corresponds to the limit value of obliquity
    over which the wind is diverted and no flow separation occurs (of course
    this is somewhat simplified.)  */
  double m_slopelimit;
  /*!  Difference in height between neighbouring sites corresponding to
    m_slopelimit.  (m_slopelimit * grid spacing) */
  double m_hdifflimit;
  static double obliquitylimit;
};


////////////////////////////////////////////////////////////
// CSepBubP3
//

/*!  Original separation bubble class which for each value of the y coordinate
  models the bubble as a third-order polynomial.  The reattachment point is
  taken at a distance of 4 times the height from the brink.  */

class CSepBubP3 : public sepbubble
{
public:
  CSepBubP3(const dunepar& P);
  virtual ~CSepBubP3();

  void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h);

protected:
  virtual void DetectFlowSep(TFktScal& stall, const TFktScal& h);

  /*!  Slope corresponding to the parameter tau.slip.  Threshold value for
    DetectSlipFace.  */
  double m_sepslope;
  double m_dkCut;
  bool m_bFilterK;
  bool m_bShapeAlg;
  double m_dLength;
  double m_dSlope;
  double m_dkx, m_dky;
  int    m_iSmooth;
  int    m_iNxFFT;
  int    m_iNyFFT;
  int m_iX0;
  int m_iY0;
  double m_dkCutY;
  int m_iParabFit;
  bool m_bPeriodicBound;

  rfftw1d_array* m_pfftH;
  rfftw1d_array* m_pfftHy;
  TFktScal* m_pS;
  TFktScal* m_pMask;

  class COpSepBubble;
  class COpFilterK;
  double CParabFit(double *y, int ndat);
  void gaussj(double **a, const int n, double **b, int m);
  void swap(double *a, double *b); 
  void fpolynom(double x,double p[], int np);
};

/*!  Separation bubble based on the transverse dune FLUENT simulation of Volker.  */
/*!  Elliptic separation bubble.  */

class sepbub_transverse : public sepbubble
{
public:
  sepbub_transverse(const dunepar& P);
  virtual ~sepbub_transverse();

  virtual void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h);

protected:
  void Interpol(const TFktScal& h, const int Xb, const int y);
  // interpolation variables
  double A, B;
  int n;

  int m_Smooth;
  double m_sepslope;
  double m_reattach_length;

  TFktScal* m_pS;
  TFktScal* m_pMask;
  /*!  Slope of the linear dependence of separation length/brink height on brink angle.  */
  static const double length_slope;
  /*!  Intercept of the linear dependence of separation length/brink height on brink angle.  */
  static const double length_intercept;
  /*!  Constant of proportionality between position of the centre of the ellipse fitting the seprating streamline and (interpolated height - slip face height).  */
  static const double shape_x0_slope;
  /*!  Distance (divided by interpolated dune height) by which the horizontal half axis is larger than the position of the reattachment point.  */
  static const double shape_b_offset;

};

/*!  Parabolic separation bubble.  */

class sepbub_parab : public sepbubble
{
public:
  sepbub_parab(const dunepar& P);
  virtual ~sepbub_parab();

  virtual void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h);

protected:
  int m_Smooth;
  bool m_Sepbub;
  double m_sepslope;
  /*!  Normalised length of the bubble in units of the brink height.  Parameter
    parabbubble.length.  */
  double m_reattach_length, m_Slope;
  
  bool m_x_periodic, m_y_periodic;

  TFktScal* m_pS;
  TFktScal* m_pMask;
};


/*!  Modification of traditional separation bubble shape: The criterion for
  determining the slip face is not the total derivative but the derivative in x
  direction.  */

class sepbub_P3derx : public CSepBubP3
{
public:
  sepbub_P3derx(const dunepar& P) : CSepBubP3(P) {}
  virtual ~sepbub_P3derx() {}

protected:
  virtual void DetectFlowSep(TFktScal& stall, const TFktScal& h);
};


/*! Big separation bubble for obtaining long deposition after a fixated
  Barchan.  This bubble in its current implementation is only smooth for
  isolated dunes on flat terrain (h=0).  */

class sepbub_tanh : public sepbubble
{
public:
  sepbub_tanh(const dunepar& p);
  virtual ~sepbub_tanh() {}
  
  virtual void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h);

private:
  /*!  Normalised length of the bubble in units of the brink height.  Parameter
    tanhbubble.length.  */
  double m_reattach_length;
};


/*!  Gaussian-shaped separation bubble in lower corner of slip face.  */

class sepbub_corner : public sepbubble
{
public:
  sepbub_corner(const dunepar& p);
  ~sepbub_corner() { delete m_fftHy; }

  virtual void Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h);

private:
  /*!  Width of the Gaussian as proportion of length (not height) of slip face.  */
  double m_sigma;
  /*!  Width of the part of the Gaussian which is added to the height profile;
    as proportion of length of slip face.  */
  double m_width;
  /*!  Height of the Gaussian as proportion of length (not height) of slip
    face.  */
  double m_height;
  /*!  Threshold slope (difference between neighbouring cells), according to
    the angle tau.slip (in degrees).  */
  double m_threshdiff;
  /*!  Cutoff for smoothing in y direction.  */
  double m_cutky;
  /*!  Size of Fourier array in y direction and offset of height profile.  */
  int m_iNyFFT, m_iY0;
  /*!  Array for Fourier transform for smoothing in y direction.  */
  rfftw1d_array *m_fftHy;
  /*!  Unit wave number in y direction after Fourier transform.  */
  double m_dky;
};


#endif

