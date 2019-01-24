#ifndef __ANALYZE_NEW_H__
#define __ANALYZE_NEW_H__

#include "globals.h"
#include "func.h"


class analyze
{
public:
  // creation
  analyze(const dunepar& P);

  void Calc(int t, double timestep, double shift_dist_x, int m_shoreline, double m_shorelinechange, int m_veget_X0,
  	double qin, double qout, double meanFlux, double meanVegetRho, const TFktScal& h, const TFktScal& m_rhoveg);

  double Center() { return m_dCenter; }

private:
  // internal helper functions
  double GetMaxPos(double f0, double f1, double f2);

private:
  double m_dTime;
  double m_dCenter;
  double m_zmin;

  bool m_bMsg;
  ofstream m_os;
};

#endif
