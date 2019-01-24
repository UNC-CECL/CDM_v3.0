/******************************************************************************
  $Id: iter_aval_new.cc,v 1.11 2005/04/14 08:31:47 duran Exp $
******************************************************************************/

#include <math.h>

#include "save.h"
#include "iter_aval_new.h"


//////////////////////////////////////////////
// CIterAval
//


CIterAvalNew::CIterAvalNew(const dunepar& P, CBoundary* boundary) :
  dunedata(P),

  m_h(0),
  m_h_nonerod(0),

  m_pactBoundary(boundary),

  m_grad_h_down(duneglobals::nx(), duneglobals::ny(), duneglobals::dx()),
  m_flux_down(duneglobals::nx(), duneglobals::ny(), duneglobals::dx())
{
  m_n_iter= P.getdefault("aval.new.maxiter", 50);

  // parameter

  double dAR_stat = P.getdefault("aval.angle_repose_stat", 34.0);
  m_tan_angle_repose_stat = tan(dAR_stat * M_PI/180.0);
  double dAR_dyn = P.getdefault("aval.angle_repose_dyn", 33.0);
  m_tan_angle_repose_dyn = tan(dAR_dyn * M_PI/180.0);

  // this is the only parameter of the model and control how fast the slope is relaxed
  m_E = P.getdefault("aval.new.relax", 0.9) * duneglobals::dx();
}

void CIterAvalNew::calc( TFktScal &h , TFktScal &h_nonerod )
{
  m_h= &h;
  m_h_nonerod= &h_nonerod;

  double iter = 0;
  double max_slope = 0;

  m_largefluxwarned= 0;
  max_slope = CalcGradDown();
  
  if(max_slope > m_tan_angle_repose_stat){
     while( iter < m_n_iter && max_slope > m_tan_angle_repose_dyn){
        max_slope = Step(max_slope);
        iter ++;
     }
     if (max_slope > m_tan_angle_repose_stat)
        cout << "CIterAvalNew::calc: Solution not converged after " << iter << " iterations, max. slope = " << max_slope << " (target: " << m_tan_angle_repose_dyn << ")\n";
    else  cout << "CIterAvalNew::calc: " << iter << " iterations, final slope " << max_slope << "\n";
  }
  else cout << "CIterAvalNew::calc: no avalanche: maximal slope " << max_slope << "\n";
}

double CIterAvalNew::CalcGradDown(){

   double grad_h2, max_grad= 0;

   for (int y=1; y < duneglobals::ny()-1; y++) {
      for (int x=1; x < duneglobals::nx()-1; x++) 
      {
        if((*m_h)(x,y) < (*m_h)(x+1,y) && (*m_h)(x,y) < (*m_h)(x-1,y))
	  m_grad_h_down(x,y)[0] = 0;
	else if( (*m_h)(x+1,y) > (*m_h)(x-1,y) )
	  m_grad_h_down(x, y)[0] = -((*m_h)(x, y) - (*m_h)(x-1,y));
	else
	  m_grad_h_down(x, y)[0] = (*m_h)(x, y) - (*m_h)(x+1,y);

	if((*m_h)(x,y) < (*m_h)(x,y+1) && (*m_h)(x,y) < (*m_h)(x,y-1))
	  m_grad_h_down(x,y)[1] = 0;
	else if( (*m_h)(x,y+1) > (*m_h)(x,y-1) )
	  m_grad_h_down(x, y)[1] = -((*m_h)(x, y) - (*m_h)(x,y-1));
	else
	  m_grad_h_down(x, y)[1] = (*m_h)(x, y) - (*m_h)(x,y+1);

	m_grad_h_down(x,y)[0]/= duneglobals::dx();
	m_grad_h_down(x,y)[1]/= duneglobals::dx();

        grad_h2= ((*m_h)(x,y) > (*m_h_nonerod)(x,y)+0.005 ? m_grad_h_down(x,y)[0] * m_grad_h_down(x,y)[0] +
	    m_grad_h_down(x,y)[1] * m_grad_h_down(x,y)[1] : 0);

	 if(max_grad <= grad_h2)
	   max_grad= grad_h2;
      }
  }
  return sqrt(max_grad);
}

double CIterAvalNew::Step(double max_slope)
{
  double grad_h, grad_h_nonerod, slope_diff, q_in, q_out;
  double surfchange, minh= 0.0, maxh= 0.0;

  // flux
  for (int y=1; y < duneglobals::ny()-1; y++) {
      for (int x=1; x < duneglobals::nx()-1; x++) {
      	grad_h = m_grad_h_down(x,y)[0] * m_grad_h_down(x,y)[0] +
	    		m_grad_h_down(x,y)[1] * m_grad_h_down(x,y)[1];
	grad_h= sqrt(grad_h);
	grad_h_nonerod = ((*m_h)(x,y) - (*m_h_nonerod)(x,y))/duneglobals::dx();

	if( (*m_h)(x, y)> maxh )
	   maxh= (*m_h)(x, y);
	if( (*m_h)(x, y)< minh )
	   minh= (*m_h)(x, y);

	if( grad_h > m_tan_angle_repose_dyn && grad_h_nonerod > 0 ) {
	   slope_diff = (grad_h_nonerod < grad_h - m_tan_angle_repose_dyn ? tanh(grad_h_nonerod) :
	   		tanh(grad_h) - tanh(0.9*m_tan_angle_repose_dyn));
           m_flux_down(x,y)[0] = slope_diff * m_grad_h_down(x,y)[0]/grad_h;
	   m_flux_down(x,y)[1] = slope_diff * m_grad_h_down(x,y)[1]/grad_h;
	 }
	 else {
	   m_flux_down(x,y)[0]= 0.0;
	   m_flux_down(x,y)[1]= 0.0;
	 }
      }
  }

  // change in h
  for (int y=1; y < duneglobals::ny()-1; y++) {
    for (int x=1; x < duneglobals::nx()-1; x++) {

	 q_out = fabs(m_flux_down(x,y)[0])+fabs(m_flux_down(x,y)[1]);
	 q_in = (m_flux_down(x-1,y)[0] > 0 ? m_flux_down(x-1,y)[0]:0)-
		(m_flux_down(x+1,y)[0] < 0 ? m_flux_down(x+1,y)[0]:0)+
		(m_flux_down(x,y-1)[1] > 0 ? m_flux_down(x,y-1)[1]:0)-
		(m_flux_down(x,y+1)[1] < 0 ? m_flux_down(x,y+1)[1]:0);
	
	grad_h = m_grad_h_down(x,y)[0] * m_grad_h_down(x,y)[0] +
	    		m_grad_h_down(x,y)[1] * m_grad_h_down(x,y)[1];
	grad_h= sqrt(grad_h);
	
	surfchange= - m_E * (q_out - q_in);
	/*if( fabs(surfchange) > 0.5*(maxh-minh) ) {
	  if( m_largefluxwarned < 3 ) {
	    cerr << "CIterAvalNew::Step: warning: excessive avalanche sand flux limited.\n";
	    ++m_largefluxwarned;
	  }
	  surfchange= (surfchange > 0? 0.5*(maxh-minh) : -0.5*(maxh-minh));
	}*/
        (*m_h)(x,y) += surfchange;
    }
  }

  m_pactBoundary->Bound(*m_h);

  // max. slope is used to stop the iteration
  return CalcGradDown();
}

void CIterAvalNew::save_arrays()
{
}

