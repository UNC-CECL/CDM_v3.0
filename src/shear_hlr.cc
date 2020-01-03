/******************************************************************************
 $Id: shear_hlr.cc,v 1.29 2005/07/01 16:26:14 duran Exp $
 ******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "shear.h"
#include "sepbubble.h"
#include "rfftw12d.h"
#include "physics_const.h"
#include "shear_hlr.h"


////////////////////////////
// Constructor
//

shearHLR::shearHLR(const dunepar& P) :
shear(P)
{
    m_dx = duneglobals::dx();
    
    m_z0 = duneglobals::z0eff();
    m_h_cut = sqrt(2.0) * P.getdefault("hlr.cut_k", 2.)/m_dx;
    
    iNx = duneglobals::nx();
    iNy = duneglobals::ny();
    
    // Implement periodic boundary conditions just by choosing an FFT array no
    // larger than the original array (since FT on a compact interval always
    // implies periodicity).  If the boundary condition is open, the array is
    // padded with zeros to form a larger one, putting the dune far enough from
    // the boundary to be unaffected by the periodicity implied by the FT.  The
    // factor 1.8 is from Sauermann.
    
    m_fftxsize = fft::GetNextPowerOf2((int)(duneglobals::periodic_x()? iNx : iNx+0*1.8));
    m_fftysize = fft::GetNextPowerOf2((int)(duneglobals::periodic_y()? iNy : iNy+0*1.8));
    m_dkx = 2.0*M_PI/(m_fftxsize*m_dx);
    m_dky = 2.0*M_PI/(m_fftysize*m_dx);
    
    m_inner_const= 2.0*physics_constants::kappa*physics_constants::kappa;
}

/*!  Currently does nothing.  */

shearHLR::~shearHLR()
{
}

/*!  Calculation of some constants needed for computing tau.  The argument \a
 lengthscale gives a typical length scale like the length of the hill.  */

/*! Calculation of the inner layer height*/
double shearHLR::innerlayer_height(double factor)
{
    // Inner layer height
    // Factor should be higher than 8.5 to ensure delta > 1
    factor*= m_inner_const;
    double delta= log(factor);
    for(int iter= 0; iter < 6; iter++)
        delta= log(factor/delta);
    return 1.0/delta;
}
/*! Calculation of the middle layer height*/
double shearHLR::middlelayer_height(double factor)
{
    // Middle layer height
    double h_m_const= log(factor);
    // Return ln(hm/z0)
    // The factor 0.08 improves the accuracy
    return h_m_const-0.5*log(h_m_const)+0.08;
}
//Auxiliar function
double shearHLR::J(double p, double p0, double dp)
{
	if(p < p0) return 1.;
	else return exp(-((p-p0)/dp)*((p-p0)/dp));
}
//////////////////////////////////
// Calc
//

double shearHLR::CalcPertTau(TFktScal& h, TFktVec& tau)
{
    // define some constants
    
    const int iX0 = (m_fftxsize-iNx);
    const int iY0 = (m_fftysize-iNy);
    
    // ---- calc HLR's shear stress pertubation ----
    
    rfftw2d_array	fft_h(m_fftxsize, m_fftysize);
    rfftw2d_array	tau_x(m_fftxsize, m_fftysize);
    rfftw2d_array	tau_y(m_fftxsize, m_fftysize);
    
    bool naninfwarned= false;
    int x, y;
    double delta, delta_h, tau_const, kx, ky, kx2, ky2, ka, kxa, kxs, px, pxa, pxs, px_function,
    taux_re, taux_im, tauy_re= 1., tauy_im= 0., taux_re_aux, taux_im_aux, tau_const_aux, L_aux;
    
    // copy height field of dune + bubble, fill up with zeros:
    
    double MinHeight= 0;
    for ( y= 0; y < iNy; y++){
        MinHeight += h(0,y);
    }MinHeight /= iNy;
    
    //  for ( x= 0; x < iNx; x++)
    //    for ( y= 0; y < iNy; y++){
    //      if(MinHeight > h(x,y))	MinHeight= h(x,y);
    //  }
    
    //cout << "!!!!!" << MinHeight << endl;
    
    for ( x= 0; x < iX0; x++)
        for ( y= 0; y < m_fftysize; y++){
            fft_h.pos(x, y) = 0*0.3*sin(M_PI*y/32.);
        }
    for ( ; x< iX0 + iNx; x++) {
        for ( y= 0; y< iY0; y++){
            fft_h.pos(x, y) = 0.;
        }
        for ( ; y< iY0+iNy; y++){
            fft_h.pos(x, y) = h(x-iX0,y-iY0) - MinHeight;
        }
        for ( ; y< m_fftysize; y++){
            fft_h.pos(x, y) = 0.;
        }
    }
    for ( ; x < m_fftxsize; x++)
        for ( y= 0; y < m_fftysize; y++){
            fft_h.pos(x, y) = 0.;
        }
    
    fft_h.transform_forw();
    
    // Computing tau
    
    for( x= 0; x< fft_h.freq_xsize(); ++x ) {
        kx= (x!=0 ? x : 1e-100)*m_dkx;
        // PARTELI TEST
        if( x> fft_h.freq_xsize()/2 ) {
            // treat upper half of frequencies as negative
            kx -= fft_h.freq_xsize()*m_dkx;
        	kxa= -kx;	// absolute value
            kxs= -1.0;	// sign
      	}
      	else {
            kxa= kx;
            kxs= 1.0;
      	}
      	kx2= kx*kx;
        
        L_aux=	(1./(kxa*m_z0) < 10.? 10. : 1./(kxa*m_z0));
        
        delta_h= middlelayer_height(L_aux);
        delta= innerlayer_height(L_aux);
        tau_const= 2.0*(delta*delta)*(delta_h*delta_h);
        
        // auxiliar constants for shear stress calculation
        pxa= 1./exp(1./delta);
        pxs= kxs;
        px= pxs*pxa;
        
        px_function= 2.*physics_constants::euler_gamma + log(pxa);
        tau_const_aux= 1./(px_function*px_function + 0.25*M_PI*M_PI);
        taux_re_aux= (1.-J(pxa,0.01,0.045))*(0.22+sqrt(pxa/2.))-J(pxa,0.01,0.045)*tau_const_aux*px_function;
        taux_im_aux= ((1.-J(pxa,0.001,0.0001))*(0.021+sqrt(pxa/2.))+J(pxa,0.001,0.0001)*0.5*M_PI*tau_const_aux)*pxs;
        
        tauy_re= (1.-0.7*M_PI/2.*pxa)*exp(-pxa/0.9);
        tauy_im= 0.84*px*log(pxa/3.)*exp(-pxa);
        
      	for( y= 0; y< fft_h.freq_ysize(); ++y )
      	{
        	ky= (y!=0 ? y : 1e-100)*m_dky;
            ky2= ky*ky;
            ka= sqrt(kx2+ky2);
            
        	if( (!finite(fft_h.freqre(x, y)) || !finite(fft_h.freqim(x, y))) && !naninfwarned ) {
          		fprintf(stderr, "CShearHLR::CalcPertTau:  NAN or INF in fft_h after FT.\n");
          		naninfwarned= true;
        	}
            
        	if( m_h_cut > 0 ) {
                double cutoff_exp= J(ka,0,m_h_cut);
          		fft_h.freqre(x, y) *= cutoff_exp;
          		fft_h.freqim(x, y) *= cutoff_exp;
     		}
            // Auxiliar terms for the new shear stress model
            tau_const_aux= 2. + delta*(1. + ky2/kx2);
	    	taux_re= -1.0+taux_re_aux*tau_const_aux/delta;
	    	taux_im= taux_im_aux*tau_const_aux/delta;
            
            double tau_const_x= tau_const*kx2/ka;
            double tau_const_y= tau_const*kx*ky/ka;
            
            tau_x.freqre(x, y)= tau_const_x *(fft_h.freqre(x, y)*taux_re-fft_h.freqim(x, y)*taux_im);
            tau_x.freqim(x, y)= tau_const_x *(fft_h.freqre(x, y)*taux_im+fft_h.freqim(x, y)*taux_re);
            tau_y.freqre(x, y)= tau_const_y *(fft_h.freqre(x, y)*tauy_re-fft_h.freqim(x, y)*tauy_im);
        	tau_y.freqim(x, y)= tau_const_y *(fft_h.freqre(x, y)*tauy_im+fft_h.freqim(x, y)*tauy_re);
            
      	}
    }
    
    tau_x.transform_back();
    tau_y.transform_back();
    
    
    for (int y=0; y<iNy; y++) {
        for (int x=0; x<iNx; x++) {
            tau(x,y)[0] = tau_x.pos(x+iX0,y+iY0);
            tau(x,y)[1] = tau_y.pos(x+iX0,y+iY0);
        }
    }
    
    // Computing <L>
    double Int_x= 0, Int= 0, L;
    for( x= 0; x< 0.25*fft_h.freq_xsize(); ++x ){
        for( y= 0; y<0.5*fft_h.freq_ysize(); ++y ){
            double fft_h_abs = sqrt(fft_h.freqre(x, y)*fft_h.freqre(x, y) + fft_h.freqim(x, y)*fft_h.freqim(x, y));
            Int_x += x * fft_h_abs;
	    	Int += fft_h_abs;
       	}
    }
    L = Int/(Int_x*m_dkx);
    
    if( !finite(L) )	L= iNx;
    
    return L;
}



