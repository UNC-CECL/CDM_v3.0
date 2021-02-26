/******************************************************************************
 $Id: sepbubble.cc,v 1.20 2005/07/05 14:11:24 duran Exp $
 ******************************************************************************/

#include <math.h>
#include <limits.h>

#include "globals.h"
#include "sepbubble.h"

//*****************************************************************************
//  class sepbubble
//

sepbubble *sepbubble::create(const dunepar& p)
{
    string type= p.getdefault<string>("sepbubble", "parabolic");
    
    if( type=="P3" )
        return new CSepBubP3(p);
    else if( type=="P3derx" )
        return new sepbub_P3derx(p);
    else if( type=="tanh" )
        return new sepbub_tanh(p);
    else if( type=="corner" )
        return new sepbub_corner(p);
    else if( type=="none" )
        return new nosepbub();
    else if( type=="small" )
        return new sepbubsmall();
    else if( type=="transverse" )
        return new sepbub_transverse(p);
    else if( type=="parabolic" )
        return new sepbub_parab(p);
    else {
        cerr << "sepbubble::create: ERROR: illegal value `" << type << "' for "
        << "parameter sepbubble.  Valid values are `P3', `P3derx' or `none'.\n";
        exit(1);
        return NULL;
    }
}


//*****************************************************************************
//  class nosepbub
//


/*!  Copies \a h to \a h_sepbub and sets \a stall to -1.  */

void nosepbub::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
{
    int x, y;
    
    for( x= 0; x< duneglobals::nx(); ++x )
        for( y= 0; y< duneglobals::ny(); ++y ) {
            stall(x, y)= -1.0;
            h_sepbub(x, y)= h(x, y);
        }
}


//*****************************************************************************
//  class sepbubsmall
//

double sepbubsmall::obliquitylimit= M_PI*70/180;

/*!  Sets m_slopelimit = cos(obliquitylimit) * tan(angle of repose), that is
 the projection of the angle of repose onto the wind direction which forms the
 angle obliquitylimit with the brink.  */

sepbubsmall::sepbubsmall()
{
    m_slopelimit= cos(obliquitylimit) * tan(duneglobals::repose_stat()*M_PI/180.0);
    m_hdifflimit= m_slopelimit * duneglobals::dx();
}

sepbubsmall::~sepbubsmall() {}

#define CURV_BACKTRACK 10

void sepbubsmall::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
{
    double curvature[CURV_BACKTRACK], avgcurv, avgslope;
    int x, y, x0, backx;
    
    for( y= 0; y< duneglobals::ny(); ++y ) {
        stall(0, y)= -1.0;
        h_sepbub(0, y)= h(0, y);
        stall(1, y)= -1.0;
        h_sepbub(1, y)= h(1, y);
    }
    
    for( y= 0; y< duneglobals::ny(); ++y )
        for( x= 2; x< duneglobals::nx(); ++x )
        {
            if( h(x-1, y)-h(x, y) > m_hdifflimit ) {
                stall(x, y)= 1.0;
                avgcurv= 0.0;
                for( backx= 0; x-backx > 1 && backx < CURV_BACKTRACK; ++backx ) {
                    curvature[backx]= h(x-backx, y) + h(x-backx-2, y) - 2*h(x-backx-1, y);
                    avgcurv += curvature[backx];
                }
                avgcurv /= backx;
                avgslope= (h(x-backx+1, y) - h(x-backx, y));
                cerr << "sepbubsmall::Calc:  average curvature: " << avgcurv << ", avg slope: " << avgslope << "\n";
                x0= x;
                for( ++x; x< duneglobals::nx(); ++x ) {
                    h_sepbub(x, y)= h(x0, y) + avgslope * (x-x0+backx) + avgcurv*(x-x0)*(x-x0);
                    if( h_sepbub(x, y) < h(x, y) ) {
                        h_sepbub(x, y) = h(x, y);
                        break;
                    }
                    stall(x, y)= 1.0;
                }
            }
            else {
                stall(x, y)= -1.0;
                h_sepbub(x, y)= h(x, y);
            }
        }
    
}


//*****************************************************************************
//  class CSepBubP3
//

////////////////////////////////////////////////////////
//  helper classes
//


class CSepBubP3::COpSepBubble
{
    double x0;
    double h0;
    double L;
    double dhdx0;
    double x1;
    double h1;
    double dhdx1;
    double a3;
    double a2;
    
public:
    COpSepBubble(double a_x0, double a_h0, double a_dhdx0, double slope) :
    x0(a_x0), h0(a_h0), dhdx0(a_dhdx0)
    {
        double a = dhdx0/slope;
        L = 1.5 * h0 / slope * (1 + a*0.25 + a*a*0.125);
        //  cout << dhdx0 << " " << a << " " << slope << endl;
        
        // cout << L << endl;
        if (L < 0.1) {
            L = 0.1;
        }
        
        a3 = (2*h0/(L*L*L) + dhdx0/(L*L));
        a2 = (-3*h0/(L*L) - 2*dhdx0/L);
        //    if( L< 0 || !finite(L) )
        //	cerr << "4-argument constructor: L<0: " << L << ", a2 " << a2 << ", a3 " << a3 << ", slope " << slope << ", dhdx0 " << dhdx0 << ", h0 " << h0 <<  "\n";
    }
    
    COpSepBubble(double a_x0, double a_h0, double a_dhdx0, double a_x1, double a_h1, double a_dhdx1) :
    x0(a_x0), h0(a_h0), dhdx0(a_dhdx0), x1(a_x1), h1(a_h1), dhdx1(a_dhdx1)
    {
        L = x1 - x0;
        a2 = (3*h1-dhdx1*L-2*dhdx0*L-3*h0)/(L*L);
        a3 = (dhdx1*L-2*h1+dhdx0*L+2*h0)/(L*L*L);
        //    if( L< 0 || !finite(L) )
        //	cerr << "6-argument constructor: L<0: " << L << ", a2 " << a2 << ", a3 " << a3 << "\n";
    }
    
    double GetL() { return L; }
    
    void apply(rfftw1d_array& arr)
    {
        for( int ind= 0; ind< arr.pos_size(); ++ind )
        {
            double x = ind * duneglobals::dx() - x0;
            if( x>0 && x<L )
                arr.pos(ind) = ((a3 * x + a2) * x + dhdx0) * x + h0;
        }
    }
};



class CSepBubP3::COpFilterK
{
    double m_L;
    double m_kCut;
    
public:
    COpFilterK(double L, double kCut) : m_L(L), m_kCut(kCut) {}
    
    inline double Gauss(double x, double s) {
        return exp(-x*x / (2*s*s));
    }
    
    void apply(rfftw1d_array & arr)
    {
        double dk= 2.0*M_PI/(arr.pos_size()*duneglobals::dx());
        
        for( int k= 0; k< arr.freq_size(); ++k )
        {
            double C = Gauss(k * dk * m_L, m_kCut);
            arr.freqre(k) *= C;
            arr.freqim(k) *= C;
        }
    }
    
};


//////////////////////////////////////////////////////
// class CSepBubP3
//


CSepBubP3::CSepBubP3(const dunepar& P)
{
    m_dkCut      = P.getdefault("sepbub.k_cut", 2.0);
    m_dSlope     = P.getdefault("sepbub.slope", 0.25);
    m_iSmooth    = P.getdefault("sepbub.smooth", 0);
    m_dkCutY     = P.getdefault("sepbub.k_cut_y", 1.0);
    m_bFilterK   = P.getdefault("sepbub.filter_k", 0);
    
    m_bPeriodicBound = P.getdefault("calc.x_periodic", false);
    
    const int iNx = duneglobals::nx();
    const int iNy = duneglobals::ny();
    const double dx = P.getdefault("dx",1.0);
    
    m_sepslope = tan(M_PI*P.getdefault("sep.angle", 30.0)/180.0);
    
    m_iNxFFT = fft::GetNextPowerOf2(iNx*2);
    m_iNyFFT = fft::GetNextPowerOf2(iNy*2);
    m_iX0 = (m_iNxFFT-iNx)/2;
    m_iY0 = (m_iNyFFT-iNy)/2;
    m_pfftH = new rfftw1d_array(m_iNxFFT);
    m_pfftHy = new rfftw1d_array(m_iNyFFT);
    m_dkx = 2.0*M_PI/(m_iNxFFT*duneglobals::dx());
    m_dky = 2.0*M_PI/(m_iNyFFT*duneglobals::dx());
    m_pS = new TFktScal(iNx, iNy, dx);
    m_pMask = new TFktScal(iNx, iNy, dx);
}

CSepBubP3::~CSepBubP3()
{
    delete m_pMask;
    delete m_pS;
    delete m_pfftH;
    delete m_pfftHy;
}


void CSepBubP3::DetectFlowSep(TFktScal& stall, const TFktScal& h)
{
    TFktVec grad_h;
    grad_h.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h.GradMid(h);
    
    TFktVec grad_h_down;
    grad_h_down.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h_down.GradUpWind(h, grad_h);
    
    for (int y=0; y< duneglobals::ny(); y++) {
        for (int x=0; x < duneglobals::nx(); x++) {
            stall(x,y) = (grad_h_down(x,y)[0] < 0 && vabs(grad_h_down(x,y)) > m_sepslope ? 1.0 : -1.0);
            if(stall(x,y)< 0 && ((stall(x,y-1)>0 && stall(x,y+1)>0) || (stall(x-1,y)>0
                                                                        && stall(x+1,y)>0))) stall(x,y)=1;
        }
    }
}


void CSepBubP3::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
{
    
    DetectFlowSep( stall, h);
    
    // define some abbreviations for readability
    
    const int iNx = h.SizeX();
    const int iNy = h.SizeY();
    const double dx = h.Delta();
    
    rfftw1d_array& fftH(*m_pfftH);
    rfftw1d_array& fftHy(*m_pfftHy);
    
    // Calculation of separation bubble with min of height as end of polynom
    double MinHeight = 1e20;
    
    // calculation of the separation bubble at the brink,
    // sliced model, (2d), y is considered as a parameter!
    for (int y=0; y<iNy; y++)
    {
        // Assume no separation
        for (int x=0; x<iNx; x++) {
            h_sepbub(x,y) = h(x,y);
            if (MinHeight > h(x,y)) MinHeight = h(x,y);
        }
        
        // Check for separation. It is assumed that separation
        // occurs at the brink of the slip-face.
        int iX0Bubble = -1;
        for (int xx=4; xx<iNx-1; xx++) {
            
            if (stall(xx,y)>0. && stall(xx-1,y)<0. && h(xx,y) > h(xx+1,y))
            {
                // x:   first slip face point
                // x-1: crest/brink point (not used due to high
                //      fluctuations while moving along the grid.
                // x-n: start of separation bubble
                iX0Bubble = xx-2;
                
                // Calculate separation bubble:
                
                // The main problem is to determine the slope near
                // the brink. Due to small fluctuations in the height
                // we get large fluctuations in the derivatives and
                // even larger ones for the shape of the separation
                // bubble, which depends on the local slope near the
                // brink.
                // Therefore, we calculate a zero order sep. bubble
                // with slope zero. The resulting shape is filtered
                // using a FFT and used in order to determined the
                // slope near the brink. Finally, we calculate the
                // sep. bubble again using this slope.
                
                double x0    = (double)(m_iX0+iX0Bubble)*dx;
                double h0    = h(iX0Bubble,y);
                
                double dhdx0 = (h(iX0Bubble,y)-h(iX0Bubble-2,y))/(2.*dx);
                // Calc sep. bubble with a zero order approximation for
                // the initial slope dhdx0.
                
                for (int x=0; x<m_iNxFFT; x++)
                    fftH.pos(x)= 0.0;
                for (int x=0; x<iNx; x++)
                    fftH.pos((int)(x+m_iX0)) = h(x,y)- MinHeight;
                
                // periodic boundary in x-direction
                if (m_bPeriodicBound) {
                    for (int x=0; x<m_iX0; x++)
                        fftH.pos(x) = h((int) (iNx-m_iX0+x-2),y)-MinHeight;
                    for (int x=iNx+m_iX0; x<m_iNxFFT; x++)
                        fftH.pos(x) = h(x-iNx-m_iX0+2,y)-MinHeight;
                }
                else {
                    for (int x=0; x<m_iX0; x++)
                        fftH.pos(x) = h(0,y)-MinHeight;
                    for (int x=iNx+m_iX0; x<m_iNxFFT; x++)
                        fftH.pos(x) = h(iNx-1,y)-MinHeight;
                }
                
                
                COpSepBubble OpSepBubble0(x0, h0-MinHeight, dhdx0, m_dSlope);
                OpSepBubble0.apply(fftH);
                
                // Filter the zero order shape to get rid of the small scale
                // fluctuations in the slope and therefore a better estimate
                // for dhdx0.
                
                if (m_bFilterK) {
                    COpFilterK op(dx, m_dkCut);
                    fftH.transform_forw();
                    op.apply(fftH);
                    fftH.transform_back();
                }
                
                
                if (m_bFilterK) {
                    dhdx0 = (fftH.pos(m_iX0+iX0Bubble) - fftH.pos(m_iX0+iX0Bubble-1))/dx;
                } else dhdx0 = (h(iX0Bubble,y)-h(iX0Bubble-1,y))/dx;
                
                // calc the sep. bubble
                
                for (int x=0; x<m_iNxFFT; x++)
                    fftH.pos(x)= 0.0;
                for (int x=0; x<iNx; x++)
                    fftH.pos(x+m_iX0) = h(x,y) - MinHeight;
                
                // periodic boundary in x-direction
                if (m_bPeriodicBound) {
                    for (int x=0; x<m_iX0; x++)
                        fftH.pos(x) = h(iNx-m_iX0+x-2,y) - MinHeight;
                    for (int x=iNx+m_iX0; x<m_iNxFFT; x++)
                        fftH.pos(x) = h(x-iNx-m_iX0+2,y) - MinHeight;
                }
                
                COpSepBubble opSepBubble1(x0, h0-MinHeight, dhdx0, m_dSlope);
                opSepBubble1.apply(fftH);
                
                // copy separation bubble
                bool cutting = false;
                int x1;
                for (int x=iX0Bubble; x<iNx; x++) {
                    if ( !cutting ) {
                        h_sepbub(x,y) = fftH.pos(x+m_iX0) + MinHeight;
                    }
                    
                    if (h_sepbub(x,y) < h(x,y) && !cutting ) {
                        if (h_sepbub(x,y) < h(x,y) && x > iX0Bubble+3) {
                            cutting = false;
                            x1 = x+1;
                            if (x1 >= iNx) cutting = false; // true?
                        }
                        
                        h_sepbub(x,y) = h(x,y);
                    }
                }
                
                // if separation bubble is cutting dune field then
                // calculating new polynom with height and slope of cutting point
                if ( cutting ) {
                    double h1 = h(x1,y);
                    double dhdx1 = (h(x1,y) - h(x1-1,y))/dx;
                    COpSepBubble opSepBubble2(x0, h0, dhdx0, (double)(m_iX0+x1)*dx, h1, dhdx1 );
                    
                    opSepBubble2.apply(fftH);
                    for (int x=iX0Bubble; x <= x1; x++) {
                        if( fftH.pos(x+m_iX0) > h(x,y) )
                            h_sepbub(x,y) = fftH.pos(x+m_iX0);
                    }
                    xx += (int)(opSepBubble2.GetL()/dx);
                    cout << "separation bubble was cutting dune field" << x1 << endl;
                } else  xx += (int)(opSepBubble1.GetL()/dx);
                
                // scan after the separation bubble for the next one
                
                // if periodic boundary and separation bubble have value over iNx
                if (xx >= iNx && m_bPeriodicBound) {
                    if (xx>=m_iNxFFT) xx = m_iNxFFT-1;
                    cutting = 0;
                    for (int x = 0; x<=xx-iNx+2; x++) {
                        if( fftH.pos(x+m_iX0+iNx-2)+MinHeight > h(x,y) )
                            h_sepbub(x,y) = fftH.pos(x+m_iX0+iNx-2) + MinHeight;
                    }
                    
                }
            }	// if (sepbubble)
        } // for xx
        if (m_bPeriodicBound) {
            h_sepbub(1,y) = h_sepbub(iNx-1,y);
            h_sepbub(0,y) = h_sepbub(iNx-2,y);
        } else {
            h_sepbub(0,y) = h_sepbub(1,y);
            h_sepbub(iNx-1,y) = h_sepbub(iNx-2,y);
        }
    } // for y
    
    
    // Smooth separation bubble in order to have a smooth object in
    // the lateral direction.
    
    if (abs(m_iSmooth) > 0) {
        // We only want to smooth the sep. bubble (h_sepbub - h), NOT the surface!
        
        for (int y=0; y<iNy; y++) {
            for (int x=0; x<iNx; x++) {
                if (h_sepbub(x,y) - h(x,y) > 1e-5) {
                    (*m_pMask)(x,y) = 1.;
                } else {
                    (*m_pMask)(x,y) = -1.;
                }
            }
        }
        
        for (int i=0; i<abs(m_iSmooth); i++) {
            m_pS->Smooth(h_sepbub);
            h_sepbub.Smooth(*m_pS);
            
            for (int y=0; y<iNy; y++) {
                for (int x=0; x<iNx; x++) {
                    if ((*m_pMask)(x,y) < 0. || h_sepbub(x,y)<h(x,y)) {
                        h_sepbub(x,y) = h(x,y);
                    }
                }
            }
            
            
        }
        
    }
    
    // making fourier filtering in y-direction
    if (m_iSmooth < 0) {
        COpFilterK op(dx, m_dkCutY);
        
        // making a fourier filtering only in y-direction (Veit)
        for (int x=0; x<iNx; x++) {
            for (int y=0; y<m_iNyFFT; y++)
                fftHy.pos(y)= 0.0;
            for (int y=0; y<iNy; y++)
                fftHy.pos(y+m_iY0) = h_sepbub(x,y)-h(x,y);
            
            fftHy.transform_forw();
            op.apply(fftHy);
            fftHy.transform_back();
            
            for (int y=0; y<iNy; y++) {
                if( fftHy.pos(y+m_iY0) > 0.0 )
                    h_sepbub(x,y) = fftHy.pos(y+m_iY0) + h(x,y);
            }
            
        }
    }
    
    for( int y= 0; y< iNy; ++y ) {
        h_sepbub(0, y)= h(0, y);
        h_sepbub(iNx-1, y)= h(iNx-1, y);
    }
}



//*****************************************************************************
//  class sepbub_tranverse
//

sepbub_transverse::sepbub_transverse(const dunepar& P)
{
    m_Smooth    = P.getdefault("sepbub.smooth", 0);
    m_reattach_length    = P.getdefault("bubble.length", 1.5);
    
    const int iNx = duneglobals::nx();
    const int iNy = duneglobals::ny();
    const double dx = duneglobals::dx();
    
    m_sepslope = tan(M_PI*P.getdefault("sep.angle", 30.0)/180.0);
    
    m_pS = new TFktScal(iNx, iNy, dx);
    m_pMask = new TFktScal(iNx, iNy, dx);
}

sepbub_transverse::~sepbub_transverse()
{
    delete m_pMask;
    delete m_pS;
}

void sepbub_transverse::Interpol(const TFktScal& h, const int Xb, const int y)
{
	// Interpolation of the maximun height and relative brink position:
	double dx= duneglobals::dx();
	double Sx= 0, Sxx= 0, Sy= 0, Sxy= 0, Delta, a, aold= 1, b, bold= -1, err= 0,
	errold= 1, Cx= 0, Cxx= 0;
	int i= 0, n0= 5, N= 20;
	while(i<N && (err<= errold || aold < 0 || bold > 0)){
		int Xi= Xb-i;
		// local slope
		double Yi= (!i ? h(Xi,y)-h(Xi-1,y): 0.5*(h(Xi+1,y)-h(Xi-1,y)));
		// local 2nd derivate
		double Ci= (!i ? 0 : 0.5*(h(Xi+1,y)+h(Xi-1,y)-2*h(Xi,y))/dx);
		Sx+= Xi;
		Sxx+= Xi*Xi;
		Sy+= Yi;
		Sxy+= Xi*Yi;
		Cx+= Ci;
		Cxx+= Ci*Ci;
		n= i+1;
		if(i>n0+1){
			errold= err;
			aold= a;
			bold= b;
		}
		if(i>n0-1){
			Delta= n*Sxx - Sx*Sx;
			a= (Sxx*Sy-Sx*Sxy)/Delta/dx;
			b= (n*Sxy-Sx*Sy)/Delta/dx/dx;
			err= (Cxx-2*b*Cx+(n-1)*b*b)/(n-1);
			if(i==n0){
				errold= err;
				aold= a;
				bold= b;
			}
		}
		i++;
	}
	A= aold;
	B= bold;
	n= (n==N? n : n-1);
}

void sepbub_transverse::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
{
    int Nx= duneglobals::nx(), Ny= duneglobals::ny();
    double dx= duneglobals::dx();
    
    // Ellipse parameters
    //center & semiaxes
    double x0, y0, ab /*(a/b)^2*/, b;
    // crest position, reattachment distance, brink pos, brink height & max heigth
    double xmax, l, xb, hb, Hmax;
    // Auxiliar const
    double dhdxb, delta, Caux, factor;
    
    // Detect flow separation
    TFktVec grad_h;
    grad_h.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h.GradMid(h);
    
    TFktVec grad_h_down;
    grad_h_down.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h_down.GradUpWind(h, grad_h);
    
    for (int y=0; y< duneglobals::ny(); y++) {
        for (int x=0; x < duneglobals::nx(); x++) {
            stall(x,y) = (grad_h_down(x,y)[0] < 0 && vabs(grad_h_down(x,y)) > m_sepslope ? 1.0 : -1.0);
            if(stall(x,y)< 0 && ((stall(x,y-1)>0 && stall(x,y+1)>0) || (stall(x-1,y)>0
                                                                        && stall(x+1,y)>0))) stall(x,y)=1;
        }
    }
    
    // calculation of the separation bubble at the brink,
    // sliced model, (2d), y is considered as a parameter!
    for (int y=0; y<Ny; y++)
    {
        Hmax= 0;
        xmax= 0;
        int x= 0;
        bool crest= true;
        while(x<Nx-1){
            // Assume no separation
            h_sepbub(x,y) = h(x,y);
            if(h(x,y) > Hmax){ Hmax= h(x,y); xmax= x*dx;}
            // Check for separation. It is assumed that separation
            // occurs at the brink of the slip-face.
            if (x>4 && stall(x,y)>0. && stall(x-1,y)<0. && h(x,y) > h(x+1,y))
            {
                // x:   first slip face point
                // x-1: crest/brink point (not used due to high
                //      fluctuations while moving along the grid.
                // x-n: start of separation bubble (x0)
                int Xb= (x-2);
                xb= Xb*dx;
                hb= h(Xb,y);
                dhdxb= 0.5*(hb-h(Xb-2,y))/dx;
                // Calculation of the maximun height:
                if(dhdxb > 0){
                    // Interpolation of the maximun or crest height and its position:
                    Interpol(h, Xb, y);
                    xmax= -A/B;	// Crest position
                    Hmax= hb-0.5*B*(xb-xmax)*(xb-xmax); // Maximum interpolated height
                    //dhdxb= B*(xb-xmax); // brink slope using the parabolic fit
                    crest= false;
                }else crest= true;
                if(Hmax < hb){ Hmax= hb; xmax= xb;}
                //Reatachment point
                if(!m_reattach_length)
                    l=hb*(dhdxb > -0.1? length_slope*dhdxb + length_intercept :
                          -0.1*(length_intercept-0.1*length_slope)/dhdxb);
                else
                    //l=hb*m_reattach_length;
                    l=hb*(dhdxb > -0.1? /*1.7*/0.5*m_reattach_length*dhdxb + m_reattach_length:
                          -0.1*(m_reattach_length-0.1*/*1.7*/0.5*m_reattach_length)/dhdxb);
                
                
                // Ellipse parameters (Volker)
                /*		x0= (crest? 0: -shape_x0_slope*(Hmax - hb)) + xmax;
                 delta= xb-x0; //auxiliar constant
                 b= (delta+l)+shape_b_offset*Hmax;
                 b= b*b; // b= b^2!
                 factor= (delta+l)*(delta+l)-delta*delta; // auxiliar constant
                 if(b<delta*delta || b<(delta+l)*(delta+l)){
                 cout << "negative sqrt(), b= " << b << ", Xb= " << Xb << ", y= " << y
                 << ", xmax= " << xmax/dx
                 << ", Hmax= " << Hmax << ", dhb= " << dhdxb << ", x0= " << x0/dx << ", l= " << l/dx
                 << ", hb= " << hb << endl;
                 exit(1);
                 }
                 Caux= (sqrt(b-delta*delta)+sqrt(b-(delta+l)*(delta+l)));
                 Caux*= Caux/(factor*factor);
                 ab= hb*hb * Caux; // ab= (a/b)^2!
                 y0= 0.5 * hb * (1-Caux*factor);
                 */
                // Ellipse parameters (matching brink slope)
                double gamma, beta, Xa, D;
                gamma= hb/l;
                beta= shape_b_offset*Hmax/l + 1;
                Xa= beta*gamma + dhdxb*(beta-1);
                D= (beta-1)*(gamma+dhdxb)*(Xa+gamma);
                y0= hb-l*beta*(Xa+sqrt(D));
                ab= gamma*(gamma+2*dhdxb)-2*(gamma+dhdxb)*y0/l;
                x0= xb + dhdxb*(hb-y0)/ab;
                b= beta*l+xb-x0;
                
                while(x< Nx-1 && x< (xb+l)/dx){
                    h_sepbub(x,y)= y0 + sqrt(ab*(b*b - (x*dx-x0)*(x*dx-x0)));
                    if(x>Xb+2 && h_sepbub(x,y)<h(x,y)){
                        h_sepbub(x,y)= h(x,y);
                        break;
                    }
                    x++;
                }
                Hmax= h(x,y);
                xmax= x*dx;
            }	// if (sepbubble)
            x++;
            // scan after the separation bubble for the next one
        } // for x
    } // for y
    
    
    // Smooth separation bubble in order to have a smooth object in
    // the lateral direction.
    
    if (abs(m_Smooth) > 0) {
        // We only want to smooth the sep. bubble (h_sepbub - h), NOT the surface!
        
        for (int y=0; y<Ny; y++) {
            for (int x=0; x<Nx; x++) {
                if (h_sepbub(x,y) - h(x,y) > 1e-5) {
                    (*m_pMask)(x,y) = 1.;
                } else {
                    (*m_pMask)(x,y) = -1.;
                }
            }
        }
        
        for (int i=0; i<abs(m_Smooth); i++) {
            m_pS->Smooth(h_sepbub);
            h_sepbub.Smooth(*m_pS);
            
            for (int y=0; y<Ny; y++) {
                for (int x=0; x<Nx; x++) {
                    if ((*m_pMask)(x,y) < 0. || h_sepbub(x,y)<h(x,y)) {
                        h_sepbub(x,y) = h(x,y);
                    }
                }
            }
        }
        
    }
    
    for( int y= 0; y< Ny; ++y ) {
        h_sepbub(0, y)= h(0, y);
        h_sepbub(Nx-1, y)= h(Nx-1, y);
    }
}

//*****************************************************************************
//  class sepbub_parabolic
//

sepbub_parab::sepbub_parab(const dunepar& P)
{
    m_Sepbub    = P.getdefault("sepbub.parabolic", false);
    m_Smooth    = P.getdefault("sepbub.smooth", 6);
    m_reattach_length    = P.getdefault("bubble.length", 0);
    m_Slope     = P.getdefault("sepbub.slope", 0.2 /*0.2*/);
        
    m_x_periodic= P.getdefault("calc.x_periodic", 0);
    m_y_periodic= P.getdefault("calc.y_periodic", 0);
    
    const int iNx = duneglobals::nx();
    const int iNy = duneglobals::ny();
    const double dx = duneglobals::dx();
    
    m_sepslope = tan(M_PI*P.getdefault("sep.angle", 20.0)/180.0);
    
    m_pS = new TFktScal(iNx, iNy, dx);
    m_pMask = new TFktScal(iNx, iNy, dx);
}

sepbub_parab::~sepbub_parab()
{
    delete m_pMask;
    delete m_pS;
}

void sepbub_parab::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
{
    int Nx= duneglobals::nx(), Ny= duneglobals::ny();
    double dx= duneglobals::dx();
    
    // Auxiliar const
    double dhdxb, hb, C, B, l; //, h_plain=1e2;
    
    // Detect flow separation
    TFktVec grad_h;
    grad_h.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h.GradMid(h);
    
    TFktVec grad_h_down;
    grad_h_down.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h_down.GradUpWind(h, grad_h);
    
    TFktScal h_plain, grad_h_x;
    h_plain.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h_x.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    
    for (int y=0; y< duneglobals::ny(); y++)
        for (int x=0; x < duneglobals::nx(); x++){
            if(x < duneglobals::nx()-1 && stall(x+1,y) < 0 && stall(x,y) > 0)	stall(x+1,y)= 0;
            (*m_pMask)(x,y) = (stall(x,y) < 0 && grad_h_down(x,y)[0] < 0 && vabs(grad_h_down(x,y)) > m_sepslope ? 1.0 : -1.0);
                        
            grad_h_x(x,y)= grad_h(x,y)[0];
            
        }
    /* Smoothing h_plain */
    for (int i=0; i<10 +0*abs(m_Smooth); i++) {
        m_pS->Smooth(grad_h_x);
        grad_h_x.Smooth(*m_pS);
    }
    for (int y=0; y< duneglobals::ny(); y++){
    	for (int x=duneglobals::nx()-1; x > -1; x--){
            if((*m_pMask)(x,y) > 0 && (*m_pMask)(x,y-1) < 0 && (*m_pMask)(x,y+1) < 0)	stall(x,y) = -1.0;
            else	stall(x,y) = (*m_pMask)(x,y);
		}
	}
    for (int y=0; y< duneglobals::ny(); y++){
    	for (int x=duneglobals::nx()-1; x > -1; x--){
            if(grad_h_x(x,y)>= 0 && grad_h_x((!x ? duneglobals::nx()-2 : x-1),y)<0)	h_plain(x,y) = h(x,y);
    		else h_plain(x,y)= h_plain((x==duneglobals::nx()-1?1:x+1),y);
            if(h_plain(x,y) < 0)	h_plain(x,y)= 0;
            if(h_plain(x,y)>=h(x,y))	h_plain(x,y)= h(x,y);
		}
        if(m_x_periodic)	h_plain(duneglobals::nx()-1,y)= h_plain(0,y);
	}
    
    // calculation of the separation bubble at the brink,
    // sliced model, (2d), y is considered as a parameter!for (int i=0; i<abs(m_Smooth); i++) {
    
    for (int y=0; y<Ny; y++)
    {
        int x= 0, x_sep= 0, x_brink= 0, x_slope= 0, x_next= 0, x_index;
        while(x < Nx){
            // Assume no separation
            h_sepbub(x,y) = h(x,y);
            // Looking for separation
            if(x<4){
                if(m_x_periodic)	x_brink = Nx-1+(x-1);
            }else	x_brink = x-1;
            if(x==Nx-1)	x_next=(m_x_periodic ? 0 : x);
            else	x_next= x+1;
            if((x>3 || m_x_periodic) && stall(x,y)>0 && stall(x_brink,y)<0 && h(x,y) > h(x_next,y))
            {
                // x:   first slip face point
                // x_brink=x-1: crest/brink point (not used due to high
                //      fluctuations while moving along the grid.
                // Xb=x-2: start of separation bubble (x0)
                if(x<4){
                    if(m_x_periodic){
                        x_sep = (x>1 ? x-2 : Nx-1+(x-2));
                        x_slope = Nx-1+(x-4);
                    }
                }else{
                    x_sep = x-2;
                    x_slope = x-4;
                }
                hb= h(x_sep,y);
                dhdxb= 0.5*(hb-h(x_slope,y))/dx;
                dhdxb= (dhdxb > 0.64 ? 0.64 : dhdxb);
                hb-= h_plain(x_sep,y);
                if(!m_reattach_length){
                    double a= dhdxb/m_Slope;
                    l=hb*1.5*(1.+0.25*a+0.125*a*a)/m_Slope;
                }else	l=hb*m_reattach_length;
                if(l<1.5*hb) l=1.5*hb;
                if(m_Sepbub)
                    C=-(hb+dhdxb*l)/(l*l);
                else{
                    B=(2*hb+dhdxb*l)/(l*l*l);
                    C=-(3*hb+2*dhdxb*l)/(l*l);
                }
                int Xb = x-2;
                while(x< Xb + (int)(l/dx)){
                    double X=(x-Xb)*dx;
                    if(x > Nx-1)
                        if(m_x_periodic)	x_index = x-Nx;
                        else break;
                        else	x_index = x;
                    h_sepbub(x_index,y)= h_plain(x_sep,y) + (m_Sepbub ? C*X*X : B*X*X*X+C*X*X)+dhdxb*X+hb;
                    if(x > Xb+2 && h_sepbub(x_index,y) < h(x_index,y)){
                        h_sepbub(x_index,y)= h(x_index,y);
                        break;
                    }
                    x++;
                }
            }	// if (sepbubble)
            x++;
            // scan after the separation bubble for the next one
        } // for x
    } // for y(!x ? duneglobals::nx()-2 : x-1)
    
    
    // Smooth separation bubble in order to have a smooth object in
    // the lateral direction.
    
    if (abs(m_Smooth) > 0) {
        // We only want to smooth the sep. bubble (h_sepbub - h), NOT the surface!
        
        for (int y=0; y<Ny; y++) {
            for (int x=0; x<Nx; x++) {
                if (h_sepbub(x,y) - h(x,y) > 1e-5) {
                    (*m_pMask)(x,y) = 1.;
                } else {
                    (*m_pMask)(x,y) = -1.;
                }
            }
        }
        
        for (int i=0; i<abs(m_Smooth); i++) {
            m_pS->Smooth(h_sepbub);
            h_sepbub.Smooth(*m_pS);
            
            for (int y=0; y<Ny; y++) {
                for (int x=0; x<Nx; x++) {
                    if ((*m_pMask)(x,y) < 0. || h_sepbub(x,y)<h(x,y)) {
                        h_sepbub(x,y) = h(x,y);
                    }
                }
            }
        }
        
    }
    
    
    for (int y=0; y<Ny; y++)
        for (int x=0; x<Nx; x++)
            stall(x,y) = grad_h(x,y)[0];
}

//void sepbub_parab::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
//{
//    int Nx= duneglobals::nx(), Ny= duneglobals::ny();
//    double dx= duneglobals::dx();
//    
//    // Auxiliar const
//    double dhdxb, hb, C, B, l; //, h_plain=1e2;
//    
//    // Detect flow separation
//    TFktVec grad_h;
//    grad_h.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
//    grad_h.GradMid(h);
//    
//    TFktVec grad_h_down;
//    grad_h_down.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
//    grad_h_down.GradUpWind(h, grad_h);
//    
//    TFktScal h_plain, grad_h_x;
//    h_plain.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
//    
//    
//    h_plain.SetAll(1e-10);
//    
//    grad_h_x.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
//    
//    for (int y=0; y< duneglobals::ny(); y++)
//        for (int x=0; x < duneglobals::nx(); x++){
//            if(x < duneglobals::nx()-1 && stall(x+1,y) < 0 && stall(x,y) > 0)	stall(x+1,y)= 0;
//            (*m_pMask)(x,y) = (stall(x,y) < 0 && grad_h_down(x,y)[0] < 0 && vabs(grad_h_down(x,y)) > m_sepslope ? 1.0 : -1.0);
//            
//            grad_h_x(x,y)= grad_h(x,y)[0];
//        }
//    
//    /* Smoothing h_plain */
//    for (int i=0; i<10 +0*abs(m_Smooth); i++) {
//        m_pS->Smooth(grad_h_x);
//        grad_h_x.Smooth(*m_pS);
//    }
//    
//    for (int y=0; y< duneglobals::ny(); y++){
//    	for (int x=duneglobals::nx()-1; x > -1; x--){
//            
//            int prevy = (y==0)?((m_y_periodic)?(duneglobals::ny()-1):0):(y-1);
//            int nexty = (y==duneglobals::ny()-1)?((m_y_periodic)?0:y):(y+1);
//            if((*m_pMask)(x,y) > 0 && (*m_pMask)(x,prevy) < 0 && (*m_pMask)(x,nexty) < 0)	stall(x,y) = -1.0;
//            else	stall(x,y) = (*m_pMask)(x,y);
//		}
//	}
//    
//    for (int y=0; y< duneglobals::ny(); y++){
//    	for (int x=duneglobals::nx()-1; x > -1; x--){
//            
//            if(grad_h_x(x,y)>= 0 && grad_h_x((!x ? duneglobals::nx()-2 : x-1),y)<0)	h_plain(x,y) = h(x,y);
//    		else h_plain(x,y)= h_plain((x==duneglobals::nx()-1?1:x+1),y);
//            if(h_plain(x,y) < 0)	h_plain(x,y)= 0;
//            if(h_plain(x,y)>=h(x,y))	h_plain(x,y)= h(x,y);
//            
//		}
//        if(m_x_periodic)	h_plain(duneglobals::nx()-1,y)= h_plain(0,y);
//	}
//    
//    
//    const double length_slope= 9.89;
//    const double length_intercept= 5.84;
//    // calculation of the separation bubble at the brink,
//    // sliced model, (2d), y is considered as a parameter!for (int i=0; i<abs(m_Smooth); i++) {
//    
//    for (int y=0; y<Ny; y++)
//    {
//        int x= 0, x_sep= 0, x_brink= 0, x_slope= 0, x_next= 0, x_index;
//        while(x < Nx){
//            // Assume no separation
//            h_sepbub(x,y) = h(x,y);
//            // Looking for separation
//            if(x<4){
//                if(m_x_periodic)	x_brink = Nx-1+(x-1);
//            }else	x_brink = x-1;
//            if(x==Nx-1)	x_next=(m_x_periodic ? 0 : x);
//            else	x_next= x+1;
//            if((x>3 || m_x_periodic) && stall(x,y)>0 && stall(x_brink,y)<0 && h(x,y) > h(x_next,y))
//            {
//                // x:   first slip face point
//                // x_brink=x-1: crest/brink point (not used due to high
//                //      fluctuations while moving along the grid.
//                // Xb=x-2: start of separation bubble (x0)
//                if(x<4){
//                    if(m_x_periodic){
//                        x_sep = (x>1 ? x-2 : Nx-1+(x-2));
//                        x_slope = Nx-1+(x-4);
//                    }
//                }else{
//                    x_sep = x-2;
//                    x_slope = x-4;
//                }
//                hb= h(x_sep,y);
//                dhdxb= 0.5*(hb-h(x_slope,y))/dx;
//                dhdxb= (dhdxb > 0.64 ? 0.64 : dhdxb);
//                hb-= h_plain(x_sep,y);
//                if(!m_reattach_length){
//                    double a= dhdxb/m_Slope;
//                    l=hb*1.5*(1.+0.25*a+0.125*a*a)/m_Slope;
//                }else	l=hb*m_reattach_length;
//                if(l<1.5*hb) l=1.5*hb;
//                if(m_Sepbub)
//                    C=-(hb+dhdxb*l)/(l*l);
//                else{
//                    B=(2*hb+dhdxb*l)/(l*l*l);
//                    C=-(3*hb+2*dhdxb*l)/(l*l);
//                }
//                int Xb = x-2;
//                while(x< Xb + (int)(l/dx)){
//                    double X=(x-Xb)*dx;
//                    if(x > Nx-1)
//                        if(m_x_periodic)	x_index = x-Nx;
//                        else break;
//                        else	x_index = x;
//                    h_sepbub(x_index,y)= h_plain(x_sep,y) + (m_Sepbub ? C*X*X : B*X*X*X+C*X*X)+dhdxb*X+hb;
//                    if(x > Xb+2 && h_sepbub(x_index,y) < h(x_index,y)){
//                        h_sepbub(x_index,y)= h(x_index,y);
//                        break;
//                    }
//                    x++;
//                }
//            }	// if (sepbubble)
//            x++;
//            // scan after the separation bubble for the next one
//        } // for x
//    } // for y(!x ? duneglobals::nx()-2 : x-1)
//    
//    
//    // Smooth separation bubble in order to have a smooth object in
//    // the lateral direction.
//    
//    if (abs(m_Smooth) > 0) {
//        // We only want to smooth the sep. bubble (h_sepbub - h), NOT the surface!
//        
//        for (int y=0; y<Ny; y++) {
//            for (int x=0; x<Nx; x++) {
//                if (h_sepbub(x,y) - h(x,y) > 1e-5) {
//                    (*m_pMask)(x,y) = 1.;
//                } else {
//                    (*m_pMask)(x,y) = -1.;
//                }
//            }
//        }
//        
//        for (int i=0; i<abs(m_Smooth); i++) {
//            m_pS->Smooth(h_sepbub);
//            h_sepbub.Smooth(*m_pS);
//            
//            for (int y=0; y<Ny; y++) {
//                for (int x=0; x<Nx; x++) {
//                    if ((*m_pMask)(x,y) < 0. || h_sepbub(x,y)<h(x,y)) {
//                        h_sepbub(x,y) = h(x,y);
//                    }
//                }
//            } 
//        } 
//        
//    }
//    
//    
//    
//}


//*****************************************************************************
//  class sepbub_P3derx
//

void sepbub_P3derx::DetectFlowSep(TFktScal& stall, const TFktScal& h)
{
    double diff_threshold= m_sepslope * duneglobals::dx();
    
    for (int y=0; y< duneglobals::ny(); y++)
        stall(0, y)= -1.0;
    
    for (int y=0; y< duneglobals::ny(); y++)
        for (int x=1; x < duneglobals::nx(); x++)
            stall(x, y)= fabs(h(x, y)-h(x-1, y)) - diff_threshold;
}



//*****************************************************************************
//  class sepbub_tanh
//

sepbub_tanh::sepbub_tanh(const dunepar& p)
{
    m_reattach_length= p.getrequired<double>("tanhbubble.length");
}


void sepbub_tanh::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
{
    int x, y, maxx, xend;
    double max, hscale, lscale;
    
    h_sepbub= h;
    for( y= 1; y < duneglobals::ny()-1; y++ )
    {
        max= -DBL_MAX;
        maxx= 0;
        for( x= 0; x < duneglobals::nx(); x++ ) {
            stall(x, y)= -1.0;
            if( h(x, y)> max ) {
                max= h(x, y);
                maxx= x;
            }
        }
        stall(maxx, y)= 1.0;
        
        x= maxx;
        hscale= h(x, y) / 1.7256217;    // h(x, y) / (tanh(3)-tanh(-1))
        lscale= 3.0 / (m_reattach_length * h(x, y)/duneglobals::dx());
        xend= x + (int)ceil(3.0/lscale);
        for( ; x< duneglobals::nx() && x< xend; ++x ) {
            //      if( hscale * (tanh(-(x-maxx)*lscale + 2.0) + .76159415) < h(x, y) )
            //        break;
            h_sepbub(x, y)= hscale * (tanh(-(x-maxx)*lscale + 2.0) + .76159415);
        }
    }
    
    for( y= 0; y < duneglobals::ny(); y++ )
        for( x= 0; x < duneglobals::nx(); x++ )
            stall(x, y)= (h_sepbub(x, y) > h(x, y)? 1.0: -1.0);
}



//*****************************************************************************
//  class sepbub_corner
//

sepbub_corner::sepbub_corner(const dunepar &p)
{
    m_sigma= p.getdefault("cornersepbub.sigma", 0.5);
    m_width= p.getdefault("cornersepbub.width", 0.5);
    m_height= p.getdefault("cornersepbub.height", 0.2);
    m_threshdiff = tan(M_PI*p.getdefault("tau.slip", 30.0)/180.0) * duneglobals::dx();
    m_cutky= p.getdefault("sepbub.k_cut_y", 1.0);
    
    m_iNyFFT = fft::GetNextPowerOf2(duneglobals::ny()*2);
    m_iY0 = (m_iNyFFT-duneglobals::ny())/2;
    m_fftHy = new rfftw1d_array(m_iNyFFT);
    m_dky = 2.0*M_PI/(m_iNyFFT*duneglobals::dx());
}


void sepbub_corner::Calc(TFktScal& h_sepbub, TFktScal& stall, const TFktScal& h)
{
    int x, y;
    
    for( y= 1; y< duneglobals::ny()-1; ++y )
        for( x= 1; x< duneglobals::nx(); ++x )
            stall(x, y)= h(x-1, y) - h(x, y) - m_threshdiff;
    
    for( y= 0; y< duneglobals::ny(); ++y )
        stall(0, y)= -1.0;
    for( x= 1; x< duneglobals::nx(); ++x )
        stall(x, 0)= stall(x, duneglobals::nx()-1)= -1.0;
    
    h_sepbub= h;
    for( y= 1; y< duneglobals::ny()-1; ++y )
    {
        double maxh= -DBL_MAX;
        int maxx= 0;
        
        for( x= 1; x< duneglobals::nx(); ++x ) {
            if( h(x, y) > maxh ) {
                maxh = h(x, y);
                maxx= x;
            }
            if( stall(x, y) > 0 )
            {
                double sigma, height, width, hshift, linslope;
                int relx, minrelx, maxgx;
                
                maxh= -DBL_MAX;
                if( x - maxx > 5 ) {
                    width= 10.5 - 0.075*(x-maxx-1);
                    //  Since the first point with stall=1 is the one after the brink, we
                    //  have to subtract 1 from the distance.
                    if( y==10 )
                        cerr << "crest-brink = " << x-maxx-1 << ", width= " << width << "\n";
                }
                else {
                    double hslope;
                    
                    if( x>=10 )
                        hslope= (h(x-1, y)-h(x-10, y)) / (10*duneglobals::dx());
                    else if( x>1 )
                        hslope= (h(x-1, y)-h(0, y)) / ((x-1)*duneglobals::dx());
                    else
                        hslope= 0.0;
                    //  Curvature -0.0025 according to Parteli et al. (cond-mat/0410178);
                    //  this was only fitted for the lee side.  Only for Lencois
                    //  Maranhenses.
                    width= 10.5 + 0.075/0.0025 * hslope;
                    if( y==10 )
                        cerr << "crest-brink = " << -hslope/0.0025 << ", width= " << width << "\n";
                }
                // We will use the width variable as the width to both sides from a centre
                // point, which is half the total width:
                width /= 2;
                maxgx= (int)floor(width);
                ++x;
                while( x< duneglobals::nx() && stall(x, y) > 0.0 )
                    ++x;
                x -= (int)floor(0.1*width+0.5);
                // centre of sep. bubble deviates slightly from end of slipface
                if( x>= duneglobals::nx() )
                    continue;
                minrelx= (x-maxgx >= 0 ? -maxgx : -x);
                sigma= width;
                height= h(x+minrelx, y) - h(x, y);
                linslope= -height/(double)(maxgx-minrelx);
                height *= 0.6;
                hshift= height * exp( -width*width/(2.0*sigma*sigma) );
                for( relx= minrelx; x+relx < duneglobals::nx() && relx< maxgx; ++relx )
                    h_sepbub(x+relx, y)= h(x+minrelx, y) + linslope*(relx-minrelx)
                    - hshift +
                    height * exp( -relx*relx/(2.0*sigma*sigma) );
            }
        }
    }
    
    if( m_cutky> 0.0 ) {
        for (int x=0; x< duneglobals::nx(); x++) {
            for( y=0; y< m_iNyFFT; y++ )
                m_fftHy->pos(y)= 0.0;
            for( y=0; y < duneglobals::ny(); y++ )
                m_fftHy->pos(y+m_iY0) = h_sepbub(x,y)-h(x,y);
            
            m_fftHy->transform_forw();
            for( int k= 0; k< m_fftHy->freq_size(); ++k )
            {
                double C = exp( (k*2*M_PI/m_iNyFFT)*(k*2*M_PI/m_iNyFFT) /
                               (2.0*m_cutky*m_cutky) );
                m_fftHy->freqre(k) *= C;
                m_fftHy->freqim(k) *= C;
            }
            m_fftHy->transform_back();
            
            for( y= 0; y < duneglobals::ny(); y++ ) {
                if( m_fftHy->pos(y+m_iY0) > 0.0 )
                    h_sepbub(x,y) = m_fftHy->pos(y+m_iY0) + h(x,y);
                else
                    h_sepbub(x, y) = h(x, y);
            }
        }
    }
}


//*****************************************************************************
//  class sepbub_transverse
//

/*!  Slope of the linear dependence of separation length/brink height on brink slope.  */
const double sepbub_transverse::length_slope= 9.89; //0.195;
/*!  Intercept of the linear dependence of separation length/brink height on brink angle.  */
const double sepbub_transverse::length_intercept= 5.84;
/*!  Constant of proportionality between position of the centre of the ellipse fitting the seprating streamline and (interpolated height - slip face height).  */
const double sepbub_transverse::shape_x0_slope= 8.0;
/*!  Distance (divided by interpolated dune height) by which the horizontal half axis is larger than the position of the reattachment point.  */
const double sepbub_transverse::shape_b_offset= 0.06;
