/******************************************************************************
 $Id: PTG_Func2dVec.cc,v 1.2 2004/09/03 11:43:57 schatz Exp $
 ******************************************************************************/


#include <math.h>

#include "PTG_Func2dVec.h"
#include "PTG_Func2dScalar.h"
#include "globals.h"

void PTG_Func2dVec::ShiftOne(int plusminus)
{
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    PTG_Func2dVec TempVec;
    TempVec.Create( SizeX(), SizeY(), Delta() );
    
    int x, y, xaux, xend, xendaux;
    
    // ---- shiftback ----
    for (y=0; y < iNy; y++) {
        if (plusminus > 0) {
            for (x=0; x < iNx-1; x++) {
                xaux = x + 1;
                TempVec(x,y)[0] = (*this)(xaux,y)[0];
                TempVec(x,y)[1] = (*this)(xaux,y)[1];
            }
        } else {
            for (x = iNx-1; x > 0; x--) {
                xaux = x - 1;
                TempVec(x,y)[0] = (*this)(xaux,y)[0];
                TempVec(x,y)[1] = (*this)(xaux,y)[1];
            }
        }
        xend = (plusminus > 0 ? iNx-1 : 0);
        xendaux = (plusminus > 0 ? iNx-2 : 1);
        
        // missing values
        TempVec(xend,y)[0] = (plusminus > 0 ? (*this)(xendaux,y)[0] : 0);
        TempVec(xend,y)[1] = (plusminus > 0 ? (*this)(xendaux,y)[1] : 0);
    }
    // output matrix (*this)(x,y)
    *this = TempVec;
}

void PTG_Func2dVec::rescale(double factor)
{
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    int x, y;
    
    for (y=0; y < iNy; y++) {
        for (x=0; x < iNx-1; x++) {
            (*this)(x,y)[0] /= factor;
        }
    }
}


void PTG_Func2dVec::GradMid(const PTG_Func2dScalar& s)
{
    const double D = 1.0/(2.0*Delta());
    const int iNx = SizeX() - 1;
    const int iNy = SizeY() - 1;
    // const int prevy0 = (duneglobals::periodic_y())?(iNy-1):iNy;
    //   const int prevx0 = (duneglobals::periodic_x())?(iNx-1):iNx;
    //   const int nextyiNy = (duneglobals::periodic_y())?1:0;
    //   const int nextxiNx = (duneglobals::periodic_x())?1:0;
    
    const int prevy0 = iNy; const int prevx0 = iNx;
    const int nextyiNy = 0; const int nextxiNx = 0;
    
    PTG_Func2dVec TempVec;
    TempVec.Create( SizeX(), SizeY(), Delta() );
    
    // ---- inside ----
    for (int y=1; y < iNy; y++) {
        for (int x=1; x < iNx; x++) {
            TempVec(x,y)[0] = opDMid( s(x-1,y), s(x+1,y), D );
            TempVec(x,y)[1] = opDMid( s(x,y-1), s(x,y+1), D );
        }
    }
    
    // ---- y boundary ----
    for (int x=1; x < iNx; x++) {
        TempVec(x,0)[0]   = opDMid( s(x-1,0), s(x+1,0), D );
        TempVec(x,0)[1]   = opDMid( s(x,prevy0), s( x ,1), D );
        
        TempVec(x,iNy)[0] = opDMid( s(x-1,iNy), s(x+1,iNy), D );
        TempVec(x,iNy)[1] = opDMid( s(x,iNy-1), s(x,nextyiNy), D );
    }
    
    // ---- x boundary ----
    for (int y=1; y < iNy; y++) {
        TempVec(0,y)[0]   = opDMid( s(prevx0,y), s(1,y), D );
        TempVec(0,y)[1]   = opDMid( s(0,y-1), s(0,y+1), D );
        
        TempVec(iNx,y)[0] = opDMid( s(iNx-1,y), s(nextxiNx,y), D );
        TempVec(iNx,y)[1] = opDMid( s(iNx,y-1), s(iNx,y+1), D );
    }
    
    // ---- edges ----
    TempVec(0,0)[0]     = opDMid( s(prevx0,0), s(1,0), D );
    TempVec(0,0)[1]     = opDMid( s(0,prevy0), s(0,1), D );
    
    TempVec(0,iNy)[0]   = opDMid( s(prevx0,iNy), s(1,iNy), D );
    TempVec(0,iNy)[1]   = opDMid( s(0,iNy-1), s(0,nextyiNy), D );
    
    TempVec(iNx,0)[0]   = opDMid( s(iNx-1,0), s(nextxiNx,0), D );
    TempVec(iNx,0)[1]   = opDMid( s(iNx,prevy0), s(iNx,1), D );
    
    TempVec(iNx,iNy)[0] = opDMid( s(iNx-1,iNy), s(nextxiNx,iNy), D );
    TempVec(iNx,iNy)[1] = opDMid( s(iNx,iNy-1), s(iNx,nextyiNy), D );
    
    // output (*this)(x,y)
    *this = TempVec;
    
}

// Laplacian (Gradient of the gradient)
void PTG_Func2dVec::DiscreteLaplacian(const PTG_Func2dScalar& s)
{
    //const double D = 1.0/Delta();
    const int iNx = SizeX();
    const int iNy = SizeY();
    const int h = Delta();

	double x_lower, x_upper, y_lower, y_upper, x_upper_tmp, y_lower_tmp, y_upper_tmp, center;
    double dx_upper, dy_lower, dy_upper;

    for (int y=0; y < iNy; y++) {
        for (int x=0; x < iNx; x++) {
		  if( (x-h)>0 ) {
			x_lower = s(x-h,y);
		  } else {
			x_lower = 0; // dirichlet boundary; set boundary condition so that vegetation density is zero at seaward boundary
		  }

		  if((x+h) < iNx) {
			x_upper = s(x+h,y);
		  } else {
			dx_upper = s(x,y)-s(x-h,y); // neumann boundary(?);
			x_upper_tmp=s(x,y)+dx_upper;
			// limit x_upper to (0,1) bounds
			if(x_upper_tmp > 1) {
				x_upper = 1;
			} else {
				if (x_upper_tmp < 0) {
					x_upper = 0;
				} else {
					x_upper = x_upper_tmp;
				}
			}
		  }

		  if((y-h)>0) {
			y_lower = s(x,y-h);
		  } else {
			/*dy_lower = s(x,y) - s(x,y+h); // neumann boundary(?);
			y_lower_tmp = s(x,y) + dy_lower;
			// limit x_upper to (0,1) bounds
			if(y_lower_tmp > 1) {
				y_lower = 1;
			} else {
				if (y_lower_tmp < 0) {
					y_lower = 0;
				} else {
					y_lower = y_lower_tmp;
				}
			}*/
			y_lower = s(x, iNy-1); //periodic boundary condition
		  }

		  if((y+h) < iNy) {
			y_upper = s(x,y+h);
		  } else {
			/*dy_upper = s(x,y) - s(x,y-h); // # neumann boundary(?);
			y_upper_tmp = s(x,y) + dy_upper;
			// limit x_upper to (0,1) bounds
			if(y_upper_tmp > 1) {
				y_upper = 1;
			} else {
				if (y_upper_tmp < 0) {
					y_upper = 0;
				} else {
					y_upper = y_upper_tmp;
				}
			}*/
			y_upper = s(x,0); //periodic boundary condition
		  }
		center = s(x,y);
		(*this)(x,y)[0] = ((x_lower + x_upper + y_lower + y_upper - 4 * center) / (h * h));
        }
    }
}

// Gradient defined my the maximum neighbor
void PTG_Func2dVec::GradMin(const PTG_Func2dScalar& s)
{
    const double D = 1.0/Delta();
    const int iNx = SizeX() - 1;
    const int iNy = SizeY() - 1;
    
    int xprev, xnext, yprev, ynext;
    
    // ---- inside ----
    for (int y=0; y <= iNy; y++) {
        for (int x=0; x <= iNx; x++) {
            
            yprev = (y==0 ? iNy : y-1);
            ynext = (y==iNy ? 0 : y+1);
            
            xprev = (x==0 ? x : x-1);
            xnext = (x==iNx ? x : x+1);
            
            if (s(x,y) > s(xprev,y) && s(x,y) > s(xnext,y)) {
                // local maxima
                (*this)(x,y)[0] = 0;
            } else {
                if (s(xnext,y) > s(xprev,y) && xnext > x) {
                    (*this)(x,y)[0] = (s(xnext,y)-s(x,y))*D;
                } else {
                    (*this)(x,y)[0] = (s(x,y)-s(xprev,y))*D;
                }
            }
            if (s(x,y) > s(x,yprev) && s(x,y) > s(x,ynext)) {
                // local maxima
                (*this)(x,y)[1] = 0;
            } else {
                if (s(x,ynext) > s(x,yprev) && ynext > y) {
                    (*this)(x,y)[1] = (s(x,ynext)-s(x,y))*D;
                } else {
                    (*this)(x,y)[1] = (s(x,y)-s(x,yprev))*D;
                }
            }
        }
    }
    
}

void PTG_Func2dVec::GradUpWind(const PTG_Func2dScalar& s, const PTG_Func2dVec& u)
{
    const double D = 1.0/Delta();
    const int iNx = SizeX() - 1;
    const int iNy = SizeY() - 1;
    // const int prevy0 = (duneglobals::periodic_y())?(iNy-1):iNy;
    //   const int prevx0 = (duneglobals::periodic_x())?(iNx-1):iNx;
    //   const int nextyiNy = (duneglobals::periodic_y())?1:0;
    //   const int nextxiNx = (duneglobals::periodic_x())?1:0;
    const int prevy0 = iNy; const int prevx0 = iNx;
    const int nextyiNy = 0; const int nextxiNx = 0;
    
    PTG_Func2dVec TempVec;
    TempVec.Create( SizeX(), SizeY(), Delta() );
    
    
    // ---- inside ----
    for (int y=1; y < iNy; y++) {
        for (int x=1; x < iNx; x++) {
            TempVec(x,y)[0] = opDUpWind( s(x-1,y), s(x,y), s(x+1,y), D, u(x,y)[0] );
            TempVec(x,y)[1] = opDUpWind( s(x,y-1), s(x,y), s(x,y+1), D, u(x,y)[1] );
        }
    }
    
    // ---- y boundary ----
    for (int x=1; x < iNx; x++) {
        TempVec(x,0)[0] = opDUpWind( s(x-1,0), s(x,0), s(x+1,0), D, u(x,0)[0] );
        TempVec(x,0)[1] = opDUpWind( s(x,prevy0), s(x,0), s( x ,1), D, u(x,0)[1] );
        
        TempVec(x,iNy)[0] = opDUpWind( s(x-1,iNy), s(x,iNy), s(x+1,iNy), D, u(x,iNy)[0] );
        TempVec(x,iNy)[1] = opDUpWind( s(x,iNy-1), s(x,iNy), s(x,nextyiNy), D, u(x,iNy)[1] );
    }
    
    // ---- x boundary ----
    for (int y=1; y < iNy; y++) {
        TempVec(0,y)[0] = opDUpWind( s(prevx0,y), s(0,y), s(1,y), D, u(0,y)[0] );
        TempVec(0,y)[1] = opDUpWind( s(0,y-1), s(0,y), s(0,y+1), D, u(0,y)[1] );
        
        TempVec(iNx,y)[0] = opDUpWind( s(iNx-1,y), s(iNx,y), s(nextxiNx,y), D, u(iNx,y)[0] );
        TempVec(iNx,y)[1] = opDUpWind( s(iNx,y-1), s(iNx,y), s(iNx,y+1), D, u(iNx,y)[1] );
    }
    
    // ---- edges ----
    TempVec(0,0)[0] = opDUpWind( s(prevx0,0), s(0,0), s(1,0), D, u(0,0)[0] );
    TempVec(0,0)[1] = opDUpWind( s(0,prevy0), s(0,0), s(0,1), D, u(0,0)[1] );
    
    TempVec(0,iNy)[0] = opDUpWind( s(prevx0,iNy), s(0,iNy), s(1,iNy), D, u(0,iNy)[0] );
    TempVec(0,iNy)[1] = opDUpWind( s(0,iNy-1), s(0,iNy), s(0,nextyiNy), D, u(0,iNy)[1] );
    
    TempVec(iNx,0)[0] = opDUpWind( s(iNx-1,0), s(iNx,0), s(nextxiNx,0), D, u(iNx,0)[0] );
    TempVec(iNx,0)[1] = opDUpWind( s(iNx,prevy0), s(iNx,0), s(iNx,1), D, u(iNx,0)[1] );
    
    TempVec(iNx,iNy)[0] = opDUpWind( s(iNx-1,iNy), s(iNx,iNy), s(nextxiNx,iNy), D, u(iNx,iNy)[0] );
    TempVec(iNx,iNy)[1] = opDUpWind( s(iNx,iNy-1), s(iNx,iNy), s(iNx,nextyiNy), D, u(iNx,iNy)[1] );
    
    // output (*this)(x,y)
    *this = TempVec;
    
}

double PTG_Func2dVec::GetMaxAbs()
{
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    double dMax2 = 0.;
    for (int y=0; y < iNy; y++) {
        for (int x=0; x < iNx; x++) {
            double dX = (*this)(x,y)[0];
            double dY = (*this)(x,y)[1];
            double dAbs2 = dX*dX + dY*dY;
            if (dMax2 < dAbs2) {
                dMax2 = dAbs2;
            }
        }
    }
    return sqrt(dMax2);
}


/*!  Copies an array of possibly different size into this one, using linear
 interpolation.  The corner values are matched exactly, and the rectangle in
 between is stretched or shrunk as necessary.  */

void PTG_Func2dVec::copyscale( const PTG_Func2dVec source )
{
    double xscale, yscale, xweight, yweight;
    int x, y, basex, basey;
    
    xscale= (double)(source.SizeX()-1)/(double)(SizeX()-1);
    yscale= (double)(source.SizeY()-1)/(double)(SizeY()-1);
    for( x=0; x < SizeX(); ++x )
        for( y=0; y < SizeY(); ++y ) {
            xweight= (double)x * xscale;
            basex= (int)trunc(xweight);
            xweight -= trunc(xweight);
            yweight= (double)y * yscale;
            basey= (int)trunc(yweight);
            yweight -= trunc(yweight);
            // We have to be careful here, because trying to use the general formula
            // for the last x or y would access an illegal location.  We ignore small
            // deviations from integer indices to be on the safe side.
            if( xweight < 1e-20 )
                if( yweight < 1e-20 )
                    (*this)(x, y)= source(basex, basey);
                else
                    (*this)(x, y)= source(basex, basey) + xweight *
                    (source(basex+1, basey) - source(basex, basey));
                else
                    if( yweight < 1e-20 )
                        (*this)(x, y)= source(basex, basey) + yweight *
                        (source(basex, basey+1) - source(basex, basey));
                    else
                        (*this)(x, y)= source(basex, basey) +
                        xweight * (source(basex+1, basey) - source(basex, basey)) +
                        yweight * (source(basex, basey+1) - source(basex, basey) +
                                   xweight * (source(basex+1, basey+1) - source(basex, basey+1)
                                              - source(basex+1, basey) + source(basex, basey)) );
        }
}



/*!  Computes the average and standard deviation of the difference of x and y
 components between neighbouring grid sites.  These values can be used to
 determine where this vecor field has discontinuities.  Zero differences are
 ignored so the standard deviation is not influenced by a constant plane
 surrounding the structure to be investigated.  */

void PTG_Func2dVec::get_diff_stat( double *avg_x, double *avg_y,
                                  double *sigma_x, double *sigma_y )
{
    double diff, du_x_sum, du_x_sqrsum, du_y_sum, du_y_sqrsum;
    int x, y, n_du_x, n_du_y;
    
    du_x_sum= du_x_sqrsum= 0.0;
    du_y_sum= du_y_sqrsum= 0.0;
    n_du_x= n_du_y= 0;
    for( x= 0; x< SizeX()-1; ++x ) {
        for( y= 0; y< SizeY()-1; ++y ) {
            diff= fabs((*this)(x, y)[0] - (*this)(x+1, y)[0]);
            if( diff ) {
                ++n_du_x;
                du_x_sum += diff;
                du_x_sqrsum += diff*diff;
            }
            diff= fabs((*this)(x, y)[0] - (*this)(x, y+1)[0]);
            if( diff ) {
                ++n_du_x;
                du_x_sum += diff;
                du_x_sqrsum += diff*diff;
            }
            diff= fabs((*this)(x, y)[1] - (*this)(x+1, y)[1]);
            if( diff ) {
                ++n_du_y;
                du_y_sum += diff;
                du_y_sqrsum += diff*diff;
            }
            diff= fabs((*this)(x, y)[1] - (*this)(x, y+1)[1]);
            if( diff ) {
                ++n_du_y;
                du_y_sum += diff;
                du_y_sqrsum += diff*diff;
            }
        }
        diff= fabs((*this)(x, y)[0] - (*this)(x+1, y)[0]);
        if( diff ) {
            ++n_du_x;
            du_x_sum += diff;
            du_x_sqrsum += diff*diff;
        }
        diff= fabs((*this)(x, y)[1] - (*this)(x+1, y)[1]);
        if( diff ) {
            ++n_du_y;
            du_y_sum += diff;
            du_y_sqrsum += diff*diff;
        }
    }
    for( y= 0; y< SizeY()-1; ++y ) {
        diff= fabs((*this)(x, y)[0] - (*this)(x, y+1)[0]);
        if( diff ) {
            ++n_du_x;
            du_x_sum += diff;
            du_x_sqrsum += diff*diff;
        }
        diff= fabs((*this)(x, y)[1] - (*this)(x, y+1)[1]);
        if( diff ) {
            ++n_du_y;
            du_y_sum += diff;
            du_y_sqrsum += diff*diff;
        }
    }
    if( n_du_x ) {
        du_x_sum /= n_du_x;
        du_x_sqrsum /= n_du_x;
        *avg_x= du_x_sum;
        *sigma_x= sqrt(du_x_sqrsum - du_x_sum*du_x_sum);
    }
    else {
        *avg_x= 0.0;
        *sigma_x= 0.0;
    }
    if( n_du_y ) {
        du_y_sum /= n_du_y;
        du_y_sqrsum /= n_du_y;
        *avg_y= du_y_sum;
        *sigma_y= sqrt(du_y_sqrsum - du_y_sum*du_y_sum);
    }
    else {
        *avg_y= 0.0;
        *sigma_y= 0.0;
    }
}



