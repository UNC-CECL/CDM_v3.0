/******************************************************************************
 $Id: PTG_Func2dScalar.cc,v 1.7 2004/10/20 14:37:39 schatz Exp $
 ******************************************************************************/

#include <iostream>
#include <math.h>
#include "globals.h"
#include "PTG_Func2dScalar.h"

void PTG_Func2dScalar::ShiftOne(int plusminus)
{
    const int iNx = SizeX();
    const int iNy = SizeY();
        
    PTG_Func2dScalar TempScal;
    TempScal.Create( SizeX(), SizeY(), Delta() );
    
    int x, y, xaux, xend, xendaux;
    
    // ---- shiftback ----
    for (y=0; y < iNy; y++) {
        if (plusminus > 0) {
            for (x=0; x < iNx-1; x++) {
                xaux = x + 1;
                TempScal(x,y) = (*this)(xaux,y);
            }
        } else {
            for (x = iNx-1; x > 0; x--) {
                xaux = x - 1;
                TempScal(x,y) = (*this)(xaux,y);
            }
        }
        xend = (plusminus > 0 ? iNx-1 : 0);
        xendaux = (plusminus > 0 ? iNx-2 : 1);
        
        // missing values
        TempScal(xend,y) = (plusminus > 0 ? (*this)(xendaux,y) : 0);
    }
    // output matrix (*this)(x,y)
    *this = TempScal;
}

void PTG_Func2dScalar::DxRight(const PTG_Func2dScalar& s)
{
    const double D = 1.0/Delta();
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    // ---- inside ----
    for (int y=0; y < iNy; y++) {
        for (int x=0; x < iNx-1; x++) {
            (*this)(x,y) = opDRight( s(x,y), s(x+1,y), D );
        }
    }
    
    // ---- x boundary ----
    for (int y=0; y < iNy; y++) {
        (*this)(iNx-1,y) = opDRight( s(iNx-1,y), s(0,y), D );
    }
}


double PTG_Func2dScalar::Integrate(int xmin) const
{
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    double dSum = 0.0;
    for (int y=0; y < iNy; y++) {
        for (int x=xmin; x < iNx; x++) {
            dSum += (*this)(x,y);
        }
    }
    dSum *= Delta()*Delta();
    return dSum;
}


double PTG_Func2dScalar::GetMax() const
{
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    double dMax = (*this)(0,0);
    for (int y=0; y < iNy; y++) {
        for (int x=0; x < iNx; x++) {
            if ((*this)(x,y) > dMax) {
                dMax = (*this)(x,y);
            }
        }
    }
    return dMax;
}

double PTG_Func2dScalar::GetFirstMax() const
{
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    int x = 2;
    bool notfound = 1;
    double hprev2, hprev, h, hnext, hnext2;
    double FirstMax;
    int Pos;
    for (int y=0; y < iNy; y++) {
        while (x < iNx-2 && notfound) {
            h = (*this)(x,y);
            hprev = (*this)(x-1,y);
            hprev2 = (*this)(x-2,y);
            hnext = (*this)(x+1,y);
            hnext2 = (*this)(x+2,y);
            
            if ( (h > hprev) && (h > hnext) && (hprev > hprev2) && (hnext > hnext2) ) {
                // first local maxima
                FirstMax = h;
                Pos = x;
                notfound = 0;
            }
            x++;
        }
    }
    return FirstMax + Pos * 10000;
}

double PTG_Func2dScalar::GetMin() const
{
    const int iNx = SizeX();
    const int iNy = SizeY();
    
    double dMin = DBL_MAX;
    for (int y=0; y < iNy; y++) {
        for (int x=0; x < iNx; x++) {
            if ((*this)(x,y) < dMin) {
                dMin = (*this)(x,y);
            }
        }
    }
    return dMin;
}


double PTG_Func2dScalar::CenterOfMassX() const
{
    double sum, xsum, centerx;
    
    sum= xsum= 0.0;
    for (int x=0; x < SizeX(); x++) {
        double linesum= 0.0;
        for (int y=0; y < SizeY(); y++)
            if( finite((*this)(x, y)) )
                linesum += (*this)(x,y);
        xsum += linesum * x;
        sum += linesum;
    }
    centerx= Delta()*xsum/sum;
    if( finite(centerx) )
        return centerx;
    else
        return Delta()*SizeX()/2.0;
}

// Laplacian
void PTG_Func2dScalar::DivGrad(const PTG_Func2dScalar& s)
{
    const double D2 = 1.0/Delta()/Delta();
    const int iNx = SizeX() - 1;
    const int iNy = SizeY() - 1;
    
    int xprev, xnext, yprev, ynext;
    
    double shere, sprevx, sprevy, snextx, snexty;
    
    // ---- inside ----
    for (int y=0; y <= iNy; y++) {
        for (int x=0; x <= iNx; x++) {
            // ---- P BC ----
            yprev = (y==0 ? iNy : y-1);
            ynext = (y==iNy ? 0 : y+1);
            // ---- O BC ----
            xprev = (x==0 ? x : x-1);
            xnext = (x==iNx ? x : x+1);
            
            shere = s(x,y);
            sprevx = s(xprev,y);
            sprevy = s(x,yprev);
            snextx = s(xnext,y);
            snexty = s(x,ynext);
            
            (*this)(x,y) = D2 * ( (snextx + sprevx - 2*shere) + (snexty + sprevy - 2*shere) );
            
        }
    }
    
}

void PTG_Func2dScalar::Smooth(const PTG_Func2dScalar& s)
{
    const int iNx = SizeX() - 1;
    const int iNy = SizeY() - 1;
    //  const int prevy0 = (duneglobals::periodic_y())?(iNy-1):iNy;
    // const int prevx0 = (duneglobals::periodic_x())?(iNx-1):iNx;
    // const int nextyiNy = (duneglobals::periodic_y())?1:0;
    // const int nextxiNx = (duneglobals::periodic_x())?1:0;
    const int prevy0 = iNy; const int prevx0 = iNx;
    const int nextyiNy = 0; const int nextxiNx = 0;
    
    PTG_Func2dScalar TempScal;
    TempScal.Create( SizeX(), SizeY(), Delta() );
    
    // ---- inside ----
    for (int y=1; y < (iNy+1)/2; y++) {
        for (int x=1; x < iNx; x++) {
            TempScal(x,y) = 0.125 * (4.*s(x,y) + s(x-1,y) + s(x+1,y) +
                                     s(x,y-1) + s(x,y+1) );
        }
    }
    for(int y=(iNy+1)/2; y < iNy; y++) {
        for (int x=1; x < iNx; x++) {
            TempScal(x,y) = 0.125 * (4.*s(x,y) + s(x-1,y) + s(x+1,y) +
                                     s(x,y+1) + s(x,y-1) );
        }
    }
    
    // ---- y boundary ----
    for (int x=1; x < iNx; x++) {
        TempScal(x,0)   = 0.125 * (4.*s(x,0) + s(x-1,0) + s(x+1,0) +
                                   s(x,prevy0) + s(x,1) );
        
        TempScal(x,iNy) = 0.125 * (4.*s(x,iNy) + s(x-1,iNy) + s(x+1,iNy) +
                                   s(x,iNy-1) + s(x,nextyiNy) );
    }
    
    // ---- x boundary ----
    for (int y=1; y < iNy; y++) {
        TempScal(0,y)   = 0.125 * (4.*s(0,y) + s(prevx0,y) + s(1,y) +
                                   s(0,y-1) + s(0,y+1) );
        
        TempScal(iNx,y) = 0.125 * (4.*s(iNx,y) + s(iNx-1,y) + s(nextxiNx,y) +
                                   s(iNx,y-1) + s(iNx,y+1) );
    }
    
    // ---- edges ----
    TempScal(0,0)     = 0.125 * (4.*s(0,0) + s(prevx0,0) + s(1,0) +
                                 s(0,prevy0) + s(0,1) );
    
    TempScal(0,iNy)   = 0.125 * (4.*s(0,iNy) + s(prevx0,iNy) + s(1,iNy) +
                                 s(0,iNy-1) + s(0,nextyiNy) );
    
    TempScal(iNx,0)   = 0.125 * (4.*s(iNx,0) + s(iNx-1,0) + s(nextxiNx,0) +
                                 s(iNx,prevy0) + s(iNx,1) );
    
    TempScal(iNx,iNy) = 0.125 * (4.*s(iNx,iNy) + s(iNx-1,iNy) + s(nextxiNx,iNy) +
                                 s(iNx,iNy-1) + s(iNx,nextyiNy) );
    
    // output matrix (*this)(x,y)
    *this = TempScal;
    
}



/*!  Copies an array of possibly different size into this one, using linear
 interpolation.  The corner values are matched exactly, and the rectangle in
 between is stretched or shrunk as necessary.  */

void PTG_Func2dScalar::copyscale( const PTG_Func2dScalar source )
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
            // We have to be careful here, because rounding errors might lead to the
            // general formula being used for the last x or y, which would mean
            // accessing an illegal location.
            if( basey+1 >= source.SizeY() )
                if( basex+1 >= source.SizeX() )
                    (*this)(x, y)= source(basex, basey);
                else
                    (*this)(x, y)= source(basex, basey) + xweight *
                    (source(basex+1, basey) - source(basex, basey));
                else
                    if( basex+1 >= source.SizeX() )
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



/*!  Rotates the scalar array by \a angle in the counter-clockwise direction
 with respect to the centre of the site (\a centre_x, \a centre_y).  \a angle
 must be in radians.  The rotation is not performed by interpolating the array
 at the new grid points, since that does not conserve the sum over all grid
 sites.  Rather, a grid rotated by the same angle in the opposite direction is
 superimposed over the old one, and contributions from each old square
 intersecting with a new one are added to it in proportion to the intersection
 area.
 
 Due to the significant errors either in the computation of the tangent or the
 constant M_PI, it is highly recommended to choose \a angle in the range from
 -pi to pi.
 
 In most cases the rotated grid will cover areas which were not defined in the
 original unrotated grid, since a rectangle or a square are not identical with
 their preimage.  This function computes the average value of all values on
 the boundary of the original array and assumes that all out-of-range sites
 have this value.  This can be problematic if the height field is not constant
 outside a maximal circle around the centre site.  The amount of sand may not
 be preserved and - even worse - there may be small discontinuities in the
 rotated array.  The only way to avoid this problem is to choose an array so
 large that it is constant on the boundary.  */

void PTG_Func2dScalar::rotate( double angle, int centre_x, int centre_y, bool rot)
{
    // Maximum difference between angle recognised as right angle and a multiple
    // of pi/2.
    const double rightangle_tolerance= 1e-10;
    
    PTG_Func2dScalar dest;
    //TFktScal* pS;
    double corner_x[4], corner_y[4], corn0_x[4], corn0_y[4], limits[16];
    double thecos, thesin, upslope, downslope, slantlarge, slantsmall;
    double leftupper, leftlower, rightupper, rightlower, meanupper, meanlower;
    double off_x, off_y, boundavg, contrib, outofrange;
    int x, src_x, y, i, lim_ind, old_x, old_y;
    bool rightangle;
    
    angle= fmod( angle, M_PI );
    if( fabs(angle) < rightangle_tolerance )
        return;
    //  If angle is so close to a multiple of a right angle that a constant array
    //  might be modified spuriously through accumulation of numerical errors, we
    //  prefer treating the angle as a multiple of pi/2 and just assign the value
    //  of the closest unrotated square to it.
    rightangle=  (fabs(angle - M_PI_2) < rightangle_tolerance ||
                  fabs(angle + M_PI_2) < rightangle_tolerance ||
                  fabs(angle - M_PI) < rightangle_tolerance ||
                  fabs(angle + M_PI) < rightangle_tolerance);
    
    //  Buffer sine and cosine values for efficiency:
    thecos= cos(angle);
    thesin= sin(angle);
    
    //  Slopes of sides of rotated square: upslope is in (0, inf], downslope in
    //  [-inf, 0].  It seems that tan(pi/2) is not actually infinite at all, just
    //  very large (of order 10^16).  To be consistent with the ordering of
    //  the corners in corn0_x/y[] (see below), upslope must be positive.  To
    //  achieve this, the fabs is necessary due to errors in M_PI_2 or the
    //  calculation of the tangent.
    if( tan(angle)> 0.0 )
        upslope= fabs(tan(angle));
    else
        upslope= fabs(tan(angle - M_PI_2));
    downslope= -1.0/upslope;
    
    //  Compute the position of corners of a unit square rotated w. r. t. its
    //  centre.  The corners are sorted in ascending x direction.  (Strictly
    //  speaking all the x and y coordinates only take the values of slantsmall
    //  and slantlarge multiplied with sqrt(0.5) with different signs.  We store
    //  them as coordinates anyway to make the code more intelligible.)
    slantlarge= fabs(cos(angle + M_PI_4));
    slantsmall= fabs(sin(angle + M_PI_4));
    if( slantsmall > slantlarge ) {
        //  Swap large and small
        slantlarge += slantsmall;
        slantsmall= slantlarge - slantsmall;
        slantlarge -= slantsmall;
    }
    corn0_x[0]= - M_SQRT1_2 * slantlarge;
    corn0_x[1]= - M_SQRT1_2 * slantsmall;
    corn0_x[2]= M_SQRT1_2 * slantsmall;
    corn0_x[3]= M_SQRT1_2 * slantlarge;
    if( sin(4.0*angle) > 0.0 || (sin(4.0*angle)==0.0 && cos(4.0*angle> 0.0)) ) {
        corn0_y[0]= corn0_x[2];
        corn0_y[1]= corn0_x[0];
    }
    else {
        corn0_y[0]= corn0_x[1];
        corn0_y[1]= corn0_x[3];
    }
    corn0_y[2]= - corn0_y[1];
    corn0_y[3]= - corn0_y[0];
    
    boundavg= 0.0;
    for( x= 0; x < SizeX(); ++x )
        boundavg += (*this)(x, 0) + (*this)(x, SizeY()-1);
    for( y= 1; y < SizeY()-1; ++y )
        boundavg += (*this)(0, y) + (*this)(SizeX()-1, y);
    boundavg /= 2.0*(SizeY() + SizeX() - 2);
    
    dest.Create( SizeX(), SizeY(), Delta() );
    //pS = new TFktScal( SizeX(), SizeY(), Delta());
    /*for( x= 0; x < SizeX(); ++x )
     for( y= 0; y < SizeY(); ++y ) {
     //  Coordinates of centre of rotated square in old coordinate system:
     off_x= centre_x + (x-centre_x)*thecos + (y-centre_y)*thesin;
     off_y= centre_y - (x-centre_x)*thesin + (y-centre_y)*thecos;
     //  Closest point in old grid:
     old_x= (int)round(off_x);
     old_y= (int)round(off_y);
     if( 1 //rightangle
     ) {
     if( old_x< 0 || old_x>=SizeX() || old_y< 0 || old_y>=SizeY() )
     dest(x, y)= boundavg;
     else
     dest(x, y)= (*this)(old_x, old_y);
     continue;
     }
     if( old_x < -1 || old_x > SizeX() || old_y < -1 || old_y > SizeY() ) {
     dest(x, y)= boundavg;
     continue;
     }
     } */
    /*for (int i=0; i<2; i++) {
     pS->Smooth(dest);
     dest.Smooth(*pS);
     }*/
    for( x= 0; x < SizeX(); ++x )
        for( y= 0; y < SizeY(); ++y ) {
            //  Coordinates of centre of rotated square in old coordinate system:
            off_x= centre_x + (x-centre_x)*thecos + (y-centre_y)*thesin;
            off_y= centre_y - (x-centre_x)*thesin + (y-centre_y)*thecos;
            //  Closest point in old grid:
            old_x= (int)round(off_x);
            old_y= (int)round(off_y);
            if( rightangle || !rot) {
                if( old_x< 0 || old_x>=SizeX() || old_y< 0 || old_y>=SizeY() )
                    dest(x, y)= boundavg;
                else
                    dest(x, y)= (*this)(old_x, old_y);
                continue;
            }
            if( old_x < -1 || old_x > SizeX() || old_y < -1 || old_y > SizeY() ) {
                dest(x, y)= boundavg;
                continue;
            }
            //  Remaining offset (modulus < 0.5) of centre of new square from centre
            //  of old square in old coordinate system:
            off_x -= old_x;
            off_y -= old_y;
            //  Corners of rotated square relative to the centre of the closest
            //  unrotated square, in the unrotated coordinate system:
            for( i= 0; i< 4; ++i ) {
                corner_x[i]= corn0_x[i] + off_x;
                corner_y[i]= corn0_y[i] + off_y;
            }
            //  Now the overlap between the old and new squares is divided up into
            //  slices.  The corresponding trapezes are then added up.  limits[]
            //  stores the x values delimiting the slices.  The slices are defined in
            //  a way that they overlap only one old square in x direction, that they
            //  are delimited in y direction by linear functions of x (not just
            //  piecewise linear) and that the number of overlapping squares in y
            //  direction (1 to 3) is a constant for each slice.  Therefore the x
            //  values -0.5 and 0.5, the x coordinates of the corners of the new
            //  square, and the x coordinates of intercection points between the
            //  lines delimiting the new and old squares, respectively, are all
            //  limits betwen slices.
            lim_ind= 0;
            limits[lim_ind++]= -0.5;
            limits[lim_ind++]= 0.5;
            for( i= 0; i< 4; ++i )
                limits[lim_ind++]= corner_x[i];
            limits[lim_ind++]= corner_x[0] + (0.5-corner_y[0])/downslope;
            limits[lim_ind++]= corner_x[0] + (0.5-corner_y[0])/upslope;
            limits[lim_ind++]= corner_x[0] + (-0.5-corner_y[0])/upslope;
            limits[lim_ind++]= corner_x[0] + (-0.5-corner_y[0])/downslope;
            limits[lim_ind++]= corner_x[3] + (0.5-corner_y[3])/downslope;
            limits[lim_ind++]= corner_x[3] + (0.5-corner_y[3])/upslope;
            limits[lim_ind++]= corner_x[3] + (-0.5-corner_y[3])/upslope;
            limits[lim_ind++]= corner_x[3] + (-0.5-corner_y[3])/downslope;
            limits[lim_ind]= limits[lim_ind++]= 100.0;
            //  Sort in ascending order:
            qsort( limits, 14, sizeof(double), cmp_dbl );
            for( lim_ind= 0; limits[lim_ind] < corner_x[0]; ++lim_ind );
            leftupper= leftlower= corner_y[0];
            if( limits[lim_ind] < -0.5 )
                src_x= old_x - 1;
            else
                src_x= old_x;
            dest(x, y)= 0.0;
            while( limits[lim_ind] < corner_x[3] )
            {
                ++lim_ind;
                while( limits[lim_ind]==limits[lim_ind-1] )
                    ++lim_ind;
                rightupper= corner_y[0] + upslope * (limits[lim_ind] - corner_x[0]);
                if( rightupper > corner_y[3] + downslope * (limits[lim_ind] - corner_x[3]) )
                    rightupper= corner_y[3] + downslope * (limits[lim_ind] - corner_x[3]);
                rightlower= corner_y[0] + downslope * (limits[lim_ind] - corner_x[0]);
                if( rightlower < corner_y[3] + upslope * (limits[lim_ind] - corner_x[3]) )
                    rightlower= corner_y[3] + upslope * (limits[lim_ind] - corner_x[3]);
                meanupper= 0.5 * (leftupper + rightupper);
                meanlower= 0.5 * (leftlower + rightlower);
                contrib= 0.0;
                outofrange= 0.0;
                if( src_x < 0 || src_x >= SizeX() )
                    outofrange= meanupper - meanlower;
                else if( meanupper > 0.5 )
                    if( meanlower > 0.5 ) {
                        if( old_y < SizeY()-1 )
                            contrib += (meanupper - meanlower) * (*this)(src_x, old_y+1);
                        else
                            outofrange += meanupper - meanlower;
                    }
                    else if( meanlower > -0.5 ) {
                        if( old_y < SizeY()-1 )
                            contrib += (meanupper - 0.5) * (*this)(src_x, old_y+1);
                        else
                            outofrange += meanupper - 0.5;
                        if( old_y >= 0 && old_y < SizeY() )
                            contrib += (0.5 - meanlower) * (*this)(src_x, old_y);
                        else
                            outofrange += 0.5 - meanlower;
                    }
                    else {
                        if( old_y < SizeY()-1 )
                            contrib += (meanupper - 0.5) * (*this)(src_x, old_y+1);
                        else
                            outofrange += meanupper - 0.5;
                        if( old_y >= 0 && old_y < SizeY() )
                            contrib += (*this)(src_x, old_y);
                        else
                            outofrange += 1.0;
                        if( old_y > 0 )
                            contrib += (-0.5 - meanlower) * (*this)(src_x, old_y-1);
                        else
                            outofrange += -0.5 - meanlower;
                    }
                    else if( meanupper > -0.5 )
                        if( meanlower > -0.5 ) {
                            if( old_y >= 0 && old_y < SizeY() )
                                contrib += (meanupper - meanlower) * (*this)(src_x, old_y);
                            else
                                outofrange += meanupper - meanlower;
                        }
                        else {
                            if( old_y >= 0 && old_y < SizeY() )
                                contrib += (meanupper + 0.5) * (*this)(src_x, old_y);
                            else
                                outofrange += meanupper + 0.5;
                            if( old_y > 0 )
                                contrib += (-0.5 - meanlower) * (*this)(src_x, old_y-1);
                            else
                                outofrange += -0.5 - meanlower;
                        }
                        else {
                            if( old_y > 0 )
                                contrib += (meanupper - meanlower) * (*this)(src_x, old_y-1);
                            else
                                outofrange += meanupper - meanlower;
                        }
                dest(x, y) += (contrib + outofrange * boundavg) * 
                (limits[lim_ind] - limits[lim_ind-1]);
                leftupper= rightupper;
                leftlower= rightlower;
                if( (limits[lim_ind]==-0.5 && limits[lim_ind-1]!=-0.5) 
                   || (limits[lim_ind]==0.5 && limits[lim_ind-1]!=0.5) )
                    ++src_x;
            }
        }
    
    *this= dest;
}


/*! Compares two double values pointed to by ]a a and \a b.  Returns -1 if the
 second argument is larger, 1 if the first is larger and 0 if they are equal.
 Used by rotate as argument to qsort().  */

int PTG_Func2dScalar::cmp_dbl( const void *a, const void *b )
{
    if( *(double*)a > *(double*)b )
        return 1;
    else if( *(double*)a < *(double*)b )
        return -1;
    else
        return 0;
}

