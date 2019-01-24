/******************************************************************************
 $Id: PTG_Func2dScalar.h,v 1.6 2004/09/30 10:05:05 schatz Exp $
 ******************************************************************************/

#ifndef __PTG_FUNC2DSCALAR_H__
#define __PTG_FUNC2DSCALAR_H__

#include <float.h>
#include "func.h"
#include "PTG_Func2d.h"

class PTG_Func2dScalar : public PTG_Func2d<double>
{
public:
    PTG_Func2dScalar() {}
    PTG_Func2dScalar(int x, int y, double dDelta) : PTG_Func2d<double>(x,y, dDelta) {}
    
    //! Construct and initialise all elements to a constant value.
    PTG_Func2dScalar(int x, int y, double dDelta, const double Value)
    : PTG_Func2d<double>(x,y, dDelta, Value) {}
    
    void ShiftOne(int plusminus);
    
    void DxRight(const PTG_Func2dScalar& s);
    void Smooth(const PTG_Func2dScalar& s);
    double Integrate(int xmin) const;
    double GetMax() const;
    double GetFirstMax() const;
    double GetMin() const;
    double GetFirstMin() const;
    double CenterOfMassX() const;
    void DivGrad(const PTG_Func2dScalar& s);
    
    void copyscale( const PTG_Func2dScalar source );
    void rotate( double angle, int centre_x, int centre_y, bool rot );
    
private:
    static int cmp_dbl( const void *a, const void *b );
};


#endif

