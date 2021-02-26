#ifndef __PTG_FUNC2DVEC_H__
#define __PTG_FUNC2DVEC_H__

#include "vec.h"
#include "PTG_Func2d.h"

class PTG_Func2dScalar;

class PTG_Func2dVec : public PTG_Func2d<vec2>
{
public:
    PTG_Func2dVec() {}
    
    PTG_Func2dVec(int x, int y, double dDelta) : PTG_Func2d<vec2>(x,y, dDelta) {}
    
    // initialize each element with a value
    PTG_Func2dVec(int x, int y, double dDelta, const vec2 & Value)
    : PTG_Func2d<vec2>(x,y, dDelta, Value) {}
    
    void ShiftOne(int plusminus);
    void rescale(double factor);
    
    void DiscreteLaplacian(const PTG_Func2dScalar& s);
    void GradMid(const PTG_Func2dScalar& s);
    void GradMin(const PTG_Func2dScalar& s);
    void GradUpWind(const PTG_Func2dScalar& s, const PTG_Func2dVec& u);
    double GetMaxAbs();
    
    void copyscale( const PTG_Func2dVec source );
    void get_diff_stat( double *avg_x, double *avg_y,
                       double *sigma_x, double *sigma_y );
};



#endif
