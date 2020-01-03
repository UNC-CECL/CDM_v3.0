/******************************************************************************
  $Id: func.h,v 1.3 2004/09/30 10:05:07 schatz Exp $

Discretised function classes.
******************************************************************************/

#ifndef FUNC_H
#define FUNC_H

#include "PTG_Func.h"

#include "vec.h"
#include "PTG_Func2dScalar.h"
#include "PTG_Func2dVec.h"

// adaption of the 1D function to 2D terminology
class CFunc1d : public PTG_Func<double>
{
public:
  CFunc1d() {}

  CFunc1d(int x, int /*y*/, double dDelta) :
    PTG_Func<double>(x, dDelta)
  {}

  void Create(int x, int /*y*/, double dDelta, const double dValue) {
    CreateAndInit(dValue, x, dDelta  );
  }
  void Create(int x, int /*y*/, double dDelta, const vec2& dValue) {
    CreateAndInit(dValue[0], x, dDelta);
  }

  void SetAll(const double dValue) {
    SetValue(dValue);
  }
  void SetAll(const vec2& Value) {
    SetValue(Value[0]);
  }

  void DivUpWind(const CFunc1d& f, const CFunc1d& u) {
    DUpWind(f,u);
  }

  void GradUpWind(const CFunc1d& f, const CFunc1d& u) {
    DUpWind(f,u);
  }
  void GradMid(const CFunc1d& f) {
    DMid(f);
  }
};


#ifdef PDE_2D
  typedef CFunc1d TFkt;
  typedef CFunc1d TFktScal;
  typedef CFunc1d TFktVec;

  typedef double TVec;
#else
  typedef PTG_Func2dScalar TFktScal;
  typedef PTG_Func2dVec    TFktVec;

  typedef vec2 TVec;
#endif

// ***********************************
// **** Boundary classes *************
// ***********************************



class CBoundary 
{

  public :
    
    CBoundary(){
  };
  
  virtual void Bound(TFktScal& f);
  virtual void Bound(TFktVec& f);
  
  virtual  void Bound(TFktScal& f, const double in);
  virtual  void Bound(TFktVec& f, const vec2 in);
  
  virtual  void Bound(TFktScal& f, const double in,  const double out);
  virtual  void Bound(TFktVec& f, const vec2 in,  const vec2 out);
  
  virtual  void Latteral(TFktScal& f);
  virtual  void Latteral(TFktVec& f);
    

};


// **** Closed boundary conditions ****

class CBoundaryClosed : public CBoundary
{ 
 public:
  
  void Bound(TFktScal& f);
  void Bound(TFktScal& f, const double in);
  void Bound(TFktScal& f, const double in,  const double out);
  void Bound(TFktVec& f);
  void Bound(TFktVec& f, const vec2 in);
  void Bound(TFktVec& f, const vec2 in,  const vec2 out);
  void Latteral(TFktScal& f);
  void Latteral(TFktVec& f);
};

class CBoundaryX : public CBoundary
{ 
 public:
  
  void Bound(TFktScal& f);
  void Bound(TFktScal& f, const double in);
  void Bound(TFktScal& f, const double in,  const double out);
  void Bound(TFktVec& f);
  void Bound(TFktVec& f, const vec2 in);
  void Bound(TFktVec& f, const vec2 in,  const vec2 out);
  void Latteral(TFktScal& f);
  void Latteral(TFktVec& f);
};

class CBoundaryY : public CBoundary
{ 
 public:
  
  void Bound(TFktScal& f);
  void Bound(TFktScal& f, const double in);
  void Bound(TFktScal& f, const double in,  const double out);
  void Bound(TFktVec& f);
  void Bound(TFktVec& f, const vec2 in);
  void Bound(TFktVec& f, const vec2 in,  const vec2 out);
  void Latteral(TFktScal& f);
  void Latteral(TFktVec& f);
};

class CBoundaryXY : public CBoundary
{ 
 public:
  
  void Bound(TFktScal& f);
  void Bound(TFktScal& f, const double in);
  void Bound(TFktScal& f, const double in,  const double out);
  void Bound(TFktVec& f);
  void Bound(TFktVec& f, const vec2 in);
  void Bound(TFktVec& f, const vec2 in,  const vec2 out);
  void Latteral(TFktScal& f);
  void Latteral(TFktVec& f);
};


#endif
