#ifndef __ROTATEMATRIX_H__
#define __ROTATEMATRIX_H__

#include "globals.h"
#include "func.h"

class CRotateMatrix 
{
public:
  CRotateMatrix();

  ~CRotateMatrix() {}

  void DoRotation(TFktScal& f, double RotAngle, bool EqualVol);
  void DoRotation(TFktVec& f, double RotAngle);

};


#endif


