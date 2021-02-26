/******************************************************************************
  $Id: func.cc,v 1.2 2004/09/03 11:43:57 schatz Exp $

Definition of methods of discretised function classes.
******************************************************************************/



#include "func.h"


// **** virtual builder functions ****
void CBoundary::Bound(TFktScal& f){}
void CBoundary::Bound(TFktVec& f){}
void CBoundary::Bound(TFktScal& f, const double in){}
void CBoundary::Bound(TFktVec& f, const vec2 in){}
void CBoundary::Bound(TFktScal& f, const double in,  const double out){}
void CBoundary::Bound(TFktVec& f, const vec2 in,  const vec2 out){}
void CBoundary::Latteral(TFktScal& f){}
void CBoundary::Latteral(TFktVec& f){}



// **** Closed boundary conditions ****

void CBoundaryClosed::Bound(TFktScal& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(1,y);
    f(iNx-1, y) = f(iNx-2, y);
  }
  Latteral(f);
}

void CBoundaryClosed::Bound(TFktScal& f, const double in)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = in;
    f(iNx-1, y) = f(iNx-2, y);
  }
      
  Latteral(f);
}

void CBoundaryClosed::Bound(TFktScal& f, const double in,  const double out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
      for (int y=0; y < iNy; y++) {
	f(0,y) = in;
	f(iNx-1, y) = out;
      }
      
      Latteral(f);
}

void CBoundaryClosed::Bound(TFktVec& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
      
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
	f(0,y) = f(1,y);
	f(iNx-1, y) = f(iNx-2, y);
  }
  
  Latteral(f);
}

void CBoundaryClosed::Bound(TFktVec& f, const vec2 in)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = in;
    f(iNx-1, y) = f(iNx-2, y);
  }
  
  Latteral(f);
}

void CBoundaryClosed::Bound(TFktVec& f, const vec2 in,  const vec2 out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
      
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = in;
    f(iNx-1, y) = out;
  }
  
  Latteral(f);
}
  
void CBoundaryClosed::Latteral(TFktScal& f)
{
  const int iNx = f.SizeX();
      const int iNy = f.SizeY();
      
      // left / right boundary (closed)
      
      for (int x=0; x < iNx; x++) {
	f(x,0) = f(x,1);
	f(x,iNy-1) = f(x,iNy-2);
      }
}

void CBoundaryClosed::Latteral(TFktVec& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
      // left / right boundary (closed)
  
  for (int x=0; x < iNx; x++) {
	f(x,0) = f(x,1);
	f(x,iNy-1) = f(x,iNy-2);
  }
}

// **** x-periodic boundary ****

void CBoundaryX::Bound(TFktScal& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
      
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1,y);
  }
      
  Latteral(f);
}
  
void CBoundaryX::Bound(TFktScal& f, const double in)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
      
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1,y);
  }
      
  Latteral(f);
}

void CBoundaryX::Bound(TFktScal& f, const double in,  const double out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
      // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1,y);
  }
  
  Latteral(f);
}

void CBoundaryX::Bound(TFktVec& f)
{
      const int iNx = f.SizeX();
      const int iNy = f.SizeY();
      
      // in / out flow boundary
      
      for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1,y);
      }
      
      Latteral(f);
}

void CBoundaryX::Bound(TFktVec& f, const vec2 in)
    {
      const int iNx = f.SizeX();
      const int iNy = f.SizeY();
      
      // in / out flow boundary
      
      for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1,y);
      }
      
      Latteral(f);
    }

void CBoundaryX::Bound(TFktVec& f, const vec2 in,  const vec2 out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1,y);
  }
  
  Latteral(f);
}

void CBoundaryX::Latteral(TFktScal& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
      // left / right boundary (closed)
  
  for (int x=0; x < iNx; x++) {
    f(x,0) = f(x,1);
    f(x,iNy-1) = f(x,iNy-2);
  }
}

void CBoundaryX::Latteral(TFktVec& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // left / right boundary (closed)
      
  for (int x=0; x < iNx; x++) {
    f(x,0) = f(x,1);
    f(x,iNy-1) = f(x,iNy-2);
  }
}

// **** x- and y-periodic boundary ****
void CBoundaryXY::Bound(TFktScal& f)
    {
      const int iNx = f.SizeX();
      const int iNy = f.SizeY();
      
      // in / out flow boundary
      
      for (int y=0; y < iNy; y++) {
	f(0,y) = f(iNx-2, y);
	f(iNx-1,y) = f(1, y);
      }
      
      Latteral(f);
    }

void CBoundaryXY::Bound(TFktScal& f, const double in)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
      
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1, y);
  }
  
  Latteral(f);
}

void CBoundaryXY::Bound(TFktScal& f, const double in,  const double out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1, y);
  }
  
  Latteral(f);
}

void CBoundaryXY::Bound(TFktVec& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1, y);
  }
  
  Latteral(f);
}

void CBoundaryXY::Bound(TFktVec& f, const vec2 in)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1, y);
  }
      
  Latteral(f);
}

void CBoundaryXY::Bound(TFktVec& f, const vec2 in,  const vec2 out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
      
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(iNx-2, y);
    f(iNx-1,y) = f(1, y);
  }
  
  Latteral(f);
}
  
void CBoundaryXY::Latteral(TFktScal& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  for (int x=0; x < iNx; x++) {
    f(x,0) = f(x,iNy-2);
    f(x,iNy-1) = f(x,1);
  }
      
}

void CBoundaryXY::Latteral(TFktVec& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();

  for (int x=0; x < iNx; x++) {
    f(x,0) = f(x,iNy-2);
    f(x,iNy-1) = f(x,1);
  }
  
}


// **** y-periodic boundary ****

void CBoundaryY::Bound(TFktScal& f)
{
  const int iNx = f.SizeX();
      const int iNy = f.SizeY();
      
      // in / out flow boundary
      
      for (int y=0; y < iNy; y++) {
	f(0,y) = f(1,y);
	f(iNx-1, y) = f(iNx-2, y);
      }
      
      Latteral(f);
}

void CBoundaryY::Bound(TFktScal& f, const double in)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = in;
    f(iNx-1, y) = f(iNx-2, y);
  }
  
  Latteral(f);
}

void CBoundaryY::Bound(TFktScal& f, const double in,  const double out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = in;
    f(iNx-1, y) = out;
  }
  
  Latteral(f);
}

void CBoundaryY::Bound(TFktVec& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = f(1,y);
    f(iNx-1, y) = f(iNx-2, y);
  }
  
  Latteral(f);
}

void CBoundaryY::Bound(TFktVec& f, const vec2 in)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = in;
    f(iNx-1, y) = f(iNx-2, y);
      }
  
  Latteral(f);
}

void CBoundaryY::Bound(TFktVec& f, const vec2 in,  const vec2 out)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  // in / out flow boundary
  
  for (int y=0; y < iNy; y++) {
    f(0,y) = in;
    f(iNx-1, y) = out;
  }
  
  Latteral(f);
}

void CBoundaryY::Latteral(TFktScal& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  for (int x=0; x < iNx; x++) {
	f(x,0) = f(x,iNy-2);
	f(x,iNy-1) = f(x,1);
  }
  
}

void CBoundaryY::Latteral(TFktVec& f)
{
  const int iNx = f.SizeX();
  const int iNy = f.SizeY();
  
  for (int x=0; x < iNx; x++) {
    f(x,0) = f(x,iNy-2);
    f(x,iNy-1) = f(x,1);
  }
  
}

