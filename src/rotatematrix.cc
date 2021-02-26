#include <math.h>

#include "rotatematrix.h"
#include "func.h"

/////////////////////////////////////////////
// Surface from MATLAB file
//

// construccion
CRotateMatrix::CRotateMatrix()  {}

void CRotateMatrix::DoRotation(TFktScal& f, double RotAngle, bool EqualVol)
{

  int Nx=f.SizeX();
  int Ny=f.SizeY();
  TFktScal tempmatrix(f);
  int CenterX = (int) floor(Nx/2.0);
  int CenterY = (int) floor(Ny/2.0);
  double X(0),Y(0);
  int x,y,x1,y1,x2,y2;

  // Calculate volumen of matrix
  double VolumenBefore = f.Integrate(0);
  
  // Rotation around center of the matrix, clockwise
  // values are interpoled within first order
  for(y=0; y < Ny; y++) {
    for(x=0; x < Nx; x++) {
      X = (x-CenterX)*cos(RotAngle)-(y-CenterY)*sin(RotAngle)+CenterX    ;
      Y = (x-CenterX)*sin(RotAngle)+(y-CenterY)*cos(RotAngle)+CenterY    ;
      
      if(X>(double) Nx-1) {
	if(Y> (double) Ny-1) {
	  f(x,y) = tempmatrix(Nx-1, Ny-1);
	} else if(Y<0) {
	  f(x,y) = tempmatrix(Nx-1, 0);
	} else {
	  f(x,y) = tempmatrix(Nx-1, (int) round(Y));
	}
      } else if(X<0) {
	if(Y> (double) Ny-1) {
	  f(x,y) = tempmatrix(0, Ny-1);
	} else if(Y<0) {
	  f(x,y) = tempmatrix(0, 0);
	} else {
	  f(x,y) = tempmatrix(0, (int) round(Y));
	}
      } else if(Y<0 || Y > (double)  Ny-1) {
	if(Y> (double) Ny-1) {
	  f(x,y) = tempmatrix((int) round(X), Ny-1);
	} else {
	  f(x,y) = tempmatrix((int) round(X), 0);
	}
      } else { 

	x1 = (int) floor(X);
	y1 = (int) floor(Y);
	x2 = (int) ceil(X);
	y2 = (int)ceil(Y);


	if(x1!=x2 || y1!=y2) {
	  double d1 = sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1));
	  double d2 = sqrt((X-x2)*(X-x2)+(Y-y1)*(Y-y1));
	  double d3 = sqrt((X-x1)*(X-x1)+(Y-y2)*(Y-y2));
	  double d4 = sqrt((X-x2)*(X-x2)+(Y-y2)*(Y-y2));

	  f(x,y) = d1*tempmatrix(x1,y1)+d2*tempmatrix(x2,y1)+d3*tempmatrix(x1,y2)+d4*tempmatrix(x2,y2);
	  f(x,y) /= d1+d2+d3+d4;
	} else {
	  f(x,y) = tempmatrix(x1,y1);
	}
	if (f(x,y) < 0.0) f(x,y) = 0.0; 
      }
    }
  }
  // calculate volumen before and adjusting new volumen to old volumen
  double VolumenAfter = f.Integrate(0);
  double VolFactor = VolumenBefore / VolumenAfter;
  
  if (EqualVol && VolFactor > 0.0) {
    for(y=0; y < Ny; y++) {
      for(x=0; x < Nx; x++) {
	f(x,y) *= VolFactor;
      }
    }
  }

}

void CRotateMatrix::DoRotation(TFktVec& f, double RotAngle)
{
  int Nx=f.SizeX();
  int Ny=f.SizeY();
  int Size=2;
  TFktVec tempmatrix(f);
  int CenterX = (int) floor(Nx/2.0);
  int CenterY = (int) floor(Ny/2.0);
  double X(0),Y(0);
  int s,x,y,x1,y1,x2,y2;
  
  // Rotation around center of the matrix, clockwise
  // values are interpoled within first order
  for(s=0; s < Size; s++) {
    for(y=0; y < Ny; y++) {
      for(x=0; x < Nx; x++) {
	X = (x-CenterX)*cos(RotAngle)-(y-CenterY)*sin(RotAngle)+CenterX;
	Y = (x-CenterX)*sin(RotAngle)+(y-CenterY)*cos(RotAngle)+CenterY;
      
	if(X>(double) Nx-1) {
	  if(Y> (double) Ny-1) {
	    f(x,y)[s] = tempmatrix(Nx-1, Ny-1)[s];
	  } else if(Y<0) {
	    f(x,y)[s] = tempmatrix(Nx-1, 0)[s];
	  } else {
	    f(x,y)[s] = tempmatrix(Nx-1, (int) round(Y))[s];
	  }
	} else if(X<0) {
	  if(Y> (double) Ny-1) {
	    f(x,y)[s] = tempmatrix(0, Ny-1)[s];
	  } else if(Y<0) {
	    f(x,y)[s] = tempmatrix(0, 0)[s];
	  } else {
	    f(x,y)[s] = tempmatrix(0, (int) round(Y))[s];
	  }
	} else if(Y<0 || Y > (double) Ny-1) {
	  if(Y> (double) Ny-1) {
	    f(x,y)[s] = tempmatrix((int) round(X), Ny-1)[s];
	  } else {
	    f(x,y)[s] = tempmatrix((int) round(X), 0)[s];
	  }
	} else { 
	  x1 = (int) floor(X);
	  y1 = (int) floor(Y);
	  x2 = (int) ceil(X);
	  y2 = (int) ceil(Y);
	  if(x1!=x2 && y1!=y2) {
	    f(x,y)[s] = tempmatrix(x1,y1)[s];
	  } else {
	    f(x,y)[s] = tempmatrix(x1,y1)[s];
	  }
	}
      }
    }
  }
}






