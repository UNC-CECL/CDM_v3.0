#include "initsurfbeach.h"
#include "globals.h"

////////////////////////////////////////
// class CInitSurfBeach
//

CInitSurfBeach::CInitSurfBeach(const dunepar& par, string prefix)
{
	m_h= par.getdefault<double>( prefix+"beach.h", duneglobals::dx());
    double m_slope= duneglobals::slope();//par.getdefault<double>( prefix+"beach.angle", 45);
    
    if (m_slope > 0) {
        m_l = m_h/(duneglobals::dx()*m_slope);
    }
}

void CInitSurfBeach::init_2d_scal(TFktScal& array)
{
	int x, y;
	
	for( x= 0; x< duneglobals::nx(); ++x )
		for( y= 0; y< duneglobals::ny(); ++y ){
           if (x < m_l) {
//                double a = 0*(x-0.5*m_l)*duneglobals::dx()/3.;
                array(x, y) = /*m_h * (1 - (1-x*1.0/m_l)*(1-x*1.0/m_l));*/ m_h / m_l * x; // + 0.2*exp(-a*a);
            } else {
                array(x, y) = m_h;
            }
 //           double a = (x-0.5*m_l)*duneglobals::dx()/3.;
 //           array(x, y) = m_h * (1 - exp(- x /m_l)) + 0.2*exp(-a*a);
        }
}

