#include <stdlib.h>

#include "vec.h"
#include "initsurfalea.h"

////////////////////////////////////////
// Plain initial surface with fluctuations
//

// THIS IS NOW ONLY FOR VEGETATION

CInitSurfAlea::CInitSurfAlea(const dunepar& par, string prefix)
{   
    m_xmin = par.getdefault("veget.xmin", 0);
    m_rho_max = par.getdefault("veget.rho.max", 1.0);
    m_nodes[0] = par.getdefault(prefix+"alea.nodes.b", 1); 
    m_nodes[1] = par.getdefault(prefix+"alea.nodes.m", 0); 
}

void CInitSurfAlea::init_2d_scal(TFktScal& array)
{
    const int iNX = array.SizeX();
    const int iNY = array.SizeY();
    
    // INIT MANUAL
    if (iNY < 8) {
//        if (iNY < 16) {
        // 2D
        int x = 0.5 * iNX;
        for( int y= 0; y< iNY; ++y ){
            array(x,y) = 0.9 * m_rho_max;
        }
    } else {
        // 3D : random seeding
        for (int i = 0; i < m_nodes[0]; i++) {
            int x = rand() % (iNX - (m_xmin + 10)) + m_xmin;
            int y = rand() % iNY;
            array(x,y) = 0.9 * m_rho_max;
        }
    }    
}

void CInitSurfAlea::init_2d_vec(TFktVec& array)
{
    const int iNX = array.SizeX();
    const int iNY = array.SizeY();
    
    // INIT MANUAL
    if (iNY < 8) {
        // 2D
        int x = 0.5 * iNX;
        int y = 0.5 * iNY;
        for (int k = 0; k < 2; k++) {
            //for( int y= 0; y< iNY; ++y ){
                array(x,y)[k] = (m_nodes[k] > 0 ? 0.9 : 0) * m_rho_max;
            //}
        }
        
    } else {
        // 3D : random seeding
        for (int k = 0; k < 2; k++) {
            for (int i = 0; i < m_nodes[k]; i++) {
                int x = rand() % (iNX - (m_xmin + 10)) + m_xmin;
                int y = rand() % iNY;
                if (k==0) {
                    //x = iNX-10;
                }
                array(x,y)[k] = 0.5 * m_rho_max;

                //cout << x << ' ' << y << ' ' << array(x,y)[k] << "#33" << endl;
            }
        }
    }    
}


