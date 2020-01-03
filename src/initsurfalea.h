/******************************************************************************
  $Id: initsurfalea.h,v 1.3 2004/08/18 10:58:10 schatz Exp $
******************************************************************************/

#ifndef __INITSURFALEA_H__
#define __INITSURFALEA_H__

#include "globals.h"
#include "initsurf.h"
#include "func.h"
#include "vec.h"


// ----  plain initial surface with aleatory fluctuations ----

class CInitSurfAlea : public arrayinit
{
public:
    CInitSurfAlea(const dunepar& P, string prefix="");
    ~CInitSurfAlea(){}

    virtual void init_2d_scal(TFktScal& array);
    virtual void init_2d_vec(TFktVec& array);

private:
    double m_rho_max;
  int m_xmin;
  int m_nodes[2];
};

#endif
