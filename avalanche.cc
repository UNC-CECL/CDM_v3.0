/******************************************************************************
 $Id: avalanche.cc,v 1.9 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#include <math.h>

#include "avalanche.h"
#include "func.h"

//*****************************************************************************
//  class avalanche

/*!  Gets static and dynamic angle of repose from duneglobals object and
 computes their tangent and the maximal height difference of adjacent points
 for use by derived classes.  */

avalanche::avalanche()
{
    m_angle_stat= duneglobals::repose_stat();
    m_tan_stat= tan(m_angle_stat * M_PI / 180.0);
    m_dh_stat= m_tan_stat * duneglobals::dx();
    m_angle_dyn= duneglobals::repose_dyn();
    m_tan_dyn= tan(m_angle_dyn * M_PI / 180.0);
    m_dh_dyn= m_tan_dyn * duneglobals::dx();
}


/*!  Returns an object of a derived class of avalanche, depending on the
 parameter "avalanche" in the parameter file.  Current valid values are legacy
 (the default, which gives a CIterAval object), legacynew and cell (cellaval).  */

avalanche *avalanche::create(const dunepar &p)
{
    avalanche *object;
    string avaltype;
    
    avaltype= p.getdefault<string>("avalanche", "flow");
    if( avaltype == "flow" )
        object= new flowaval(p);
    else if( avaltype == "cell" )
        object= new cellaval(p);
    else if( avaltype == "dummy" )
        object= new dummyaval();
    else {
        cerr << "avalanche::create: ERROR: illegal value `" << avaltype << "' for "
        << "parameter avalanche.  Valid values are `flow', "
        "`cell' or `dummy'.\n";
        exit(1);
    }
    return object;
}

//*****************************************************************************
//  class flowaval

flowaval::flowaval(const dunepar& p) :
m_h(0),
m_h_nonerod(0),
m_grad_h_down(duneglobals::nx(), duneglobals::ny(), duneglobals::dx()),
m_flux_down(duneglobals::nx(), duneglobals::ny(), duneglobals::dx())
{
    m_n_iter= p.getdefault("aval.new.maxiter", 20);
    m_x_periodic= duneglobals::periodic_x();
    m_y_periodic= duneglobals::periodic_y();
    
    // this is the only parameter of the model and control how fast the slope is relaxed
    m_E = p.getdefault("aval.new.relax", 0.9) * duneglobals::dx();
}

void flowaval::calc( TFktScal &h , TFktScal &h_nonerod )
{
    m_h= &h;
    m_h_nonerod= &h_nonerod;
    
    double iter = 0;
    double max_slope = 0;
    
    max_slope = CalcGradDown();
    
    if(max_slope > m_tan_stat){
        while( iter < m_n_iter && max_slope > m_tan_dyn){
            max_slope = Step(max_slope);
            iter ++;
        }
        if (max_slope > m_tan_stat)
            cout << "CIterAvalNew::calc: Solution not converged after " << iter << " iterations, max. slope = " << max_slope << " (target: " << m_tan_dyn << ")\n";
        else  cout << "CIterAvalNew::calc: " << iter << " iterations, final slope " << max_slope << "\n";
    }
    else cout << "CIterAvalNew::calc: no avalanche: maximal slope " << max_slope << "\n";
}

double flowaval::CalcGradDown(){
    
    double grad_h2, max_grad= 0;
    int prevy, nexty;
    
    for (int y=0; y < duneglobals::ny(); y++) {
        for (int x=0; x < duneglobals::nx(); x++)
        {
            m_grad_h_down(x,y)[0] = m_grad_h_down(x,y)[1] = 0;
            
            if(x==0 && m_x_periodic){
                if((*m_h)(0,y) < (*m_h)(1,y) && (*m_h)(0,y) < (*m_h)(duneglobals::nx()-1,y))
                    m_grad_h_down(0,y)[0] = 0;
                else if( (*m_h)(1,y) > (*m_h)(duneglobals::nx()-1,y) )
                    m_grad_h_down(0, y)[0] = -((*m_h)(0, y) - (*m_h)(duneglobals::nx()-1,y));
                else
                    m_grad_h_down(0, y)[0] = (*m_h)(0, y) - (*m_h)(1,y);
            }else if(x==0 && (*m_h)(0,y) > (*m_h)(1,y)){
                m_grad_h_down(0, y)[0] = (*m_h)(0, y) - (*m_h)(1,y);
            }else if(x==duneglobals::nx()-1 && m_x_periodic){
                if((*m_h)(x,y) < (*m_h)(0,y) && (*m_h)(x,y) < (*m_h)(x-1,y))
                    m_grad_h_down(x,y)[0] = 0;
                else if( (*m_h)(0,y) > (*m_h)(x-1,y) )
                    m_grad_h_down(x, y)[0] = -((*m_h)(x, y) - (*m_h)(x-1,y));
                else
                    m_grad_h_down(x, y)[0] = (*m_h)(x, y) - (*m_h)(0,y);
            }else if(x==duneglobals::nx()-1 && (*m_h)(x,y) < (*m_h)(x-1,y)){
                m_grad_h_down(x, y)[0] = -((*m_h)(x, y) - (*m_h)(x-1,y));
            }else if(x > 0 && x < duneglobals::nx()-1)
                if((*m_h)(x,y) < (*m_h)(x+1,y) && (*m_h)(x,y) < (*m_h)(x-1,y))
                    m_grad_h_down(x,y)[0] = 0;
                else if( (*m_h)(x+1,y) > (*m_h)(x-1,y) )
                    m_grad_h_down(x, y)[0] = -((*m_h)(x, y) - (*m_h)(x-1,y));
                else
                    m_grad_h_down(x, y)[0] = (*m_h)(x, y) - (*m_h)(x+1,y);
            
            prevy = (y==0)?(duneglobals::ny()-1):(y-1);
            nexty = (y==duneglobals::ny()-1)?0:(y+1);
            
            
            if((*m_h)(x,y) < (*m_h)(x,nexty) && (*m_h)(x,y) < (*m_h)(x,prevy))
                m_grad_h_down(x,y)[1] = 0;
            else if( (*m_h)(x,nexty) > (*m_h)(x,prevy) )
                m_grad_h_down(x,y)[1] = -((*m_h)(x, y) - (*m_h)(x,prevy));
            else
                m_grad_h_down(x,y)[1] = (*m_h)(x, y) - (*m_h)(x,nexty);
            
            
            
            m_grad_h_down(x,y)[0]/= duneglobals::dx();
            m_grad_h_down(x,y)[1]/= duneglobals::dx();
            
            grad_h2= ((*m_h)(x,y) > (*m_h_nonerod)(x,y)+0.005 ? m_grad_h_down(x,y)[0] * m_grad_h_down(x,y)[0] +
                      m_grad_h_down(x,y)[1] * m_grad_h_down(x,y)[1] : 0);
            
            if(max_grad <= grad_h2)
                max_grad= grad_h2;
            
        }
    }
    return sqrt(max_grad);
    
}

double flowaval::Step(double max_slope)
{
    double grad_h, grad_h_nonerod, slope_diff, q_in, q_out;
    int nextx, prevx, nexty, prevy, flux_prevx, flux_nextx, flux_prevy, flux_nexty;
    
    // flux
    for (int y=0; y < duneglobals::ny(); y++) {
        for (int x=0; x < duneglobals::nx(); x++) {
            grad_h = m_grad_h_down(x,y)[0] * m_grad_h_down(x,y)[0] +
            m_grad_h_down(x,y)[1] * m_grad_h_down(x,y)[1];
            grad_h= sqrt(grad_h);
            grad_h_nonerod = ((*m_h)(x,y) - (*m_h_nonerod)(x,y))/duneglobals::dx();
            
            if( grad_h > m_tan_dyn && grad_h_nonerod > 0 ) {
                slope_diff = (grad_h_nonerod < grad_h - m_tan_dyn ? tanh(grad_h_nonerod) :
                              tanh(grad_h) - tanh(0.9*m_tan_dyn));
                m_flux_down(x,y)[0] = slope_diff * m_grad_h_down(x,y)[0]/grad_h;
                m_flux_down(x,y)[1] = slope_diff * m_grad_h_down(x,y)[1]/grad_h;
            }
            else {
                m_flux_down(x,y)[0]= 0.0;
                m_flux_down(x,y)[1]= 0.0;
            }
        }
    }
    
    // change in h
    
    
    for (int y=0; y < duneglobals::ny(); y++) {
        for (int x=0; x < duneglobals::nx(); x++) {
            q_out = fabs(m_flux_down(x,y)[0])+fabs(m_flux_down(x,y)[1]);
            
            prevx = (x==0 ? duneglobals::nx()-1 : x-1);
            nextx = (x==duneglobals::nx()-1 ? 0 : x+1);
            prevy = (y==0)?(duneglobals::ny()-1):(y-1);
            nexty = (y==duneglobals::ny()-1)?0:(y+1);
            
            
            flux_prevx = (x==0 && !m_x_periodic ? 0 : 1);
            flux_nextx = (x==duneglobals::nx()-1 && !m_x_periodic ? 0 : 1);
            flux_prevy = (y==0 && !m_y_periodic ? 0 : 1);
            flux_nexty = (y==duneglobals::ny()-1 && !m_y_periodic ? 0 : 1);
            
            q_in = flux_prevx * (m_flux_down(prevx,y)[0] > 0 ? m_flux_down(prevx,y)[0] : 0)-
            flux_nextx * (m_flux_down(nextx,y)[0] < 0 ? m_flux_down(nextx,y)[0] : 0)+
            flux_prevy * (m_flux_down(x,prevy)[1] > 0 ? m_flux_down(x,prevy)[1] : 0)-
            
            flux_nexty * (m_flux_down(x,nexty)[1] < 0 ? m_flux_down(x,nexty)[1] : 0);
            
            (*m_h)(x,y) += m_E * (q_in - q_out);
        }
    }
    
    // max. slope is used to stop the iteration
    return CalcGradDown();
}

//*****************************************************************************
//  class cellaval


/*!  Pointer to h field used by calc.  For cmp_coordpair, which has to be
 able to access it.  */
const TFktScal *cellaval::hpointer;


/*!  Allocates memory for m_done and m_sliporder and writes coordinate pairs in
 ascending order to m_sliporder.  */

cellaval::cellaval(const dunepar &p)
{
    int x, y, i;
    
    m_maxiter= p.getrequired<int>( "aval.cell.maxiter" );
    m_done= new bool[duneglobals::nx() * duneglobals::ny()];
    m_sliporder= new int[2 * duneglobals::nx() * duneglobals::ny()];
    i= 0;
    for( x= 0; x< duneglobals::nx(); ++x )
        for( y= 0; y< duneglobals::ny(); ++y ) {
            m_sliporder[i++]= x;
            m_sliporder[i++]= y;
        }
}


/*!  Frees the memory for m_done.  */

cellaval::~cellaval( )
{
    delete[] m_sliporder;
    delete[] m_done;
}


/*!  This function recursively relaxes avalanches on a sand surface.
 
 The grid is searched for sites whose difference to its nearest neighbours
 exceed the value corresponding to the static angle of repose.  If so, the
 amount in excess of the _dynamic_ angle of repose is distributed among its
 neighbours.  This is done recursively until all sites differ less from its four
 neighbours than the difference corresponding to the static angle of repose.
 
 The issue of how to distribute the sand from the centre site to its neighbours
 is non-trivial because an amount of sand transferred to one neighbour also
 lowers the surplus slope with a different neighbour.  Here, it is assumed that
 the velocity of the flowing sand is proportional to the height difference
 between neighbouring sites.  Taking into account the total amount of sand
 transferred away from a given site, one gets that the time until the angle of
 repose is reached is proportional to
 
 sand excess compared to one neighbour / (total sand excess + excess wrt the same neighbour),
 
 where the sand excess is the height difference minus the one corresponding to
 the dynamic angle of repose.  The sand transported during that time to this
 neighbour is then proportional to the same expression with the sand transfer in
 the numerator squared.  All that remains to be done now is to find the
 proportionality constant.  The difference in height after the sand slide must
 be such as not to exceed the dynamical angle of repose.  The difference between
 the centre site and a neighbour after the sand transfer is
 
 initial difference - transfer to this site - total transfer.
 
 Taking into account the dependence of the transfer on the old difference, one
 finds that this is largest for the site with the largest difference.  Since
 sand transport stops only after all differences are relaxed and the relaxation
 time is largest for the neighbour with the largest initial difference, the
 largest final difference must be the one corresponding to the dynamical angle
 of repose.  In other words, the final excess difference of the site with the
 largest initial difference is zero.  Knowing the ratios between the sand
 transfer and the total sand transfer and the initial height difference, we can
 compute the normalisation constant after the sand transfer.  */

void cellaval::calc( TFktScal &h, TFktScal &h_nonerod )
{
    double surplus[4], flux[4];
    double surplus_mx, surplus_px, surplus_my, surplus_py, total, maxsurplus;
    double totalflux, maxflux, fluxnorm;
    const double statdyndiff = m_dh_stat - m_dh_dyn;
    int iteration, coordind, maxcoordind, x, y, exceeding, neighbour;
    bool all_done;
    
    set_sliporder( h );
    
    for( x= 0; x< duneglobals::nx(); ++x )
        for( y= 0; y< duneglobals::ny(); ++y )
            done(x, y)= false;
    
    maxcoordind= 2*duneglobals::nx()*duneglobals::ny();
    
    for( all_done= false, iteration= 0; !all_done && iteration< m_maxiter; ++iteration )
    {
        //cout << "cellaval::calc: iteration " << iteration << endl;
        all_done= true;
        for( coordind= 0; coordind< maxcoordind; coordind += 2 )
        {
            x= m_sliporder[coordind];
            y= m_sliporder[coordind+1];
            
            if( done(x, y) )
                continue;
            exceeding= 0;
            total= 0.0;
            //	surplus_mx= surplus_px = surplus_my= surplus_py= 0.0;
            if( x && (surplus_mx= h(x, y) - h(x-1, y) - m_dh_stat) > 0.0 ) {
                surplus[exceeding]= surplus_mx + statdyndiff;
                total += surplus[exceeding];
                ++exceeding;
            }
            if( x< duneglobals::nx()-1 &&
               (surplus_px= h(x, y) - h(x+1, y) - m_dh_stat) > 0.0 ) {
                surplus[exceeding]= surplus_px + statdyndiff;
                total += surplus[exceeding];
                ++exceeding;
            }
            if( y && (surplus_my= h(x, y) - h(x, y-1) - m_dh_stat) > 0.0 ) {
                surplus[exceeding]= surplus_my + statdyndiff;
                total += surplus[exceeding];
                ++exceeding;
            }
            if( y< duneglobals::ny()-1 &&
               (surplus_py= h(x, y) - h(x, y+1) - m_dh_stat) > 0.0 ) {
                surplus[exceeding]= surplus_py + statdyndiff;
                total += surplus[exceeding];
                ++exceeding;
            }
            done(x, y)= true;
            if( exceeding )
            {
                all_done= false;
                totalflux= 0.0;
                maxsurplus= 0.0;
                maxflux= 0.0;
                for( neighbour= 0; neighbour< exceeding; ++neighbour ) {
                    flux[neighbour]= surplus[neighbour] * surplus[neighbour] /
                    (total + surplus[neighbour]);
                    totalflux += flux[neighbour];
                    if( surplus[neighbour] > maxsurplus ) {
                        maxsurplus= surplus[neighbour];
                        maxflux= flux[neighbour];
                    }
                }
                // Normalisation factor= real total flux / unnormalised total flux
                fluxnorm= maxsurplus / (1 + maxflux/totalflux) / totalflux;
                if( m_maxiter > 5 )
                    //  aval.cell.maxiter should be chosen large enough
                    fluxnorm /= 0.2 * m_maxiter;
                h(x, y) -= fluxnorm * totalflux;
                neighbour= 0;
                if( x ) {
                    if( surplus_mx > 0.0 ) {
                        h(x-1, y) += fluxnorm * flux[neighbour++];
                        done(x-1, y)= false;
                    }
                    else
                        done(x-1, y)= done(x-1, y) && !(h(x-1, y) - h(x, y) > m_dh_stat);
                }
                if( x< duneglobals::nx()-1 ) {
                    if( surplus_px > 0.0 ) {
                        h(x+1, y) += fluxnorm * flux[neighbour++];
                        done(x+1, y)= false;
                    }
                    else
                        done(x+1, y)= done(x+1, y) && !(h(x+1, y) - h(x, y) > m_dh_stat);
                }
                if( y ) {
                    if( surplus_my > 0.0 ) {
                        h(x, y-1) += fluxnorm * flux[neighbour++];
                        done(x, y-1)= false;
                    }
                    else
                        done(x, y-1)= done(x, y-1) && !(h(x, y-1) - h(x, y) > m_dh_stat);
                }
                if( y< duneglobals::ny()-1 ) {
                    if( surplus_py > 0.0 ) {
                        h(x, y+1) += fluxnorm * flux[neighbour++];
                        done(x, y+1)= false;
                    }
                    else
                        done(x, y+1)= done(x, y+1) && !(h(x, y+1) - h(x, y) > m_dh_stat);
                }
            }
        }
    }
    
    cout << "cellaval took " << iteration << " iterations\n";
    
    if( !all_done )
        //    cerr << "cellaval::calc: giving up after " << iteration << " iterations.\n";
        cout << "cellaval::calc: giving up after " << iteration << " iterations.\n";
}


/*!  This function handles the call of the qsort function for sorting the
 coordinate pairs in m_sliporder.  They are ordered so that the highest points
 are processed first.  */

void cellaval::set_sliporder( TFktScal &h )
{
    hpointer= &h;
    qsort( m_sliporder, duneglobals::nx()*duneglobals::ny(), 2*sizeof(int), 
          cmp_coordpair );
}


/*!  Comparison function for sorting the coordinate pairs of grid sites in
 descending order.  */

int cellaval::cmp_coordpair(const void *a, const void *b)
{
    double vala= (*hpointer)(*(int*)a, ((int*)a)[1]);
    double valb= (*hpointer)(*(int*)b, ((int*)b)[1]);
    
    if( vala > valb )
        return -1;
    else if( valb > vala )
        return 1;
    else
        return 0;
}



