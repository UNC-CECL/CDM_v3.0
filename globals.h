/******************************************************************************
 $Id: globals.h,v 1.13 2005/02/01 12:40:21 schatz Exp $
 
 Contains three globally used classes: dunepar, which handles the parameter
 file, duneglobals, the singleton containing global parameters, and dunedata,
 the parent class of all classes storing variable fields.
 ******************************************************************************/

#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "func.h"

///*!
//  Exception to be thrown in the case of missing or illegal parameter values.
// */
//class para_error
//{
//public:
//  para_error() {}
//  ~para_error() {}
//};

using namespace std;

/*!  Parameter handling class which reads the parameter file.  This class is
 too slow to be used within a calculation.  It should be called in the
 constructor of some computing object to initialise member variables, which
 are subsequently used.  */

class dunepar
{
public:
    dunepar();
    ~dunepar();
    
    void scan( int argc, char **argv );
    
    template<class T> bool get( const char *s, T& retval ) const;
    
    bool exists(const char *s) const { return 0<=find(s); }
    bool exists(const string& s) const { return 0<=find(s.c_str()); }
    
    /*!  Get a parameter and exit the program if it does not exist.  The return
     type must be given explicitly when calling this function since g++ does not
     recognise it from the l-value (this may be standard).  */
    template<class T> T getrequired( const char *s, bool verbose= true ) const
    {
        T retval;
        if( !get( s, retval ) ) {
            cerr << "ERROR: Required parameter `" << s << "' not found!" << endl;
            exit(1);
        }
        else if( verbose )
            cout << "dunepar: Read parameter " << s << " = " << retval << "\n";
        return retval;
    }
    
    /*!  Same as the eponymous function which takes a char * as its argument.  */
    template<class T> T getrequired( string s, bool verbose= true ) const
    { return getrequired<T>(s.c_str(), verbose); }
    
    /*!  Get a parameter or return default if it does not exist.    */
    template<class T> T getdefault( const char *s, const T& defaultval,
                                   bool verbose= true ) const
    {
        T retval;
        if( !get( s, retval ) ) {
            if( verbose )
                cout << "dunepar: Using default parameter " << s << " = " << defaultval << "\n";
            return defaultval;
        }
        else {
            if( verbose )
                cout << "dunepar: Read parameter " << s << " = " << retval << "\n";
            return retval;
        }
    }
    
    /*!  Same as the eponymous function which takes a char * as argument.  */
    template<class T> T getdefault( string s, const T& defaultval,
                                   bool verbose= true ) const
    { return getdefault<T>(s.c_str(), defaultval, verbose); }
    
private:
    int find( const char *key ) const;
    
    /*!  Storage space for contents of parameter file.  */
    char *m_parfile;
    /*!  Pointers to keys (parameter names).  */
    char **m_key;
    /*!  Pointers to values of parameters.  */
    char **m_value;
};


/*!  Singleton class for global parameters and quantities derived from them.
 See the member variables for a description.  The eponymous static functions
 just return the members.
 
 initialise has to be called with a valid dunepar object before anything
 else, otherwise a segmentation fault will occur.  It was decided against
 checking the existence of the instance every time for reasons of efficiency.  */

class duneglobals
{
public:
    /*! Creates the unique private instance of this class from the
     parameters. If called again, the old instance is destroyed first.  */
    static void initialise( const dunepar& parameters )
	{ if( instance ) delete instance;
        instance= new duneglobals( parameters ); }
    
    //  template<type T> static T checkreturn( const T& member );
    
    static int startstep() { return instance->m_startstep; }
    static double starttime() { return instance->m_starttime; }
    static bool sim3d() { return instance->m_3d; }
    static double length() { return instance->m_length; }
    static double width() { return instance->m_width; }
    static int nx() { return instance->m_nx; }
    static int ny() { return instance->m_ny; }
    static double dx() { return instance->m_dx; }
    static double z0eff() { return instance->m_z0eff; }
    static double repose_stat() { return instance->m_repose_stat; }
    static double repose_dyn() { return instance->m_repose_dyn; }
    static bool periodic_x() { return instance->m_periodic_x; }
    static bool periodic_y() { return instance->m_periodic_y; }
    static bool messages() { return instance->m_messages; }
    static bool timing() { return instance->m_timing; }
    static double rho_fluid() { return instance->m_rho_fluid; }
    static double rho_grains() { return instance->m_rho_grains; }
    static double rho_sand() { return instance->m_rho_sand; }
    /*! parameters for flux calculation*/
    static double d_grain() { return instance->m_d_grain; }
    static double M() { return instance->m_M; }
    static double log_z1z0() { return instance->m_log_z1z0; }
    static double b_z1() { return instance->m_b_z1; }
    static double DeltaU() { return instance->m_DeltaU; }
    static double Cd_dx() { return instance->m_Cd_dx; }
    static double alpha() { return instance->m_alpha; }
    static double alpha2_g() { return instance->m_2alpha_g; }
    static double u_star_t() { return instance->m_u_star_t; }
    static double u_star_ft() { return instance->m_u_star_ft; }
    static double tau_t() { return instance->m_tau_t; }
    static double tau_ft() { return instance->m_tau_ft; }
    static double gamma() { return instance->m_gamma; }
    static double gamma_2alpha_g2() { return instance->m_gamma_2alpha_g2; }
    static double beta_2alpha_g_tau_ft() { return instance->m_beta_2alpha_g_tau_ft; }
    
    static string datadir() { return instance->m_datadir; }
    // Units
    static double secday() { return instance->m_secday; }
    static double secmonth() { return instance->m_secmonth; }
    static double secyear() { return instance->m_secyear; }
    // Time conversion
    static double timefrac() { return instance->m_rtime; }
    // Shore param
    static double HMWL() { return instance->m_shore_HMWL; }
    static double MSL() { return instance->m_shore_watertable; }
    static double angle() { return instance->m_beach_angle; }
    static double slope() { return instance->m_beach_slope; }
        
    static double C_drag(double Re);
    
    // PARTELI BEGIN
    static bool check_error() { return instance->m_check_error; }
    static bool high_precision() { return instance->m_high_precision; }
    // PARTELI END
    
private:
    static class duneglobals *instance;
    
    duneglobals( const dunepar& parameters );
    ~duneglobals() {}
    
    /*!  Number of iteration steps already done (for continuing simulations).  */
    int m_startstep;
    /*!  Time which has already passed (for continuing simulations).  */
    double m_starttime;
    /*! True for 3D simulation, currently signaled only by ny>3.  */
    bool m_3d;
    /*! Grid size in wind direction (nx) and perpendicular to wind direction
     (ny).  */
    int m_nx, m_ny;
    /*!  Length (in x direction) and width (in y direction) of the simulation
     region.  */
    double m_length, m_width;
    /*! Grid spacing in metres.  */
    double m_dx;
    /*! Roughness length for fluid shear stress calculations.  This is an
     effective roughness length taking into account sand entrainment.  Therefore
     use of this value implies a stationary situation.  */
    double m_z0eff;
    /*! Static/dynamic angle of repose in degrees.  */
    double m_repose_stat, m_repose_dyn;
    /*! Periodic boundary conditions in x and y directions.  Bear in mind that
     since we use a discrete Fourier transform in several places, even for open
     boundary conditions the simulated dune(s) will not be truly alone on an
     infinite plane.  (A discrete Fourier transform always describes a periodic
     function.)  For open boundary conditions the Fourier-transformed arrays
     will be extended by a factor 1.8, which is thought sufficient for the
     dune not to be influenced by its imaginary neighbours.  */
    bool m_periodic_x, m_periodic_y;
    
    /*! Verbose message output.  */
    bool m_messages;
    /*! Time computation steps if this is true. Useless without m_messages.  */
    bool m_timing;
    /*! Density of fluid driving the grains (normally air) in kg/m^3.  */
    double m_rho_fluid;
    /*! Density of the material the grains consist of (normally quartz for
     sand) in kg/m^3. */
    double m_rho_grains;
    /*! Density of the sand including the gaps between the grains (= m_rho_grains
     * packing density) in kg/m^3.  */
    double m_rho_sand;
    
    /*! Parameters for flux calculation*/
    double m_M, m_d_grain, m_alpha, m_log_z1z0, m_b_z1, m_DeltaU, m_Cd_dx, m_2alpha_g, m_u_star_t, m_u_star_ft, m_tau_t, m_tau_ft, m_gamma, m_gamma_2alpha_g2, m_beta_2alpha_g_tau_ft;
    
    /*!  Directory all data are to be saved to.  */
    string m_datadir;
    
    /* Real Time / Sim Time */
    double m_rtime;
    /*! Units change*/
    double m_secday, m_secmonth, m_secyear;
    /*Shore param*/
    double m_shore_HMWL, m_shore_watertable, m_beach_angle, m_beach_slope;
    
    //  PARTELI BEGIN
    bool m_check_error, m_high_precision;
    // PARTELI END
    
};


/*!  Base class of all objects which contain array data in this dune software.
 Objects derived from this class are registered automatically, and all their
 data can be saved with one function call.  A virtual function for saving
 arrays has to be implemented by subclasses.  */

class dunedata
{
public:
    dunedata(const dunepar &par);
    virtual ~dunedata();
    
    /*!  Save arrays (discretised physical quantities) to disk.  */
    virtual void save_arrays()=0;
    
    static void save_all_data();
    
protected:
    ofstream *open_writeout( string basename );
    void save_2d_scalarray( string basename, const TFktScal& data );
    void save_2d_vecarray( string basename, const TFktVec& data );
    void save_1d_scalarray( string basename, const CFunc1d& data );
    
private:
    /*!  True if x is to be the column index in the saved files (varying along a
     line).  */
    bool m_xline;
    /*!  Reference to parameter file object for checking whether a variable field
     should be saved.  */
    const dunepar& m_par;
    /*!  Next object in list.  */
    dunedata *m_next;
    
    static dunedata *first_object;
};




#endif // GLOBALS_H


