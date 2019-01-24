/******************************************************************************
 $Id: globals.cc,v 1.22 2005/04/13 15:00:18 parteli Exp $
 ******************************************************************************/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

#include <cstring>

#include <sstream>
#include <iomanip>

#include "globals.h"
#include "evolution.h"
#include "rfftw12d.h"

//*****************************************************************************
//  class dunepar


/*!  Initialises m_parfile, m_key and m_value to NULL.  */

dunepar::dunepar()
{
    m_parfile= NULL;
    m_key= NULL;
    m_value= NULL;
}


/*!  Frees the space occupied by m_parfile, m_key and m_value.  */

dunepar::~dunepar()
{
    free(m_parfile);
    delete[] m_key;
    delete[] m_value;
}


/*!  Reads the parameters from the command line and from a parameter file,
 which may be given in the command line or will be assumed to be default.par.
 Command line parameters override those from the file.  A line in the parameter
 file can be empty or a key, followed by optional white space, followed by an
 equals sign, optional white space and the value corresponding to the key.  Keys
 and values must not contain white space or equals signs.  Command-line
 arguments have the same format, except that there must be no white spaces.  In
 the parameter file, everything after a "#" character is consideren a comment
 and ignored.  */

void dunepar::scan( int argc, char **argv )
{
    FILE *filep;
    char *fname, *newparfile, *scan, *line;
    long size, nlines;
    int argind, parind, parind2, parerr, fileerr;
    
    fname= NULL;
    //  If the first command-line argument does not contain an =, we interpret it
    //  as the name of the parameter file
    if( argc > 1 ) {
        for( scan= argv[1]; *scan && *scan!='='; ++scan );
        if( !*scan )
            fname= argv[1];
    }
    for( argind= 1; argind< argc; ++argind )
        if( !strncmp(argv[argind], "parafile=", 9) ) {
            if( !fname )
                fname= argv[argind]+9;
            else
                cerr << "dunepar::scan: Name of parameter file given more than once"
                " on the command line.  Ignoring `" << argv[argind] << "\n";
        }
    
    if( !fname )
        fname= "default.par";
    
    if( !strcmp(fname, "-") )
    {
        char *write;
        
        cout << "dunepar::scan: Reading parameters from standard input...\n";
        size= 1024;
        //  Use malloc to be able to use realloc later
        newparfile= (char*)malloc(size);
        write= newparfile;
        while(true)
        {
            cin.read(write, newparfile + size - write);
            if( write + cin.gcount() < newparfile + size )
                break;
            newparfile= (char*)realloc(newparfile, 2*size);
            write= newparfile + size;
            size *= 2;
        }
        write[cin.gcount()]= 0;
    }
    else {
        cout << "dunepar::scan: Reading parameter file `" << fname << "'...\n";
        filep= fopen( fname, "r" );
        if( !filep ) {
            cerr << "dunepar::scan: Error opening parameter file `" << fname << "'.  Aborting.\n";
            exit(1);
        }
        fseek( filep, 0L, SEEK_END );
        size= ftell( filep );
        fseek( filep, 0L, SEEK_SET );
        newparfile= (char*)malloc(size+1);
        fileerr= (size_t)size!=fread( newparfile, 1L, size, filep );
        newparfile[size]= 0;
        fclose( filep );
        if( fileerr ) {
            cerr << "dunepar::scan: Error reading parameter file `" <<
            fname << "'. Aborting.\n";
            //delete[] newparfile;
            //exit(-1);
        }
    }
    
    free(m_parfile);
    delete[] m_key;
    delete[] m_value;
    m_parfile= newparfile;
    nlines= 0;
    for( scan= m_parfile; *scan; ++scan )
        if( *scan=='\n' )
            ++nlines;
    m_key= new char*[nlines + argc];
    m_value= new char*[nlines + argc];
    for( parind= nlines + argc - 1; parind>=0; --parind )
        m_key[parind]= m_value[parind]= NULL;
    line= m_parfile;
    nlines= 0;
    parind= 0;
    parerr= 0;
    while( true )
    {
        ++nlines;
        while( *line==' ' || *line=='\t' )
            ++line;
        for( scan= line; *scan && *scan!='\n' &&
		    *scan!='#' && *scan!='='; ++scan );
        if( *scan=='#' )
            while( *scan && *scan!='\n' )
                ++scan;
        if( !*scan )
            break;
        if( *scan=='\n' ) {  // This line contained at most a comment
            line= scan+1;
            continue;
        }
        //  Now scan the line more carefully:
        if( *line == '=' ) {
            cerr << "dunepar::scan: Error in parameter file: empty key in line " << nlines << ".\n";
            parerr= 1;
            line= scan + 1;
            continue;
        }
        scan= line;
        while( *scan!=' ' && *scan!='\t' && *scan!='=' )
            ++scan;
        *scan++= 0;
        if( 0 <= find(line) ) {
            cerr << "dunepar::scan: Error in parameter file: key in line " << nlines << " has been defined before.\n";
            parerr= 1;
            for( line= scan; *line && *line!='\n'; ++line );
            if( *line )
                ++line;
            continue;
        }
        while( *scan==' ' || *scan=='\t' || *scan=='=' )
            ++scan;
        if( *scan=='\n' || *scan=='#' || !*scan ) {
            cerr << "dunepar::scan: Error in parameter file: empty value in line " << nlines << ".\n";
            parerr= 1;
            for( line= scan; *line && *line!='\n'; ++line );
            if( *line )
                ++line;
            continue;
        }
        m_key[parind]= line;
        m_value[parind]= scan;
        ++parind;
        while( *scan!=' ' && *scan!='\t' && *scan!='\n' && *scan && *scan!='#' )
            ++scan;
        if( *scan=='\n' ) {
            *scan++= 0;
            line= scan;
        }
        else if( !*scan )
            line= scan;
        else {
            *scan++ = 0;
            while( *scan && *scan!='\n' )
                ++scan;
            if( *scan=='\n' )
                ++scan;
            line= scan;
        }
    }
    
    // Now scan the command line again:
    for( argind= 1; argind< argc; ++argind ) {
        for( scan= argv[argind]; *scan && *scan!='='; ++scan );
        if( !*scan ) {
            if( argind!=1 ) {
                cerr << "dunepar::scan: Error: Command line argument " << argind <<
                " does not have the form <key>=<value>.\n";
                parerr= 1;
            }
            continue;
        }
        *scan++= 0;
        if( !*scan || scan==argv[argind]+1 ) {
            cerr << "dunepar::scan: Error: Empty key or value in command line "
            "argument " << argind << ".\n";
            parerr= 1;
            continue;
        }
        if( (parind2= find(argv[argind])) >= 0 )
            m_value[parind2]= scan;
        else {
            m_key[parind]= argv[argind];
            m_value[parind]= scan;
            ++parind;
        }
    }
    
    if( parerr )
        exit(1);
    
    //  Convert "true" to 1 and "false" to 0 so we can use them for boolean
    //  parameters:
    for( parind= 0; m_value[parind]; ++parind )
        if( !strcmp(m_value[parind], "true") ) {
            *m_value[parind]= '1';
            m_value[parind][1]= 0;
        }
        else if( !strcmp(m_value[parind], "false") ) {
            *m_value[parind]= '0';
            m_value[parind][1]= 0;
        }
}


/*!  Writes the value corresponding to the key \a s to \a retval if it can
 find() it and convert it to the desired type.  Returns true on success.  */

template<class T> bool dunepar::get( const char *s, T& retval ) const
{
    int parind= find(s);
    
    if( parind< 0 )
        return false;
    
    string value(m_value[parind]);
    std::istringstream strm(value);
    
    strm >> retval;
    if( strm.fail() )
        cerr << "dunepar::get: Could not convert parameter `" << s << "' to the right  type.\n";
    return !strm.fail();
}


//  Instantiate explicitly here to save space - otherwise it might be inlined
template bool dunepar::get<bool>( const char *s, bool& retval ) const;
template bool dunepar::get<int>( const char *s, int& retval ) const;
template bool dunepar::get<double>( const char *s, double& retval ) const;
template bool dunepar::get<std::string>( const char *s, std::string& retval ) const;


/*!  Searches for the key by comparing the \a key with each key from m_key.
 This is pretty slow but it will do for now.  */

int dunepar::find(const char *key) const
{
    int ind;
    
    for( ind= 0; m_key[ind]; ++ind )
        if( !strcmp(key, m_key[ind]) )
            return ind;
    
    return -1;
}


//*****************************************************************************
//  class duneglobals


//! Pointer to single instance of duneglobals class.
class duneglobals *duneglobals::instance= NULL;


///*!  Guards against non-existence of instance before returning a member
//	variable. */
//
//template<type T> static T duneglobals::checkreturn( const T& member )
//{
//	if( !instance ) {
//	  cerr << "ERROR: duneglobals asked for value before instance exists!" << endl;
//	  return T( 0L );
//	}
//	else {
//	  return member;
//	}
//}


/*! The constructor extracts global parameters from the parameter file via
 dunepar.  It performs some sanity checks and aborts the program if vital
 parameters were not found or had illegal values.  */

duneglobals::duneglobals( const dunepar& parameters )
{
    struct stat savedirstat;
    
    m_startstep= parameters.getdefault("Nt0", 0);
    m_starttime= parameters.getdefault("time0", 0.0);
    m_dx= parameters.getrequired<double>("dx");
    if( parameters.exists("NX") && parameters.exists("length") ) {
        cerr << "duneglobals constructor: Parameters `NX' and `length' must not both be given!  Aborting...\n";
        exit(1);
    }
    else if( parameters.exists("length") ) {
        m_length= parameters.getrequired<double>("length");
        m_nx= (int)ceil(m_length/m_dx);
        m_length= m_nx*m_dx;
    }
    else {
        m_nx= parameters.getrequired<int>( "NX" );
        m_length= m_nx*m_dx;
    }
    if( m_nx<=3 ) {
        cerr << "ERROR: duneglobals constructor found NX<=3" << endl;
        exit(1);
    }
    if( parameters.exists("NY") && parameters.exists("width") ) {
        cerr << "duneglobals constructor: Parameters `NY' and `width' must not both be given!  Aborting...\n";
        exit(1);
    }
    else if( parameters.exists("width") ) {
        m_width= parameters.getrequired<double>("width");
        m_ny= (int)ceil(m_width/m_dx);
        m_width= m_ny*m_dx;
    }
    else {
        m_ny= parameters.getrequired<int>( "NY" );
        m_width= m_ny*m_dx;
    }
    
    m_3d= m_ny > 1;
    m_z0eff= parameters.getdefault("hlr.z0", 1e-3/*1e-3*/);
    
    m_repose_stat= parameters.getdefault("aval.angle_repose_stat",34);
    m_repose_dyn= parameters.getdefault("aval.angle_repose_dyn",33);
    
    m_periodic_x= parameters.getdefault( "calc.x_periodic", false );
    if( m_periodic_x && m_nx!=fft::GetNextPowerOf2(m_nx) ) {
        cerr << "duneglobals constructor:  Periodic boundary condition in x direction requires NX to be power of 2!  Aborting...\n";
        exit(1);
    }
    m_periodic_y= parameters.getdefault( "calc.y_periodic", true );
    if( m_periodic_y && m_ny!=fft::GetNextPowerOf2(m_ny) ) {
        cerr << "duneglobals constructor:  Periodic boundary condition in y direction requires NY to be power of 2!  Aborting...\n";
        exit(1);
    }
    
    m_messages= parameters.getdefault( "msg", false );
    m_timing= parameters.getdefault( "timing", m_messages );
    
    m_rho_fluid= parameters.getdefault( "fluid_density", 1.225);
    m_rho_grains= parameters.getdefault( "grain_density", 2650);
    m_rho_sand= m_rho_grains * parameters.getdefault( "packing", 0.6226); //1650
    
    /* Shore parameters*/
    m_shore_HMWL = parameters.getdefault("shore.MHWL", 0.0);
    m_shore_watertable = parameters.getdefault("shore.sealevel", 0.0);
    m_beach_angle = parameters.getdefault("beach.angle", 1.0);
    m_beach_slope = tan(m_beach_angle * M_PI / 180.);

    /* Real Time / Sim Time */
    m_rtime = parameters.getdefault( "wind.fraction", 1.0);
    /* Units conversion*/
    m_secday = 60*60*24; // seconds in a day
    m_secmonth = 60*60*24*30; // seconds in a month (30 days)
    m_secyear = 60*60*24*365; // seconds in a year (365 days)
    
    /*!Parameters for flux calculation*/
    m_d_grain= parameters.getdefault( "salt.d", 250e-6);
    double m_fluid_viscosity= parameters.getdefault( "fluid_viscosity", 1.8e-5);
    double fluid_viscosity_kin = m_fluid_viscosity / m_rho_fluid;
    double m_g= parameters.getdefault( "salt.g", 9.8 );
    
    double s = m_rho_grains / m_rho_fluid;
    double Aft = parameters.getdefault( "Bagnold-shields_parameter", 0.11);
    double At = 0.8 * Aft;
    m_u_star_ft = Aft * sqrt(m_g*m_d_grain*(s-1.));
    m_u_star_t = (At/Aft) * m_u_star_ft;
    cout << "duneglobals constructor: u*ft = " << m_u_star_ft << "m/s has been calculated from Bagnold-Shields Parameter\n";
    cout << "duneglobals constructor: u*t = " << m_u_star_t << "m/s (= ratio_u*t/u*ft multiplied by u*ft)\n";
        
    double r= parameters.getdefault( "salt.Z/zm", 0.);
    
    /*! Auxiliar quantities */
    double t_v = exp(1./3. * log(fluid_viscosity_kin / (m_g * m_g)));	//characteristic time
    double l_v = exp(1./3. * log(fluid_viscosity_kin * fluid_viscosity_kin / (At*At * m_g*(s-1))));	//characteristic lenght
    double S = 0.25 * m_d_grain / fluid_viscosity_kin * sqrt(m_g*m_d_grain*(s-1.));
    
    // Parameters
    double m_C_diff= parameters.getdefault( "salt.D", 0.0 );
    double m_beta= parameters.getdefault( "beta", 5.7e-4 );
    m_gamma= parameters.getdefault( "gamma", 0.1);
        
    // derived quantities:
    double m_z0= m_d_grain / 20.;
    cout << "duneglobals constructor: z0 = " << m_z0 << endl;
    double m_zm= 14./(1.+1.4*r) * m_u_star_t * t_v;
    cout << "duneglobals constructor: zm = " << m_zm << endl;
    double m_z1= 35. * l_v;
    cout << "duneglobals constructor: z1 = " << m_z1 << endl;
    m_alpha= 0.18 * (m_d_grain / l_v);
    cout << "duneglobals constructor: alpha = " << m_alpha << endl;
    
    double Ad = 0.95, Bd = 5.12; //natural sand
    double m_C_drag = 4./3. * (Ad + sqrt(2 * m_alpha) * Bd / S)*(Ad + sqrt(2 * m_alpha) * Bd / S);
    cout << "duneglobals constructor: Cd = " << m_C_drag << endl;
    
    m_DeltaU = sqrt(4.0/3.0 * m_g * (s-1) * m_d_grain / (m_C_drag * 2 * m_alpha));
    
    double z0zm= m_z0/m_zm;
    double z1zm= m_z1/m_zm;
    m_log_z1z0= log(m_z1/m_z0);
    
    if( r < z0zm )	r = z0zm;
    double log_r = log(r/z0zm);    
    double b_inf = -0.5772 - log(z0zm) + z0zm;
    double b_Z = (r < 2. ? log_r-(r - z0zm)+0.25*(r*r - z0zm*z0zm)-0.042*(r*r*r - z0zm*z0zm*z0zm) : b_inf - exp(-r)/r);
    m_M = (log_r > 0 ? log_r/b_Z : 1.);
    cout << "duneglobals constructor: M = " << m_M << endl;
    
    m_b_z1 = (z1zm < 2. ? m_log_z1z0-(z1zm - z0zm)+0.25*(z1zm*z1zm - z0zm*z0zm)-0.042*(z1zm*z1zm*z1zm - z0zm*z0zm*z0zm) : b_inf - exp(-z1zm)/z1zm);
    m_Cd_dx= m_C_diff / (m_dx*m_dx);
    m_tau_ft= m_rho_fluid * m_u_star_ft * m_u_star_ft;
    m_tau_t= m_rho_fluid * m_u_star_t * m_u_star_t;
    m_2alpha_g= 2.0*m_alpha / m_g;
    m_gamma_2alpha_g2= m_gamma / (m_2alpha_g * m_2alpha_g);
    m_beta_2alpha_g_tau_ft = m_beta / (m_2alpha_g * m_tau_ft);
    
    // PARTELI BEGIN
    m_check_error = parameters.getdefault( "check.error", false );
    m_high_precision = parameters.getdefault( "high.precision", false );
    // PARTELI END
    
    m_datadir= parameters.getdefault( "save.dir", string("DAT") );
    if( stat( m_datadir.c_str(), &savedirstat ) ) {
        cerr << "duneglobals constructor:  Cannot stat `" << m_datadir <<
        "'.  Data directory doesn't exist or permission denied.  Aborting.\n";
        exit(1);
    }
    if( !S_ISDIR(savedirstat.st_mode) ) {
        cerr << "duneglobals constructor:  Data directory `" << m_datadir <<
        "' is not a directory.  Aborting.\n";
        exit(1);
    }
}




/*!  Returns the drag coefficient of a sphere depending on the (grain-scale)
 Reynolds number.  */

double duneglobals::C_drag(double Re)
{
    static double drag_lut[]= { 0.2, 100, 0.5, 45, 1, 30, 2, 14.5, 5, 7.5,
	    10, 6.2, 20, 4.2, 50, 1.5, 100, 1, 200, 0.75, 500, 0.55,
	    1000, 0.45, 2000, 0.3, 5000, 0.37, 1e4, 0.4, 2e4, 0.47, 5e4, 0.5,
	    1e5, 0.47, 2e5, 0.425, 2.4e5, 0.125, 4e5, 0.175, 6e5, 0.2, -1, -1 };
    static int lutsize, init_done= 0;
    
    int left, right, mid;
    
    if( !init_done ) {
        int i;
        
        for( lutsize= 0; drag_lut[lutsize]> 0; ++lutsize );
        lutsize /= 2;
        for( i= 0; i< 2*lutsize; ++i )
            drag_lut[i]= log(drag_lut[i]);
        init_done= 1;
    }
    
    Re= log(Re);
    if( Re< *drag_lut )
        return exp(drag_lut[1]);
    if( Re> drag_lut[2*lutsize-2] )
        return exp(drag_lut[2*lutsize-1]);
    left= 0;
    right= lutsize-1;
    while( left < right-1 )
    {
        mid= (left + right)/2;
        if( Re> drag_lut[2*mid] )
            left= mid;
        else
            right= mid;
    }
    return exp( drag_lut[2*left+1] + (Re - drag_lut[2*left]) *
               (drag_lut[2*right+1] - drag_lut[2*left+1]) /
               (drag_lut[2*right] - drag_lut[2*left]) );
}




//*****************************************************************************
//  class dunedata


/*!  First object in list.  */
dunedata *dunedata::first_object= NULL;

/*!  Appends itself to the list starting with first_object.  */

dunedata::dunedata(const dunepar &par) :
m_par(par)
{
    dunedata *p;
    
    m_next= NULL;
    if( !first_object )
        first_object= this;
    else {
        for( p= first_object; p->m_next; p= p->m_next );
        p->m_next= this;
    }
    m_xline= par.getdefault("save.x-line", true);
}


/*!  Removes itself from the record list.  */

dunedata::~dunedata()
{
    dunedata *p;
    
    for( p= first_object; p; p= p->m_next )
        if( p->m_next==this )
        {
            p->m_next= m_next;
            break;
        }
}


/*!  Calls save_arrays function of all dunedata objects.  */

void dunedata::save_all_data()
{
    dunedata *p;
    
    cout << "dunedata::save_all_data: saving discretised functions (grid spacing " << duneglobals::dx() << " metres) ... ";
    for( p= first_object; p; p= p->m_next )
        p->save_arrays();
    cout << "done" << endl;
}


/*!  Opens a file stream for the output file \a basename "." evolution::steps()
 ".dat" in the directory duneglobals::datadir().  The pointer to a stream
 object is returned, which has to be deleted by the caller.  (Due to some copy
 constructors of component classes being private, it seems to be impossible to
 return a stream by value.)  */

ofstream *dunedata::open_writeout( string basename )
{
    std::ostringstream namestream;
    
    namestream << duneglobals::datadir() << '/' << basename << '.'
    << setw(5) << setfill('0') << evolution::steps() << ".dat";
    return new ofstream( namestream.str().c_str() );
}


/*!  Saves two-dimensional scalar array.  The x index becomes the column
 index, y the line index.  */

void dunedata::save_2d_scalarray( string basename, const TFktScal& data )
{
    ofstream *os;
    int x, y;
    
    if( m_par.getdefault( "dontsave." + basename, false, false ) )
        return;
    os= open_writeout( basename );
    
    // PARTELI BEGIN
    if( duneglobals::high_precision()) (*os) << fixed << setprecision(16);
    // PARTELI END
    
    if( m_xline )
        for( y= 0; y< data.SizeY(); ++y ) {
            for( x= 0; x< data.SizeX(); ++x )
                (*os) << data(x, y) << ' ';
            (*os) << endl;
        }
    else
        for( x= 0; x< data.SizeX(); ++x ) {
            for( y= 0; y< data.SizeY(); ++y )
                (*os) << data(x, y) << ' ';
            (*os) << endl;
        }
    delete os;
}


/*!  Saves two-dimensional vector array to two different files with basenames
 \a basename "_x" and \a basename "_y".  x is the column index, y the line
 index.  */

void dunedata::save_2d_vecarray( string basename, const TFktVec& data )
{
    ofstream *os;
    int x, y;
    
    if( m_par.getdefault( "dontsave." + basename, false, false ) )
        return;
    {
        os= open_writeout( basename+"_x" );
        //  PARTELI BEGIN
        if( duneglobals::high_precision()) (*os) << fixed << setprecision(16);
        //  PARTELI END
        if( m_xline )
            for( y= 0; y< data.SizeY(); ++y ) {
                for( x= 0; x< data.SizeX(); ++x )
                    (*os) << data(x, y)[0] << ' ';
                (*os) << endl;
            }
        else
            for( x= 0; x< data.SizeX(); ++x ) {
                for( y= 0; y< data.SizeY(); ++y )
                    (*os) << data(x, y)[0] << ' ';
                (*os) << endl;
            }
        delete os;
    }
    
    if( !m_par.getdefault( "dontsave." + basename + "_y", false, false ) )
    {
        os= open_writeout( basename+"_y" );
        //  PARTELI BEGIN
        if( duneglobals::high_precision()) (*os) << fixed << setprecision(16);
        //  PARTELI END
        if( m_xline )
            for( y= 0; y< data.SizeY(); ++y ) {
                for( x= 0; x< data.SizeX(); ++x )
                    (*os) << data(x, y)[1] << ' ';
                (*os) << endl;
            }
        else
            for( x= 0; x< data.SizeX(); ++x ) {
                for( y= 0; y< data.SizeY(); ++y )
                    (*os) << data(x, y)[1] << ' ';
                (*os) << endl;
            }
        delete os;
    }
}


/*!  Saves one-dimensional scalar array.  The x index becomes the column
 index, y the line index.  */

void dunedata::save_1d_scalarray( string basename, const CFunc1d& data )
{
    ofstream *os;
    int x;
    
    if( m_par.getdefault( "dontsave." + basename, false, false ) )
        return;
    os= open_writeout( basename );
    for( x= 0; x< data.GetSize(); ++x )
        (*os) << data[x] << endl;
    delete os;
}

