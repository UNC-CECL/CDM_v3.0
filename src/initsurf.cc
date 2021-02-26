/******************************************************************************
  $Id: initsurf.cc,v 1.14 2005/01/21 14:36:08 schatz Exp $
******************************************************************************/

#include <fstream>
#include <ctype.h>

#include "func.h"
#include "initsurf.h"
#include "initsurfbeach.h"
#include "initsurfalea.h"

//*****************************************************************************
// class arrayinit

/*!  Performs the initialisation of a one-dimensional scalar array.  The
  default implementation aborts the program with an error message: a call to
  a not reimplemented initialise function is interpreted as an unsuitable
  argument type.  */
void arrayinit::init_1d_scal(CFunc1d& array)
{
  array_mismatch("init_1d_scal");
}


/*!  Performs the initialisation of a two-dimensional scalar array.  The
  default implementation aborts the program with an error message: a call to
  a not reimplemented initialise function is interpreted as an unsuitable
  argument type.  */
void arrayinit::init_2d_scal(TFktScal& array)
{
  array_mismatch("init_2d_scal");
}


/*!  Performs the initialisation of a two-dimensional vector-valued array.
  The default implementation aborts the program with an error message: a call
  to a not reimplemented initialise function is interpreted as an unsuitable
  argument type.  */
void arrayinit::init_2d_vec(TFktVec& array)
{
  array_mismatch("init_2d_vec");
}


/*!  Creates an arrayinit object of the type requested by the parameter file.
  If \a prefix is given, it is prepended to all parameter names.  This allows
  creating objects for the initialisation of different variable fields.  */

arrayinit *arrayinit::create(const dunepar& par, string prefix)
{
  arrayinit *object;
  
  string strType = par.getrequired<string>(prefix + "Init-Surf");
  if (strType == "alea") {
    object = new CInitSurfAlea(par, prefix);
  } else if (strType == "init_h") {
    object = new arrayinit_ascii(par, prefix);
  } else if (strType == "plain") {
    object = new CInitSurfPlain(par, prefix);
  } else if (strType == "beach") {
    object = new CInitSurfBeach(par, prefix);
  } else {
    cout << "arrayinit::create: FATAL ERROR: Unknown value `" << strType << "\' for parameter `" << prefix << "Init-Surf\' !" << endl;
    cout << "  Valid values are: gauss, ridge, paraboloid, matlab, alea, init_h, plain" << endl;
    exit(1);
  }
  return object;
}


/*!  Prints an error message in the name of the function \a function and exits
  program.  Called by all default implemetations of initialisation functions.  */
void arrayinit::array_mismatch(const char *function)
{
  cerr << "Fatal error: arrayinit::" << function 
	<< ":  Initialisation of this array type"
	<< " not reimplemented in the selected initialisation class!" << endl;
  exit(1);
}


//*****************************************************************************
// class CInitSurfPlain

/*!  Writes the constant m_const to all positions of the array.  */
void CInitSurfPlain::init_1d_scal(CFunc1d& array)
{
  int x;
  
  for( x= 0; x< array.GetSize(); ++x )
    array[x]= m_const;
}


/*!  Writes the constant m_const to all positions of the array.  */
void CInitSurfPlain::init_2d_scal(TFktScal& array)
{
  int x, y;
  
  for( x= 0; x< array.SizeX(); ++x )
    for( y= 0; y< array.SizeY(); ++y )
      array(x, y)= m_const;
}
void CInitSurfPlain::init_2d_vec(TFktVec& array)
{
    int x, y;
    
    for( x= 0; x< array.SizeX(); ++x )
        for( y= 0; y< array.SizeY(); ++y ){
            array(x, y)[0]= m_const;
            array(x, y)[1]= m_const_aux;
        }

}


//*****************************************************************************
// class arrayinit_ascii


/*!  Reads parameters init_h.file and init_h.x-line (optional) from parameter
  file.  If init_h.x-line is false, x is assumed to be the line index in the
  data table in the file, otherwise it is the column index (default).  */

arrayinit_ascii::arrayinit_ascii( const dunepar &par, string prefix)
{
  m_filename= par.getrequired<string>( prefix+"init_h.file" );
  m_filename_aux= par.getdefault<string>( prefix+"init_h.file_aux", "ABC");
  m_xline= par.getdefault( prefix+"init_h.x-line", false );
  m_zscale= par.getdefault( prefix+"init_h.zscale", 1.0 );
  m_hplain= par.getdefault( prefix+"init_h.plain.height", 0.0 );
}


/*!  Alternative constructor used for re-initialisation when a simulation is
    continued.  The parameters are given as arguments instead of being read
    from the parameter file.  Except for storing them in member variables, the
    constructor does nothing.  */

arrayinit_ascii::arrayinit_ascii(string filename, bool xline, double zscale) :
  m_xline(xline), m_zscale(zscale), m_filename(filename)
{}


/*!  Opens the file given by the parameter init_h.file.  If init_h.x-line is
  false, x is assumed the line index in the data table in the file, otherwise
  it is the column index (default).  */

void arrayinit_ascii::init_2d_scal(TFktScal& array)
{
  ifstream is( m_filename.c_str() );
  TFktScal filearray;
  int file_ncol, file_nline, col, line;
  char ch;
  
  if( !is ) {
    cerr << "ERROR: arrayinit_ascii::init_2d_scal():  Cannot open initialisation file `" << m_filename << "'!" << endl;
    exit(1);
  }
  
  file_ncol= 0;
  is.get(ch);
  while( is )
  {
    while( ch!='\n' && isspace(ch) && is.get(ch) );
    if( ch=='\n' )
      break;
    ++file_ncol;
    while( !isspace(ch) && is.get(ch) );
  }
  
  is.clear();
  is.seekg( 0 );
  file_nline= 0;
  while( is.get(ch) )
    if( ch=='\n' )
      ++file_nline;
  
  if( !file_ncol || !file_nline ) {
    cerr << "arrayinit_ascii::init_2d_scal(): File seems to have zero number of columns or lines. Aborting." << endl;
    is.close();
    exit(1);
  }
  
  is.clear();
  is.seekg( 0 );
  if( m_xline ) {
    if( file_ncol==array.SizeY() && file_nline==array.SizeX() )
      cerr << "arrayinit_ascii::init_2d_scal:  Warning: Size of array in file is size of simulation arrays with x and y swapped.  Did you set init_h.x-line correctly?" << endl;
    filearray.Create(file_ncol, file_nline, 1.0);
    for( line= 0; line < file_nline; ++line )
      for( col= 0; col < file_ncol; ++col )
	is >> filearray(col, line);
  }
  else {
    if( file_ncol==array.SizeX() && file_nline==array.SizeY() )
      cerr << "arrayinit_ascii::init_2d_scal:  Warning: Size of array in file is size of simulation arrays with x and y swapped.  Did you set init_h.x-line correctly?" << endl;
    filearray.Create(file_nline, file_ncol, 1.0);
    for( line= 0; line < file_nline; ++line )
      for( col= 0; col < file_ncol; ++col )
	is >> filearray(line, col);
  }
  
  if( !is ) {
    cerr << "arrayinit_ascii::init_2d_scal():  Error while reading initialisation file `" << m_filename << "'!" << endl;
    is.close();
    exit(1);
  }
  
  array.copyscale( filearray );
  
  for( int x= 0; x< array.SizeX(); ++x )
    for( int y= 0; y< array.SizeY(); ++y )
    { array(x, y) += m_hplain;
      array(x, y) *= m_zscale;
    }
}

void arrayinit_ascii::init_2d_vec(TFktVec& array)
{
    TFktScal auxfield;
    auxfield.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    
    init_2d_scal(auxfield);
    
    int x, y;
	for( x= 0; x< duneglobals::nx(); ++x ){
		for( y= 0; y< duneglobals::ny(); ++y ){
            array(x, y)[0] = auxfield(x, y);
        }
    }
    
    m_filename = m_filename_aux;
    
    init_2d_scal(auxfield);
    
	for( x= 0; x< duneglobals::nx(); ++x ){
		for( y= 0; y< duneglobals::ny(); ++y ){
            array(x, y)[1] = auxfield(x, y);
        }
    }
    
}

//*****************************************************************************
// class arrayinit_cube


arrayinit_cube::arrayinit_cube(const dunepar &p, string prefix)
{
  m_height= p.getrequired<double>(prefix + "cube.height");
  m_x_frac= p.getrequired<double>(prefix + "cube.x_frac");
  if( m_x_frac<= 0.0 || m_x_frac > 1.0 ) {
    cerr << "arrayinit_cube::arrayinit_cube: illegal value `" << m_x_frac 
	    << "' for x size fraction of cube!\n";
    exit(1);
  }
  m_y_frac= p.getrequired<double>(prefix + "cube.y_frac");
  if( m_y_frac<= 0.0 || m_y_frac > 1.0 ) {
    cerr << "arrayinit_cube::arrayinit_cube: illegal value `" << m_y_frac 
	    << "' for y size fraction of cube!\n";
    exit(1);
  }
  m_x0= p.getdefault<double>(prefix+"cube.x_0_percent",0.5);
  m_y0= p.getdefault<double>(prefix+"cube.y_0_percent",0.5);
}

void arrayinit_cube::init_2d_scal(TFktScal& array)
{
  int x, y, minx, miny, maxxp1, maxyp1;
  
  minx= (int)round(array.SizeX() * (m_x0 - 0.5 * m_x_frac));
  maxxp1= (int)round(array.SizeX() * 2. * m_x0) - minx;
  miny= (int)round(array.SizeY() * (m_y0 - 0.5 * m_y_frac));
  maxyp1= (int)round(array.SizeY() * 2. * m_y0) - miny;
  
  for( x= 0; x< minx; ++x )
    for( y= 0; y< array.SizeY(); ++y )
      array(x, y)= 0.0;
  for( ; x< maxxp1; ++x ) {
    for( y= 0; y< miny; ++y )
      array(x, y) = 0.0;
    for( ; y< maxyp1; ++y )
      array(x, y) = m_height;
    for( ; y< array.SizeY(); ++y )
      array(x, y) = 0.0;
  }
  for( ; x< array.SizeX(); ++x )
    for( y= 0; y< array.SizeY(); ++y )
      array(x, y)= 0.0;
}


//*****************************************************************************
//  class arrayinit_paraboloid

/*!  Gets the length, width, height and position of the parboloid from the
  parameter file.  */

arrayinit_paraboloid::arrayinit_paraboloid(const dunepar &p, string prefix)
{
  m_l= p.getrequired<double>( prefix+"paraboloid.l" );
  m_w= p.getrequired<double>( prefix+"paraboloid.w" );
  m_h= p.getrequired<double>( prefix+"paraboloid.h" );
  m_l2= 0.25*m_l*m_l/(duneglobals::dx()*duneglobals::dx());
  m_w2= 0.25*m_w*m_w/(duneglobals::dx()*duneglobals::dx());
  m_x0= p.getdefault<double>( prefix+"paraboloid.x0", duneglobals::nx()/2.0 )
	     / duneglobals::dx();;
  m_y0= p.getdefault<double>( prefix+"paraboloid.y0", duneglobals::ny()/2.0 )
	     / duneglobals::dx();;
}


/*!  Puts the upper half of a paraboloid with length m_l, width m_w and height
  m_h into the array.  Its centre is placed at m_x0, m_y0.  */

void arrayinit_paraboloid::init_2d_scal(TFktScal& array)
{
  double dist;
  int x, y;
  
  for( x= 0; x< duneglobals::nx(); ++x )
    for( y= 0; y< duneglobals::ny(); ++y )
      if( 1.0 > (dist= (x-m_x0)*(x-m_x0)/m_l2 + (y-m_y0)*(y-m_y0)/m_w2) )
	array(x, y)= m_h*sqrt(1.0-dist);
      else
        array(x, y)= 0.0;
}


//*****************************************************************************
//  class arrayinit_ridge

arrayinit_ridge::arrayinit_ridge(const dunepar &p, string prefix)
{
  m_ypos= p.getdefault( prefix+"ridge.ypos", 0.5*(double)duneglobals::ny()*duneglobals::dx() );
  m_sigma= p.getrequired<double>(prefix+"ridge.sigma");
  m_height= p.getrequired<double>(prefix+"ridge.height");
}

void arrayinit_ridge::init_2d_scal(TFktScal& array)
{
  double normarg, value;
  int x, y;

  for( y= 0; y< duneglobals::ny(); ++y ) {
    normarg= ((double)y*duneglobals::dx() - m_ypos)/(M_SQRT2*m_sigma);
    value= m_height * exp( -normarg*normarg );
    for( x= 0; x< duneglobals::nx(); ++x )
      array(x, y) = value;
  }
}
