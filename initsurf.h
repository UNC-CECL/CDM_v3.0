#ifndef __INITSURF_H__
#define __INITSURF_H__

#include <string>
using std::string;

#include "func.h"
#include "globals.h"

/*!  Base class for intialisation of one- or two-dimensional arrays, usually
 the dune height profile h.  */

class arrayinit
{
public:
    arrayinit() {}
    virtual ~arrayinit() {}
    
    virtual void init_1d_scal(CFunc1d& array);
    virtual void init_2d_scal(TFktScal& array);
    virtual void init_2d_vec(TFktVec& array);
    
    static arrayinit *create(const dunepar& par, string prefix= "");
    
private:
    void array_mismatch(const char *function);
};


/*!  Constant initialisation class for scalar arrays (both one- and
 two-dimensional).  */

class CInitSurfPlain : public arrayinit
{
public:
    CInitSurfPlain(const dunepar& P, string prefix= "") {
        m_const= P.getdefault(prefix+"plain.Height", 0.0);
        m_const_aux= P.getdefault(prefix+"plain.Height_aux", 0.0);
    }
    
    CInitSurfPlain(double initval) { m_const= initval; }
    
    virtual ~CInitSurfPlain() {}
    
    virtual void init_1d_scal(CFunc1d& array);
    virtual void init_2d_scal(TFktScal& array);
    virtual void init_2d_vec(TFktVec& array);
    
private:
    double m_const, m_const_aux;
};


/*!  Initialisation from an ASCII file.  For instance for continuing
 interrupted simulations.  */

class arrayinit_ascii : public arrayinit
{
public:
    arrayinit_ascii(const dunepar& par, string prefix= "");
    arrayinit_ascii(string filename, bool xline, double zscale);
    ~arrayinit_ascii() {}
    
    virtual void init_2d_scal(TFktScal& array);
    virtual void init_2d_vec(TFktVec& array);

private:
    /*!  True if x is the column index in the file (that is, x varies along each
     line).  "init_h.x-line" in parameter file.  */
    bool m_xline;
    /*!  Scale factor with which every value from the file is multiplied.
     "init_h.zscale" in parameter file.  */
    double m_zscale, m_hplain;
    /*!  Name of file to read.  */
    string m_filename, m_filename_aux;
};


/*!  Puts a cube in the middle of the height profile.  This was created just
 for testing the array rotation.  The cube's heigth is given by the parameter
 cube.height, and the fraction of the array to be taken up is given by
 cube.x_frac and cube.y_frac.  */

class arrayinit_cube : public arrayinit
{
public:
    arrayinit_cube(const dunepar &p, string prefix= "");
    virtual ~arrayinit_cube() {}
    
    virtual void init_2d_scal(TFktScal& array);
    
private:
    /*!  Height of the cube  */
    double m_height;
    /*!  Fraction of the array to be taken up by the centred cube in x and y
     direction.  */
    double m_x_frac, m_y_frac;
    /*!  Fraction of the array to center the cube in x and y
     direction.  */
    double m_x0, m_y0;
};


/*!  Puts the upper half of a paraboloid into the array.  */

class arrayinit_paraboloid : public arrayinit
{
public:
    arrayinit_paraboloid(const dunepar &p, string prefix= "");
    virtual ~arrayinit_paraboloid() {}
    
    virtual void init_2d_scal(TFktScal& array);
    
private:
    /*!  Length (in x direction), width and height of the paraboloid.  (In
     metres.)  */
    double m_l, m_w, m_h;
    /*!  Square of half of length and width.  In units of the grid spacing.  */
    double m_l2, m_w2;
    /*!  Position of the centre of the parboloid in units of the grid spacing.  */
    double m_x0, m_y0;
};


/*!  Initialisation of height profile with a ridge running in wind direction at
 a specified y position and a specified width.  */

class arrayinit_ridge : public arrayinit
{
public:
    arrayinit_ridge(const dunepar &p, string prefix= "");
    virtual ~arrayinit_ridge() {}
    
    virtual void init_2d_scal(TFktScal& array);
    
private:
    /*!  Transverse (y) position of the ridge.  */
    double m_ypos;
    /*!  RMS of Gaussian ridge profile in y direction.  */
    double m_sigma;
    /*!  Height of the crest of the ridge.  */
    double m_height;
};


#endif


