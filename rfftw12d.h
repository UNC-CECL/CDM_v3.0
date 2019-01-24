#include <fftw3.h>

#ifndef RFFTW12D_H
#define RFFTW12D_H


/*!  Contains just one function to obtain the next largest power of two (for
  the right array size for efficient FFT).  */

namespace fft
{
  inline int GetNextPowerOf2(int iN)
  {
    int iFFT = 1;
    --iN;
    while( iN> 0 ) {
      iFFT <<= 1;
      iN >>= 1;
    }
    return iFFT;
  }
}


/*!  \brief Array for in-place 1-dimensional real fftw transforms, double
  precision.  */

class rfftw1d_array 
{
public:
  rfftw1d_array( int size );
  ~rfftw1d_array();

  //! Access array element in position space (no range check!)
  double &pos( int x ) { return m_buf[x]; }
  
  //! Access array element in frequency space (no range check!)
  fftw_complex &freq( int fx ) { return m_cbuf[fx]; }
  
  //! Real part of array element in frequency space (no range check!)
  double &freqre( int fx ) { return m_cbuf[fx][0]; }
  
  //! Imaginary part of array element in frequency space (no range check!)
  double &freqim( int fx ) { return m_cbuf[fx][1]; }
  
  //! Size of position array 
  int pos_size() const { return m_size; }

  //! Size of frequency array 
  int freq_size() const { return m_padsize/2; }

  //! Data buffer for use with fftw functions (position space type: real*)
  double *buf() { return m_buf; }

  //! Data buffer for use with fftw functions (frequency space type: complex*)
  fftw_complex *cbuf() { return m_cbuf; }
  
  //! Forward Fourier transformation: real position space to complex frequency space.
  void transform_forw() { fftw_execute(m_plan_forw); }
  
  //! Backward Fourier transformation: complex frequency space to real position space.
  void transform_back() { fftw_execute(m_plan_back);  renormalize(); }
  
  void dumppos(void);
  void dumpfreq(void);
  
private:

  void renormalize(void);
  
  int m_size, m_padsize;
  double *m_buf;
  fftw_complex *m_cbuf;
  fftw_plan m_plan_forw, m_plan_back;
};


/*!  \brief Array for in-place 2-dimensional real fftw transforms, double
  precision.  */

class rfftw2d_array 
{
public:
  rfftw2d_array( int xsize, int ysize );
  ~rfftw2d_array();

  //! Access array element in position space (no range check!)
  double &pos( int x, int y ) { return m_buf[x*m_lineoff + y]; }
  
  //! Access array element in frequency space (no range check!)
  fftw_complex &freq( int fx, int fy ) { return m_cbuf[fx*m_clineoff + fy]; }
  
  //! Real part of array element in frequency space (no range check!)
  double &freqre( int fx, int fy ) { return m_cbuf[fx*m_clineoff + fy][0]; }
  
  //! Imaginary part of array element in frequency space (no range check!)
  double &freqim( int fx, int fy ) { return m_cbuf[fx*m_clineoff + fy][1]; }
  
  //! Size of position array in x direction
  int pos_xsize() const { return m_xsize; }

  //! Size of position array in y direction
  int pos_ysize() const { return m_ysize; }

  //! Size of frequency array in x direction
  int freq_xsize() const { return m_xsize; }

  //! Size of frequency array in y direction
  int freq_ysize() const { return m_clineoff; }

  //! Data buffer for use with fftw functions (position space type: real*)
  double *buf() { return m_buf; }

  //! Data buffer for use with fftw functions (frequency space type: complex*)
  fftw_complex *cbuf() { return m_cbuf; }
  
  //! Forward Fourier transformation: real position space to complex frequency space.
  void transform_forw() { fftw_execute(m_plan_forw); }
  
  //! Backward Fourier transformation: complex frequency space to real position space.
  void transform_back() { fftw_execute(m_plan_back);  renormalize(); }
  
  void dumppos(void);
  void dumpfreq(void);
  
private:

  void renormalize(void);
  
  int m_xsize, m_ysize, m_lineoff, m_clineoff;
  double *m_buf;
  fftw_complex *m_cbuf;
  fftw_plan m_plan_forw, m_plan_back;
};




#endif //  RFFTW12D_H

