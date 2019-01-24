#include "rfftw12d.h"

/******************************************************************************
    Class rfftw2d_array
******************************************************************************/

rfftw1d_array::rfftw1d_array( int size )
{
  m_size= size;
  m_padsize= 2*(size/2+1);
  m_buf= (double*)fftw_malloc(sizeof(double)*m_padsize);
  m_cbuf= (fftw_complex*)m_buf;
  
  m_plan_forw= fftw_plan_dft_r2c_1d( m_size, m_buf, m_cbuf, FFTW_MEASURE );
  m_plan_back= fftw_plan_dft_c2r_1d( m_size, m_cbuf, m_buf, FFTW_MEASURE );
}

rfftw1d_array::~rfftw1d_array()
{
  fftw_destroy_plan( m_plan_forw );
  fftw_destroy_plan( m_plan_back );
  fftw_free(m_buf);
}


/*!  Divides all values by the factor xsize*ysize by which they are too large
  after forward+backward transformation by fftw.  This is called by
  transform_back() automatically.  */

void rfftw1d_array::renormalize()
{
  int x;
  double scale= 1.0/(double)(m_size);
  
  for( x= 0; x< m_size; ++x )
    pos(x) *= scale;
}


/******************************************************************************
    Class rfftw2d_array
******************************************************************************/

rfftw2d_array::rfftw2d_array( int xsize, int ysize )
{
  m_xsize= xsize;
  m_ysize= ysize;
  m_clineoff= m_ysize/2 + 1;
  m_lineoff= 2*m_clineoff;
  m_buf= (double*)fftw_malloc(sizeof(double)*m_xsize*m_lineoff);
  m_cbuf= (fftw_complex*)m_buf;
  
  m_plan_forw= fftw_plan_dft_r2c_2d( m_xsize, m_ysize, m_buf, m_cbuf,
										FFTW_MEASURE );
  m_plan_back= fftw_plan_dft_c2r_2d( m_xsize, m_ysize, m_cbuf, m_buf,
										FFTW_MEASURE );
}

rfftw2d_array::~rfftw2d_array()
{
  fftw_destroy_plan( m_plan_forw );
  fftw_destroy_plan( m_plan_back );
  fftw_free(m_buf);
}


/*!
  Prints out all position-space values. y determines the column, x the row.
*/
void rfftw2d_array::dumppos(void)
{
  int x, y;
  
  for( y= 0; y< pos_ysize(); ++y ) {
    for( x= 0; x< pos_xsize(); ++x )
	  printf("%g  ", pos(x, y));
	printf("\n");
  }
}


/*!
  Prints out all frequency-space values as (real, imaginary) value pairs. 
  y determines the column, x the row.
*/
void rfftw2d_array::dumpfreq(void)
{
  int x, y;
  
  for( y= 0; y< freq_ysize(); ++y ) {
    for( x= 0; x< freq_xsize(); ++x )
	  printf( "(%g, %g)  ", freqre(x, y), freqim(x, y) );
	printf("\n");;
  }
}


/*!
  Divides all values by the factor xsize*ysize by which they are too large 
  after forward+backward transformation by fftw.  This is called by 
  transform_back() automatically.
*/
void rfftw2d_array::renormalize()
{
  int x, y;
  double scale= 1.0/(double)(m_xsize*m_ysize);
  
  for( x= 0; x< m_xsize; ++x )
    for( y= 0; y< m_lineoff; ++y )
      pos(x, y) *= scale;
}

//void rfftw2d_array::foreachpos( 
//	  void (*func)(double *elem, int x, int y) )
//{
//	int x, y;
//	
//	for( x= 0; x< m_xsize; ++x )
//	  for( y= 0; y< m_ysize; ++y )
//		func( &pos(x, y), x, y);
//}
//
//void rfftw2d_array::foreachfreq( 
//	  void (*func)(fftw_complex *elem, int fx, int fy) )
//{
//	int x, y;
//	
//	for( x= 0; x< m_xsize; ++x )
//	  for( y= 0; y< m_clineoff; ++y )
//		func( &freq(x, y), x, y);
//}


