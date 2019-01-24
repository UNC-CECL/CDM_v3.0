#if !defined(__PTG_Func_h__)
#define __PTG_Func_h__

#include <stdlib.h>
#include <math.h>
#include <assert.h>

/**
* PTG_Func is a class providing basic numeric operations for 
* 1 dim. functions. The class also takes care of the memory,
* which is used for the data points. It is templated on the "value" 
* type T of the stored data. The data points are discrete, 
* their index represented by integer values. Where the interface refers to 
* values from the possibly continuous range of ordinate (x) values of the 
* function, a second template INDEX_T is used (by default T and INDEX_T are 
* the same type)
*
* IMPORTANT: All operations assume periodic boundary conditions.
*
* Check out PTG_Func and start make to compile the sample program.
*
* For drawing functions with the XT-Library see also XTFunction1d.
* 
*
* @short Basic operations for 1 dim. functions.
* @version $Id: PTG_Func.h,v 1.1 2004/09/30 10:05:05 schatz Exp $
* @author gerd
*/

template <class T, class INDEX_T=T>
class PTG_Func  
{
public:
  // ---- construction, destruction ----

  /**
    * Creates a copy of the function given as the argument. 
    */
  PTG_Func( const PTG_Func & f );

  /**
    * Creates an instance of a function without initializing anything.
    * To complete the creation of the function you must call
    * Create() or CreateAndInit() before using the function.
    *
    * This way of construction is useful if the function is an instance of a
    * class and at creation time of this class the size of the function
    * is not known.
    *
    * @see PTG_Func#CreateAndInit
    * @see PTG_Func#Create
    */
  PTG_Func();

  /**
    * Allocates memory for the function data points and sets the grid
    * size of the function. Call this function only if
    * the instance was created with the default constructor PTG_Func().
    *
    * @param iCount  Number of data points to allocate for the function.
    * @param tDx     Distance between two function data points. (grid size)
    * @param tMin    Value of the first grid point. (Not used, yet.)
    *
    * @see PTG_Func#PTG_Func
    * @see PTG_Func#CreateAndInit
    */
  inline void Create(int iCount, INDEX_T tDx=1, INDEX_T tMin=0); 

  /**
    * This function works like Create(), but you can specify a value
    * to initialize the function data points.
    *
    * tInit: Initial value of the function data points.
    *
    * @see PTG_Func#PTG_Func
    * @see PTG_Func#Create
    */
  inline void CreateAndInit(T tInit, int iCount, INDEX_T tDx=1, INDEX_T tMin=0); 

  /**
    * This contructor is equivalent to the default constructor
    * and a call of the Create() methode.
    *
    * @see PTG_Func#PTG_Func
    * @see PTG_Func#Create
    */
  explicit PTG_Func(int iCount, INDEX_T tDx=1, INDEX_T tMin=0); 
 
  /**
    * Frees the memory, which was allocated for the data points.
    */
  virtual ~PTG_Func();

  // attributes
  inline T* GetBuffer() const { return m_ptF; }
  inline int GetSize() const { return m_iCount; }

  inline const INDEX_T& GetDelta() const { return m_tDx; }
  inline void SetDelta(const INDEX_T& t) { m_tDx = t; }

  // function value access
  inline const T& operator[]( int i ) const {
    assert( i>=0 && i<m_iCount); return m_ptF[i];
  }
  inline T& operator[]( int i ) {
    assert( i>=0 && i<m_iCount); return m_ptF[i];
  }

  // interpolation
  inline T GetInterpolatedValue( const INDEX_T& x );

  // operators
  inline const PTG_Func& operator=( const PTG_Func& f );

  inline const PTG_Func& operator+=( const PTG_Func& f );
  inline PTG_Func operator+( const PTG_Func& f ) const;
  inline const PTG_Func& operator+=( const T& f );
  inline PTG_Func operator+( const T& f ) const;

  inline const PTG_Func& operator-=( const PTG_Func& f );
  inline PTG_Func operator-( const PTG_Func& f ) const;
  inline const PTG_Func& operator-=( const T& f );
  inline PTG_Func operator-( const T& f ) const;

  inline const PTG_Func& operator*=( const PTG_Func& f );
  inline PTG_Func operator*( const PTG_Func& f ) const;
  inline const PTG_Func& operator*=( const T& f );
  inline PTG_Func operator*( const T& f ) const;

  inline const PTG_Func& operator/=( const PTG_Func& f );
  inline PTG_Func operator/( const PTG_Func& f ) const;
  inline const PTG_Func& operator/=( const T& f );
  inline PTG_Func operator/( const T& f ) const;

  // operations
  inline void SetValue(T tValue);

  // extrema
  inline void GetExtrema( int& iXMin, T& tMin, 
			  int& iXMax, T& tMax ) const;

  inline void GetMinMax(T &tMin, T &tMax) const 
    { int i,n; GetExtrema(i, tMin, n, tMax); }

  inline void GetMax(int& i, T& tMax) const
    { int n; T t; GetExtrema(n,t, i,tMax); }
  inline int GetMaxPos() const 
    { T t; int i; GetMax(i,t); return i; }
  inline T GetMax() const 
    { T t; int i; GetMax(i,t); return t; }

  inline void GetMin(int& i, T& tMin) const
    { int n; T t; GetExtrema(i,tMin, n,t); }
  inline int GetMinPos() const 
    { T t; int i; GetMin(i,t); return i; }
  inline T GetMin() const 
    { T t; int i; GetMin(i,t); return t; }

  // find positions
  inline int FindPosRightGreaterThan( int iStart, const T& t ) const;
  inline int FindPosRightSmallerThan( int iStart, const T& t ) const;

  inline int FindPosLeftGreaterThan( int iStart, const T& t ) const;
  inline int FindPosLeftSmallerThan( int iStart, const T& t ) const;

  // shift / rotation
  inline void ShiftLeft(int iEntries = 1);
  inline void RotLeft(int iEntries = 1);

  // integration
  inline T Integrate() const { return Integrate(0, m_iCount-1); }// entire function
  inline T Integrate( int iMin, int iMax ) const;
  inline void Integrate( const PTG_Func &f );

  // derivatives 1. order
  inline T DLeft( int i ) const;
  inline T DRight( int i ) const;
  inline T DMid( int i ) const;
  inline T DUpWind( int i, double v ) const;
  inline void DLeft( const PTG_Func &f );
  inline void DRight( const PTG_Func &f );
  inline void DMid( const PTG_Func &f );
  inline void DUpWind( const PTG_Func &f, const PTG_Func &v );

  // derivatives 2. order
  inline T D2Mid( int i ) const;
  inline void D2Mid( const PTG_Func &f );

  // convective term: (u * d/dx) u
  inline void Convect( const PTG_Func &u );

  // averages
  inline T Average() const { return Average(0, m_iCount-1); }// entire function
  inline T Average( int iMin, int iMax ) const;

  inline void AverageRect( const PTG_Func &f, int i );
  // i=1: (f[i-1]+f[i]+f[i+1])/3, i=2: (f[i-2]+ .. +f[i+2])/5, i=3: (f[i-3]+ ..)/7, etc.

  // correlation
  inline T Corr( const PTG_Func &f, int i ) const;
  inline void Corr( const PTG_Func &f, const PTG_Func &g );

  //virtual void FFT
  //virtual void DFT

protected:
  T *m_ptF;
  int m_iCount;
  T m_tDx;
  T m_tMin;
};



////////////////////////////////////
// construction destruction
//

template <class T, class INDEX_T>
PTG_Func<T,INDEX_T>::PTG_Func( const PTG_Func<T,INDEX_T> & f ) :
  m_ptF(NULL)
{
  m_iCount = f.m_iCount;

  if (m_iCount > 0) {
    Create( f.m_iCount, f.m_tDx, f.m_tMin );
    for( int i=0; i < m_iCount; i++) {
      m_ptF[i] = f.m_ptF[i];
    }
  }
}

template <class T,class INDEX_T>
PTG_Func<T,INDEX_T>::PTG_Func() :
  m_ptF(NULL),
  m_iCount(0),
  m_tDx(1),
  m_tMin(0)
{
}

template <class T,class INDEX_T>
PTG_Func<T,INDEX_T>::PTG_Func( int iCount, INDEX_T tDx, INDEX_T tMin) :
  m_ptF(NULL)
{
  assert( iCount > 0 );
  Create( iCount, tDx, tMin);
}

template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::CreateAndInit(T tInit, int iCount, INDEX_T tDx, INDEX_T tMin)
{
  assert( iCount > 0 );
  Create( iCount, tDx, tMin );
  SetValue( tInit );
}

template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::Create( int iCount, INDEX_T tDx, INDEX_T tMin )
{
  m_iCount = iCount;
  m_tDx = tDx;
  m_tMin = tMin;

  // delete previously allocated memory
  if (m_ptF)
    delete [] m_ptF;

  assert(m_iCount>0);
  m_ptF = new T[ m_iCount ];
  if (m_ptF == NULL) {
    // can´t allocate memory
    m_iCount = 0;
    assert( false );
  }
}

template <class T, class INDEX_T>
PTG_Func<T,INDEX_T>::~PTG_Func()
{
  if (m_ptF)
    delete [] m_ptF;
}


///////////////////////////////////
// operators
//
template <class T, class INDEX_T>
inline T PTG_Func<T,INDEX_T>::GetInterpolatedValue( const INDEX_T& x )
{
  int i1 = (int) floor( (x - m_tMin)/m_tDx );
  assert( i1>=0 && i1<m_iCount );
  int i2 = (i1+1)%m_iCount;

  T dx = (x - m_tMin) - i1*m_tDx;
  
  return (m_ptF[i2]-m_ptF[i1])*dx + m_ptF[i1];
}



///////////////////////////////////
// operators
//

template <class T,class INDEX_T>
const PTG_Func<T,INDEX_T>& PTG_Func<T,INDEX_T>::operator=( const PTG_Func<T,INDEX_T>& f )
{
  assert(m_iCount == f.m_iCount);

  m_tDx = f.m_tDx;
  m_tMin = f.m_tMin;

  for( int i=0; i < m_iCount; i++) {
    m_ptF[i] = f.m_ptF[i];
  }

  return *this;
}

///////////////////////////

#define IMPL_OP_ASSIGN(OP) \
template <class T,class INDEX_T>\
inline const PTG_Func<T,INDEX_T>& PTG_Func<T,INDEX_T>::operator OP ( const PTG_Func<T,INDEX_T>& f )\
{\
  assert(m_iCount == f.m_iCount);\
  for( int i=0; i < m_iCount; i++) {\
    m_ptF[i] OP f.m_ptF[i];\
  }\
  return *this;\
}

IMPL_OP_ASSIGN(+=)
IMPL_OP_ASSIGN(-=)
IMPL_OP_ASSIGN(*=)
IMPL_OP_ASSIGN(/=)

#undef IMPL_OP_ASSIGN

////////////////////////////

#define IMPL_OP(OP) \
template <class T, class INDEX_T>\
inline PTG_Func<T,INDEX_T> PTG_Func<T,INDEX_T>::operator OP ( const PTG_Func<T,INDEX_T>& f ) const\
{\
  assert(m_iCount == f.m_iCount);\
  PTG_Func<T,INDEX_T> ftemp( m_iCount, m_tDx, m_tMin );\
  for( int i=0; i < m_iCount; i++) {\
    ftemp.m_ptF[i] = m_ptF[i] OP f.m_ptF[i];\
  }\
  return ftemp;\
}

IMPL_OP(+)
IMPL_OP(-)
IMPL_OP(*)
IMPL_OP(/)

#undef IMPL_OP

///////////////////////////

#define IMPL_OP_ASSIGN_VALUE(OP) \
template <class T,class INDEX_T>\
inline const PTG_Func<T,INDEX_T>& PTG_Func<T,INDEX_T>::operator OP ( const T& f )\
{\
  for( int i=0; i < m_iCount; i++) {\
    m_ptF[i] OP f;\
  }\
  return *this;\
}

IMPL_OP_ASSIGN_VALUE(+=)
IMPL_OP_ASSIGN_VALUE(-=)
IMPL_OP_ASSIGN_VALUE(*=)

#undef IMPL_OP_ASSIGN_VALUE

template <class T,class INDEX_T>
inline const PTG_Func<T,INDEX_T>& PTG_Func<T,INDEX_T>::operator/=( const T& f )
{
  return operator*=(1/f);
}


///////////////////////////

#define IMPL_OP_VALUE(OP) \
template <class T,class INDEX_T>\
inline PTG_Func<T,INDEX_T> PTG_Func<T,INDEX_T>::operator OP ( const T& f ) const\
{\
  PTG_Func<T,INDEX_T> ftemp( m_iCount, m_tDx, m_tMin );\
  for( int i=0; i < m_iCount; i++) {\
    ftemp.m_ptF[i] = m_ptF[i] OP f;\
  }\
  return ftemp;\
}

IMPL_OP_VALUE(+)
IMPL_OP_VALUE(-)
IMPL_OP_VALUE(*)

#undef IMPL_OP_VALUE

template <class T,class INDEX_T>
inline PTG_Func<T,INDEX_T> PTG_Func<T,INDEX_T>::operator/( const T& f ) const
{
  return operator*(1/f);
}

///////////////////////////


//////////////////////////////////
// operations
//

template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::SetValue(T tValue)
{
  assert( m_iCount > 0 );
  for( int i=0; i < m_iCount; i++) {
    m_ptF[i] = tValue;
  }
}


// ---- extreme values ----

template <class T,class INDEX_T>
inline void PTG_Func<T,INDEX_T>::GetExtrema( int& iXMin, T& tMin, int& iXMax, T& tMax ) const
{
  assert( m_iCount > 0 );

  tMin = tMax = m_ptF[0];
  iXMin = iXMax = 0;

  for( int i=1; i < m_iCount; i++) {
    if (m_ptF[i] < tMin) {
      tMin = m_ptF[i];
      iXMin = i;
    } else if (m_ptF[i] > tMax) {
      tMax = m_ptF[i];
      iXMax = i;
    }
  }
}


// ---- find positions ----
template <class T,class INDEX_T>
inline int PTG_Func<T,INDEX_T>::FindPosRightGreaterThan( int iStart, const T& t ) const
{
  assert( iStart>=0 && iStart < m_iCount );

  int i=iStart;
  for(; i < m_iCount; i++ ) {
    if (m_ptF[i] > t)
      return i;
  }
  for( i=0; i < iStart; i++ ) {
    if (m_ptF[i] > t)
      return i;
  }
}

template <class T,class INDEX_T>
inline int PTG_Func<T,INDEX_T>::FindPosRightSmallerThan( int iStart,  const T& t ) const
{
  assert( iStart>=0 && iStart < m_iCount );

  int i=iStart;
  for(; i < m_iCount; i++ ) {
    if (m_ptF[i] < t)
      return i;
  }
  for( i=0; i < iStart; i++ ) {
    if (m_ptF[i] < t)
      return i;
  }
}

template <class T,class INDEX_T>
inline int PTG_Func<T,INDEX_T>::FindPosLeftGreaterThan( int iStart,  const T& t ) const
{
  assert( iStart>=0 && iStart < m_iCount );

  int i=iStart;
  for(; i >= 0; i-- ) {
    if (m_ptF[i] > t)
      return i;
  }
  for( i=m_iCount-1; i>iStart; i-- ) {
    if (m_ptF[i] > t)
      return i;
  }
}

template <class T,class INDEX_T>
inline int PTG_Func<T,INDEX_T>::FindPosLeftSmallerThan( int iStart,  const T& t ) const
{
  assert( iStart>=0 && iStart < m_iCount );

  int i=iStart;
  for(; i >= 0; i-- ) {
    if (m_ptF[i] < t)
      return i;
  }
  for( i=m_iCount-1; i>iStart; i-- ) {
    if (m_ptF[i] < t)
      return i;
  }
}


template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::ShiftLeft(int iEntries)
{
  assert( iEntries < m_iCount );
  int i=0;
  for(; i < m_iCount-iEntries; i++ ) {
    m_ptF[i] = m_ptF[i+iEntries];
  }
  for( i = m_iCount-iEntries; i < m_iCount; i++ ) {
    m_ptF[i] = 0;
  }
}

template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::RotLeft(int iEntries)
{
  assert( iEntries < m_iCount );
  
  PTG_Func<T,INDEX_T> temp(*this);

  int i=0;
  for(; i < m_iCount-iEntries; i++ ) {
    m_ptF[i] = temp[i+iEntries];
  }
  for( i = m_iCount-iEntries; i < m_iCount; i++ ) {
    m_ptF[i] = temp[i+iEntries-m_iCount];
  }
}


template <class T,class INDEX_T>
T PTG_Func<T,INDEX_T>::Integrate( int iMin, int iMax ) const
{
  T tSum = 0;
  for( int i=iMin; i <= iMax; i++) {
    tSum += m_ptF[i];
  }
  return tSum * m_tDx;
}


template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::Integrate( const PTG_Func<T,INDEX_T> &f )
{
  T tSum = 0;
  for( int i=0; i < m_iCount; i++) {
    tSum += f.m_ptF[i];
    m_ptF[i] = tSum * m_tDx;
  }
}


template <class T,class INDEX_T>
T PTG_Func<T,INDEX_T>::DLeft( int i ) const
{
  assert( i>=0 && i<m_iCount);

  if (i==0)
    return (m_ptF[0]-m_ptF[m_iCount-1]) / m_tDx;
  else
    return (m_ptF[i]-m_ptF[i-1]) / m_tDx;
}

template <class T,class INDEX_T>
T PTG_Func<T,INDEX_T>::DRight( int i ) const
{
  assert( i>=0 && i<m_iCount);

  if (i==m_iCount-1)
    return (m_ptF[0]-m_ptF[m_iCount-1]) / m_tDx;
  else
    return (m_ptF[i+1]-m_ptF[i]) / m_tDx; 
}

template <class T,class INDEX_T>
T PTG_Func<T,INDEX_T>::DMid( int i ) const
{
  assert( i>=0 && i<m_iCount);

  if (i==m_iCount-1)
    return (m_ptF[0]-m_ptF[m_iCount-2]) / (2*m_tDx);
  else if (i==0)
    return (m_ptF[1]-m_ptF[m_iCount-1]) / (2*m_tDx);
  else
    return (m_ptF[i+1]-m_ptF[i-1]) / (2*m_tDx); 
}

template <class T,class INDEX_T>
T PTG_Func<T,INDEX_T>::DUpWind( int i, double v ) const
{
  assert( i>=0 && i<m_iCount);

  /*if (u[x] > 0)
    return DLeft(x);
  else
    return DRight(x);*/
}


template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::DLeft( const PTG_Func<T,INDEX_T> &f )
{
  assert(m_iCount == f.m_iCount);
  double d = 1/m_tDx;
  m_ptF[0] = (f.m_ptF[0] - f.m_ptF[m_iCount-1]) * d;
  for (int i=1; i<m_iCount; i++)
    m_ptF[i] = (f.m_ptF[i] - f.m_ptF[i-1]) * d;
}

template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::DRight( const PTG_Func<T,INDEX_T> &f )
{
  assert(m_iCount == f.m_iCount);
  double d = 1/m_tDx;
  m_ptF[m_iCount-1] = (f.m_ptF[0] - f.m_ptF[m_iCount-1]) * d;
  for (int i=0; i<m_iCount-1; i++)
    m_ptF[i] = (f.m_ptF[i+1] - f.m_ptF[i]) * d;
}

template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::DMid( const PTG_Func<T,INDEX_T> &f )
{
  assert(m_iCount == f.m_iCount);
  double d = 1/(2*m_tDx);
  m_ptF[m_iCount-1] = (f.m_ptF[0] - f.m_ptF[m_iCount-2]) * d;
  m_ptF[0]          = (f.m_ptF[1] - f.m_ptF[m_iCount-1]) * d;
  for (int i=1; i<m_iCount-1; i++)
    m_ptF[i] = (f.m_ptF[i+1] - f.m_ptF[i-1]) * d;
}

template <class T,class INDEX_T>
T PTG_Func<T,INDEX_T>::D2Mid( int i ) const
{
  assert( i>=0 && i<m_iCount);

  if (i==m_iCount-1)
    return ( (m_ptF[0]          - m_ptF[m_iCount-1])
            -(m_ptF[m_iCount-1] - m_ptF[m_iCount-2])
	   ) / (m_tDx*m_tDx);
  else if (i==0)
    return ( (m_ptF[1]          - m_ptF[0])
            -(m_ptF[0] - m_ptF[m_iCount-1])
	   ) / (m_tDx*m_tDx);
  else
    return ( (m_ptF[i+1] - m_ptF[i])
            -(m_ptF[i]   - m_ptF[i-1])
	   ) / (m_tDx*m_tDx);
}

template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::DUpWind( const PTG_Func<T,INDEX_T> &f, const PTG_Func<T,INDEX_T> &v )
{
  for (int i=0; i<m_iCount; i++) {
    if (v[i] > 0)
      m_ptF[i] = f.DLeft(i);
    else
      m_ptF[i] = f.DRight(i);
  }
}


template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::D2Mid( const PTG_Func<T,INDEX_T> &f )
{
  assert(m_iCount == f.m_iCount);
  double d2 = 1 / (m_tDx*m_tDx);

  m_ptF[m_iCount-1] = ( (f.m_ptF[0]          - f.m_ptF[m_iCount-1])
                       -(f.m_ptF[m_iCount-1] - f.m_ptF[m_iCount-2])
	              ) * d2;

  m_ptF[0] = ( (f.m_ptF[1] - f.m_ptF[0])
              -(f.m_ptF[0] - f.m_ptF[m_iCount-1])
	     ) * d2;

  for (int i=1; i<m_iCount-1; i++)
    m_ptF[i] = ( (f.m_ptF[i+1] - f.m_ptF[i])
                -(f.m_ptF[i]   - f.m_ptF[i-1])
	       ) * d2;
}


template <class T,class INDEX_T>
void PTG_Func<T,INDEX_T>::Convect( const PTG_Func<T,INDEX_T> &u )
{
  DUpWind(u,u);
  (*this) *= u;
}


// ---- Average ----

template <class T,class INDEX_T>
T PTG_Func<T,INDEX_T>::Average( int iMin, int iMax ) const
{
  assert( iMin >= iMax );
  T tSum = 0;
  for( int i = iMin; i <= iMax; i++) {
    tSum += m_ptF[i];
  }
  return tSum / (iMax-iMin+1);
}

template <class T, class INDEX_T>
void PTG_Func<T,INDEX_T>::AverageRect( const PTG_Func<T,INDEX_T> &f, int i )
{
  double d;
  double dNorm = 1/(double)(2*i+1);
  for (int x=0; x<m_iCount; x++) {
    d = 0;
    for (int n=-i; n<=i; n++) {
      d += f[(x+n+m_iCount)%m_iCount];
    }
    m_ptF[x] = d * dNorm;
  }
}


// ---- Correlation ----

template <class T,class INDEX_T>
inline T PTG_Func<T,INDEX_T>::Corr( const PTG_Func<T,INDEX_T> &f, int i ) const
{
  T tSum = 0;
  int n = m_iCount-1;
  for (; n >= m_iCount-i; n--) {
    tSum += m_ptF[n] * f.m_ptF[n+i-m_iCount];
  }
  for (; n >= 0; n--) {
    tSum += m_ptF[n] * f.m_ptF[n+i];
  }
  return tSum*m_tDx;
}



template <class T,class INDEX_T>
inline void PTG_Func<T,INDEX_T>::Corr( const PTG_Func<T,INDEX_T> &f, const PTG_Func<T,INDEX_T> &g )
{
  for (int i=m_iCount-1; i>=0; i--)
    m_ptF[i] = f.Corr(g, i);
}



#endif
