#ifndef __PTG_FUNC2D_H__
#define __PTG_FUNC2D_H__

//#include <iostream>

template<class T>
class PTG_Func2d
{
public:
  // ---- construction of a new array ----

  PTG_Func2d(int x, int y, double dDelta);

  // initialize each element with value
  PTG_Func2d(int x, int y, double dDelta, const T& Value);

  PTG_Func2d();
  void Create(int x, int y, double dDelta);
  void Create(int x, int y, double dDelta, const T& Value);

  // ---- destruction ----

  ~PTG_Func2d();

  // ---- copy-constructor ----

  PTG_Func2d(const PTG_Func2d& Source);

  // ---- assignement operator ----

  const PTG_Func2d& operator=(const PTG_Func2d& Source);

  const PTG_Func2d& operator*=(const PTG_Func2d& Source);
  const PTG_Func2d& operator*=(const double Value);

  const PTG_Func2d& operator+=(const PTG_Func2d& Source);
  const PTG_Func2d& operator-=(const PTG_Func2d& Source);

  // ---- f(x,y) => fast access ----

  const T& operator()(int x, int y) const {
    return m_pArr[ x*m_iNY + y ];
  }
  T& operator()(int x, int y) {
//    if( x<0 || y<0 || x>=SizeX() || y>=SizeY() )
//	std::cerr << "PTG_Func2d::operator(): index out of range: (" << x << ", " << y << "\n";
    return m_pArr[ x*m_iNY + y ];
  }

  // ---- f[x][y] => C array style compatibility ----

  class YIndex {
  public:
    YIndex( T* pYArr ) : m_pYArr( pYArr ) {}
    inline const T& operator[]( int y ) const { return m_pYArr[ y ]; }
    inline T& operator[]( int y ) { return m_pYArr[ y ]; }
  private:
    T* m_pYArr;
  };

  const YIndex operator[]( int x ) const {
    return YIndex(m_pArr + x*m_iNY);
  }

  YIndex operator[](int x) {
    return YIndex(m_pArr + x*m_iNY);
  }

  // ---- attributes ----
  int SizeX() const { return m_iNX; }
  int SizeY() const { return m_iNY; }
  
  double Delta() const { return m_dDelta; }

  // ---- operations ----
  void SetAll(const T& Value);
//  void WriteAscii(std::ostream& os) const;

private:

  // ---- internal helper functions ----
  void Copy(const PTG_Func2d& Source);

protected:

  // ---- internal operators ----

  inline double opDUpWind(const double l, const double m, const double r,
			  const double d_rez, const double u)
  {
    return (u>0.0) ? (m-l)*d_rez : (r-m)*d_rez;
  }

  inline double opDMid(const double l, const double r, const double d_rez)
  {
    return (r-l)*d_rez;
  }

  inline double opDRight(const double l, const double r, const double d_rez)
  {
    return (r-l)*d_rez;
  }

  inline double opDLeft(const double l, const double r, const double d_rez)
  {
    return (r-l)*d_rez;
  }

  inline double opD2Mid(const double l, const double m, 
			const double r, const double d_rez2)
  {
    return (r-2*m+l)*d_rez2;
  }

private:
  int m_iNX;
  int m_iNY;
  T* m_pArr;
  double m_dDelta;
};


//template<class T>
//std::ostream& operator<<(std::ostream& os, const PTG_Func2d<T>& f) {
//  f.WriteAscii(os);
//  return os;
//}



///////////////////////////////////////////////////
// implementation
//

// ---- construction ----

template<class T>
PTG_Func2d<T>::PTG_Func2d(int x, int y, double dDelta) :
  m_iNX(x), m_iNY(y), m_dDelta(dDelta)
{
  m_pArr = new T[m_iNX * m_iNY];
}

template<class T>
PTG_Func2d<T>::PTG_Func2d(int x, int y, double dDelta, const T& Value) :
  m_iNX(x), m_iNY(y), m_dDelta(dDelta)
{
  m_pArr = new T[m_iNX * m_iNY];
  SetAll( Value );
}


template<class T>
PTG_Func2d<T>::PTG_Func2d() :
  m_iNX(0), m_iNY(0), m_pArr(0), m_dDelta(0.0)
{
}


template<class T>
void PTG_Func2d<T>::Create(int x, int y, double dDelta)
{
  m_iNX = x; 
  m_iNY = y;
  m_dDelta = dDelta;
  m_pArr = new T[m_iNX * m_iNY];
}


template<class T>
void PTG_Func2d<T>::Create(int x, int y, double dDelta, const T& Value)
{
  Create(x,y, dDelta);
  SetAll( Value );
}


// ---- copy construction ----

template<class T>
PTG_Func2d<T>::PTG_Func2d(const PTG_Func2d& Source) :
  m_iNX(Source.m_iNX), m_iNY(Source.m_iNY), m_dDelta(Source.m_dDelta)
{ 
  m_pArr = new T[m_iNX * m_iNY];
  Copy(Source);
}


// ---- assignement operator ----

template<class T>
const PTG_Func2d<T>& PTG_Func2d<T>::operator=(const PTG_Func2d<T>& Source) { 
  Copy(Source);
  return *this;
}


// ---- destruction ----

template<class T>
PTG_Func2d<T>::~PTG_Func2d()
{
  delete[] m_pArr;
}


// ---- operators *= ----

template<class T>
const PTG_Func2d<T>& PTG_Func2d<T>::operator*=(const PTG_Func2d& Source) {
  for (int i = m_iNX*m_iNY - 1; i >= 0; i--) {
    m_pArr[i] *= Source.m_pArr[i];
  }
  return *this;
}

template<class T>
const PTG_Func2d<T>& PTG_Func2d<T>::operator*=(const double Value) {
  for (int i = m_iNX*m_iNY - 1; i >= 0; i--) {
    m_pArr[i] *= Value;
  }
  return *this;
}


// ---- operators += ----

template<class T>
const PTG_Func2d<T>& PTG_Func2d<T>::operator+=(const PTG_Func2d& Source) {
  for (int i = m_iNX*m_iNY - 1; i >= 0; i--) {
    m_pArr[i] += Source.m_pArr[i];
  }
  return *this;
}


// ---- operators -= ----

template<class T>
const PTG_Func2d<T>& PTG_Func2d<T>::operator-=(const PTG_Func2d& Source) {
  for (int i = m_iNX*m_iNY - 1; i >= 0; i--) {
    m_pArr[i] -= Source.m_pArr[i];
  }
  return *this;
}



// ---- operations ----

template<class T>
void PTG_Func2d<T>::SetAll(const T& Value) {
  for (int i = m_iNX*m_iNY - 1; i >= 0; i--) {
    m_pArr[i] = Value;
  }
}

template<class T>
void PTG_Func2d<T>::Copy(const PTG_Func2d<T>& Source) {
  for (int i = m_iNX*m_iNY - 1; i >= 0; i--) {
    m_pArr[i] = Source.m_pArr[i];
  }
}


#endif



