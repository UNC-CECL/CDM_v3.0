#ifndef __VEC_H__
#define __VEC_H__

#include <math.h>
#include <iostream>

class vec2
{
private:
  double m_data[2];

public:
  vec2() {}
  vec2(double x) { m_data[0]= x; m_data[1]= x; }
  vec2(double x1, double x2) { m_data[0]= x1; m_data[1]= x2; }

  inline double& operator[](const int ind) { return ind? m_data[1]: m_data[0]; }
  inline const double& operator[](const int ind) const { return ind? m_data[1]: m_data[0]; }
  inline vec2 operator-() const { vec2 v(-m_data[0], -m_data[1]); return v; }
  inline vec2 operator+(const vec2 rhs) const { vec2 v(rhs[0]+m_data[0], rhs[1]+m_data[1]); return v; }
  inline vec2 operator-(const vec2 rhs) const { vec2 v(rhs[0]-m_data[0], rhs[1]-m_data[1]); return v; }
  inline vec2 operator*(const double rhs) const { vec2 v(m_data[0]*rhs, m_data[1]*rhs); return v; }
  inline vec2 operator/(const double rhs) const { vec2 v(m_data[0]/rhs, m_data[1]/rhs); return v; }
  inline vec2& operator=(const vec2 rhs) { m_data[0]=rhs[0]; m_data[1]=rhs[1]; return *this; }
  inline vec2& operator+=(const vec2 rhs) { m_data[0]+=rhs[0]; m_data[1]+=rhs[1]; return *this; }
  inline vec2& operator-=(const vec2 rhs) { m_data[0]-=rhs[0]; m_data[1]-=rhs[1]; return *this; }
  inline vec2& operator*=(const double rhs) { m_data[0]*=rhs; m_data[1]*=rhs; return *this; }
  inline vec2& operator/=(const double rhs) { m_data[0]/=rhs; m_data[1]/=rhs; return *this; }
  inline bool operator==(const vec2 rhs) const { return (*this)[0]==rhs[0] && (*this)[1]==rhs[1]; }
  inline bool operator!=(const vec2 rhs) const { return (*this)[0]!=rhs[0] || (*this)[1]!=rhs[1]; }
};


inline vec2 operator*(const double fac, const vec2 vec)
{
  vec2 v(fac*vec[0], fac*vec[1]);
  return v;
}

inline double vabs(vec2 v)
{
  return sqrt( v[0]*v[0] + v[1]*v[1] );
}

inline std::ostream& operator<<(std::ostream& os, const vec2& v)
{
  return os << "(" << v[0] << ", " << v[1] << ")";
}


#endif




