#include <string>

#include "analyze_new.h"

using std::string;

//////////////////////////////////////
// helper classes

struct CMinMax
{
  int m_iMin;
  int m_iMax;
  double m_dMin;
    double m_dMax;
  double m_dh;

  CMinMax(int iMax, double dh) :
	m_iMin(iMax), m_iMax(0),
	m_dMin(iMax), m_dMax(0),
	m_dh(dh)
  {}


  double Interpol(double h0, double h1, double href) {
	if (h1-h0 > 1e-10) {
	  double retval;
	  retval= (href-h0)/(h1-h0);
	  if( retval < 0 )
		retval= 0;
	  return (retval> 1.0? 1.0 : retval);
	}
	else
	  return 0.5;
  }


  void operator()(int i, double dl, double d, double dr)
  {
	if (i<=m_iMin) {
	  m_iMin = i;
	  if( m_dMin > i-1+Interpol(dl,d,m_dh) )
		m_dMin = i-1+Interpol(dl,d,m_dh);
	}

	if (i>=m_iMax) {
	  m_iMax = i;
	  if( m_iMax < i+Interpol(d,dr,m_dh) )
		m_dMax = i+Interpol(d,dr,m_dh);
	}
  }

  void Check(const char* pszName)
  {
	if (m_dMin > m_dMax) {
	  double d = m_dMin;
	  m_dMin = m_dMax;
	  m_dMax = d;
	  cout << pszName << ": Min/max detection failed!" << endl;
	}
  }

  int Center() {
	return (m_iMin+m_iMax)/2;
  }

};


//////////////////////////////////////
// analyze
//

analyze::analyze(const dunepar& P)
{
  m_zmin = P.getdefault("veget.zmin", 0.0);

  std::string strFile = P.getdefault<string>("save.dir", "./");
  if (strFile[strFile.length()-1] != '/') {
	strFile += '/';
  }
  strFile += "time.dat";
  m_os.open(strFile.c_str());

  if (!m_os) {
	cout << "Open file \"" << strFile << "\" failed!" << endl;
	exit(1);
  }

  m_os 	<< "%# 1: iterations \n"
		<< "%# 2: time in yr \n"
		<< "%# 3: maximum height \n"
		<< "%# 4: maximum cover	\n"
		<< "%# 5: volume / mass of sand \n"
		<< "%# 6: distance traveled by the dune in X \n"
		<< "%# 7: dune in flux \n"
		<< "%# 8: dune out flux \n"
		<< "%# 9: surge above MHWL \n"
	   << endl;
}


void analyze::Calc(int t, double time, double shift_dist_x, int m_shoreline, double m_shorelinechange, int m_veget_X0,
	double qin, double qout, double m_ustar0, double surge, const TFktScal& m_h, const TFktScal& m_rhoveg)
{
	//m_dTime += timestep;

	// calc dune properties
	const double dHMax = m_h.GetMax();
	const double RhoMax = m_rhoveg.GetMax();
    
    double H1X1 = m_h.GetFirstMax();
    double X1 = floor(H1X1/10000);
    double H1 = H1X1 - X1*10000;

    double H0 = 0;
    double veget0 = 0;
    for (int y=0; y<m_h.SizeY(); y++) {
        H0 += 0.5*(m_h(m_shoreline+m_veget_X0+4,y)+m_h(m_shoreline+m_veget_X0+5,y));
        veget0 += m_rhoveg(m_shoreline+m_veget_X0+1,y);
    }
	H0 /= m_h.SizeY();
	veget0 /= m_h.SizeY();

	// ---- volume / mass ----

	const double dVol= m_h.Integrate(0);

	// ----- max h ----
	double hmax = 0, hmaxc = 0, hmaxd = 0;
	int yc = 0;
	double hmaxY = 0;
	for (int y=0; y<m_h.SizeY(); y++) {
		hmaxY = 0;
		for (int x=0; x<m_h.SizeX(); x++) {
			if (hmaxY < m_h(x,y)) {
				hmaxY = m_h(x,y);
			}
		}
		// include only regions < Hc
		if (hmaxY < duneglobals::HMWL() + m_zmin)
		{
			hmaxc += hmaxY;
			yc++;
		} else {
			hmaxd += hmaxY;
		}
		hmax+=hmaxY;
	}
	hmax /= m_h.SizeY();
	hmaxc /= (yc > 0 ? yc : 1);
	hmaxd /= m_h.SizeY()-yc;
    
    // --- TIME ----
    double realtime = time/duneglobals::secyear()/duneglobals::timefrac(); // years
    // --- Shoreline rate of change ---
    double SR = m_shorelinechange; // / realtime;

	// ---- write data ----

	m_os << t << " "                    								// 1: iterations
		<< time/duneglobals::secyear()/duneglobals::timefrac() << " " 	// 2: time in yr
		<< dHMax-duneglobals::HMWL() << " "                				// 3: maximum height
		<< RhoMax << " " 												// 4: maximum cover	
		<< dVol << " "                 									// 5: volume / mass of sand
		<< SR << " "	      											// 6: shoreline change (m)
		<< qin << " "            										// 7: dune in flux
		<< qout << " "             										// 8: dune out flux
		<< endl;
}


double analyze::GetMaxPos(double f0, double f1, double f2)
{
  // assert that: f0 < f1 and f1 > f2
  if (f0 > f1 || f2 > f1)
	return 0;

  // calc derivatives: d1 > 0, d2 < 0
  double d1 = f1-f0;
  double d2 = f2-f1;

  // calc zero pos.  d1 + epsilon * (d2-d1) = 0
  double epsilon;
  if (d1 != d2) {
	epsilon = d1 / (d1-d2);
  } else {
	epsilon = 0.5;
  }
  return epsilon - 0.5;
}
