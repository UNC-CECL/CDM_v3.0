#include "save.h"



CSaveFields::CSaveFields(const dunepar& P) :
  m_FileName(0,5)
{
  m_bOff = P.getdefault("save.off", false);
  m_iEveryT = P.getdefault("save.every", 100);
  m_iT = 0;
  m_bXLine = P.getdefault("save.x-line", false);

  m_strDir = P.getdefault<string>("save.dir", "./");
  if (m_strDir[m_strDir.length()-1] != '/') {
    m_strDir += '/';
  }
}



bool CSaveFields::Save(const TFktScal& f, const string& strName)
{
 
  return DoSave(f, strName);
}


class CFieldComp
{
  const TFktVec& m_f;
  int m_i;
public:
  CFieldComp(const TFktVec& f, int i) :
    m_f(f), m_i(i)
  {}

  double operator()(int x, int y) const { return m_f(x,y)[m_i]; }
  double SizeX() const { return m_f.SizeX(); }
  double SizeY() const { return m_f.SizeY(); }
};


bool CSaveFields::Save(const TFktVec& f, const string& strName)
{
  return DoSave(CFieldComp(f,0), strName+"_x")
    && DoSave(CFieldComp(f,1), strName+"_y");
}


template<class T>
bool CSaveFields::DoSave(const T& f, const string& strName)
{
  bool bSpecial = false;
  for (int i=0; i<(int) m_veciSpecialT.size(); i++) {
    if (m_iT == m_veciSpecialT[i]) {
      bSpecial = true;
      break;
    } 
  }

  if (!(m_bForce || bSpecial)) {
    if (m_bOff || m_iT % m_iEveryT != 0)
      return true;
  }

  string strFileName = m_FileName.GetName(m_strDir + strName + ".", m_iT);

  ofstream os( strFileName.c_str() );
  if (! os) {
    cerr << "CSaveFields::DoSave: Creating file " << strFileName << " failed!\n";
    return false;
  }

  cout << " Save(" << strFileName << ")" << flush;

  const int iNx = (int) f.SizeX();
  const int iNy = (int) f.SizeY();



  if (m_bXLine) {
    for (int y=0; y < iNy; y++) {
      for (int x=0; x < iNx-1; x++) {
	os << f(x,y) << " ";
      }
      os << f(iNx-1,y) << "\n";
    }  
  } else {
    for (int x=0; x < iNx; x++) {
      for (int y=0; y < iNy-1; y++) {
	os << f(x,y) << " ";
      }
      os << f(x,iNy-1) << "\n";
    }  
  }

  os.close();
  return true;
}

