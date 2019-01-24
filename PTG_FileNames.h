#ifndef __PTG_FILENAMES_H__
#define __PTG_FILENAMES_H__

#include <iostream>
#include <sstream>
#include <iomanip>

namespace PTG {

using std::string;

class FileNameCounter
{
public:
  // ---- construction ----
  FileNameCounter(int iInit = 0,
		  int iDigits = 3,
		  const string & strPost = ".dat",
		  const string & strPre = "unknown.") :
    m_iValue(iInit),
    m_iDigits(iDigits),
    m_strPost(strPost),
    m_strPre(strPre)
  {}

  // ---- GetName ----
  string GetName(const string & strPre, int iValue, const string & strPost) {
    std::ostringstream sFileName;
    sFileName << strPre 
	      << std::setw(m_iDigits) << std::setfill('0') << iValue
	      << strPost << std::ends;
    return sFileName.str();
  }

  string GetName() {
    return GetName(m_strPre, m_iValue, m_strPost);
  }
  string GetName(const string & strPre) {
    return GetName(strPre, m_iValue, m_strPost);
  }
  string GetName(const string & strPre, int iValue) {
    return GetName(strPre, iValue, m_strPost);
  }
  string GetName(const string & strPre, const string & strPost) {
    return GetName(strPre, m_iValue, strPost);
  }

  // ---- attributes ----
  void SetValue(int iValue) { m_iValue = iValue; }
  void SetPre(const string & strPre) { m_strPre = strPre; }
  void SetPost(const string & strPost) { m_strPost = strPost; }
  void SetDigits(int iDigits) { m_iDigits = iDigits; }

  int GetValue() { return m_iValue; }
  const string & GetPre() { return m_strPre; }
  const string & GetPost() { return m_strPost; }
  int GetDigits() { return m_iDigits; }

  // ---- operations ----
  void Next() { m_iValue++; }

private:
  int m_iValue;
  int m_iDigits;
  string m_strPost;
  string m_strPre;
};

} // namespace PTG


#endif
