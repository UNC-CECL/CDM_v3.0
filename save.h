#ifndef __SAVE_H__
#define __SAVE_H__

#include <vector>

#include "PTG_FileNames.h"

#include "globals.h"
#include "vec.h"
#include "func.h"


class CSaveFields 
{
  bool m_bOff;
  bool m_bForce;
  int m_iEveryT;
  vector<double> m_veciSpecialT;
  int m_iT;
  bool m_bXLine;
  string m_strDir;
  PTG::FileNameCounter m_FileName;

public:
  CSaveFields(const dunepar& P);

  void SetTime(int t) { m_iT = t; }
  void IncTime() { m_iT++; }

  void SetForce(bool b) { m_bForce = b; }

  bool Save(const TFktScal& f, const string& strName);
  bool Save(const TFktVec& f, const string& strName);

private:
  template<class T>
  bool DoSave(const T& f, const string& strName);
};


#endif
