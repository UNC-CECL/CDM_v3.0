/******************************************************************************
    $Id: main.cc,v 1.21 2005/04/13 15:06:11 parteli Exp $
******************************************************************************/

#include <time.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <unistd.h>

#include "globals.h"
#include "dune_evolution.h"

//! Application class.

class CApp {
public:
  CApp(int argc, char **argv);
  virtual ~CApp();

  //! Execute the simulation
  virtual int Run();

private:
  /*! Class handling parameter file.  */
  dunepar m_para;

  /*!  New top-level simulation object  */
  evolution *m_evol;
};


/*!  Constructor scans command line and parameter file for global parameters.
The name of the parameter file can be given on the command line by
"parafile=..." and defaults to "default.par".  The unique instance of the
duneglobals class is created via duneglobals::initialise().  Depending on the
parameter file version, a CPde object (for backward compatibility; version <
2.0) or a dune_evolution object is created.  */

CApp::CApp(int argc, char **argv) : m_evol(0)

{
  struct utsname u;
  time_t t= time(NULL);
  
  if( !uname(&u) )
    cout << "Running on " << u.nodename << ", machine type " << u.machine << ", with process id " << getpid() << "\n";
  cout << ctime( &t ) << "\n";

  // find name for parameter file
  // Format for command line arguments: <token>=<value>, for instance 
  // parafile=myparms.par
  m_para.scan(argc, argv);

  // Initialise global parameter object:
  duneglobals::initialise( m_para );
  
  m_evol= new dune_evol_3d(m_para);
}


/*!  Destructor deletes all objects created in the constructor.  */

CApp::~CApp()
{
  delete m_evol;
}


int CApp::Run()
{
  int iNt = m_para.getrequired<int>("Nt");
  int iNt0 = m_para.getdefault<int>("Nt0", 0);
  int saveinterval = m_para.getrequired<int>( "save.every" );

  int step, savecounter;
    
  savecounter= saveinterval - 1;
  for (step= iNt; step > iNt0; --step, --savecounter) {
      m_evol->step();
      cout << "Main program: step " << evolution::steps() << " done, time is " << evolution::time() << " seconds\n";
      if( savecounter<=0 ) {
        dunedata::save_all_data();
        savecounter= saveinterval;
      }
  }
  if( savecounter!=saveinterval-1 || iNt==0 )
      dunedata::save_all_data();
  
  return 0;
}


////////////////
// main
//

int main(int argc, char **argv)
{
  CApp App(argc, argv);

  App.Run();

  return 0;
}

