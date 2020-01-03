
#include <iostream>

#include "evolution.h"

//*****************************************************************************
//  class evolution

/*!  Pointer to what ist meant to be the only instance of this class.  See the
    constructor.  */
evolution *evolution::instance= NULL;


/*! The constructor checks that there is no other evolution object already.  If
  there is, the program is terminated with an error message.  Otherwise, this
  instance is written to the static variable instance, and the member variables
  are initialised.  */
evolution::evolution(const dunepar &par) : dunedata(par)
{
  if( instance ) {
    cerr << "evolution constructor: There can be only one instance of an evolution class!  Aborting..." << endl;
    exit(1);
  }
  instance= this;
  m_started= false;
  m_steps= duneglobals::startstep();
  m_time= duneglobals::starttime();
}

/*!  The destructor resets the static variable instance.  */
evolution::~evolution()
{
  instance = NULL;
}

