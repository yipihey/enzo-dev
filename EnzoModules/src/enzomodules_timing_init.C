/***********************************************************************
/
/  ENZOMODULES: initialize Enzo's global performance timer.
/
/  Some libenzo routines (e.g. grid::SolveHydroEquations) use the TIMER_START
/  macros, which dereference the global `enzo_timer`.  That object is normally
/  created in enzo.C's main(); a bridge that calls those routines directly must
/  create it first.  This is isolated in its own translation unit -- including
/  only EnzoTiming.h and the standard headers it pulls in, with NO
/  macros_and_parameters.h -- to avoid the min/max/math macro collisions that
/  occur when EnzoTiming.h is included after the Enzo macros.
/
************************************************************************/

#define EXTERN extern
#include "EnzoTiming.h"

extern "C" void enzomodules_init_timer()
{
  if (enzo_timer == NULL)
    enzo_timer = new enzo_timing::enzo_timer();
}
