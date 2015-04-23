//timer.cpp

#include <sys/time.h>
#include "timer.h"

  void timer::start()
  {
    gettimeofday(&tvi,0); 
  }
  
  void timer::stop()
  {
    gettimeofday(&tvf,0);
  }
  
  double timer::how_long()
  {
    double initial = tvi.tv_sec + (tvi.tv_usec/1000000.);
    double final = tvf.tv_sec + (tvf.tv_usec/1000000.);
    
    return (final - initial);
  }
