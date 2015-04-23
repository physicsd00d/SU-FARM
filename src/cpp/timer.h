//timer.h

#include <sys/time.h>

class timer
{
  private:
  struct timeval tvi;
  struct timeval tvf;


  public:
  void start();
  void stop();
  double how_long();



};
