#include "f2c.h"
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

int cpu_time__(doublereal *r)
{
   struct rusage usage;
   int err;
   err = getrusage (RUSAGE_SELF, &usage);

   *r = usage.ru_utime.tv_sec;
   return err;

}
