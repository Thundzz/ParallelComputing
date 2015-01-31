#include "cdt_bords.h"
#include "tools.h"

#include <math.h>

double f1( double posx, double posy, double t)
{
  UNUSED(t);
  double function;
  function = sin(posx) + cos(posy);
  return(function);
}

double f2 ( double posx, double posy, double t)
{
  UNUSED(t);
  return 2*(posy - posy*posy + posx - posx*posx);
}

double func_zero(double posx, double posy, double t)
{
  UNUSED(t);
  UNUSED(posx);
  UNUSED(posy); 
  return 0;
}

double func_one(double posx, double posy, double t)
{
  UNUSED(t);
  UNUSED(posx);
  UNUSED(posy);
  return 1;
}
