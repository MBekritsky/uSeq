/*
Timer.cc

ABSTRACT:
A timer class to track time elapsed for a program or
function.

CREATION DATE:
08.08.2010

LAST REVISION:
26.08.2010

AUTHOR:
Michael Schatz
Mitchell Bekritsky (modified elapsed function)
*/

#include "Timer.h"
#include <cstring>
#include <cstdio>

using namespace std;

Timer_t::Timer_t()
{
  start ();
  memset (&tEnd_, 0, sizeof(struct timeval) );
}

void Timer_t::start()
{
  gettimeofday (&tStart_, NULL);
}

void Timer_t::end()
{
  gettimeofday (&tEnd_, NULL);
}

double Timer_t::length()
{
  end();
  return ( ((tEnd_ . tv_sec - tStart_ . tv_sec) * 1000000) + (tEnd_ . tv_usec -	tStart_ . tv_usec) );
}

string Timer_t::elapsed()
{
  double r = length();
  char elapsedString[1024];

  int hour   = (int) (r / USECTOHOUR);
  int minute = (int) ((r - (hour * USECTOHOUR)) / USECTOMIN);
  int second = (int) ((r - (hour * USECTOHOUR) - (minute * USECTOMIN) ) / USECTOSEC);
  int millisecond = (int) ((r - (hour * USECTOHOUR) - (minute * USECTOMIN) - (second * USECTOSEC) ) / USECTOMSEC);
  int microsecond = (int) (r - (hour * USECTOHOUR) - (minute * USECTOMIN) - (second * USECTOSEC) - (millisecond * USECTOMSEC));

  sprintf(elapsedString, "%02d:%02d:%02d.%03d.%03d", hour, minute, second, millisecond, microsecond);
  return elapsedString;
}
