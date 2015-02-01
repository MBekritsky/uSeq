/*
Timer.h

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

#ifndef __TIMER_H__
#define __TIMER_H__

#include <string>
#include <ctime>
#include <iostream>
#include <sys/time.h>

using namespace std;

class Timer_t
{
 private:
  struct timeval tStart_;
  struct timeval tEnd_;
  static const long int USECTOHOUR = 60L * 60L * 1000000L;
  static const long int USECTOMIN  = 60 * 1000000;
  static const long int USECTOSEC  = 1000000;
  static const long int USECTOMSEC = 1000;

 public:
  Timer_t();
  void start();
  void end();
  double length();
  string elapsed();
};

#endif /* Timer.h */
