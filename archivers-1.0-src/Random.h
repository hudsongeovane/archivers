/***********************************************-*- mode: c++ -*-********

 Random Number Generation Class (C++)

 ---------------------------------------------------------------------

                       Copyright (c) 2010
          Manuel Lopez-Ibanez  <manuel.lopez-ibanez@ulb.ac.be>

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

  IMPORTANT NOTE: Please be aware that the fact that this program is
  released as Free Software does not excuse you from scientific
  propriety, which obligates you to give appropriate credit! If you
  write a scientific paper describing research that made substantive
  use of this program, it is your obligation as a scientist to
  acknowledge its use.  Moreover, as a personal note, I would
  appreciate it if you would email manuel.lopez-ibanez@ulb.ac.be with
  citations of papers referencing this work so I can mention them to
  my funding agent and tenure committee.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at http://www.gnu.org/copyleft/gpl.html

*************************************************************************/
#ifndef RANDOM_H
#define RANDOM_H

#include "debug.h"
#include "gsl_rng.h"

#include <vector>
#include <cassert>
#include <cstdlib>

class Random {

private:
  gsl_rng * rng;

public:
  Random(unsigned long int arg) : seed(arg) 
  {
    rng = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(rng, seed);
  }

  unsigned long int seed;

  /* next() and rand_01() generate a pseudo-random real number in the
     range [0,1). */
  double rand_01() 
  { 
    assert(rng != NULL);
    return(gsl_rng_uniform(rng));
  }

  double next() { return rand_01(); }

  double rand_double(double low, double high)
  {
    assert(rng != NULL);
    return(low + (gsl_rng_uniform(rng)*high));
  }

  // Generate a pseudo-random integer in the range [0, HIGH - 1]
  int rand_int (int high)
  {
    assert(high > 0);
    return int(rand_01 () * high);
  }

  // Generate a pseudo-random integer in the range [LOW, HIGH]
  int rand_int (int low, int high)
  {
    int tmp = high - low;
    assert(rng != NULL);
    assert(tmp >= 0);

    if (tmp <= 0) {
      return high;
    } else if (tmp == 1) {
      if(this->next() >= 0.5) return high;
      else return low;
    } else {
#if DEBUG > 0
      assert( ((unsigned long)tmp + 1UL) < (gsl_rng_max (rng) - gsl_rng_min(rng)));
      assert( ((unsigned long)tmp + 1UL) > 0);
#endif
      return(low + (int) gsl_rng_uniform_int(rng, (unsigned long)tmp + 1UL));
    }
  }

  ptrdiff_t operator() (ptrdiff_t high)
  {
    assert (high > 0);
    return ptrdiff_t(rand_01 () * high);
  }
};
#endif
// Local Variables:
// mode: c++
// End:
