/***********************************************-*- mode: c++ -*-*********

 Sequential Online Archiving of Objective Vectors

 ---------------------------------------------------------------------

                          Copyright (c) 2011
         Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
             Joshua Knowles <j.knowles@manchester.ac.uk>
                 Marco Laumanns <mlm@zurich.ibm.com>

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at http://www.gnu.org/copyleft/gpl.html

 IMPORTANT NOTE: Please be aware that the fact that this program is
 released as Free Software does not excuse you from scientific
 propriety, which obligates you to give appropriate credit! If you
 write a scientific paper describing research that made substantive use
 of this program, it is your obligation as a scientist to (a) mention
 the fashion in which this software was used in the Methods section;
 (b) mention the algorithm in the References section. The appropriate
 citation is:
 
     M. Lopez-Ibanez, J. Knowles, and M. Laumanns. On Sequential
     Online Archiving of Objective Vectors. In R. Takahashi et al., 
     editors, Evolutionary Multi-criterion Optimization (EMO 2011),
     volume 6576 of Lecture Notes in Computer Science, pages 46-60.
     Springer, Heidelberg, Germany, 2011.
 
 Moreover, as a personal note, I would appreciate it if you would email
 manuel.lopez-ibanez@ulb.ac.be with citations of papers referencing this
 work so I can mention them to my funding agent and tenure committee.

 ----------------------------------------------------------------------

 Relevant literature:

 [1] M. Lopez-Ibanez, J. Knowles, and M. Laumanns. On Sequential
     Online Archiving of Objective Vectors. In EMO 2011,
     LNCS. Springer, 2011.

*************************************************************************/
#ifndef _eAPPROX_ARCHIVE_H_
#define _eAPPROX_ARCHIVE_H_

#include "Archive.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

template<class T>
class eApproxArchive : public BaseArchive<T> {

public:
 
  typedef typename BaseArchive<T>::element_type element_type;
  typedef typename BaseArchive<T>::iterator iterator;
  typedef typename BaseArchive<T>::const_iterator const_iterator;
  using BaseArchive<T>::ubound;
  using BaseArchive<T>::lbound;
  using BaseArchive<T>::dimension;

  class eApproxArchiveElementData : public BaseArchiveElementData
  {
  public:
    ~eApproxArchiveElementData() { }

    unsigned int _seq_num;
  };

  eApproxArchive (unsigned int maxsize, unsigned int dim, double eps = 0)
    : BaseArchive<T>(maxsize, dim),
      seq_max(0), epsilon(eps)
  { 
    assert(eps >= 0);
  }

  static eApproxArchiveElementData & getData(const element_type &s) {
    return *(static_cast<eApproxArchiveElementData*>(s.data));
  }

  static unsigned int& seq_num (const element_type &s) {
    return getData(s)._seq_num;
  }
  
  dominance_t add (const T & a)
  {
    element_type s (a);

    DEBUG2_FUNPRINT ("Adding point %d: ", this->sequence);
    DEBUG2 (vector_fprintf(stderr, "%.6f", a.o);
            fprintf(stderr, ", archive size = %d\n", this->size()));

    // lines 2-3
    if (this->eDominates(s)) {
      return IS_DOMINATED_BY;
    }
    
    // line 4-6
    dominance_t result = this->update(s);
    assert (result == DOMINATES || result == NONDOMINATED);
    this->push_back(s);

    // initialize epsilon
    if (this->epsilon == 0.0 && this->size() > 1)
      this->update_epsilon();
  
    // enlarge epsilon if archive still too large
    if (this->size() > this->max_size()) {
      this->update_epsilon();
      this->truncate();
      assert (this->size() <= this->max_size());
    }
    return result;
  }

private:
  unsigned int seq_max;
  double epsilon;

  bool
  eDominates(const element_type &a,
             const element_type &b) const
  {
    int dim = a.num_objs();
    double epsilon = this->epsilon;
    bool a_leq_b = true;

    for (int k = 0; a_leq_b && k < dim; k++) {
      //      a_leq_b = a_leq_b && (a.o[k] <= (1. + epsilon) * b.o[k]); // multiplicative version
      // additive approximation:
      a_leq_b = a_leq_b && (a.o[k] <= epsilon + b.o[k]); 
    }
    return a_leq_b;
  }

  bool
  eDominates(const element_type &s) const
  {
    const_iterator iter = this->begin();
    while (iter != this->end()) {
      const element_type &a = **iter;
      if (this->eDominates(a, s))
        return true;
      ++iter;
    }
    return false;
  }

  void 
  push_back (const element_type &s)
  {
    element_type * a = new element_type(s);
    a->data = new eApproxArchiveElementData();
    seq_num(*a) = this->seq_max++;
    this->BaseArchive<T>::push_back(a);
  }

  void
  update_epsilon(void)
  {
    int k;
    int dim = this->dimension();
    double K = 0;

    // We need the Pareto front bounds, but since these are unknown,
    // we use the current archive's bounds.
    this->calculate_bounds();

    for (k = 0; k < dim; k++) {
      // K = max (K, this->ubound[k] / this->lbound[k]); // multiplicative
      K = max (K, this->ubound[k] - this->lbound[k]); // additive
    }
    assert(K >= 0);

    // use floor to have a safe upper bound
    // this->epsilon = pow (double(K), 1. / floor(pow (double(this->max_size()), 1. / dim)) ) - 1.;     // multiplicative
    this->epsilon = K / floor(pow (double(this->max_size()), 1. / dim)); // additive  

    // enlarging epsilon a bit to avoid rounding errors
    this->epsilon = this->epsilon + this->epsilon / 65536.0;

    DEBUG2(this->debug ());
    assert (this->epsilon > 0.0);
  }

  void
  truncate(void)
  {
    for (iterator j = this->begin(); j != this->end(); ) {
      for (iterator k  = j + 1; k != this->end(); ) {
        const element_type * aj = *j;
        const element_type * ak = *k;
	if (j < k) assert(seq_num(*aj) < seq_num(*ak));
        if (seq_num(*aj) < seq_num(*ak)) {
	  assert (j < k);
          if (eDominates(*aj, *ak)) {
            this->erase(k);
            continue;
          } 
        } else { 
          assert(1 < 0); // we should never be in this branch!
          if (eDominates(*ak, *aj)) {
            j = this->erase(j);
            k = j + 1;
            continue;
          }
        }
        ++k;
      }
      ++j;
    }
  }

  static void print_element (double i) { fprintf (stderr, " %f", i); }
  
  void debug(void) {
    fprintf (stderr, "size = %d, e = %f", this->size(), epsilon);
    fprintf (stderr, ", ubound =");
    for_each (ubound.begin(), ubound.end(), print_element);
    fprintf (stderr, ", lbound =");
    for_each (lbound.begin(), lbound.end(), print_element);
    fprintf (stderr, "\n");
  }

};

#endif

