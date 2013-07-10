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
#ifndef _ePARETO_ARCHIVE_H_
#define _ePARETO_ARCHIVE_H_

#include "Archive.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

template<class T>
class eParetoArchive : public BaseArchive<T> {

public:

  typedef typename BaseArchive<T>::element_type element_type;
  typedef typename BaseArchive<T>::iterator iterator;
  typedef typename BaseArchive<T>::const_iterator const_iterator;
  using BaseArchive<T>::ubound;
  using BaseArchive<T>::lbound;

  class eParetoArchiveElementData : public BaseArchiveElementData
  {
  private: 
    eParetoArchiveElementData () {}

  public:
    ~eParetoArchiveElementData() { }

    eParetoArchiveElementData (const T &s) 
      : _box(s.num_objs(), 0)
    {
    }

    void calculate_box (const T &s, double epsilon)
    {
      assert(epsilon > 0);
      int dim = _box.size();
      for (int k = 0; k < dim; k++) {
        // _box[k] = int (floor (log10(double(s.o[k])) / log10 (1.0 + epsilon))); // multiplicative
        _box[k] = (int) floor(s.o[k] / epsilon); // additive
      }
    }

    unsigned int _seq_num;
    vector<int> _box;

  };

  eParetoArchive (unsigned int maxsize, unsigned int dim, double eps = 0)
    : BaseArchive<T>(maxsize, dim),
      seq_max(0), epsilon(eps)
  { 
    assert(eps >= 0);
  }

  static eParetoArchiveElementData & getData (const element_type &s)
  {
    return *(static_cast<eParetoArchiveElementData*>(s.data));
  }

  static unsigned int& seq_num (const element_type &s) {
    return getData(s)._seq_num;
  }

  static vector<int>& box (const element_type &s) {
    return getData(s)._box;
  }

  static void calculate_box (const element_type &s, double epsilon)
  {
    getData(s).calculate_box (s, epsilon);
  }

  dominance_t box_dominance (const element_type &a,
                             const element_type &b) const
  {
    int dim = box(a).size();
    bool a_leq_b, b_leq_a;
    a_leq_b = b_leq_a = true;
    for (int d = 0; d < dim; d++) {
      a_leq_b = a_leq_b && (box(a)[d] <= box(b)[d]);
      b_leq_a = b_leq_a && (box(b)[d] <= box(a)[d]);
    }
    // FIXME: This could use bit-wise operations to be much faster.
    if (!a_leq_b && !b_leq_a) {
      return NONDOMINATED;
    } else if (!a_leq_b) {
      return IS_DOMINATED_BY;
    } else if (!b_leq_a) {
      return DOMINATES;
    } else {
      return EQUALS;
    }
  }
  
  dominance_t add (const T & a)
  {
    element_type s (a);
    s.data = new eParetoArchiveElementData(s);

    DEBUG2_FUNPRINT ("Adding point %d: ", this->sequence);
    DEBUG2 (vector_fprintf(stderr, "%.6f", a.o);
            fprintf(stderr, ", archive size = %d\n", this->size()));
  
    dominance_t result;

    // optional initial phase until first epsilon is set
    if (this->epsilon == 0.0) {
      result = this->update(s);
      if (result == DOMINATES || result == NONDOMINATED) {
        this->push_back(s);
      }
      // calculate initial epsilon when first two nondom solutions found
      if (this->size() > 1) {
	this->update_epsilon();
	// compute box values
	for (iterator j = this->begin(); j != this->end(); ++j) {
	  calculate_box (**j, this->epsilon);
	}
	this->truncate();
      }
    }
    else {
      // Here we assume there is a current epsilon > 0,
      // and all archive members have up-to-date box values.
      assert(epsilon > 0);
      assert (this->size() <= this->max_size());

      calculate_box (s, epsilon);

      result = this->box_update (s);
      if (!(result == DOMINATES || result == NONDOMINATED))
	return result;

      this->push_back(s);
    }
  
    // it might take several iterations of joining boxes
    while (this->size() > this->max_size()) {
      assert(epsilon > 0);
      // enlarge epsilon as joining adjacent boxes ('doubling boxes')
      // epsilon = (1. + epsilon) * (1. + epsilon) - 1.; // multiplicative
      epsilon = 2 * epsilon; // additive
      // recompute box values
      for (iterator j = this->begin(); j != this->end(); ++j) {
        calculate_box (**j, this->epsilon);
      }
      this->truncate();
    }
    assert (this->size() <= this->max_size());
    return result;
  }

private:
  unsigned int seq_max;
  double epsilon;

  dominance_t
  box_update(const element_type &s)
  {
    bool dominates_p = false;
    iterator iter = this->begin();
    while (iter != this->end()) {
      const element_type &a = **iter;
      switch (box_dominance(s, a))
        {
        case IS_DOMINATED_BY:
          return IS_DOMINATED_BY;
        
        case EQUALS:
          if (s.dominance(a) == DOMINATES) {
            iter = this->erase(iter);
            dominates_p = true;
          } else
            return EQUALS;
          break;

        case DOMINATES:
          iter = this->erase(iter);
          dominates_p = true;
          break;
          
        case NONDOMINATED:
          ++iter;
          break;
        }
    }
    return (dominates_p) ? DOMINATES : NONDOMINATED;
  }

  void 
  push_back (const element_type &s)
  {
    element_type * a = new element_type(s);
    assert (s.data);
    a->data = new eParetoArchiveElementData(getData(s));
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

    // This epsilon calculation is different from eApproxArchive. See
    // Theorems 5 and 6 of Laumanns's PhD thesis.  The eParetoArchive
    // needs less points than eApproxArchive for the same epsilon.  So
    // for a given maximum number of points, a smaller epsilon is
    // necessary.
    // use floor to have a safe upper bound
    // this->epsilon = pow (double(K), 1. / floor(pow (double(this->max_size()), 1. / (dim - 1))) ) - 1; // multiplicative
    this->epsilon = K / floor(pow (double(this->max_size()), 1. / (dim-1))); // additive  

    // enlarging epsilon a bit to avoid rounding errors
    this->epsilon = this->epsilon + this->epsilon / 65536.0;

    DEBUG2 (this->debug());
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
	  dominance_t dom = box_dominance(*aj, *ak);
          if (dom == EQUALS || dom == DOMINATES) {
            this->erase(k);
            continue;
          } 
          if (dom == IS_DOMINATED_BY) {
            this->erase(j);
	    j--;
            break;
          } 
        } else { 
	  assert(1 < 0); // we should never be in this branch!
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

