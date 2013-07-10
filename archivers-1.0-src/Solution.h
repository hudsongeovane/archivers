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
#ifndef SOLUTION_H
#define SOLUTION_H

#include <cstdio>
#include <vector>

#define point_printf_format "%-16.15g"

void
vector_fprintf(FILE *stream, const char * format,
               const std::vector<double> & vec, const char *sep = " "); //const

using namespace std;

enum dominance_t { IS_DOMINATED_BY, DOMINATES, NONDOMINATED, EQUALS };


class Solution {
  
public:

  vector<double> o;

  static void Initialise (unsigned d) {
    Solution::_num_objs = d;
  };

  Solution (void) 
    : o ()
  { 
    o.reserve(_num_objs); 
  };
  
  void setObjectives(vector<double> src)
  {
    this->o = src;
  };

  unsigned int num_objs(void) const { return Solution::_num_objs; };

  Solution * clone(void) { return new Solution(*this); };

  void print(FILE *stream=stdout) const
  {
    vector<double>::const_iterator iter = o.begin();
    fprintf (stream, point_printf_format, *iter);
    for (++iter; iter != o.end(); ++iter) {         
      fprintf (stream, "\t" point_printf_format, *iter);
    }
    fprintf (stream, "\n");
  };

  dominance_t 
  dominance (const Solution &b) const
  {
    bool a_leq_b, b_leq_a;
    a_leq_b = b_leq_a = true;
    unsigned n_objs = this->o.size();
    for (unsigned int d = 0; d < n_objs; d++) {
      a_leq_b = a_leq_b && (this->o[d] <= b.o[d]);
      b_leq_a = b_leq_a && (b.o[d] <= this->o[d]);
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

  bool is_equal (const Solution &b) const
  {
    bool equal = true;
    unsigned i = 0;
    unsigned n_objs = this->o.size();
    while (equal && i < n_objs) {
      // FIXME: This should use some epsilon difference.
      equal = (this->o[i] == b.o[i]);
      ++i;
    }
    return equal;
  }

  bool dominates (const Solution &b) const
  {
    dominance_t result = this->dominance(b);
    return result == DOMINATES;
  }

  bool weaklydominates (const Solution &b) const
  {
    dominance_t result = this->dominance(b);
    return (result == DOMINATES || result == EQUALS);
  }

  
private:
  static unsigned int _num_objs;
};

#endif
