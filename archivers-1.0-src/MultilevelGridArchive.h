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
#ifndef _MGA_ARCHIVE_H_
#define _MGA_ARCHIVE_H_

#include "Archive.h"
#include <iostream>

template<class T>
class MultilevelGridArchive : public BaseArchive<T>
{
public:

  typedef typename BaseArchive<T>::size_type size_type;
  typedef typename BaseArchive<T>::iterator iterator;
  typedef typename BaseArchive<T>::element_type element_type;

  MultilevelGridArchive (unsigned int maxsize, unsigned int dim)
    : BaseArchive<T>(maxsize, dim)
  { }

  dominance_t add(const T &s)
  {
    // line 1-3
    dominance_t result = this->update(s);
    if (result == IS_DOMINATED_BY || result == EQUALS)
      return result;

    // line 4
    this->push_back(new element_type(s));

    // line 5-7
    if (this->size() <= this->max_size())
      return result;

    // line 8
    int b = compute_b_bar();
    int index = find_box_dominated(b);
    
    // line 9-10
    if (index == -1) {
      this->pop_back(); // return original A, without new solution
      return result;
    }

    // line 11-13
    unsigned int index_to_delete = 0;
    while (index >= 0) {
      b--;
      index_to_delete = (unsigned int) index; // remember box-dominated
      DEBUG2 (
              cout << "--------------------" << endl;
              cout << "Comparing for b=" << b << ":" << endl << endl;
              );
      index = find_box_dominated(b);
      DEBUG2 (cout << "... dominated index: " << index << endl << endl);
    }

    // line 14-19
    
    this->erase(this->begin() + index_to_delete);

    assert (this->size() <= this->max_size());
    return result;
  }
  
private:

  int compute_b_bar()
  {
    double fabs_max = 0.0;
    iterator iter;
    for (iter = this->begin(); iter < this->end(); ++iter) {
      vector<double>::iterator f;
      for (f = (*iter)->o.begin(); f < (*iter)->o.end(); ++f) {
	double my_fabs = fabs(*f);
	if (my_fabs > fabs_max)
	  fabs_max = my_fabs;
      }
    }
    return (int)floor(log2(fabs_max)) + 1;
  }

  void box_index_vector(vector<double>& y, vector<double>& r, int b)
  { 
    for (unsigned int i = 0; i < y.size(); i++) 
      r[i] = floor(y[i] / pow(2.0, b));
  }
  
  int find_box_dominated(int b)
  {
    Solution ra;
    Solution rb;
    vector<double> f(this->dimension(),0);
    ra.setObjectives(f);
    rb.setObjectives(f);

    int i = (int)this->size() - 1;

    for (; i >= 0; i--) {
      box_index_vector(this->at(i)->o, ra.o, b);
      for (int j = ((int)this->size()) - 1; j >= 0; j--) {
	if (i == j)
	  continue;
        DEBUG2 (
                cout << "Element " << i << ":    ";
                this->at(i)->print();
                cout << "boxvec  " << i << ":    ";
                ra.print();
                );
	box_index_vector(this->at(j)->o, rb.o, b);
        DEBUG2 (
                cout << "Element " << j << ":    ";
                this->at(j)->print();
                cout << "boxvec  " << j << ":    ";
                rb.print();
                cout << endl;
                );
	dominance_t dom = ra.dominance(rb);
	if (dom == IS_DOMINATED_BY || dom == EQUALS)
	  return i;
      }
    }
    assert (i < (int)this->size());
    return i;
  }
};

#endif
