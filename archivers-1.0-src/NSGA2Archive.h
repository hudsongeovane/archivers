/***********************************************-*- mode: c++ -*-********

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

 Parts of this code used under the terms of PISA_LICENSE.txt and from:

  ========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
  ========================================================================
  NSGA2

  Implementation in C for the selector side.
  
  Implements Petri net.
  
  file: nsga2.c
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
  ========================================================================

*************************************************************************/
#ifndef _ARCHIVE_NSGA2_H_
#define _ARCHIVE_NSGA2_H_

#include "Archive.h"
#include "Random.h"

template<class T>
class NSGA2Archive : public BaseArchive<T> {

public:
  typedef typename BaseArchive<T>::element_type element_type;
  typedef typename BaseArchive<T>::size_type size_type;
  typedef typename BaseArchive<T>::iterator iterator;
  typedef typename BaseArchive<T>::const_iterator const_iterator;
  using BaseArchive<T>::ubound;
  using BaseArchive<T>::lbound;
  using BaseArchive<T>::dimension;

  class NSGA2ArchiveElementData : public BaseArchiveElementData
  {
  public:
    ~NSGA2ArchiveElementData() { }
    NSGA2ArchiveElementData() : index(-1), fitness(-1) { }

    unsigned int index;
    double fitness;
  };

  NSGA2Archive (unsigned int maxsize, unsigned int dim, Random &rng)
    : BaseArchive<T>(maxsize + 1, dim),
      copies (maxsize + 1, 1),
      front (maxsize + 1, vector<int> (maxsize + 1)),
      dist (maxsize + 1),
  //  double  *f_max;
  //  double  *f_min;
  //  double  *f_norm;
      rng(rng)
  { // We reserve maxsize + 1 but then set max_size to the real value
    this->max_size(maxsize);
  }

  static NSGA2ArchiveElementData & getData(const element_type &s) {
    return *(static_cast<NSGA2ArchiveElementData*>(s.data));
  }

  dominance_t add (const T & a)
  {
    element_type s (a);
    dominance_t result = NONDOMINATED;

    this->push_back(s);
    if (!this->overfull())
      return result;

    /* Calculates NSGA2 fitness values for all individuals */
    nondominatedSort();

    /* Calculates distance cuboids */
    calcDistances();

    /* Performs environmental selection
       (truncates 'pp_all' to size 'alpha') */
    environmentalSelection();

    assert (!this->overfull());
    return result;
  }

private:
  // NSGA2 internal global variables
  vector<int>  copies;
  vector< vector<int> > front;
  vector<double> dist;
  //  double  *f_max;
  //  double  *f_min;
  //  double  *f_norm;
  Random &rng;

  void 
  push_back (const element_type &s)
  {
    element_type * a = new element_type(s);
    a->data = new NSGA2ArchiveElementData();
    this->BaseArchive<T>::push_back(a);
  }

  void
  nondominatedSort()
  {
    int i, j, l;
    int size;

    size = this->size();

    vector<bool> d (size, true);
    vector<bool> f (size, true);
    for (i = 0; i < size; i++)      {
      element_type & s = *((*this)[i]);
      getData(s).fitness = 0;
      copies[i] = 0;
    }

    int num = size;
    for (l = 0; l < size; l++)    {
      /* find next front */
      for (i = 0; i < size; i++)	{
        element_type * si = (*this)[i];
        d[i] = false;
        if (f[i] != false) {
          for (j = 0; j < i && d[i] == false; j++)
            if (f[j] != false) {
              element_type * sj = (*this)[j];
              if (sj->dominates(*si))
                d[i] = true;
            }
          for(j = i + 1; j < size && d[i] == false; j++)
            if (f[j] != false) {
              element_type * sj = (*this)[j];
              if (sj->dominates(*si))
                d[i] = true;
            }
        }
      }
      /* extract front */
      for (i = 0; i < size; i++) {
        element_type * si = (*this)[i];
        if (f[i] != false && d[i] == false) {
          getData(*si).fitness = l;
          f[i] = false;
          num--;
          front[l][copies[l]] = i;
          copies[l] += 1;
        }
      }
      if (num == 0)
        break;
    }
    
    DEBUG2_FUNPRINT ("fronts:\n");
    DEBUG2 (
            for (i = 0; i < size; i++) {
              for (int k = 0; k < copies[i]; k++) {
                fprintf (stderr, "%d[%d] : ", i, k);
                vector_fprintf(stderr, "%f", (*this)[front[i][k]]->o);
                fprintf (stderr, "\n");
              }
            });
    DEBUG1 (
            int sum = 0;
            int last_empty = -1;
            for (l = 0; l < size; l++) {
              sum += copies[l];
              if (copies[l] == 0) {
                if (last_empty < 0)
                  last_empty = l;
              } else if (last_empty > 0) {
                fprintf (stderr, "error: assert failed: found non-empty front %d after last empty one %d!\n",
                         l, last_empty);
                assert (0);
              }
            }
            assert (sum == size);
            for (l = 0; l < size; l++) {
              for (i = 0; i < copies[l]; i++) {
                for (j = i + 1; j < copies[l]; j++) {
                  dominance_t result = (*this)[front[l][i]]->dominance(*((*this)[front[l][j]]));
                  // NSGA2 does keep weakly dominated individuals.
                  if (result != NONDOMINATED && result != EQUALS) {
                    fprintf (stderr, "error: assert failed: %d[%d] : ", l, i);
                    vector_fprintf(stderr, "%f", (*this)[front[l][i]]->o);
                    fprintf (stderr, " !nondominated with ");
                    fprintf (stderr, "%d[%d] : ", l, j);
                    vector_fprintf(stderr, "%f", (*this)[front[l][j]]->o);
                    fprintf (stderr, "\n");
                    assert(0);
                  }
                }
              }
            });
  }


  double calcDistance(const element_type &a, const element_type &b)
  {
    int i;
    int dim = dimension();
    double distance = 0;

    if (a.is_equal (b))
      return 0;

    for (i = 0; i < dim; i++) {
      double tmp_double = a.o[i] - b.o[i];
      distance += tmp_double * tmp_double;
    }

    return sqrt(distance);
  }

  void
  calcDistances()
  {
    int i, j, l, d;
    int size = this->size();
    int dim = dimension();

    // ??? Why divide?
    double dmax = HUGE_VAL / (dim + 1);

    // initialize copies[] vector and NN[][] matrix
    for (i = 0; i < size; i++) {
      dist[i] = 1;
    }

    for (l = 0; l < size; l++)    {
      for (d = 0; d < dim; d++)	{
        /* sort accorting to d-th objective */
        for (i = 0; i < copies[l]; i++)	    {
          int min_index = -1;
          int min = i;
          for (j = i + 1; j < copies[l]; j++) {
            if ( (*this)[front[l][j]]->o[d] < (*this)[front[l][min]]->o[d])
              min = j;
          }
          min_index = front[l][min];
          front[l][min] = front[l][i];
          front[l][i] = min_index;
        }

        /* add distances */
        for (i = 0; i < copies[l]; i++) {
          if (i == 0 || i == copies[l] - 1)
            dist[front[l][i]] += dmax;
          else {
            dist[front[l][i]] +=
              (*this)[front[l][i+1]]->o[d] - (*this)[front[l][i-1]]->o[d];
          }
        }
      }
    }
  }

  void environmentalSelection()
  {
    size_type i, j;
    size_type size = this->size();
    
    for (i = 0; i < size; i++)
      {
	getData(*((*this)[i])).fitness += 1.0 / dist[i];
      }

    for (i = 0; i < this->max_size(); i++)
      {
	int min = i;
	for (j = i + 1; j < size; j++)
	{
          if (getData(*((*this)[j])).fitness < getData(*((*this)[min])).fitness)
            min = j;
	}
	element_type * p_min = (*this)[min];
	(*this)[min] = (*this)[i];
	(*this)[i] = p_min;
    }
    
    for (i = this->max_size(); i < size; i++)
      {
        this->erase(this->end() - 1);
      }

    assert (this->size() <= this->max_size());
  }

};

#endif

