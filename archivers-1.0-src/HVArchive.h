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

 Relevant literature:

 [1] M. Lopez-Ibanez, J. Knowles, and M. Laumanns. On Sequential
     Online Archiving of Objective Vectors. In EMO 2011,
     LNCS. Springer, 2011.

*************************************************************************/
#ifndef _HYPERVOLUME_ARCHIVE_H_
#define _HYPERVOLUME_ARCHIVE_H_

#include "Archive.h"
#include "Random.h"
#include <limits.h>
#include "hv.h"

#include "debug.h"

template<class T>
class HVArchive : public BaseArchive<T> {

public:
  typedef typename BaseArchive<T>::element_type size_type;
  typedef typename BaseArchive<T>::element_type element_type;
  typedef typename BaseArchive<T>::iterator iterator;
  typedef typename BaseArchive<T>::const_iterator const_iterator;
  using BaseArchive<T>::ubound;
  using BaseArchive<T>::lbound;
  using BaseArchive<T>::dimension;
  using BaseArchive<T>::uev;

  class HVArchiveElementData : public BaseArchiveElementData
  {
  public:
    ~HVArchiveElementData() { }
    HVArchiveElementData() : index(-1), fitness(-1) { }

    unsigned int index;
    double fitness;
  };

  HVArchive (unsigned int maxsize, unsigned int dim, Random &rng, bool only_nondominated_p = true)
    : BaseArchive<T>(maxsize + 1, dim),
      _only_nondominated_p (only_nondominated_p),
      copies (maxsize + 1, 1),
      front (maxsize + 1, vector<int> (maxsize + 1)),
      dist (maxsize + 1),
      dominated (maxsize + 1, 0),
  //  double  *f_max;
  //  double  *f_min;
  //  double  *f_norm;
      rng(rng)
  { 
    // We reserve maxsize + 1 but then set max_size to the real value
    this->max_size(maxsize);
    // The original SMS says that uevs are only kept for d == 2.
    _keep_uevs = (dim == 2);
    //_keep_uevs = true;
  }

  static HVArchiveElementData & getData(const element_type &s) {
    return *(static_cast<HVArchiveElementData*>(s.data));
  }

  dominance_t add (const T & a)
  {
    element_type s (a);
    dominance_t result = NONDOMINATED;
    int pos;

    this->sequence++;
    DEBUG2 (
            fprintf(stderr, "Adding point %d: ", this->sequence);
            vector_fprintf(stderr, "%.6f", a.o);
            fprintf(stderr, ", archive size = %d\n", this->size())
            );
    /*
    if (_only_nondominated_p) {
      result = this->update(s);
      if (result == IS_DOMINATED_BY || result == EQUALS)
        return result;
    }
    */
    this->push_back(s);
    if (!this->overfull())
      return result;

    this->calculate_bounds(); // Calculate ranges of grid NOT INCLUDING the new point.

    // Calculates SMS fitness values for all individuals. 
    // FIXME: It is not clear in the original SMS whether fronts may
    // contain equal objective vector.
    nondominatedSort();

    if (copies[1] > 0) {
      pos = find_most_dominated();
    } else {
      pos = find_least_hv_contributor();
    }

    this->erase(this->begin() + pos);

    assert (!this->overfull());
    return result;
  }

private:
  // SMS internal global variables
  bool _only_nondominated_p; // Keep only non-weakly-dominated objective vectors.
  bool _keep_uevs;
  vector<int>  copies;
  vector< vector<int> > front;
  vector<double> dist;
  vector<int> dominated;
  //  double  *f_max;
  //  double  *f_min;
  //  double  *f_norm;
  Random &rng;

  void 
  push_back (const element_type &s)
  {
    element_type * a = new element_type(s);
    a->data = new HVArchiveElementData();
    this->BaseArchive<T>::push_back(a);
  }

  // FIXME: How to avoid GCC inlining this, so it is visible in GDB?
  void
  debug_fronts ()
  {
    for (int i = 0; i < int(this->size()); i++) {
      for (int k = 0; k < copies[i]; k++) {
        fprintf (stderr, "%d[%d] : ", i, k);
        vector_fprintf(stderr, point_printf_format, (*this)[front[i][k]]->o);
        fprintf (stderr, "\n");
      }
    }
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
    DEBUG2 (debug_fronts());
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
                    vector_fprintf(stderr, point_printf_format,
                                   (*this)[front[l][j]]->o);
                    fprintf (stderr, "\n");
                    assert(0);
                  }
                }
              }
            });
  }

  void calc_number_of_dominating ()
  {
    int i,j;

    int size = this->size();
    
    for (i = 0; i < size; i++)  {
      int sum = 0;
      element_type * si = (*this)[i];
      for (j = 0; j < size; j++)  {
          element_type * sj = (*this)[j];
          if (sj->dominates(*si)) {
            sum++;
          }
        }
      this->dominated[i] = sum;
    }
  }

  int find_least_hv_contributor()
  {
    int worst_front = calculate_worst_front ();
    int worst_front_size = copies[worst_front];
    vector<double> hv (worst_front_size, 0);
    vector<int> position (worst_front_size  + 1, 0);
    vector<double> ref = ubound;
    const int dim = this->dimension();

    // If keep_uevs, then use ref == ubound; otherwise, ref == bound + 1.0
    if (!_keep_uevs) {
      for (int i = 0; i < int(ref.size()); i++) ref[i] += 1.0;
    }
    // The original SMS says that uevs are only kept for d == 2.
    this->calculate_uev();
    
    double * data = (double *)malloc (sizeof(double) * (worst_front_size+1) * dimension());
    
    for (int i = 0; i < worst_front_size; i++) {
      int pos = 0;
      const double * pointi = &((*this)[front[worst_front][i]]->o[0]);
      for (int j = 0; j < worst_front_size; j++) {
        const double * pointj = &((*this)[front[worst_front][j]]->o[0]);
        if (i != j) {
          memcpy (data + pos, pointj, sizeof (double) * dimension());
          pos += dim;
        }
      }
      // For d = 2, keep uev solutions.
      if (_keep_uevs && this->uev[front[worst_front][i]]) {
        assert (pointi[0] == ubound[0] || pointi[1] == ubound[1]);
        hv[i] = -1.0;
      }
      else
        hv[i] = fpli_hv(data, dimension(), worst_front_size - 1, &ref[0]);
    }
    free (data);

    DEBUG2 (debug_fronts());
    DEBUG2 (
            for (int i = 0; i < worst_front_size; i++) {
              vector_fprintf (stderr, point_printf_format, (*this)[front[worst_front][i]]->o);
              fprintf (stderr, ": hv = "point_printf_format"\n", hv[i]);
            });
    DEBUG1 (
            for (int i = 0; i < worst_front_size; i++) {
              for (int j = 0; j < worst_front_size; j++) {
                dominance_t result = (*this)[front[worst_front][i]]->dominance(*((*this)[front[worst_front][j]]));
                switch (result) {
                case IS_DOMINATED_BY:
                  assert (hv[i] > hv[j]);
                  break;
                case EQUALS:
                  if (hv[i] < 0.0)
                    assert (uev[i]);
                  else if (hv[j] < 0.0)
                    assert (uev[j]);
                  else
                    assert (fabs(hv[i] - hv[j]) < 1e-15);
                  break;
                case DOMINATES:
                  assert (hv[i] < hv[j]);
                  break;
                case NONDOMINATED:
                  // Anything here.
                  break;
                }
              }
            });

    // For d = 2, hv(ubound) < 0
    double max_hv = 0.0;
    int index_max_hv = -1;
    for (int i = 0; i < worst_front_size; i++) {
      if (max_hv < hv[i]) {
        max_hv = hv[i];
        index_max_hv = i;
      }
    }
    assert (max_hv > 0);
    assert (index_max_hv >= 0);
    return front[worst_front][index_max_hv];
  }

  int calculate_worst_front ()
  {
    int size = this->size();
    int worst_front = 1;
    assert(copies[0] > 0);
    while (worst_front < size && copies[worst_front] > 0)
      worst_front++;
    worst_front--;
    return worst_front;
  }

  int find_most_dominated ()
  {
    int worst_front = calculate_worst_front();

    calc_number_of_dominating();
    int dominated_max = front[worst_front][0];
    for (int i = 1; i < copies[worst_front]; i++) {
      if (dominated[dominated_max] < dominated[front[worst_front][i]]) {
        dominated_max = front[worst_front][i];
      }
    }
    assert (dominated[dominated_max] >= worst_front);
    return dominated_max;
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
