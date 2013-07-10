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
  SPEA2 - Strength Pareto EA 2

  Implementation in C for the selector side.
  
  Implements Petri net.
  
  file: spea2.c
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by: Stefan Bleuler, bleuler@tik.ee.ethz.ch
  ========================================================================

*************************************************************************/
#ifndef _ARCHIVE_SPEA2_H_
#define _ARCHIVE_SPEA2_H_

#include "Archive.h"
#include "Random.h"

template<class T>
class SPEA2Archive : public BaseArchive<T> {

public:
  typedef typename BaseArchive<T>::element_type element_type;
  typedef typename BaseArchive<T>::iterator iterator;
  typedef typename BaseArchive<T>::const_iterator const_iterator;

  enum ind_valid { IND_NOT_VALID = 0, IND_VALID };

  class SPEA2ArchiveElementData : public BaseArchiveElementData
  {
  public:
    ~SPEA2ArchiveElementData() { }
    SPEA2ArchiveElementData() : flag(IND_VALID), fitness(0) { }

    unsigned int flag;
    double fitness;
  };

  SPEA2Archive (unsigned int maxsize, unsigned int dim, Random &rng)
    : BaseArchive<T>(maxsize + 1, dim),
      fitness_bucket ((maxsize + 1) * (maxsize + 1)),
      fitness_bucket_mod (maxsize + 1),
      copies (maxsize + 1, 1),
      NN (maxsize + 1, vector<int> (maxsize + 1)),
      dist (maxsize + 1, vector<double> (maxsize + 1)),
  //  double  *f_max;
  //  double  *f_min;
  //  double  *f_norm;
      strength (maxsize + 1),
      rng(rng)
  { // We reserve maxsize + 1 but then set max_size the real value
    this->max_size(maxsize);
  }

  static SPEA2ArchiveElementData & getData(const element_type &s) {
    return *(static_cast<SPEA2ArchiveElementData*>(s.data));
  }

  dominance_t add (const T & a)
  {
    element_type s (a);
    dominance_t result = NONDOMINATED;

    this->push_back(s);

    calcFitnesses(); /* Calculates SPEA2 fitness values for all
                        individuals */

    calcDistances(); /* Calculates distance matrix dist[][] */

    environmentalSelection(); /* Performs environmental selection
                                 (truncates 'pp_all' to size 'alpha')*/

    DEBUG2_FUNPRINT("archive size = %5d\n", this->size());

    return result;
  }

private:
  // SPEA2 internal global variables
  vector<int>  fitness_bucket;
  vector<int>  fitness_bucket_mod;
  vector<int>  copies;
  vector<int>  old_index;
  vector< vector<int> > NN;
  vector< vector<double> > dist;
  //  double  *f_max;
  //  double  *f_min;
  //  double  *f_norm;
  vector<int> strength;
  Random &rng;

  void 
  push_back (const element_type &s)
  {
    element_type * a = new element_type(s);
    a->data = new SPEA2ArchiveElementData();
    this->BaseArchive<T>::push_back(a);
  }

  void
  calcFitnesses()
  {
    int i, j;
    int size;

    size = this->size();

    // initialize fitness and strength values
    for (i = 0; i < size; i++)      {
      element_type & s = *((*this)[i]);
      getData(s).fitness = 0;
      strength[i] = 0;
      fitness_bucket[i] = 0;
      fitness_bucket_mod[i] = 0;
      for (j = 0; j < size; j++)
        {
          fitness_bucket[i * size + j] = 0;
        }
    }
    
    // calculate strength values
    for (i = 0; i < size; i++)  {
      element_type * si = (*this)[i];
      for (j = 0; j < size; j++)
        {
          element_type * sj = (*this)[j];
          if (si->dominates(*sj))
            {
              strength[i]++;
            }
        }
    }

    // Fitness values =  sum of strength values of dominators
    for (i = 0; i < size; i++)  {
      int sum = 0;
      element_type * si = (*this)[i];
      for (j = 0; j < size; j++)
        {
          element_type * sj = (*this)[j];
          if (sj->dominates(*si))
            {
              sum += strength[j];
            }
        }
      getData(*si).fitness = sum;
      fitness_bucket[sum]++;
      fitness_bucket_mod[(sum / size)]++;
#if DEBUG > 2
      fprintf(stderr, "%.3f\t", pp_all->ind_array[i]->fitness);
      for (j = 0; j < NUM_OBJ; j++)
        fprintf(stderr, "%.2f ", pp_all->ind_array[i]->f[j]);
      SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
      fprintf(stderr, "\n");
#endif
    }
  }


  double calcDistance(const element_type &a, const element_type &b)
  {
    int i;
    int dim = this->dimension();
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
    int i, j;
    int size = this->size();

    // initialize copies[] vector and NN[][] matrix
    for (i = 0; i < size; i++) {
      copies[i] = 1;
      for (j = 0; j < size; j++) {
        NN[i][j] = -1;
      }
    }

    // calculate distances
    for (i = 0; i < size; i++)  {
      NN[i][0] = i;
      for (j = i + 1; j < size; j++) {
        dist[i][j] = calcDistance(*((*this)[i]), *((*this)[j]));
        assert(dist[i][j] < HUGE_VAL);
        dist[j][i] = dist[i][j];
        if (dist[i][j] == 0) {
          NN[i][copies[i]] = j;
          NN[j][copies[j]] = i;
          copies[i]++;
          copies[j]++;
        }
      }
      dist[i][i] = 0;
    }
  }

  int getNN(int index, int k)
  /* lazy evaluation of the k-th nearest neighbor
     pre-condition: (k-1)-th nearest neigbor is known already */
  {
    assert(index >= 0);
    assert(k >= 0);
    assert(copies[index] > 0);
    int size = this->size();

    if (NN[index][k] < 0)
      {
        int i;
        double min_dist = HUGE_VAL;
        int min_index = -1;
        int prev_min_index = NN[index][k-1];
        double prev_min_dist = dist[index][prev_min_index];
        assert(prev_min_dist >= 0);

        for (i = 0; i < size; i++) {
          double my_dist = dist[index][i];

          if (my_dist < min_dist && index != i)  {
                if (my_dist > prev_min_dist 
                    || (my_dist == prev_min_dist && i > prev_min_index))
                  {
                    min_dist = my_dist;
                    min_index = i;
                }
            }
        }
        NN[index][k] = min_index;
    }
    return NN[index][k];
}

  double getNNd(int index, int k)
  /* Returns the distance to the k-th nearest neigbor
     if this individual is still in the population.
     For for already deleted individuals, returns -1 */
  {
    int neighbor_index = getNN(index, k);

    if (copies[neighbor_index] == 0)
      return (-1);
    else
      return (dist[index][neighbor_index]);
  }

  void truncate_nondominated()
  /* truncate from nondominated individuals (if too many) */
  {
    int i;
    int size = this->size();
    int max_size = this->max_size();
    /* delete all dominated individuals */
    /* ??? There should not be any dominated individuals?  */
    for (int i = 0; i < size; ++i) {
      element_type &s = *((*this)[i]);
      if (getData(s).fitness > 0) {
        getData(s).flag = IND_NOT_VALID;
        copies[i] = 0;
        fprintf(stderr, "dominated: ");
        s.print(stderr);
        fprintf(stderr, "\n");
      }
    }

    // truncate from non-dominated individuals
    while (fitness_bucket[0] > max_size)
    {
      int max_copies = 0;
      int count = 0;
      int delete_index;
      vector<int> marked;
      marked.reserve(size);

      // compute inds with maximal copies
      for (i = 0; i < size; i++) {
        if (copies[i] > max_copies)
          {
            count = 0;
            max_copies = copies[i];
          }
        if (copies[i] == max_copies)
          {
            marked[count] = i;
            count++;
          }
      }
             
      assert(count >= max_copies);

      if (count > max_copies)
        {
          vector<int> neighbor (count, 1);
          while (count > max_copies)
            {
              double min_dist = HUGE_VAL;
              int count2 = 0;
              for (i = 0; i < count; i++)
                {
                  double my_dist = -1;
                  while (my_dist == -1 && neighbor[i] < size)
                    {
                      my_dist = getNNd(marked[i],neighbor[i]);
                      neighbor[i]++;
                    }
                  
                  if (my_dist < min_dist)
                    {
                      count2 = 0;
                      min_dist = my_dist;
                    }
                  if (my_dist == min_dist)
                    {
                      marked[count2] = marked[i];
                      neighbor[count2] = neighbor[i];
                      count2++;
                    }
                }
              count = count2;
              if (min_dist == -1) // all have equal distances
                {
                  break;
                }
            }
        }

      // remove individual from population
      int temp_rand = rng.rand_int (0, count - 1);
      assert((temp_rand >= 0) && (temp_rand < size));
      delete_index = marked[temp_rand];

      getData(*((*this)[delete_index])).flag = IND_NOT_VALID;

#if DEBUG >= 3
      fprintf(stderr, "removed: ");
      (*this)[delete_index]->print(stderr);
      fprintf(stderr, "\n");
#endif

      for (i = 0; i < count; i++)
        {
          if (dist[delete_index][marked[i]] == 0)
            {
              copies[marked[i]]--;
            }
        }
      copies[delete_index] = 0; // Indicates that this index is empty
      fitness_bucket[0]--;
      fitness_bucket_mod[0]--;
    }
  }

  void truncate_dominated()
  /* truncate from dominated individuals */
  {
    int i, j;
    int num = 0;
    int size = this->size();
    int max_size = this->max_size();

    i = -1;
    while (num < max_size)
      {
        i++;
        num += fitness_bucket_mod[i];
      }

    j = i * size;
    num = num - fitness_bucket_mod[i] + fitness_bucket[j];
    while (num < max_size)
      {
        j++;
        num += fitness_bucket[j];
      }

    if (num == max_size) {
      for (i = 0; i < size; i++)        {
        element_type * si = (*this)[i];
        if (getData(*si).fitness > j) {
          getData(*si).flag = IND_NOT_VALID;
#if DEBUG > 2
          fprintf(stderr, "dominated2: ");
          SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
          fprintf(stderr, "\n");
#endif
        }
      }
    } else {// if not all fit into the next generation
      int k;
      int fill_level = 0;

      int free_spaces = max_size - (num - fitness_bucket[j]);
      vector<int> best (free_spaces, 0);

      for (i = 0; i < size; i++)       {
        element_type * si = (*this)[i];
        if (getData(*si).fitness > j)          {
          getData(*si).flag = IND_NOT_VALID;
#if DEBUG > 2
          fprintf(stderr, "removed2: ");
          SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
          fprintf(stderr, "\n");
#endif
        } else if (getData(*si).fitness == j) {
          if (fill_level < free_spaces) {
            best[fill_level] = i;
            fill_level++;
            for (k = fill_level - 1; k > 0; k--)                  {
              int temp;
              if (getNNd(best[k], 1) <= getNNd(best[k - 1], 1))                 
                break;
              temp = best[k];
              best[k] = best[k-1];
              best[k-1] = temp;
            }
          } else   {
            if (getNNd(i, 1) <= getNNd(best[free_spaces - 1], 1)) {
              getData(*si).flag = IND_NOT_VALID;
#if DEBUG > 2
              fprintf(stderr, "removed3: ");
              SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
              fprintf(stderr, "\n");
#endif
            } else {
              element_type * b = (*this)[best[free_spaces - 1]];
              getData(*b).flag=IND_NOT_VALID;
#if DEBUG > 2
              fprintf(stderr, "removed4: ");
              SolPrintOneLine(stderr, pp_all->ind_array[i]->solution);
              fprintf(stderr, "\n");
#endif
              best[free_spaces - 1] = i;
              for (k = fill_level - 1; k > 0; k--) {
                int temp;
                if (getNNd(best[k], 1) <= getNNd(best[k - 1], 1)) 
                  break;
                temp = best[k];
                best[k] = best[k-1];
                best[k-1] = temp;
              }
            }
          }
        }
      }
    }
  }

  void
  environmentalSelection()
  {
    if (fitness_bucket[0] > (int) this->max_size())
    {
      truncate_nondominated();
    }
    else if (this->size() > this->max_size())
    {
      truncate_dominated();
    }

    // Move remaining individuals to top of archive
    // FIXME: A smarter erase would move the deleted ones to the bottom.

    int size = this->size();
    for (int i = 0; i < size;) {
      element_type * tmp = (*this)[i];
      if (getData(*tmp).flag != IND_NOT_VALID) {
        assert(copies[i] > 0);
        ++i;
      } else {
        (*this)[i] = this->back();
        this->back() = tmp;
        this->erase(this->end() - 1);
        --size;
        copies[i] = copies[size];
      }
    }
    assert(this->size() <= this->max_size());
  }

};

#endif

