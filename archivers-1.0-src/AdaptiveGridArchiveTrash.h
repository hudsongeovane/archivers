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

 Based on code from:

   =======================================
   AGA - Adaptive grid archiving algorithm
   2/12/2002 (C) Joshua Knowles
   
   For the PISA interface tool
   =======================================   
 ----------------------------------------------------------------------

*************************************************************************/
#ifndef _AGA_ARCHIVE_TRASH_H_
#define _AGA_ARCHIVE_TRASH_H_

#include "Archive.h"
#include "Random.h"
#include <limits.h>

template<class T>
class AdaptiveGridArchiveTrash : public BaseArchive<T> {

public:
 
  typedef typename BaseArchive<T>::element_type element_type;
  typedef typename BaseArchive<T>::iterator iterator;
  typedef typename BaseArchive<T>::const_iterator const_iterator;
  using BaseArchive<T>::ubound;
  using BaseArchive<T>::lbound;
  using BaseArchive<T>::dimension;
  using BaseArchive<T>::uev;

  class Region {
  public:
    Region() : index(0), pop(0) {}
    int index; // FIXME: Unused ???
    int pop;
  }; 

  Random &rng;

  int num_cr;  // number of crowded regions
  int grid_levels; // quadtree `levels' of the grid
  vector<Region> grid;
  vector<Region> f; // FIXME: Not used?

  vector<int> crowded_region;
  vector<double> old_ub; // FIXME: Unused?
  vector<double> old_lb; // FIXME: Unused?
  int fall; //FIXME: Unused?
  int grid_updates;

  class AdaptiveGridArchiveElementData : public BaseArchiveElementData
  {
  public:

    AdaptiveGridArchiveElementData() : index(-1), population(-1), region(-1) { }
    ~AdaptiveGridArchiveElementData() { }

    int index; // FIXME: Never used?
    int population; // FIXME: Never used?
    int region; // FIXME: Never used?
  };

  AdaptiveGridArchiveTrash (unsigned int maxsize, unsigned int dimension, Random & rng, int grid_levels = 5)
    : BaseArchive<T>(maxsize + 1, dimension),
      rng (rng),
      grid_levels (grid_levels),
      grid (int(pow(2.0, (double)grid_levels * dimension)), Region()),
      crowded_region (maxsize+1, 0),
      old_ub (dimension, 0), old_lb(dimension, 0),
      fall (0), grid_updates (0)
  { 
    if (this->max_size() <= 2 * dimension) {
      fprintf(stderr,  "warning: the archive capacity N should be greater than 2 * dimension to ensure that all uniquely extremal vectors can be accommodated. You may continue but algorithm could fail to operate as expected.\n");
    }
    assert (dimension > 1);
    assert (grid_levels * dimension <= 32);
    // We reserve above one more than the maxsize, but set the limit correctly here.
    this->max_size(maxsize);
  }

  // See [Knowles2002PhD].
  static int default_grid_levels (unsigned int maxsize, unsigned int dimension)
  {
    assert (maxsize > 2.0);
    assert (dimension > 1);
    int levels;
    if (dimension == 2)
      levels = int(floor (log2(maxsize/2.0)));
    else if (dimension == 3) {
      if (maxsize <= 50)
        levels = 2;
      else if (maxsize > 50 && maxsize <= 200)
        levels = 3;
      else if (maxsize > 200 && maxsize <= 800)
        levels = 4;
      else /*if (maxsize > 800)*/
        levels = 5;
    } else if (dimension == 4) {
      if (maxsize < 50)
        levels = 1;
      else if (maxsize >= 50 && maxsize <= 350)
        levels = 2;
      else if (maxsize > 350)
        levels = 3;
    } else /*(dimension > 4)*/
      levels = 1;

    DEBUG2_FUNPRINT("grid_levels = %d\n", levels);
    assert (dimension * levels <= 32);
    return levels;
  }

  static AdaptiveGridArchiveElementData & getData(const element_type &s) {
    return *(static_cast<AdaptiveGridArchiveElementData*>(s.data));
  }

  dominance_t add (const T & a)
  {
    element_type s (a);

    bool extends = false;

    //1. if archive is full update_ranges and check if new point extends
    //  2. if (result != -1) i.e. new point not dominated
    //     3. if (result == 1) i.e. if new point dominates
    //        4. archive new point and remove all dominated
    //     5. else if (archive is full)
    //        6. if (new point extends the ranges)
    //           7. add the new point and remove a non-extremal from the most crowded
    //	8. else if (new point not in a most crowded region) 
    //           9. add the new point and remove a non-extremal from the most crowded
    //     10 else
    //          11. archive the new point
    
    if (this->size() >= this->max_size())
      {
        this->calculate_bounds(); //calculate ranges of grid NOT INCLUDING the new point
        // printf("ranges calculated!\n");
        
        // Now, set a flag if the new point would extend the range of the grid in any dimension
        extends = extends_ranges_p (s);
      }

    // Call BaseArchive update procedure: removes dominated solutions if found.
    dominance_t result = this->update(s);
    if (result != NONDOMINATED) {
      if (result == DOMINATES)
        this->push_back(s); // Add the new point to the archive.
      else
        assert (result == IS_DOMINATED_BY || result == EQUALS);
      assert (this->size() <= this->max_size());
      return result;
    }

    if (this->size() < this->max_size()) {  // if arcsize < N
      //  printf("And I should be added\n");
      this->push_back (s);
      assert (this->size() <= this->max_size());
      return result;
    }

    // Archive capacity reached.

    grid_updates++;
    //	  printf("grid updates needed=%d\n",updates);
    this->update_bounds (s);
    int old_size = this->size();
    this->push_back (s); // add the new point (just temporarily to update the grid)
    update_grid();

    if (extends) { // if the new point extends the range of the objective space
      //	      printf("and I extend..\n");
      // getchar();
      int pos = get_index_nue(old_size); // then select a single non-uniquely extremal vector
      
      trash.push_back(* new ArchiveElement<T>(* this->at(pos)));
      
      this->erase(this->begin() + pos); // and remove it
      // for(i=0;i<arcsize+1;i++)
      //	printf("%d",kill[i]);
      // printf("\n");
      // Since the point was already archived, we are done.
      assert (this->size() <= this->max_size());
      return result;
    } else {
      //  printf("$$$$$$$$$$$$$$$$$$$$$$$$$ num_cr=%d\n", num_cr);
      // printf("$$$$$$$$$$$$$$$$$$$$$$$$$ location = %d, crowded region0=%d\n", find_loc(newpoint),crowded_region[0]);
      if (in_crowded_p (s)) {
        //  printf("and I am not in a most crowded region\n");
        // getchar();
        int pos = get_index_nue(old_size);
	
	
	trash.push_back(* new ArchiveElement<T>(* this->at(pos)));
	
        this->erase(this->begin() + pos); // and remove it
        // Since the point was already archived, we are done.
        assert (this->size() <= this->max_size());
        return result;
      } else {
        //	  printf("I suppose I must be in a most crowded region!!!!\n");
        // don't archive the point
        // Since the point was temporarily archived, delete it.
        this->pop_back();
        assert (this->size() <= this->max_size());
        return result;
      }
    }
  }

private:

  void 
  push_back (const element_type &s)
  {
    element_type * a = new element_type(s);
    a->data = new AdaptiveGridArchiveElementData();
    this->BaseArchive<T>::push_back(a);
  }

  bool extends_ranges_p (const element_type &a)
  {
    int dim = this->dimension();
    for(int j = 0; j < dim; j++) {
      if(a.o[j] > ubound[j])
        return true;
      if(a.o[j] < lbound[j])
        return true;
    }
    return false;
  }

  int find_loc(const element_type &a)
  {
    const int div = (int) pow(2.0, grid_levels);
    const int num_grids = int(pow(2.0, (int) dimension() * grid_levels));
    int loc = 0;
    
    for (int i = 0; i < int(dimension()); i++) {
      //printf("%g, %g, %g\n", a.o[i], ubound[i], lbound[i]);
      int s = int((a.o[i] - lbound[i]) / (ubound[i] - lbound[i]) * (div - 1));
      // printf("s=%d div=%d\n", s, div);
      loc += s * (int)pow(double(div), i);
    }
    DEBUG2_FUNPRINT("square = %d\n", loc);
    assert (loc >= 0);
    // FIXME: This doesn't pass.
    assert (loc < num_grids);
    if (loc >= num_grids)
      return num_grids - 1;

    return loc;

    /*
    // find the grid location of a solution given a vector of its objective values

    int loc = 0;
    int d;
    int n = 1;  
    int i;
    int div;  
    int inc[MAX_OBJ];
    double width[MAX_OBJ];
    double offset[MAX_OBJ];
    
    div = (int)pow(2.0,l);
    
    for (i = 0; i < k; i++)
    {
      offset[i]=lbound[i]-(1.0/(2*div))*(ubound[i]-lbound[i]);
      inc[i] = n;
      n *=2;
      width[i] = (1.0+(1.0/(2*div)))*(ubound[i]-lbound[i]);
    }
	    
    for (d = 1; d <= l; d++)
    {
      for (i = 0; i < k; i++)
	{
	  if(eval[i] < width[i]/2+offset[i])
	    loc += inc[i];
	  else
	    offset[i] += width[i]/2;
	}
      for (i = 0; i < k; i++)
	{
	  inc[i] *= (k *2);
	  width[i] /= 2;
	}
    }
    return(loc);
    */
  } 

  bool in_crowded_p (const element_type &p)
  {
    int loc = find_loc(p);
    for (int i = 0; i < this->num_cr; ++i)
      if (loc == crowded_region[i])
        return true;
    return false;
  }

  int get_index_nue(int arcsize)
  {
    // get any point from a crowded region except if it is uniquely extremal
    // uses the global variable, num_cr
    int j = 0;
    int idx_cr[arcsize];
    int num = 0;
    int nuev = 0;
    for (j = 0; j < arcsize; j++) {
      if (!in_crowded_p (*((*this)[j])))
        continue;
      if (this->uev[j] == 1) {
        nuev++;
        continue;
      }
      idx_cr[num] = j;
      num++;
    }
    assert (num > 0);
    int rand_pos = rng.rand_int (num);
    assert (rand_pos >= 0);
    assert (rand_pos < num);

    j = idx_cr[rand_pos];

    assert (this->uev[j] != 1);
    assert (in_crowded_p (*((*this)[j])));
    return j;
  }

  void update_grid()
  {
    // this function uses the following global variables:
    // int num_cr
    int a, i;
    int square;
    //double offset[MAX_OBJ];
    int dim = dimension();
    vector<double> ub (dim);
    vector<double> lb (dim);
    int cf_max; // maximum crowding factor
    
    int div = (int) pow(2.0, this->grid_levels);
    int arcsize = this->size();
    //  printf("l=%d k=%d div=%d\n",l, k, div);
    //  printf("maxima:\n");
    // print_double_vector(ubound,k);
    // printf("minima:\n");
    // print_double_vector(lbound,k);
  
    for (a = 0; a < dim; a++) {
      ub[a] = ubound[a] + (1.0 / (2 * div)) * (ubound[a] - lbound[a]);
      lb[a] = lbound[a] - (1.0 / (2 * div)) * (ubound[a] - lbound[a]);
    }
    
    // if((old_ub != ub || old_lb != lb))
    //   {
    //      fprintf(stderr,"%d ",iter);
    //      for(a=0;a<k;a++)     
    //	fprintf(stderr,"%lf %lf ",lbound[a], ubound[a]);   
    //      fprintf(stderr,"\n");
    //     }

    this->calculate_uev();

    int num_grids = int(pow(2.0, dim * grid_levels));
    assert (num_grids == (int)grid.size());
    for (i = 0; i < num_grids; i++) {
      grid[i].pop = 0;
    }
  
    for (a = 0; a < arcsize; a++)  {
      square = find_loc(*((*this)[a]));
      //  printf("square=%d\n", square);
      // FIXME: this is never used so commented out
      // getData((*this)[a]).region = square;
      if (uev[a] != 1)
        grid[square].pop++;
    }

    cf_max = -1;
    for (a = 0; a < num_grids; a++)  {
      //      printf("gridpop[%d]=%d\n", a, grid[a].pop);
      if (grid[a].pop > cf_max)
        cf_max = grid[a].pop;
    }
    
    num_cr = 0;
    for (a = 0; a < (int)this->max_size(); a++)
      crowded_region[a] = -1;
    
    for (a = 0; a < num_grids; a++)   {
      if (grid[a].pop == cf_max)	{
        crowded_region[num_cr] = a;
        num_cr++;
      }
    }
    
    old_lb = lb;
    old_ub = ub;
    //  printf("Number of crowded regions = %d\n", num_cr);
    // for(i=0;i<num_cr;i++)
    // printf("%d\n", crowded_region[i]);
  }
  
  
  
  dominance_t addNovo(const T &s) {
    dominance_t result = this->update(s);
    if ((result == DOMINATES || result == NONDOMINATED)) {
      this->push_back(s);
    }
    return result;
  }

  void finish() {
    typename vector<T>::iterator iter = trash.begin();
    while (iter != trash.end()) {
      this->addNovo(*iter);
      iter++;
    }
  }
  vector<T> trash;

};

#endif

