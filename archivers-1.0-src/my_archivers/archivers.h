#ifndef ARCHIVERS_H
#define ARCHIVERS_H

#include <vector>
#include <iostream>
#include "Solution.h"
#include <cmath>


class IdealArchive {
public:
  vector<Solution> solutions;
  unsigned int max_archiver_size;
  unsigned int dimension;
  
  Solution ideal;
  
  IdealArchive(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
  }
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < this->dimension; ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  dominance_t add(Solution nova) {
    vector<Solution>::iterator it = solutions.begin();
    while (it != solutions.end()) {
      bool dominates = false;
      switch (nova.dominance(*it)) {
	case IS_DOMINATED_BY:
	  return IS_DOMINATED_BY;
	case DOMINATES:
	  dominates = true;
	  solutions.erase(it);
	  break;
	case NONDOMINATED:
	  it++;
	  break;
	case EQUALS:
	  return EQUALS;
      }
    }
    solutions.push_back(nova);
    
    
    if (solutions.size() == 1) {
      ideal = nova;
    }
    else {
      for(int k = 0; k < dimension; k++) {
	ideal.o[k] = min(ideal.o[k],nova.o[k]);
      }
    }
    
    //If i'm here, the nova solution is non dominated by vector content and was added to it. But how is the capacity?
    
    //So we have to delete anyone. Who?
    if (solutions.size() > max_archiver_size) filter();
    
  }
  void filter() {
    vector<Solution>::iterator to_remove;
    double max_distance = -1.;
    for(vector<Solution>::iterator it = solutions.begin(); it != solutions.end(); it++) {
      double d = distance(it->o,ideal.o);
      if (d > max_distance) {
	to_remove = it;
	max_distance = d;
      }
    }
    solutions.erase(to_remove);
  }
};

class IdealArchiveTrash {
public:
  vector<Solution> solutions;
  unsigned int max_archiver_size;
  unsigned int dimension;
  
  Solution ideal;
  
  IdealArchiveTrash(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
  }
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < this->dimension; ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  dominance_t add(Solution nova) {
    vector<Solution>::iterator it = solutions.begin();
    while (it != solutions.end()) {
      bool dominates = false;
      switch (nova.dominance(*it)) {
	case IS_DOMINATED_BY:
	  return IS_DOMINATED_BY;
	case DOMINATES:
	  dominates = true;
	  solutions.erase(it);
	  break;
	case NONDOMINATED:
	  it++;
	  break;
	case EQUALS:
	  return EQUALS;
      }
    }
    solutions.push_back(nova);
    
    
    if (solutions.size() == 1) {
      ideal = nova;
    }
    else {
      for(int k = 0; k < dimension; k++) {
	ideal.o[k] = min(ideal.o[k],nova.o[k]);
      }
    }
    
    //If i'm here, the nova solution is non dominated by vector content and was added to it. But how is the capacity?
    
    //So we have to delete anyone. Who?
    if (solutions.size() > max_archiver_size) filter();
    
  }
  void filter() {
    vector<Solution>::iterator to_remove;
    double max_distance = -1.;
    for(vector<Solution>::iterator it = solutions.begin(); it != solutions.end(); it++) {
      double d = distance(it->o,ideal.o);
      if (d > max_distance) {
	to_remove = it;
	max_distance = d;
      }
    }
    trash.push_back(* to_remove);
    solutions.erase(to_remove);
  }
  vector<Solution> trash;
  
  dominance_t addNovo(Solution nova) {
    vector<Solution>::iterator it = solutions.begin();
    while (it != solutions.end()) {
      bool dominates = false;
      switch (nova.dominance(*it)) {
	case IS_DOMINATED_BY:
	  return IS_DOMINATED_BY;
	case DOMINATES:
	  dominates = true;
	  solutions.erase(it);
	  break;
	case NONDOMINATED:
	  it++;
	  break;
	case EQUALS:
	  return EQUALS;
      }
    }
    solutions.push_back(nova);
  }
  
  void finish() {
    vector<Solution>::iterator iter = trash.begin();
    while (iter != trash.end()) {
      this->addNovo(*iter);
      iter++;
    }
  }
  
};

#endif