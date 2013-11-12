#ifndef ARCHIVERS_H
#define ARCHIVERS_H

#include <set>
#include <vector>
#include <iostream>
#include "Solution.h"
#include <cmath>

#define PI 3.14159265

class Archiver {
public:
  virtual dominance_t add(Solution nova) = 0;
  virtual void finish() = 0;
  virtual vector<Solution> getSolutions() = 0;
};

class IdealArchive: public Archiver {
private:
  unsigned int max_archiver_size;
  unsigned int dimension;
  Solution ideal;
  vector<Solution> solutions;
public:
  vector<Solution> getSolutions() {
    return solutions;
  }
  IdealArchive(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
  }
  void finish() {
    return;
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
    calculate_ideal_point();

    //If i'm here, the nova solution is non dominated by vector content and was added to it. But how is the capacity?
    
    //So we have to delete anyone. Who?
    if (solutions.size() > max_archiver_size) filter();
    
  }
private:
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
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < p.size(); ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  
  void calculate_ideal_point() {
    vector<double> idealp(dimension,99999999.0);
    for(int d = 0; d < dimension; d++) {
      for(int i = 0; i < solutions.size(); i++) {
	idealp[d] = min(idealp[d],solutions[i].o[d]);
      }
    }
    ideal.setObjectives(idealp);
  }
};

class IdealArchiveTrash: public Archiver {
private:
  unsigned int max_archiver_size;
  unsigned int dimension;
  
  Solution ideal;
  
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < p.size(); ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  vector<Solution> solutions;
public:  
  vector<Solution> getSolutions() {
    return solutions;
  }
  IdealArchiveTrash(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
  }
  dominance_t add(Solution nova) {
    vector<Solution>::iterator it = solutions.begin();
    while (it != solutions.end()) {
      bool dominates = false;
      switch (nova.dominance(*it)) {
	case IS_DOMINATED_BY:
	  trash.push_back(* (it->clone()));
	  return IS_DOMINATED_BY;
	case DOMINATES:
	  dominates = true;
	  trash.push_back(* (it->clone()));
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
    calculate_ideal_point();
    
    
    //If i'm here, the nova solution is non dominated by vector content and was added to it. But how is the capacity?
    
    //So we have to delete anyone. Who?
    if (solutions.size() > max_archiver_size) filter();
    
  }
  void finish() {
    vector<Solution>::iterator iter = trash.begin();
    while (iter != trash.end()) {
      this->addNovo(*iter);
      iter++;
    }
    calculate_ideal_point();
  }
private:
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
  void calculate_ideal_point() {
    vector<double> idealp(dimension,99999999.0);
    for(int d = 0; d < dimension; d++) {
      for(int i = 0; i < solutions.size(); i++) {
	idealp[d] = min(idealp[d],solutions[i].o[d]);
      }
    }
    ideal.setObjectives(idealp);
  } 
};

/*
 * Distributed Archiver, by BRITTO e POZO
 * 
 */

class DistributedArchiver: public Archiver {
private:
  unsigned int max_archiver_size;
  unsigned int dimension;
  
  Solution ideal;
  vector<Solution> reference_points;
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < p.size(); ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  vector<Solution> solutions;
public:

  vector<Solution> getSolutions() {
    return solutions;
  }
  DistributedArchiver(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
    reference_points.reserve(d+1);
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
    calculate_reference_points();    
    
    if (solutions.size() > max_archiver_size) filter();
  }
  void finish() {
    return;
  }
private:
  
  void filter() {
    int solution_per_region = solutions.size() / (dimension+1);
    vector<Solution> copy = solutions;
    vector<Solution> result;
    
    while (copy.size() != 1) {
      for(int i = 0; i < dimension+1; i++) {
	if (copy.size() == 1) break;
	
	int closest = 0;
	for(int it = 1; it < copy.size(); it++) {
	  if (distance(copy[it].o,reference_points[i].o) < distance(copy[closest].o,reference_points[i].o))
	    closest = it;
	}
	result.push_back(copy[closest]);
	copy.erase(copy.begin() + closest);
      }
    }
    solutions = result;
  }
  
  
  void calculate_ideal_point() {
    vector<double> idealp(dimension,99999999.0);
    for(int d = 0; d < dimension; d++) {
      for(int i = 0; i < solutions.size(); i++) {
	idealp[d] = min(idealp[d],solutions[i].o[d]);
      }
    }
    ideal.setObjectives(idealp);
  }
  void calculate_reference_points() {
    reference_points.clear();
    for(int d = 0; d < dimension; d++) {
      int best_in_dimension = 0;
      for(int i = 1; i < solutions.size(); i++) {
	if (solutions[i].o[d] < solutions[best_in_dimension].o[d]) {
	  best_in_dimension = i;
	}
      }
      reference_points.push_back(solutions[best_in_dimension]);
    }
    
    calculate_ideal_point();
    reference_points.push_back(ideal);
  }
};

class DistributedArchiverTrash: public Archiver {
private:
  unsigned int max_archiver_size;
  unsigned int dimension;
  
  Solution ideal;
  vector<Solution> reference_points;
  
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < p.size(); ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  vector<Solution> solutions;
public:
  vector<Solution> getSolutions() {
    return solutions;
  }
  DistributedArchiverTrash(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
    reference_points.reserve(d+1);
  }
  dominance_t add(Solution nova) {
    vector<Solution>::iterator it = solutions.begin();
    while (it != solutions.end()) {
      bool dominates = false;
      switch (nova.dominance(*it)) {
	case IS_DOMINATED_BY:
	  trash.push_back(* (nova.clone()));
	  return IS_DOMINATED_BY;
	case DOMINATES:
	  dominates = true;
	  trash.push_back(* (it->clone()));
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
    
    calculate_reference_points();    
    
    if (solutions.size() > max_archiver_size) filter();
  }
  void finish() {
    vector<Solution>::iterator iter = trash.begin();
    while (iter != trash.end()) {
      this->addNovo(*iter);
      iter++;
    }
  }
private:
  void filter() {
    int solution_per_region = solutions.size() / (dimension+1);
    set<int> regions;
    while (regions.size() != solutions.size()-1) {
      for (int i = 0; i < dimension+1; i++) {
	int min_distance = -1;
	//Selecionar o primeiro ponto ainda não incluido em regions.
	for(int k = 0; k < solutions.size(); k++) {
	  if (regions.count(k) == 0) {
	    min_distance = k;
	    break;
	  }
	}
	
	if (regions.size() == max_archiver_size) break;
	
	for(int k = min_distance+1; k < solutions.size(); k++) {
	  if (distance(solutions[k].o,reference_points[i].o) < distance(solutions[min_distance].o,reference_points[i].o) && regions.count(k) == 0) {
	    min_distance = k;
	  }
	}
	regions.insert(min_distance);
      }
    }
    for(int i = 0; i < solutions.size(); i++) {
      if (regions.count(i) == 0) {
	trash.push_back(* solutions[i].clone());
	solutions.erase(solutions.begin() + i);
      }
    }
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
  void calculate_ideal_point() {
    vector<double> idealp(dimension,99999999.0);
    for(int d = 0; d < dimension; d++) {
      for(int i = 0; i < solutions.size(); i++) {
	idealp[d] = min(idealp[d],solutions[i].o[d]);
      }
    }
    ideal.setObjectives(idealp);
  }
  void calculate_reference_points() {
    reference_points.clear();
    for(int d = 0; d < dimension; d++) {
      int best_in_dimension = 0;
      for(int i = 1; i < solutions.size(); i++) {
	if (solutions[i].o[d] < solutions[best_in_dimension].o[d]) {
	  best_in_dimension = i;
	}
      }
      reference_points.push_back(solutions[best_in_dimension]);
    }
    calculate_ideal_point();
    reference_points.push_back(ideal);
  }
};

/*
 * DistanceArchiver
 * 
 */

class DistanceArchiver: public Archiver {
private:
  unsigned int max_archiver_size;
  unsigned int dimension;
  
  Solution ideal;
  vector<Solution> reference_points;
  
  
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < p.size(); ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  vector<Solution> solutions;
public:
  vector<Solution> getSolutions() {
    return solutions;
  }
  DistanceArchiver(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
    reference_points.reserve(d+1);
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
    calculate_reference_points();
    
    if (solutions.size() > max_archiver_size) filter();
  }
  void finish() {
    return;
  }
private:
  void filter() {
    double smaller_distance[solutions.size()];
    for(int i = 0; i < solutions.size(); i++)
      smaller_distance[i] = 999999.;
    
    for(int i = 0; i < solutions.size(); i++) {
      for(int k = 0; k < reference_points.size(); k++) {
	smaller_distance[i] = min(smaller_distance[i],distance(solutions[i].o,reference_points[k].o));
      }
    }
    int max_distance_id = 0;
    for(int i = 1; i < solutions.size(); i++) {
      if (smaller_distance[max_distance_id] < smaller_distance[i]) max_distance_id = i;
    }
    
    solutions.erase(solutions.begin() + max_distance_id);
  }
  void calculate_ideal_point() {
    vector<double> idealp(dimension,99999999.0);
    for(int d = 0; d < dimension; d++) {
      for(int i = 0; i < solutions.size(); i++) {
	idealp[d] = min(idealp[d],solutions[i].o[d]);
      }
    }
    ideal.setObjectives(idealp);
  }
  void calculate_reference_points() {
    reference_points.clear();
    for(int d = 0; d < dimension; d++) {
      int best_in_dimension = 0;
      for(int i = 1; i < solutions.size(); i++) {
	if (solutions[i].o[d] < solutions[best_in_dimension].o[d]) {
	  best_in_dimension = i;
	}
      }
      reference_points.push_back(solutions[best_in_dimension]);
    }
    
    calculate_ideal_point();
    reference_points.push_back(ideal);
  }
};

class DistanceArchiverTrash: public Archiver {
private:
  unsigned int max_archiver_size;
  unsigned int dimension;
  vector<Solution> trash;
  Solution ideal;
  vector<Solution> reference_points;
  
  double distance(vector<double> p, vector<double> q) {
    double sum = 0.;
    for(unsigned int i=0; i < p.size(); ++i)
      sum += pow(p[i]-q[i],2);
    return sqrt(sum);
  }
  vector<Solution> solutions;
public:
  vector<Solution> getSolutions() {
    return solutions;
  }

  DistanceArchiverTrash(unsigned int ms, unsigned int d) {
    this->max_archiver_size = ms;
    this->dimension = d;
    solutions.reserve(ms+1);
    reference_points.reserve(d+1);
    trash.reserve(10);
  }
  dominance_t add(Solution nova) {
    vector<Solution>::iterator it = solutions.begin();
    while (it != solutions.end()) {
      bool dominates = false;
      switch (nova.dominance(*it)) {
	case IS_DOMINATED_BY:
	  trash.push_back(* (nova.clone()));
	  return IS_DOMINATED_BY;
	case DOMINATES:
	  dominates = true;
	  trash.push_back(* (it->clone()));
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
    calculate_reference_points();
    if (solutions.size() > max_archiver_size) filter();
  }
  void finish() {
    vector<Solution>::iterator iter = trash.begin();
    while (iter != trash.end()) {
      addNovo(*iter);
      iter++;
    }
    calculate_ideal_point();
  }
private:
  void filter() {
    double smaller_distance[solutions.size()];
    for(int i = 0; i < solutions.size(); i++)
      smaller_distance[i] = 999999.;
    
    for(int i = 0; i < solutions.size(); i++) {
      for(int k = 0; k < reference_points.size(); k++) {
	smaller_distance[i] = min(smaller_distance[i],distance(solutions[i].o,reference_points[k].o));
      }
    }
    int max_distance_id = 0;
    for(int i = 1; i < solutions.size(); i++) {
      if (smaller_distance[max_distance_id] < smaller_distance[i]) max_distance_id = i;
    }
    
    trash.push_back(* (solutions[max_distance_id].clone()));
    solutions.erase(solutions.begin() + max_distance_id);
  }
  void calculate_ideal_point() {
    vector<double> idealp(dimension,99999999.0);
    for(int d = 0; d < dimension; d++) {
      for(int i = 0; i < solutions.size(); i++) {
	idealp[d] = min(idealp[d],solutions[i].o[d]);
      }
    }
    ideal.setObjectives(idealp);
  }
  void calculate_reference_points() {
    reference_points.clear();
    for(int d = 0; d < dimension; d++) {
      int best_in_dimension = 0;
      for(int i = 1; i < solutions.size(); i++) {
	if (solutions[i].o[d] < solutions[best_in_dimension].o[d]) {
	  best_in_dimension = i;
	}
      }
      reference_points.push_back(solutions[best_in_dimension]);
    }
    
    calculate_ideal_point();
    reference_points.push_back(ideal);
  }
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
};

/*
 * ARA - By Jin & Wong.
 * 
 */

class ARArchive: public Archiver {
private:
  unsigned int dimension;
  Solution a_min;
  Solution a_max;
  
  vector<double> E;
  
  vector<Solution> references; //Called in Jin & Wong 2010 by A^(min)
  vector<Solution> next_archiver;
  vector<Solution> solutions;
public:
  vector<Solution> getSolutions() {
    vector<Solution> v = solutions;
    for(int i = 0; i < references.size(); i++) {
      bool equal = false;
      for(int j = 0; j < solutions.size() && !equal; j++) {
	if (solutions[j].dominance(references[i]) == EQUALS)
	  equal = true;
      }
      if (!equal) v.push_back(references[i]);
    }
    
    return v;
  }
  
  ARArchive(unsigned int d) {
    this->dimension = d;
    
    vector<double> bad_solution(dimension,999999.0);
    vector<double> good_solution(dimension,-1.0);
    Solution infinite;
    infinite.setObjectives(bad_solution);
    
    
    for(int i = 0; i < d; i++) {
      references.push_back(* (infinite.clone()));
      E.push_back(PI/1000.);
    }
    a_min = * (infinite.clone());
    a_max.setObjectives(good_solution);
    
  }
  dominance_t add(Solution y) {      
    if (y_dominates_Amin(y) || y.dominance(a_min) == NONDOMINATED) { //Line 1
      for(int i = 0; i < dimension; i++) { //Line 2
	if (y.o[i] < references[i].o[i]) //Line 3
	  references[i] = y; //Line 4 - RECEDES
	else if (y.dominance(references[i]) == DOMINATES)
	  references[i] = y; //Line 6 - DOMINATES
      }
//       cout << "REFORMS!" << endl;
      next_archiver.clear(); //Line 9
      vector<Solution> copy = solutions;
      for(int i = 0; i < copy.size(); i++) {
	if (! set_dominates_solution(references,solutions[i])) {
	  insert_rect(solutions[i],next_archiver);
	}
      }
      copy = references;
      for(int i = 0; i < copy.size(); i++) {
	if (! set_dominates_solution(references,references[i])) {
	  insert_rect(references[i],next_archiver);
	}
      }
    }
    else if (!set_dominates_solution(references,y)) {
      next_archiver = solutions;
      insert_rect(y,next_archiver);
    }
    solutions = next_archiver;
  }
  
  void finish() {
    solutions = getSolutions();
    return;
  }
private:
  //test if Y dominates anyone from A^min
  //if-then-else at Line 1
  bool y_dominates_Amin(Solution y) {
    for(int i = 0; i < dimension; i++) {
      if (y.dominance(references.at(i)) == DOMINATES)
	return true;
    }
    return false;
  }
  
  //test if a set has one solution that dominates Y
  bool set_dominates_solution(vector<Solution> a, Solution y) {
    for(int i = 0; i < a.size(); i++) {
      if (a[i].dominance(y) == DOMINATES) return true;
    }
    return false;
  }
  
  
  
  dominance_t compare_rectangle(vector<double> x1, vector<double> x2) {
    Solution a, b;
    a.setObjectives(x1);
    b.setObjectives(x2);
    
    return a.dominance(b);
  }
  void insert_rect(Solution x, vector<Solution> & aarc) {
    
    vector<double> where = rect(x,references);
    
    vector<Solution> d;
    for(int i = 0; i < aarc.size(); i++) {
      if (compare_rectangle(where,rect(aarc[i],references)) == DOMINATES)
	d.push_back(aarc[i]);
    }
    
    if (!d.empty()) {
//       cout << "INTERRECTDOM: Inserting " << x.o[0] << " " << x.o[1] << " " << x.o[2] << " in rect: " << where[0] << " " << where[1] << " " << where[2] << endl;
      
      aarc.push_back(x);
      aarc = remover_solucoes(aarc,d);
    }
    else {
      for(int i = 0; i < aarc.size(); i++) {
	if (compare_rectangle(rect(aarc[i],references),where) == EQUALS && x.dominance(aarc[i]) == DOMINATES) {
	  aarc[i] = x;
// 	  cout << "INTRARECTDOM: Inserting " << x.o[0] << " " << x.o[1] << " " << x.o[2] << "  in rect: " << where[0] << " " << where[1] << " " << where[2] << endl;
	  return;
	}
      }
      
      for(int i = 0; i < aarc.size(); i++) {
	if (! (compare_rectangle(rect(aarc[i],references),where) == NONDOMINATED) ){
	  //cout << "STEADYSTATE" << endl;
	  return;
	}
      }
//       cout << "OCCUPIES: Inserting " << x.o[0] << " " << x.o[1] << " " << x.o[2] << " in rect: " << where[0] << " " << where[1] << " " << where[2] << endl;
      aarc.push_back(x);
    }
  }
  vector<Solution> remover_solucoes(vector<Solution> solucoes, vector<Solution> remover) {
    for(int i = 0; i < remover.size(); i++) {
      vector<Solution>::iterator it = solucoes.begin();
      while (it != solucoes.end()) {
	if (it->dominance(remover[i]) == EQUALS)
	  solucoes.erase(it);
	else
	  it++;
      }
    }
    return solucoes;
  }
  void update_min_and_max(Solution one_more) {
    vector<Solution> total = getSolutions();
    total.push_back(one_more);
    for(int j = 0; j < dimension; j++) {
      for(int i = 0; i < total.size(); i++) {
	if (total[i].o[j] < a_min.o[j])
	  a_min.o[j] = total[i].o[j];
	
	if (total[i].o[j] > a_max.o[j])
	  a_max.o[j] = total[i].o[j];
      }
    }
  }
  vector<double> rect(Solution x, vector<Solution> a) {
    update_min_and_max(x);
    vector<double> r(dimension,-1.);
    for(int i = 0; i < dimension; i++) {
      double scale = tan(PI/2. - E[i]) / (a_max.o[i] - a_min.o[i]);
      double alfa = atan((x.o[i] - a_min.o[i]) * scale);
      r[i] = 1 + ceil(alfa/E[i]);
    }
    return r; 
  }
};



class ARArchiveTrash: public Archiver
{
private:
  unsigned int dimension;
  Solution a_min;
  Solution a_max;
  
  vector<double> E;
  
  vector<Solution> references; //Called in Jin & Wong 2010 by A^(min)
  vector<Solution> next_archiver;
  vector<Solution> solutions;
  vector<Solution> trash;
  
public:
  vector<Solution> getSolutions() {
    vector<Solution> v = solutions;
    for(int i = 0; i < references.size(); i++) {
      bool equal = false;
      for(int j = 0; j < solutions.size() && !equal; j++) {
	if (solutions[j].dominance(references[i]) == EQUALS)
	  equal = true;
      }
      if (!equal) v.push_back(references[i]);
    }
    
    return v;
  }
  
  ARArchiveTrash(unsigned int d) {
    this->dimension = d;
    
    vector<double> bad_solution(dimension,999999.0);
    vector<double> good_solution(dimension,-1.0);
    Solution infinite;
    infinite.setObjectives(bad_solution);
    
    
    for(int i = 0; i < d; i++) {
      references.push_back(* (infinite.clone()));
      E.push_back(PI/1000.);
    }
    a_min = * (infinite.clone());
    a_max.setObjectives(good_solution);
  }
  dominance_t add(Solution y) {
    trash.push_back(y);
    if (y_dominates_Amin(y) || y.dominance(a_min) == NONDOMINATED) { //Line 1
      for(int i = 0; i < dimension; i++) { //Line 2
	if (y.o[i] < references[i].o[i]) //Line 3
	  references[i] = y; //Line 4 - RECEDES
	else if (y.dominance(references[i]) == DOMINATES)
	  references[i] = y; //Line 6 - DOMINATES
      }
      //cout << "REFORMS!" << endl;
      next_archiver.clear(); //Line 9
      vector<Solution> copy = solutions;
      for(int i = 0; i < copy.size(); i++) {
	if (! set_dominates_solution(references,solutions[i])) {
	  insert_rect(solutions[i],next_archiver);
	}
      }
      copy = references;
      for(int i = 0; i < copy.size(); i++) {
	if (! set_dominates_solution(references,references[i])) {
	  insert_rect(references[i],next_archiver);
	}
      }
    }
    else if (!set_dominates_solution(references,y)) {
      next_archiver = solutions;
      insert_rect(y,next_archiver);
    }
    solutions = next_archiver;
  } 
 
  void finish() {
    for(int i = 0; i < references.size(); i++)
      addNovo(references[i]);
    
    vector<Solution>::iterator iter = trash.begin();
    while (iter != trash.end()) {
      addNovo(*iter);
      iter++;
    }
  }
  
private:
  //test if Y dominates anyone from A^min
  //if-then-else at Line 1
  bool y_dominates_Amin(Solution y) {
    for(int i = 0; i < dimension; i++) {
      if (y.dominance(references.at(i)) == DOMINATES)
	return true;
    }
    return false;
  }
  
  //test if a set has one solution that dominates Y
  bool set_dominates_solution(vector<Solution> a, Solution y) {
    for(int i = 0; i < a.size(); i++) {
      if (a[i].dominance(y) == DOMINATES) return true;
    }
    return false;
  }
  
  
  
  dominance_t compare_rectangle(vector<double> x1, vector<double> x2) {
    Solution a, b;
    a.setObjectives(x1);
    b.setObjectives(x2);
    
    return a.dominance(b);
  }
  void insert_rect(Solution x, vector<Solution> & aarc) {
    
    vector<Solution> d;
    for(int i = 0; i < aarc.size(); i++) {
      if (compare_rectangle(rect(x,references),rect(aarc[i],references)) == DOMINATES)
	d.push_back(aarc[i]);
    }
    
    if (!d.empty()) {
//       cout << "INTERRECTDOM: Inserting " << x.o[0] << " " << x.o[1] << " " << x.o[2] << " in rect: " << rect(x,references).at(0) << " " << rect(x,references).at(1) << endl;
      
      aarc.push_back(x);
      aarc = remover_solucoes(aarc,d);
    }
    else {
      for(int i = 0; i < aarc.size(); i++) {
	if (compare_rectangle(rect(aarc[i],references),rect(x,references)) == EQUALS && x.dominance(aarc[i]) == DOMINATES) {
	  aarc[i] = x;
// 	  cout << "INTRARECTDOM: Inserting " << x.o[0] << " " << x.o[1] << " " << x.o[2] << "  in rect: " << rect(x,references).at(0) << " " << rect(x,references).at(1) << endl;
	  return;
	}
      }
      bool get_out = false;
      for(int i = 0; i < aarc.size(); i++) {
	if (! (compare_rectangle(rect(aarc[i],references),rect(x,references)) == NONDOMINATED) )
	  return;
      }
//       cout << "OCCUPIES: Inserting " << x.o[0] << " " << x.o[1] << " " << x.o[2] << " in rect: " << rect(x,references).at(0) << " " << rect(x,references).at(1) << endl;
      aarc.push_back(x);
    }
  }
  vector<Solution> remover_solucoes(vector<Solution> solucoes, vector<Solution> remover) {
    for(int i = 0; i < remover.size(); i++) {
      vector<Solution>::iterator it = solucoes.begin();
      while (it != solucoes.end()) {
	if (it->dominance(remover[i]) == EQUALS) {
	  solucoes.erase(it);
	}
	else
	  it++;
      }
    }
    return solucoes;
  }  
  void update_min_and_max(Solution one_more) {
    vector<Solution> total = getSolutions();
    total.push_back(one_more);
    for(int j = 0; j < dimension; j++) {
      for(int i = 0; i < total.size(); i++) {
	if (total[i].o[j] < a_min.o[j])
	  a_min.o[j] = total[i].o[j];
	
	if (total[i].o[j] > a_max.o[j])
	  a_max.o[j] = total[i].o[j];
      }
    }
  }
  vector<double> rect(Solution x, vector<Solution> a) {
    update_min_and_max(x);
    vector<double> r(dimension,-1.);
    for(int i = 0; i < dimension; i++) {
      double scale = tan(PI/2. - E[i]) / (a_max.o[i] - a_min.o[i]);
      double alfa = atan((x.o[i] - a_min.o[i]) * scale);
      r[i] = 1 + ceil(alfa/E[i]);
    }
    return r;
  }
  dominance_t addNovo(Solution nova) {
    vector<Solution>::iterator it = solutions.begin();
    while (it != solutions.end()) {
      bool dominates = false;
      switch (nova.dominance(*it)) {
	case IS_DOMINATED_BY:
	  return IS_DOMINATED_BY;
	case DOMINATES:
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
};

#endif
