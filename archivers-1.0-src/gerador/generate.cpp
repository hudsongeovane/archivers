#include <iostream>
#include <ctime>
#include <cstdlib>
#include "Solution.h"
#include <algorithm>
#include <vector>
#include <cstring>
#include <cmath>
using std::cin;
using std::cout;
using std::endl;
using std::vector;


double rand0to1() {
  return ((double)rand()/(double)RAND_MAX);
}
double adistance(const vector<double> a, const vector<double> b) {
	double sum = 0;
	for(int i = 0; i < a.size(); i++) {
		sum += (b[i] - a[i])*(b[i] - a[i]);
	}
	return sqrt(sum);
}
//Use: ./generate <number_of_solutions> <dimension>
int main(int argc, char ** argv) {
  srand(time(0));
  
  if (argc < 4) { cout << "Use: ./generate <number_of_solutions> <dimension> <output file> [-trash <number_of_dominated_solutions>]" << endl; exit(0); }
  
  unsigned int number_of_solutions = atoi(argv[1]), dimension = atoi(argv[2]);
  
  FILE * file = fopen(argv[3],"w");
  int number_of_dominated_solutions = 0;
  if (argc > 4) {
    if (!strcmp(argv[4],"-trash")) number_of_dominated_solutions = atoi(argv[5]);
  }
  
  Solution::Initialise(dimension);
  vector<Solution> solutions;
  solutions.reserve(number_of_solutions);
  
	vector<double> pointA(3,0.4);
	vector<double> pointB;
	vector<double> pointC;
	pointB.push_back(.2); pointB.push_back(0.6); pointB.push_back(0.4);
	pointC.push_back(0.6); pointC.push_back(0.8); pointC.push_back(0.2) ;


  while (solutions.size() != number_of_solutions) {
    vector<double> o;
    for(int i = 0; i < dimension; i++) {
      o.push_back(rand0to1());
    }
    
    	cout << "add " << adistance(o,pointA) << endl;
	if ((adistance(o,pointA) < 0.1) || (adistance(o,pointB) < 0.1) || (adistance(o,pointC) < 0.1)) {
		    Solution to_add;
    	to_add.setObjectives(o);

    	solutions.push_back(to_add);
	}
  }
  /*
  vector<unsigned int> dominated;
  int count = 0;
  for(vector<Solution>::iterator it = solutions.begin(); it != solutions.end(); it++) {
    bool isdominated = false;
    for(int i = 0; i < solutions.size() && !isdominated; i++) {
      if (it->dominance(solutions[i]) == IS_DOMINATED_BY) { dominated.push_back(count); isdominated = true; }
    }
    count++;
  }
  
  for(int i = 0; i < dominated.size(); i++) {
    bool isdominated_or_domines = true;
    Solution to_replace;
    while (isdominated_or_domines) {
      vector<double> o;
      for(int j = 0; j < dimension; j++) {
		o.push_back(rand0to1());
      }
      to_replace.setObjectives(o);
      
      for(int k = 0; k < solutions.size() && !isdominated_or_domines; k++) {
		if (to_replace.dominance(solutions[k]) != NONDOMINATED) isdominated_or_domines = true;
      }
    }
    solutions[dominated[i]] = to_replace;
    cout << "Replacing " << dominated[i] << endl;
  }
  
  //! Adding dominated solutions
  for(int i = 0; i < number_of_dominated_solutions; i++) {
    vector<double> o;
    for(int i = 0; i < dimension; i++) {
      o.push_back(rand0to1());
    }
    Solution to_add;
    to_add.setObjectives(o);
    
    
    for(int k = 0; k < solutions.size(); k++) {
		if (to_add.dominance(solutions[k]) == IS_DOMINATED_BY) { solutions.push_back(to_add); break; }
    }
  } */
  random_shuffle(solutions.begin(),solutions.end());
  
  for(int i = 0; i < solutions.size(); i++) {
    for(int j = 0; j < dimension; j++) {
      fprintf(file,"%f",solutions[i].o[j]);
      if (j != dimension-1) fprintf(file," ");
    }
    fprintf(file,"\n");
  }
  
}
