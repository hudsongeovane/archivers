#include <iostream>
#include <cstdlib>
#include "archivers.h"
#include <string.h>
#include <ctime>

bool strtovector_double (char * str, vector<double> &p) {
  Solution s;
  char * endp = str;
  char * cursor = NULL;

  p.reserve (20);

  do {
    cursor = endp;
    double value = strtod (cursor, &endp);
    if (cursor == endp) break;
    p.push_back (value);
  } while (1);

  // not end of string: error
  while (*cursor != '\0') {
    if (!isspace(*cursor)) { return false; }
    cursor++;
  }
  
  if (p.size() == 0) return false;
  return true;
}

bool readSolution (FILE *stream, Solution &s) {
  char buffer[500];
  vector<double> v;
  static bool first_time = true;

  if (!fgets (buffer, 499, stream))
    return false;

  if (!strtovector_double(buffer, v))
    return false;

  if (first_time) {
    /* Set the number of objectives.  */
    Solution::Initialise (v.size());
    first_time = false;
  }

  s.setObjectives(v);
  return true;
}
void printError() {
  cout << "Error. Use ./archive <input file> <archiver> [options]" << endl;
  cout << "\n<archiver> can be:" << endl;
  cout << "ideal -> Ideal Archiver" << endl;
  cout << "distance -> Distance Archiver" << endl;
  cout << "distributed -> Distributed Archiver" << endl;
  cout << "ara -> Adaptive Rectangle Archiver" << endl;
  cout << "\nOptions:" << endl;
  cout << "-trash -> Enables recycling removed objects at the end of reading input file" << endl;
  cout << "-n <MAX_SIZE> -> Limits the maximum size of the archiver. Value 100 by default" << endl;
  cout << "-hide -> The final output is not shown in sdtout" << endl;
  cout << "-check -> Check deterioration at the end of process" << endl;
  cout << "-random <number of points> -> Does not use input file, but generate N points randomly" << endl;
  cout << "-max_points <number of points> -> Limit the number of points in the input to N\n" << endl;
  exit(0);
}

double rand0to1() {
  return ((double)rand()/(double)RAND_MAX);
}
bool deteriorated(vector<Solution> setone, vector<Solution> settwo) {
	for(int i = 0; i < setone.size(); i++) {
		for(int j = 0; j < settwo.size(); j++) {
			if (setone[i].dominance(settwo[j]) == DOMINATES) return true;
		}
	}
}
bool isBetterThan(vector<Solution> setone, vector<Solution> settwo) {
  int dominated_elements = 0;
  int equal_elements = 0;
  int non_dominated = 0;
  for(int i = 0; i < settwo.size(); i++) {
    bool get_out = false;
    bool found_one_equal = false;	
    bool found_nondominated = false;	
    dominance_t d;
    for(int j = 0; j < setone.size() && !get_out; j++) {
      d = setone[j].dominance(settwo.at(i));
      if (d == DOMINATES) {
	get_out = true;
	cout << "DETERIORATE!" << endl;
	dominated_elements++;
      }
      else if (d == EQUALS)
	found_one_equal = true;
      else if (d == IS_DOMINATED_BY)
	return false;
      else if (d == NONDOMINATED)
	found_nondominated = true;
    }
    if (!get_out && found_one_equal) equal_elements++;
    else if (!get_out && found_nondominated) return false;
  }
  return (dominated_elements > 0);
}


int main(int argc, char * argv[]) {
  clock_t t;
  t = clock();
  
  srand(time(0));
  if (argc < 3) {
    printError();
  }
  int max_size = 100, dimension = 0;
  FILE * stream = fopen(argv[1], "r");
  
  bool trash = false, hide = false, paren = false, randomize = false, checkDeteriorate = false;
  /* ideal = 0
   * distance = 1
   * distributed = 2
   * ara = 3
   */
  int n_archiver = -1;
  int random_number = -1;
  int max_points = 9999999;
  const char* archivers[] = {"ideal","distance","distributed","ara"};
  for(int i = 0; i < 4 && n_archiver == -1; i++) {
    if (!strcmp(archivers[i],argv[2])) {
      n_archiver = i;
    }
  }
  if (n_archiver == -1) {
    printError();
  }
  for(int i = 3; i < argc; i++) {
    if (!strcmp("-trash",argv[i])) 
      trash = true;
    else if (!strcmp("-n",argv[i])) {
      max_size = atoi(argv[++i]);
    }
    else if (!strcmp("-hide",argv[i])) {
      hide = true;
    }
    else if (!strcmp("-p",argv[i])) {
      paren = true;
    }
    else if (!strcmp("-random",argv[i])) {
      randomize = true;
      random_number = atoi(argv[++i]);
    }
    else if (!strcmp("-check",argv[i])) {
      checkDeteriorate = true;
    }
    else if (!strcmp("-max_points",argv[i])) {
      max_points = atoi(argv[++i]);
    }
    else printError();
  }
  
  
  Solution s;
  if (!readSolution(stream, s)) {
    fprintf(stderr, "error reading input file %s\n", argv[1]);
    exit (1);
  }
  if (dimension == 0) {
    dimension = s.num_objs();
  }
  Archiver * archiver;
  if (n_archiver == 0) {
    if (trash) archiver = new IdealArchiveTrash(max_size,dimension);
    else archiver = new IdealArchive(max_size,dimension);
  }
  else if (n_archiver == 1) {
    if (trash) archiver = new DistanceArchiverTrash(max_size,dimension);
    else archiver = new DistanceArchiver(max_size,dimension);
  }
  else if (n_archiver == 2) {
    if (trash) archiver = new DistributedArchiverTrash(max_size,dimension);
    else archiver = new DistributedArchiver(max_size,dimension);
  }
  else if (n_archiver == 3) {
    if (trash) archiver = new ARArchiveTrash(dimension);
    else archiver = new ARArchive(dimension);
  }
  
  vector< vector<Solution> > conjuntos;
  
  
  if (randomize) {
    FILE * entrada = fopen("entrada.txt","w");
    int counter = 0;
    Solution p;
    while(counter++ < random_number) {
      vector<double> sol;
      for(int i = 0; i<dimension; i++)
        sol.push_back(rand0to1());
      
      p.setObjectives(sol);
      archiver->add(p);
      p.print(entrada);
       if (true) {
        conjuntos.push_back(archiver->getSolutions());
      }
    }
  }
  else {
    int inputed = 0;
    do {
      archiver->add(s);
      conjuntos.push_back(archiver->getSolutions());
    } while (readSolution (stream, s) && ++inputed < max_points);
  }

  if(checkDeteriorate) {
    for(int i = 0; i < conjuntos.size()-1; i++) {
      for(int j = i+1; j < conjuntos.size(); j++) {
        if (isBetterThan(conjuntos[i],conjuntos[j])) {
  	  cout << ">-DETERIORATE" << endl;
  	
	  cout << "Set " << j << endl;
	  for(int k = 0; k < conjuntos[j].size(); k++) {
	    conjuntos[j].at(k).printLatex();
	  }
	
	  cout << "Set " << i << endl;
	  for(int k = 0; k < conjuntos[i].size(); k++) {
	    conjuntos[i].at(k).printLatex();
	  }
	
        }
	else if (deteriorated(conjuntos[i],conjuntos[j])) {
	  cout << "In set after " << i << "th element, there is a point dominated by one point of the set after " << j << "th element." << endl;
	}
      }
    }
  }
  
  
  archiver->finish();
  if (!hide) {
    for(int i = 0; i < archiver->getSolutions().size(); i++) {
      if (paren) archiver->getSolutions().at(i).printLatex();
      else archiver->getSolutions().at(i).print();
    }
  }
  cout << "Size: " << archiver->getSolutions().size() << endl;
  t = clock() - t;
  FILE * append = fopen("../Resultados/results.txt","a");
  fprintf(append,"%s Archiver %s\n",archivers[n_archiver],trash ? "+ Trash" : "");
  fprintf(append,"Max size: %d\n",max_size);
  fprintf(append,"Input: %s\n",argv[1]);
  fprintf(append,"Solutions in output: %d\n",(int) archiver->getSolutions().size());
  fprintf(append,"Time running: %f seconds\n\n",((float)t)/CLOCKS_PER_SEC);
  fclose(append);
}
