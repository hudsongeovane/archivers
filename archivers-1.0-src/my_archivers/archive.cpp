#include <iostream>
#include <cstdlib>
#include "archivers.h"
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


int main(int argc, char * argv[]) {
  if (argc < 3) {
    cout << "Error. Use ./archive file max_size" << endl;
    return 0;
  }
  FILE * stream = fopen(argv[1], "r");
  int max_size = atoi(argv[2]), dimension = 0;
  Solution s;
  if (!readSolution(stream, s)) {
    fprintf(stderr, "error reading input file %s\n", argv[1]);
    exit (1);
  }
  if (dimension == 0) {
    dimension = s.num_objs();
  }
  IdealArchiveTrash archiver(max_size,dimension);
  //archiver.add(s);
  //cout << archiver.solutions.size() << endl;
  do {
    archiver.add(s);
  } while (readSolution (stream, s));
  
  archiver.finish();
  for(int i = 0; i < archiver.solutions.size(); i++) {
    archiver.solutions[i].print();
  }
  
  cout << "Ideal Point: ";
  archiver.ideal.print();
  cout << "Size: " << archiver.solutions.size() << endl;
}
