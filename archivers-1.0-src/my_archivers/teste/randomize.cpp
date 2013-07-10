#include <stdio.h>
#include <cstdlib>
#include <time.h>

using namespace std;

int main() {
  srand(time(NULL));
  FILE * p;
  p = fopen("testRandom.txt","w");
  for (int i = 0; i < 100000; i++) {
    int x1 = rand() % 10000 + 1;
    int x2 = rand() % 10000 + 1;
    int x3 = rand() % 10000 + 1;
    int x4 = rand() % 10000 + 1;
    int x5 = rand() % 10000 + 1;
    fprintf(p,"%d %d %d %d %d\n",x1,x2,x3,x4,x5);
  }
  fclose(p);
}
    