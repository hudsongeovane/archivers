/*************************************************************************

 archiver.cpp

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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include "Solution.h"
#include "eApproxArchive.h"
#include "eParetoArchive.h"
#include "SPEA2Archive.h"
#include "NSGA2Archive.h"
#include "MultilevelGridArchive.h"
#include "AdaptiveGridArchive.h"
#include "HVArchive.h"


#include "NSGA2ArchiveTRASH.h"
#include "MultilevelGridArchiveTrash.h"
#include "SPEA2ArchiveTrash.h"
#include "AdaptiveGridArchiveTrash.h"
#include "HVArchiveTrash.h"
//for program running time
#include <time.h>

enum archive_types { UNBOUND_ARCHIVE = 0, DOMINATING_ARCHIVE, ePARETO_ARCHIVE, eAPPROX_ARCHIVE, SPEA2_ARCHIVE, NSGA2_ARCHIVE, AGA_ARCHIVE, HV_ARCHIVE, MGA_ARCHIVE, TRASH_ARCHIVE, NSGA2_ARCHIVE_TRASH, MGA_ARCHIVE_TRASH, SPEA2_ARCHIVE_TRASH, AGA_ARCHIVE_TRASH, HV_ARCHIVE_TRASH, UNDEFINED_ARCHIVE };
static const char * const archive_names[] =
  { "Unbound Archiver", "Dominating Archiver", "ePareto Archiver", "e-approx Archiver", "SPEA2 Archiver", "NSGA2 Archiver", "Adaptive Grid Archiver (AGA)", "Hypervolume Archiver (AA_S)", "Multilevel Grid Archiver (MGA)", "Trash Archive", "NSGA2 Archiver+Trash", "MGA+TRASH", "SPEA2 Archiver+Trash", "AGA+TRASH", "HyperVolumeArchive+Trash", "UNDEFINED" };

static  long seed = 0;
static  unsigned max_size = 100;
static  unsigned dimension = 0;
static  int grid_levels = 0;
static  int seq_length = -1;
static  double epsilon = 0.0001;
static  const char * fileprefix = NULL;
static  const char * seq_filename = NULL;
static  enum archive_types archive_type = UNDEFINED_ARCHIVE;


bool 
strtovector_double (char * str, vector<double> &p)
{
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

bool
readSolution (FILE *stream, Solution &s)
{
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

static void version(void)
{
  printf("%s version " VERSION 
#ifdef ARCH
           " (optimised for " ARCH ")"
#endif
           "\n\n", "archiver");

    printf(
"Copyright (C) 2011"
"\nManuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>"
"\nJoshua Knowles <j.knowles@manchester.ac.uk>"
"\nMarco Laumanns <mlm@zurich.ibm.com>"
"\n"
"This is free software, and you are welcome to redistribute it under certain\n"
"conditions.  See the GNU General Public License for details. There is NO   \n"
"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
"\n"        
"\nIMPORTANT NOTE: Please be aware that the fact that this program is         "
"\nreleased as Free Software does not excuse you from scientific              "
"\npropriety, which obligates you to give appropriate credit! If you          "
"\nwrite a scientific paper describing research that made substantive use     "
"\nof this program, it is your obligation as a scientist to (a) mention       "
"\nthe fashion in which this software was used in the Methods section;        "
"\n(b) mention the algorithm in the References section. The appropriate       "
"\ncitation is:                                                               "
"\n                                                                           "
"\n    M. Lopez-Ibanez, J. Knowles, and M. Laumanns. On Sequential            "
"\n    Online Archiving of Objective Vectors. In EMO 2011,                    "
"\n    LNCS. Springer, 2011.                                                  "
"\n                                                                           "
"\nMoreover, as a personal note, I would appreciate it if you would email     "
"\nmanuel.lopez-ibanez@ulb.ac.be with citations of papers referencing this    "
"\nwork so I can mention them to my funding agent and tenure committee.       "
"\n"
           );
}

void usage (void)
{
    printf("\n"
           "Usage: %s [OPTIONS]\n\n", "archiver");

    printf(
"Implementation of various multi-objective archiving algorithms.\n\n"
"Options:\n"
"Other optional parameters are:\n"
"-t integer : archive type\n");
  for (int i = 0; i < UNDEFINED_ARCHIVE; ++i)  {
    printf("\t%2d\t%s\n", i, archive_names[i]);
  }
  printf(
"-f character string : file name of sequence data\n"
"-N positive integer : capacity of the archive\n"
//"-k positive integer : number of objectives\n"
"-len positive integer : length of the input sequence\n"
"-s positive long : random seed\n"
"-o character string: output filename for sequence output, otherwise, print only the final result to stdout.\n"
"-g positive integer : number of levels of the adaptive grid; #grid regions=2^(l*k)\n"
"-e positive float : epsilon value for epsilon archivers\n"
"-v                : print version and copyright information\n"
"\n");
  printf ("Default values are:\n"
          "N = %d, len %d%s, seed = time(NULL), g = floor(log_2(N/2)), e = %f\n",
          max_size,
          seq_length, (seq_length < 0) ? " (all points)" : "",
          epsilon);
}

void
print_archive(BaseArchive<Solution> *archive, const char *fileprefix, int iteration)
{
  char name[1000];
  FILE *fp;
  
  assert (strlen(fileprefix) < 950);
  sprintf(name, "%s.%d", fileprefix, iteration);
  puts(name);
  if((fp = fopen(name, "w")))
    {
      archive->print(fp);
      fclose(fp);
    }
  else
    {
      perror(name);
      exit(1);
    }
}

int main(int argc, char *argv[])
{
  clock_t t;
  t = clock();
  if (argc < 2) {
    usage();
    exit(1);
  }

  seed =  time(NULL);
  // printf("argc = %d\n", argc);
  for (int i = 1; i < argc; i += 2)    {
    if(strcmp("-h", argv[i]) == 0) {
	  usage();
	  exit(0);
    }
    else if(strcmp("-v", argv[i]) == 0) {
      version();
      exit(0);
    }
    else if(strcmp("-t", argv[i])==0)
      archive_type = archive_types(atoi(argv[i+1]));
    else if(strcmp("-N", argv[i])==0)
      max_size = atoi(argv[i+1]);
    else if(strcmp("-k", argv[i])==0)
      dimension = atoi(argv[i+1]);
    else if(strcmp("-g", argv[i])==0)
      grid_levels = atoi(argv[i+1]);
    else if(strcmp("-len", argv[i])==0)
      seq_length = atoi(argv[i+1]);
    else if(strcmp("-f", argv[i])==0)
      seq_filename = argv[i+1];
    else if(strcmp("-s", argv[i])==0)
      seed = atol(argv[i+1]);
    else if(strcmp("-e", argv[i])==0)
      epsilon = atof(argv[i+1]);
    else if(strcmp("-o", argv[i])==0)
      fileprefix = argv[i+1];
    else
      {
        fprintf(stderr, "Undefined command line parameter entered. "
                "Do \"./archiver -h\" for help with parameters.\n");
        exit(1);
      }
  }

  if (seq_filename == NULL) {
    fprintf(stderr, "error: not input file (use -f file)\n");
    exit (1);
  }

  FILE *fich = fopen(seq_filename, "r");
  if (!fich) {
    perror (seq_filename);
    exit(1);
  }

  Solution s;
  
  if (!readSolution(fich, s)) {
    fprintf(stderr, "error reading input file %s\n", seq_filename);
    exit (1);
  }

  if (dimension == 0) {
    dimension = s.num_objs();
  }

  Random rng (seed);
  BaseArchive<Solution> * archive;

  switch (archive_type) {
  case UNBOUND_ARCHIVE:
    archive = new UnboundedArchive<Solution> (s.num_objs());
    break;
  case DOMINATING_ARCHIVE:
    archive = new DominatingArchive<Solution> (max_size, s.num_objs());
    break;
  case ePARETO_ARCHIVE:
    archive = new eParetoArchive<Solution> (max_size, s.num_objs(), epsilon);
    break;
  case eAPPROX_ARCHIVE:
    archive = new eApproxArchive<Solution> (max_size, s.num_objs(), epsilon);
    break;
  case SPEA2_ARCHIVE:
    archive = new SPEA2Archive<Solution> (max_size, s.num_objs(), rng);
    break;
  case NSGA2_ARCHIVE:
    archive = new NSGA2Archive<Solution> (max_size, s.num_objs(), rng);
    break;
  case AGA_ARCHIVE:
    if (grid_levels == 0) {
      grid_levels = AdaptiveGridArchive<Solution>::default_grid_levels (max_size, s.num_objs());
    }
    archive = new AdaptiveGridArchive<Solution> (max_size, s.num_objs(), rng, grid_levels);
    break;
  case HV_ARCHIVE:
    archive = new HVArchive<Solution> (max_size, s.num_objs(), rng);
    break;
  case MGA_ARCHIVE:
    archive = new MultilevelGridArchive<Solution> (max_size, s.num_objs());
    break;
  case TRASH_ARCHIVE:
    archive = new TrashArchive<Solution> (max_size, s.num_objs());
    break;
  case NSGA2_ARCHIVE_TRASH:
    archive = new NSGA2ArchiveTrash<Solution> (max_size, s.num_objs(), rng);
    break;
  case MGA_ARCHIVE_TRASH:
    archive = new MultilevelGridArchiveTrash<Solution> (max_size, s.num_objs());
    break;
  case SPEA2_ARCHIVE_TRASH:
    archive = new SPEA2ArchiveTrash<Solution> (max_size, s.num_objs(), rng);
    break;
  case AGA_ARCHIVE_TRASH:
    if (grid_levels == 0) {
      grid_levels = AdaptiveGridArchiveTrash<Solution>::default_grid_levels (max_size, s.num_objs());
    }
    archive = new AdaptiveGridArchiveTrash<Solution> (max_size, s.num_objs(), rng, grid_levels);
    break;
  case HV_ARCHIVE_TRASH:
    archive = new HVArchiveTrash<Solution>(max_size, s.num_objs(), rng);
    break;
    
  default:
    printf ("error: undefined archive type: %d\n", archive_type);
    exit(1);
  }
  int iteration = 0;
  do {
    iteration++;
    archive->add(s);
    if (fileprefix)
      print_archive(archive, fileprefix, iteration);
  } while (readSolution (fich, s) && iteration != seq_length);

 
  if (archive_type == TRASH_ARCHIVE || archive_type == MGA_ARCHIVE_TRASH || archive_type == NSGA2_ARCHIVE_TRASH || archive_type == SPEA2_ARCHIVE_TRASH || archive_type == AGA_ARCHIVE_TRASH || archive_type == HV_ARCHIVE_TRASH) {
    archive->finish();
  }
  //if (!fileprefix)
    //archive->print();

  printf("Tamanho do arquivo: %d\n",(int) archive->size());
  t = clock() - t;
  printf("Tempo rodando: %f\n",((float)t)/CLOCKS_PER_SEC);
 
  FILE * append = fopen("Resultados/results.txt","a");
  fprintf(append,"%s\n",archive_names[archive_type]);
  fprintf(append,"Max size: %d\n",max_size);
  fprintf(append,"Input: %s\n",seq_filename);
  fprintf(append,"Solutions in output: %d\n",(int) archive->size());
  fprintf(append,"Time running: %f seconds\n\n",((float)t)/CLOCKS_PER_SEC);
  fclose(append);

  return 0;
}

