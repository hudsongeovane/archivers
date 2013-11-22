#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iomanip>
using namespace std;

int main(int argc, char * argv[])
{
	//Usage: ./calc_result <input file> <archiver name>
	const char *archivers[] = {"Dominating Archiver", "Trash Archive", "NSGA2 Archiver", "NSGA2 Archiver+Trash", "SPEA2 Archiver", "SPEA2 Archiver+Trash", "Adaptive Grid Archiver (AGA)", "AGA+TRASH", "Hypervolume Archiver (AA_S)", "HyperVolumeArchive+Trash", "Multilevel Grid Archiver (MGA)", "MGA+TRASH", "ideal Archiver ", "ideal Archiver + Trash", "distance Archiver ", "distance Archiver + Trash", "distributed Archiver ", "distributed Archiver + Trash", "ara Archiver ", "ara Archiver + Trash"};
	vector<double> allresult;
	for(int i = 0; i < 20; i++) {
		ifstream input(argv[1]);
		string read, tmp;
		int number = 0;
		double resultado = 0., temp_d;
		while (getline(input,read)) {
			temp_d = 0.0;
			if (!read.compare(archivers[i])) {
				getline(input,tmp);
				getline(input,tmp);
				getline(input,tmp);
				input >> tmp >> tmp >> temp_d >> tmp;
				number++;
				cout << temp_d << endl;
			}
			resultado += temp_d;
		}
		cout << "Average: " << resultado/((double) number) << " seconds" << endl;
		allresult.push_back(resultado/((double) number));
	}
	for(int i = 0; i < allresult.size(); i++) {
		cout << setiosflags (ios::fixed) << setprecision(2) << allresult[i] << "s & ";
	}
}
