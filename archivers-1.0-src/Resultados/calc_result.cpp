#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>

using namespace std;

int main(int argc, char * argv[])
{
	//Usage: ./calc_result <input file> <archiver name>
	ifstream input(argv[1]);
	string read, tmp;
	int number = 0;
	double resultado = 0., temp_d;
	while (getline(input,read)) {
		temp_d = 0.0;
		if (!read.compare(argv[2])) {
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
}
