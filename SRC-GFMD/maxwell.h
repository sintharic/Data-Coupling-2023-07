#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <complex>
using namespace std;
typedef complex<double> Complex;

// switches
unsigned short int fRunRelax = 2;
unsigned short int fUpdateParams = 2;

// parameters for relaxation curves
const double RELAX_TIME_FAC = 100;
const int RELAX_MODE = 0; 
// 0: displ. relax. from F/stiffHigh to F/stiffness, 1: from F/stiffness to 0


// auxiliary functions
vector<string> readLines(string);
bool find_maxwellFile();
void assert_exist(string);
double vmin(vector<double>);
double vmax(vector<double>);


// default params
double stiffness = 0.5;
double stiffHigh = 500;
int nMaxwell = 1;
double exponent  = 1;
double timeFast = 1.*stiffness/stiffHigh;
double timeScalFac = 6.;
double timeRelGFMD = 0.5;


// Maxwell params
vector<double> stiffMw, etaMw, tauMw;
string maxwellFile = "maxwell.in";


// GFMD params
double massGFMD, dampGlobal, dTime;
vector<double> stiffFacMw, invTauMw;
const double TIME_STEPS_PER_PERIOD = 15;
string paramsFile = "params.in";
int ID;


double vmax(vector<double> v) {
  double result = v.at(0);
  for (double val : v) {if (val>result) result=val;}
  return result;
}

double vmin(vector<double> v) {
  double result = v.at(0);
  for (double val : v) {if (val<result) result=val;}
  return result;
}


void assert_exist(string fileName) {
  ifstream test(fileName);
  if (!test.is_open()) {
    test.close();
    cerr << "ERROR: file " << fileName << " does not exist.\n";
    exit(1);
  }
  test.close();
}


bool find_maxwellFile() {
  ifstream test;
  for (int i = 0; i < 10; ++i) {
    string testFile = "maxwell" + to_string(i) +".in";
    test.open(testFile);
    if (test.is_open()) {
      maxwellFile = testFile;
      test.close();
      break;
    }
  }
  ID = (int) (maxwellFile[7] - '0');
  //cout << ID << " ID\n";//DEBUG
  if (maxwellFile == "maxwell.in") return 0;
  else return 1;
}


vector<string> readLines(string fileName) {
  vector<string> result;
  ifstream input(fileName);
  while (!input.eof()) {
    string line;
    getline(input, line);
    result.push_back(line);
  }
  input.close();
  return result;
}