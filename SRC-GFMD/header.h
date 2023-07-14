#ifndef HEADER_H
#define HEADER_H
//endif at EOF

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <cmath>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <bitset>

#include <vector>
#include <stdlib.h>
#include <cstring>
#include <limits>
const double INF = std::numeric_limits<double>::infinity();

#include<ctime>

using namespace std;

typedef ptrdiff_t Lint;
typedef std::complex<double> Complex;

#ifndef _PI_
  #define _PI_
  const double PI = 4.*atan(1.), TWOPI = 2*PI;
#endif

template <class T> int sign(T a) {return (a>0) - (a<0);}

class fieldReal {
public:
  unsigned int nx, ny;
  Lint size;
  double *data;
  fieldReal(unsigned int x, unsigned int y) : nx(x), ny(y) {
    size = nx*ny;
    Lint nReal = nx*(ny+2);
    data = (double *) fftw_malloc( nReal * sizeof(double) );
  }
};

class fieldComplex {
public:
  unsigned int nx, ny;
  Lint size;
  Complex *data;
  fieldComplex(unsigned int nx, unsigned int ny) : nx(nx), ny(ny) {
    size = nx*(ny/2+1);
    data = (Complex *) fftw_malloc( size * sizeof(double) );
  }
};

class fieldManager {
private:
  vector <fieldReal> fieldsReal;
  vector <fieldComplex> fieldsComplex;
public:
  fieldReal getReal(unsigned int nx, unsigned int ny) {
    for (fieldReal f : fieldsReal) {if ((f.nx==nx) && (f.ny==ny)) return f;}
    fieldsReal.push_back(fieldReal(nx,ny));
    return fieldsReal.back();
  }
  fieldComplex getComplex(unsigned int nx, unsigned int ny) {
    for (fieldComplex f : fieldsComplex) {if ((f.nx==nx) && (f.ny==ny)) return f;}
    fieldsComplex.push_back(fieldComplex(nx,ny));
    return fieldsComplex.back();
  }
  void add(fieldReal field) {fieldsReal.push_back(field);};
  void add(fieldComplex field) {fieldsComplex.push_back(field);};
};

#endif