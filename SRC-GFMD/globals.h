#ifndef GLOBALS_H
#define GLOBALS_H
//endif at EOF

#include "header.h"

// global variables from contMech

extern Lint nxGlobal, nyGlobal;
extern double lengthX, lengthY, areaXY;

extern Lint iTime, nTime, nRelax;
extern double dTime, dTime2, dampGlobal, mdTime;

extern int fFire, fFireOn;
extern double fireRedrct, fireIncrmt, fireDecrmt, fireRedrctInit;

extern int fLangevin;
extern double temp, tempInit,  tempFinal;

extern bool fExternalDriving;

extern int nSheet, nInter;

// global functions from auxiliary

extern double rrand();

extern void termination(const string&);
extern void termination(const string&, int);
extern void termination(const string&, double);

extern void dumpRealCS(vector<double*>, vector<string>, Lint, Lint, int, int, const string&);
extern void dumpReal(vector<double*>, Lint, Lint, int, int, const string&);
extern void dumpFour(vector<Complex*>, Lint, Lint, const string&);
extern int  readReal(double*, Lint, Lint, const string&, int);
extern int  readFour(Complex*, Lint, Lint, const string&, int);

extern void spectralAnalysis(Complex *, vector <double> *, Lint nx, Lint ny);
extern void realSpaceAnalysis(double *, vector<double> *, Lint, Lint);

extern double mix64(void);
extern void changeSeedMix64(int);

#endif