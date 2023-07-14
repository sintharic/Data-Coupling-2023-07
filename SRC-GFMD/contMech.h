#include "header.h"

// functions called in main

void initParams();
void initSystem();
void writeParams();
void initMeasure();
void propagate();
void constrain();
void getForces();
void fire(), fireRedirecT(Complex *, Complex *, Complex *, Complex *, double *, Lint);
void measure();
bool finished();
void outMeasure();
void outSystem();

// additional functions

void initParamsDefault();
bool initParamsFile();
void writeParamsDefault();

// global parameters

double lengthX, lengthY, areaXY;
Lint nxGlobal, nyGlobal;

Lint iTime, nTime, nRelax;
double dTime, dTimeInit, dTime2, dampGlobal, mdTime;

int fFire, fFireOn, iTimeLFS; // LFS = last fire stop
double fireRedrct, fireIncrmt, fireDecrmt, fireRedrctInit;
bool fFireConstraint;

int fLangevin;
double temp, tempInit,  tempFinal;
bool fExternalDriving;

int nSheet, nInter, nAtomL;

ifstream readParams;
ofstream moni;

int iFrame, frameInterval, fLogMeasure, incrmtMeasure = 1;

double timeTotal;

// parameters only needed in condMat

int randSeed;

double tKinGlobal, vPotGlobal;
double tKinOld, vPotOld;

clock_t tPropagate, tGetForces, tConstrain, tFire, tMeasure;
