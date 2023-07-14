#include "header.h"
#include "contMech.h"
#include "gfmdSheet.h"
#include "interSheet.h"
#include "atomicSheet.h"

const int nSheetMax=3, nInterMax=3, nAtomLMax=2;
vector <gfmdSheet> sheet(nSheetMax);
vector <interSheet> inter(nInterMax);
vector <atomicSheet> atomL(nAtomLMax);

extern void termination(const string&);
extern double rrand();
extern void changeSeedMix64(int);

int main(){

  initParams();
  writeParams();
  initSystem();
  initMeasure();

  //if (iTime < 3) cerr << "---\n";//FLOW
  while (!finished()){
    iTime += 1;
    if (iTime) propagate();
    constrain();
    getForces();
    if (fFire) fire();
    measure();
  }
  //if (iTime < 3) cerr << "---\n";//FLOW

  outMeasure();
  outSystem();

}

void initParams(){

  // setting parameters
  initParamsDefault();
  if (!initParamsFile()) cerr << "\n# Using global default parameters\n";
  else cerr << "\n# Reading params.in\n";
  
  // further initializations
  iTime = -1;
  temp = tempInit;
  srand(randSeed);
  changeSeedMix64(randSeed);

  // sanity checks and parameter post processing

  // http://www.graphics.stanford.edu/~seander/bithacks.html
  if ( ((nxGlobal&(nxGlobal-1))!=0) || ((nyGlobal&(nyGlobal-1))!=0) )
  cerr << "\n### nxGlobal and/or nyGlobal is not a power of 2. Caution!!!\n\n";

  if ( (lengthX<=0) || (lengthY<=0) )
  termination("lengthX and/or lengthY not positive.\n");

  if (fFire&&fLangevin) termination("fFire + fLangevin used simultanseously.\n");

  areaXY = lengthX*lengthY;
  dTimeInit = dTime;
  dTime2 = dTime*dTime;

  vPotGlobal = tKinGlobal = 0; //CM-change: if we ever read in old velocities...

  // initialize object parameters
  for (int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].initParams(iSheet);
  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].initParams(iInter);
  for (int iAtomL = 0; iAtomL < nAtomL; ++iAtomL) atomL[iAtomL].initParams(iAtomL);

}

void initSystem(){
  for (int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].initSystem();
  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].initSystem();
  for (int iAtomL = 0; iAtomL < nAtomL; ++iAtomL) atomL[iAtomL].initSystem();
}

void writeParams(){

  writeParamsDefault();
  ofstream output("params.out");

#ifdef VERSION
  string version = to_string(VERSION);
  string msg = version.substr(0,4) + "-" + version.substr(4,2) + "-" + version.substr(6,2);
  msg += " " + version.substr(8,2) + ":" + version.substr(10,2);
  output << 0 << "\t# compiled " << msg << "\n";
#endif

  output << lengthX     << "\t\t# lengthX #" << "\n";
  if (lengthY!=lengthX) output << lengthY     << "\t\t# lengthY #" << "\n\n";

  output << nxGlobal    << "\t\t# nxGlobal #" << "\n";
  if (nyGlobal!=nxGlobal) output << nyGlobal    << "\t\t# nyGlobal #" << "\n\n";

  output << "\n"<<nRelax << "\t\t# nRelax #" << "\n";
  output << nTime << "\t\t# nTime #" << "\n";
  output << dTime << "\t\t# dTime #" << "\n\n";

  output << dampGlobal  << "\t\t# dampGlobal #" << "\n";
  output << randSeed    << "\t\t# randSeed #" << "\n\n";

  if (fLogMeasure) output << fLogMeasure << "\t\t# fLogMeasure #" << "\n";
  if (frameInterval) output << frameInterval << "\t\t# frameInterval #" << "\n";

  if (fFire) {
    output << fFire       << "\t\t# fFire #" << "\n";
    output << fireRedrct  << "\t\t# fireRedrct #" << "\n";
    output << fireIncrmt  << "\t\t# fireIncrmt #" << "\n";
    output << fireDecrmt  << "\t\t# fireDecrmt #" << "\n\n";
  }

  if (fLangevin) {
    output << fLangevin << "\t\t# fLangevin #" << "\n";
    output << tempInit  << "\t\t# tempInit #"  << "\n";
    if (tempInit!=tempFinal) output << tempFinal  << "\t\t# tempFinal #" << "\n\n";
  }

  output << nSheet      << "\t\t# nSheet #" << "\n";
  output << nInter      << "\t\t# nInter #" << "\n";
  output << nAtomL      << "\t\t# nAtomL #" << "\n\n";
  output << 0           << "\t\t# end global parameters \n\n";
  output.close();

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].writeParams();
  for (int iInter=0; iInter<nInter; ++iInter) inter[iInter].writeParams();
  for (int iAtomL=0; iAtomL<nAtomL; ++iAtomL) atomL[iAtomL].writeParams();

} // writeParams

void initMeasure(){
  //if (iTime < 3) cerr << " " << iTime << " initMeasure()\n";//FLOW

  mdTime = -dTime*nRelax;

  moni.open("gMoni.dat", ofstream::out);
  timeTotal = clock();

  tPropagate = tGetForces = tConstrain = tMeasure = 0;

  for (int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].initMeasure();
  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].initMeasure();
  for (int iAtomL = 0; iAtomL < nAtomL; ++iAtomL) atomL[iAtomL].initMeasure();

  if (frameInterval) {
    cerr << "# Animating " << nTime/frameInterval << " frames.\n";
    system("mkdir -p frames");
    system("rm -f frames/*");
  }

}

void propagate(){
  //if (iTime < 3) cerr << " " << iTime << " propagate()\n";//FLOW

  clock_t timerLocal = clock();

  if (fFire) {
    if (vPotGlobal>vPotOld) {
      dTime *= fireDecrmt;
      for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].fireHalt();
      
    } else {
      dTime *= fireIncrmt;
      if (fFireConstraint && (dTime>2*dTimeInit) ) dTime = 2*dTimeInit;
      if (dTime>2*dTimeInit) dTime = 2*dTimeInit; // MM change
    }
  }

  mdTime += dTime;

  temp = tempInit;
  if (iTime>nRelax) temp += (iTime-nRelax)*(tempFinal-tempInit)/nTime;

  tKinGlobal = 0;

  for (int iSheet=0; iSheet<nSheet; ++iSheet) 
  tKinGlobal += sheet[iSheet].propagate();

  for (int iAtomL=0; iAtomL<nAtomL; ++iAtomL) 
  tKinGlobal += atomL[iAtomL].propagate();

  tPropagate += clock() - timerLocal;

}

void constrain(){
  //if (iTime < 3) cerr << " " << iTime << " constrain()\n";//FLOW

  clock_t timerLocal = clock();

  for (int iInter = 0; iInter < nInter; ++iInter)
  if (inter[iInter].fConstraint) inter[iInter].constraint();

  tConstrain += clock() - timerLocal;
}

void getForces(){ 
  //if (iTime < 3) cerr << " " << iTime << " getForces()\n";//FLOW

  clock_t timerLocal = clock();

  if (fFire) fFireOn=1;
  vPotGlobal = 0;

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].zeroStress();

  for (int iInter=0; iInter<nInter; ++iInter) 
    vPotGlobal += inter[iInter].addStress();

  for (int iSheet=0; iSheet<nSheet; ++iSheet) 
    vPotGlobal += sheet[iSheet].addStress();

  tGetForces += clock() - timerLocal;

}

void fire(){
  //if (iTime < 3) cerr << " " << iTime << " fire()\n";//FLOW

  // recompute kinetic energy on half steps
  tKinGlobal = 0;
  for (int iSheet = 0; iSheet<nSheet; ++iSheet) 
    tKinGlobal += sheet[iSheet].tKinHalfStep();

  if (tKinOld>0.999*tKinGlobal) {
    dTime *= fireDecrmt;
    for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].fireHalt();
    tKinGlobal = 0;
    vPotGlobal = vPotOld;
  }
}

void fireRedirecT(Complex *force, Complex *rN, Complex *rC, 
      Complex *rO, double *mass, Lint nR) {
  // cannot use template for fireRedirect, as complex and real squares differ
  // rN = rNew, rC = rCopy (of new), rO = rOld

  if(iTime==1) cerr << "# fireRedirecT not well tested" << endl;

  Complex wvDotWV = 0, fDotF = 0; // wv: weighted velocity
  for(Lint k = 0; k < nR; ++ k){
    Complex velocity = rN[k] - rO[k];
    wvDotWV += mass[k] * 1. * velocity   * conj(velocity);
    fDotF += mass[k] * 1. * force[k] * conj(force[k]);
  }

  double scaleV = (1. - fireRedrct);
  double scaleF = fireRedrct * sqrt( wvDotWV.real() / fDotF.real() );

  for(Lint k = 0; k < nR; ++ k){
    rC[k] = rN[k];
    rN[k]  = rO[k] + scaleV * (rN[k] - rO[k]);
    rN[k] += scaleF * force[k];
  }
}

void measure(){
  //if (iTime < 3) cerr << " " << iTime << " measure()\n";//FLOW

  clock_t timerLocal = clock();

  if ( (frameInterval) && (iTime%frameInterval == 0) ) {
    for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].dumpFrame(iFrame);
    for (int iInter=0; iInter<nInter; ++iInter) inter[iInter].dumpFrame(iFrame);
    ++iFrame;
  }
  
  if (fLogMeasure) {
    if (iTime == 64*incrmtMeasure) incrmtMeasure *= 2;
    if (iTime%incrmtMeasure && !(iTime==nTime+nRelax-1)) {tMeasure += clock() - timerLocal; return;}
  } 

  if (iTime==0) moni << "# mdTime\ttKinGlobal\tvPotGlobal\n";
  moni << mdTime << "\t" << tKinGlobal << "\t" << vPotGlobal;
  if (fFire) moni << "\t" << dTime;
  moni << endl;

  vPotOld = vPotGlobal;
  tKinOld = tKinGlobal;

  if (iTime<nRelax) {tMeasure += clock() - timerLocal; return;}

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].measure();

  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].measure();

  tMeasure += clock() - timerLocal;
  
}

bool finished(){ 

  // planned stop
  if (iTime==(nTime+nRelax)) return(true);

  // fFire stop
  if ( fFire && (dTime<1e-10*dTimeInit) ) {
    cerr << "# Fire: dTime small at iTime = " << iTime << "\n";
    getForces();
    return(true);
  } else if ( fFire && (dTime>1e+100*dTimeInit) ) {
    cerr << "# Fire: dTime large at iTime = " << iTime << "\n";
    getForces();
    return(true);
  }

  // hard stop
  ifstream testStop1("stop");
  if (testStop1.is_open()) {
    system("rm stop");
    testStop1.close();
    termination("stop requested.");
  } testStop1.close();

  // early stop
  ifstream testStop2("finish");
  if (testStop2.is_open()) {
    system("rm finish");
    testStop2.close();
    cerr << "\n# finish requested at iTime = " << iTime << "\n\n";
    return(true);
  } testStop2.close();

  return(false);
}

void outMeasure(){
  //cerr << " " << iTime << " outMeasure()\n";//FLOW

  timeTotal = clock() - timeTotal;
  clock_t timeKnown = tPropagate + tGetForces + tConstrain + tFire + tMeasure; 

  ofstream timing("params.out",ofstream::app);
  timing << (double) timeTotal/CLOCKS_PER_SEC/60 << "\t# absolute computing time (min)\n";
  timing << (double) timeKnown/CLOCKS_PER_SEC/60 << "\t# measured computing time (min)\n\n";

  timing.precision(3);
  timing << "0 \t# relative computing times:\n" << fixed;

  timing << (double) tPropagate/timeKnown << "\t\t# tPropagate\n";
  timing << (double) tGetForces/timeKnown << "\t\t# tGetForces\n";

  if (tConstrain>0.001) timing << (double) tConstrain/timeKnown << "\t\t# tConstrain\n";
  if (tFire>0.001) timing << (double) tFire/timeKnown << "\t\t# tFire\n";
  if (tMeasure>0.001) timing << (double) tMeasure/timeKnown << "\t\t# tMeasure\n";

  timing << 0 << "\t# relative computing times:\n\n";
  timing.close();

  // up-date of real-space stress field required
  getForces();  //MM2XY: Why is the elastic stress (sometimes) destroyed in konfig-D.dat?

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].outMeasure();
  for (int iInter=0; iInter<nInter; ++iInter) inter[iInter].outMeasure();

} // outMeasure

void outSystem(){
  //cerr << " " << iTime << " outSystem()\n";//FLOW

  for (int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].outSystem(); 
  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].outSystem(); 

  ofstream output("params.out",ofstream::app);
  output << 1 << "\t\t# D O N E";
  output.close(); 

  moni.close();
}

void initParamsDefault(){

  lengthX = 1.; lengthY = 1.;
  nxGlobal = 128; nyGlobal = 128;

  nTime = 100; dTime = dTimeInit = 0.25; dampGlobal = 1.5;
  nRelax = 0;
  randSeed = 4711;

  fFire = 0;
  fireRedrct = 0.0;
  fireIncrmt = 1.2;
  fireDecrmt = 0.5;
  fFireConstraint = false;

  fLangevin = 0;
  tempInit  = 0.01;
  tempFinal = tempInit;

  frameInterval = 0;
  iFrame = 0;

  fExternalDriving = false;

  nSheet = 2;
  nInter = 1;
  nAtomL = 0;

}

bool initParamsFile(){
  ifstream readParams;
  readParams.open("params.in");
  if ( !readParams.is_open() ) return 0;

  double param;
  std::string ROL; // rest of line
  std::size_t NIS = std::string::npos; // NIS == Not In String
  int iCount = 0; 
  while ( !readParams.eof() ) {
    readParams >> param; getline(readParams,ROL);
    if (++iCount>128) break;
    if      (ROL.find("# lengthX #")  !=NIS){lengthX = param;  lengthY = lengthX;}
    else if (ROL.find("# lengthY #")  !=NIS) lengthY = param;  
    else if (ROL.find("# nxGlobal #") !=NIS){nxGlobal = param; nyGlobal = nxGlobal;}
    else if (ROL.find("# nyGlobal #") !=NIS) nyGlobal = param; 
    else if (ROL.find("# nRelax #")   !=NIS) nRelax = param;
    else if (ROL.find("# nTime #")    !=NIS) nTime = param;
    else if (ROL.find("# dTime #")    !=NIS) dTime = param;
    else if (ROL.find("# dampGlobal #") !=NIS) dampGlobal = param;
    else if (ROL.find("# randSeed #") !=NIS) randSeed = param;
    else if (ROL.find("# fFire #")    !=NIS) fFire = param;
    else if (ROL.find("# fireRedrct #") !=NIS) fireRedrct = param;
    else if (ROL.find("# fireIncrmt #") !=NIS) fireIncrmt = param;
    else if (ROL.find("# fireDecrmt #") !=NIS) fireDecrmt = param;
    else if (ROL.find("# fLangevin #")!=NIS) fLangevin = param;
    else if (ROL.find("# tempInit #") !=NIS) {tempInit = param; tempFinal = tempInit;}
    else if (ROL.find("# tempFinal #")!=NIS) tempFinal = param;
    else if (ROL.find("# fLogMeasure #")!=NIS) fLogMeasure = param;
    else if (ROL.find("# frameInterval #")!=NIS) frameInterval = param;
    else if (ROL.find("# freqFrame #")!=NIS) {
      frameInterval = param;
      cerr << "# WARNING: freqFrame is deprecated. Use frameInterval instead.\n";
    }
    else if (ROL.find("# nSheet #")   !=NIS) nSheet = param;
    else if (ROL.find("# nInter #")   !=NIS) nInter = param;
    else if (ROL.find("# nAtomL #")   !=NIS) nAtomL = param;
    else if (ROL.find("# end global") !=NIS) break;
  }
  readParams.close();
  return(1);
}

void writeParamsDefault(){

  ofstream output("params.def");
  output << "! ! ! ! ! ! ! ! ! global (default) values\n\n";
  
  output << 1 << "\t\t# dampGlobal #" << "\n";
  output << 4711 << "\t\t# randSeed #" << "\n\n";

  output << 1 << "\t\t# fFire #" << "\n";
  output << 1e-3  << "\t\t# fireRedrct #" << "\n";
  output << 1.2  << "\t\t# fireIncrmt #" << "\n";
  output << 0.5  << "\t\t# fireDecrmt #" << "\n\n";

  output << 1 << "\t\t# fLangevin #" << "\n";
  output << 0.01  << "\t\t# tempInit #" << "\n";
  output << 0.01 << "\t\t# tempFinal #" << "\n\n";

  output << 1 << "\t\t# fLogMeasure #" << "\n";
  output << 100 << "\t\t# frameInterval #" << "\n\n";

  output.close();

}
