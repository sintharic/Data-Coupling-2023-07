#include "maxwell.h"

void printHelp(string);

void readMaxwell();
void initMaxwell();
void initGFMD();
void updateMaxwell();
void updateParams();

void runRelax();
void freqMod();
void stressAna(double, string);
void dispAna(double, string);
void dispGFMD(double, string);


int main(int argc, const char *argv[]) {

  // process inputs
  vector<string> args(argc);
  for (int iarg = 0; iarg < argc; ++iarg) {
    args[iarg] = argv[iarg];
    if (args[iarg] == "--help" || args[iarg] == "-h") printHelp(args[0]);
  }
  if (argc > 1) {
    maxwellFile = args[1]; 
    assert_exist(maxwellFile); 
    ID = (int) (maxwellFile[7] - '0'); 
    readMaxwell();
  }
  else if (find_maxwellFile()) readMaxwell();
  if (argc > 2) {
    paramsFile = args[2]; 
    assert_exist(paramsFile); 
    fUpdateParams = 1;
  }
  if (argc > 3) {fRunRelax = stoi(args[3]);}

  initMaxwell();
  freqMod();
  initGFMD();
  updateMaxwell();
  updateParams();
  runRelax();
}


void printHelp(string name) {
  printf("Usage: %s [maxwellFile] [paramsFile] [fRunRelax]\n\n", name.c_str());
  printf("Options and arguments:\n");
  printf("<maxwellFile> input file containing Maxwell-specific parameters.\n");
  printf("              This file is updated to be readable by contMech.exe.\n");
  printf("              If no file is specified, the current directory is searched\n");
  printf("              for a compatible one, otherwise a new one is generated.\n");
  printf("<paramsFile>  contMech params.in file, which will be updated to be \n");
  printf("              compatible with <maxwellFile>. If no file is specified,\n");
  printf("              the current directory is searched for a compatible one.\n");
  printf("<fRunRelax>   do (<fRunRelax>=1) or do not (<fRunRelax>=0) calculate\n");
  printf("              and export displacement and stress relaxation curves.\n");
  exit(0);
}

void readMaxwell() {
  //cerr << "  readMaxwell()" << endl;//FLOW

  ifstream input(maxwellFile);
  if (!input.is_open()) return;
  cerr << "Reading " << maxwellFile << endl;

  std::string ROL; //REST OF LINE
  std::size_t NIS = std::string::npos; //NIS = NOT IN STRING
  
  while (input.peek()!=EOF){

    // skip empty lines and lines starting with '#'
    if (input.peek()=='#' || input.peek()=='\n') {
      getline(input, ROL);
      //cout << "skipping: _" << ROL << "_\n";//DEBUG
      continue;
    }

    // read stiffness, stiffHigh, nMaxwell
    double param;
    input >> param;
    getline(input, ROL);
    if (ROL.find("# stiffness #") !=NIS) stiffness = param;
    if (ROL.find("# stiffHigh #") !=NIS) stiffHigh = param;
    if (ROL.find("# nMaxwell #") !=NIS) nMaxwell = param;
    if (ROL.find("# exponent #") !=NIS) exponent = param;
    if (ROL.find("# timeFast #") !=NIS) timeFast = param;
    if (ROL.find("# timeScalFac #") !=NIS) timeScalFac = param;
    if (ROL.find("# timeRelGFMD #") !=NIS) timeRelGFMD = param;

  }
  input.close();
}

void initMaxwell() {
  //cerr << "  initMaxwell()" << endl;//FLOW

  // calculate prony series
  double stiffScal  = 1. / pow(timeScalFac, exponent );
  stiffMw.resize(nMaxwell);
  tauMw.resize(nMaxwell);
  etaMw.resize(nMaxwell);

  double kHighTest = stiffMw[0] = 1;
  for (int iMw=1; iMw<nMaxwell; ++iMw) {
    double tau = timeFast*pow(timeScalFac,1.*iMw);
    stiffMw[iMw] = stiffScal*stiffMw[iMw-1];
    kHighTest += stiffMw[iMw];
  }

  for (int iMw=0; iMw<nMaxwell; ++iMw) {
    tauMw[iMw] = timeFast*pow(timeScalFac,1.*iMw);
    stiffMw[iMw] *= (stiffHigh-stiffness)/kHighTest;
    etaMw[iMw] = tauMw[iMw]*stiffMw[iMw];
  } 
}


void initGFMD() {
  //cerr << "  initGFMD()" << endl;//FLOW

  // calculate optimal GFMD params
  massGFMD = 1.*stiffHigh*pow(timeRelGFMD*timeFast/(2*M_PI),2);
  dampGlobal = 2*sqrt(stiffHigh/massGFMD);
  double minPeriod = 2*M_PI*sqrt(massGFMD/stiffHigh);
  dTime = minPeriod/TIME_STEPS_PER_PERIOD;

  stiffFacMw.resize(nMaxwell);
  invTauMw.resize(nMaxwell);
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    stiffFacMw[iMw] = stiffMw[iMw]/stiffHigh;
    invTauMw[iMw] = 1./tauMw[iMw];
  }

  /*//DEBUG
  cout << "# Check for critical damping condition:\n";
  cout << "# 4 * stiffness/m  = " << 4*stiffness/massGFMD << "\n";
  cout << "# 4 * stiffHigh/m = " << 4*stiffHigh/massGFMD << "\n";
  cout << "# damping^2   = " << dampGlobal*dampGlobal << "\n";//*/
}

void updateParams() {
  //cerr << "  updateParams()" << endl;//FLOW

  if (fUpdateParams==2) {
    cout << "Do you want to update " << paramsFile << " now? (y/n): ";
    string cmd; cin >> cmd;
    if (cmd != "y") return;
  }
  else if (fUpdateParams==0) return;

  assert_exist(paramsFile);
  cerr << "Reading " << maxwellFile << endl;
  vector<string> oldContent = readLines(paramsFile);
  std::size_t NIS = std::string::npos; // NIS == Not In String
  vector<string> newContent;
  int currentID = -1;
  for (string line : oldContent) {
    //cout << ID << "_" << currentID << "_" << line << "\n";//DEBUG

    // global params
    if (line.find("# dampGlobal #") != NIS) {
      ostringstream val; val << dampGlobal;
      newContent.push_back(val.str() + "\t\t# dampGlobal #");
      continue;
    }
    else if (line.find("# sheet start") != NIS) currentID = (int) (line[0] - '0');
    else if (line.find("# sheet end") != NIS) currentID = -1;
    else if (currentID != ID) {newContent.push_back(line); continue;}

    // sheet params
    else if ( (line.find("# nElast #")!=NIS) && (line[0] != '0') ) {
      ostringstream val; 
      newContent.push_back("1\t\t# nElast #");
      val.str(string()); val << stiffness; newContent.push_back(val.str() + "\t\t# stiffness0 #");
      val.str(string()); val << stiffHigh; newContent.push_back(val.str() + "\t\t# stiffHigh0 #");
      newContent.push_back("1\t\t# fMaxwell #");
      val.str(string()); val << nMaxwell; newContent.push_back(val.str() + "\t\t# nMaxwell #");
      newContent.push_back("1\t\t# fMassWeightg #");
      val.str(string()); val << massGFMD; newContent.push_back(val.str() + "\t\t# massGFMD #");
      newContent.push_back("");
      continue;
    }
    // avoid repeating the same params
    else if (line.find("# stiffness0 #") != NIS) continue;
    else if (line.find("# stiffHigh0 #") != NIS) continue;
    else if (line.find("# fMaxwell #") != NIS) continue;
    else if (line.find("# nMaxwell #") != NIS) continue;
    else if (line.find("# fMassWeightg #") != NIS) continue;
    else if (line.find("# massGFMD #") != NIS) continue;
    // avoid incompatible settings
    else if (line.find("# stiffness1 #") != NIS) continue;
    else if (line.find("# stiffness2 #") != NIS) continue;
    else if (line.find("# stiffness3 #") != NIS) continue;
    else if (line.find("# zeroModeMass #") != NIS) continue;
    else if (line.find("# fKelvinVoigt #") != NIS) {
      newContent.push_back("0\t\t# fKelvinVoigt #");
      continue;
    }

    newContent.push_back(line);
  }

  // update params file
  ofstream params(paramsFile);
  for (string line : newContent) {
    params << line << "\n";
  }
  params.close();
  cout << "Please adjust time-related parameters manually:\n";
  cout << dTime << "\t\t# dTime # (suggested)\n\n";
}


void updateMaxwell() {
  //cerr << "  updateMaxwell()" << endl;//FLOW

  // generate input file
  ofstream input(maxwellFile);
  input << stiffness << "\t\t# stiffness #\n";
  input << stiffHigh << "\t\t# stiffHigh #\n";
  input << nMaxwell << "\t\t# nMaxwell #\n";
  input << exponent << "\t\t# exponent #\n";
  input << timeFast << "\t\t# timeFast #\n";
  input << timeScalFac << "\t\t# timeScalFac #\n";
  input << timeRelGFMD << "\t\t# timeRelGFMD #\n\n";

  input << "# suggested GFMD parameters:\n";
  input << dTime << "\t\t# dTime #\n";
  input << massGFMD << "\t\t# massGFMD #\n";
  input << dampGlobal << "\t\t# dampGlobal #\n";

  /*// old version
  input << "\n#iMw\tkMw\tetaMw\n";
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    input << iMw << "\t" << stiffMw[iMw] << "\t" << etaMw[iMw] << "\t# MWelement #\n";
  }//*/

  input << "\n#iMw\tstiffMw\ttauMw\n";
  for (int iMw = 0; iMw < nMaxwell; ++iMw)
    input << iMw << "\t" << stiffMw[iMw] << "\t" << tauMw[iMw] << "\t# MWelement #\n";

  input.close();
}

void freqMod() {
  //cerr << "  freqMod()" << endl;//FLOW

  ofstream output("freqMod"+to_string(ID)+".dat");
  output << "#frequency\tRe(E)\tIm(E)\tabs(E)\n";
  double freqMin = (2*M_PI/vmax(tauMw))*(stiffness/stiffHigh)/RELAX_TIME_FAC;
  double freqMax = RELAX_TIME_FAC*M_PI/vmin(tauMw);
  for(double freq = freqMin; freq<freqMax; freq *= 1.1) {
    Complex modulus = stiffness;
    for (int iMw=0; iMw<nMaxwell; ++iMw) {
      Complex dampg = Complex(0,1.)*etaMw[iMw]*freq;
      modulus += stiffMw[iMw]*dampg/(stiffMw[iMw]+dampg);
    }
    output << freq << "\t" << modulus.real() << "\t" << modulus.imag()
    << "\t" << abs(modulus) << endl;
  }
  output.close();
}

void dispAna(double F0, string filename) {
  //cerr << "  dispAna()" << endl;//FLOW

  double deltaTime = vmin(tauMw)/RELAX_TIME_FAC;
  unsigned int nTime = RELAX_TIME_FAC*vmax(tauMw)*(stiffHigh/stiffness)/deltaTime;
  vector<double> dispNow(nMaxwell);
  double dispGlobal, forceGlobal;

  // initial conditions: relax from dispGlobal=forceGlobal/stiffness
  if (RELAX_MODE==1) {
    for (int iMw = 0; iMw < nMaxwell; ++iMw) dispNow[iMw] = F0/stiffness;
    dispGlobal = F0*(1./stiffness - 1./stiffHigh);
    forceGlobal = 0;
  }
  // initial conditions: relax toward dispGlobal=forceGlobal/stiffness
  else {
    for (int iMw = 0; iMw < nMaxwell; ++iMw) dispNow[iMw] = 0;
    dispGlobal = F0/stiffHigh;
    forceGlobal = F0;
  }
  
  // start time loop
  ofstream output;
  output.open(filename);
  unsigned int increment = 1;
  output << 0 << "\t" << dispGlobal << "\n";
  for (int iTime = 1; iTime < nTime; ++iTime) {
    double forceMw;

    // dump state in quasi-log steps
    if (iTime == 64*increment) increment *= 2;
    if (iTime%increment==0) output << (iTime-0.5)*deltaTime << "\t" << dispGlobal << "\n";

    // update Maxwell elements
    for (int iMw = 0; iMw < nMaxwell; ++iMw) {
      double du_dt = invTauMw[iMw]*(dispGlobal-dispNow[iMw]);
      dispNow[iMw] += du_dt*deltaTime;
    }

    // update dispGlobal
    forceMw = 0;
    for (int iMw = 0; iMw < nMaxwell; ++iMw) forceMw += stiffFacMw[iMw]*stiffHigh*dispNow[iMw];
    dispGlobal = (forceGlobal + forceMw)/stiffHigh;
  }
  output.close();
}

void stressAna(double u0, string filename) {
  //cerr << "  stressAna()" << endl;//FLOW

  double deltaTime = vmin(tauMw)/RELAX_TIME_FAC;
  unsigned int nTime = RELAX_TIME_FAC*vmax(tauMw)*(stiffHigh/stiffness)/deltaTime;
  vector<double> dispNow(nMaxwell);
  double dispGlobal, forceGlobal;

  // initial conditions: relax toward forceGlobal=dispGlobal*stiffness
  for (int iMw = 0; iMw < nMaxwell; ++iMw) dispNow[iMw] = 0;
  dispGlobal = u0;
  forceGlobal = u0*stiffHigh;
  
  // start time loop
  ofstream output;
  output.open(filename);
  unsigned int increment = 1;
  for (int iTime = 0; iTime < nTime; ++iTime) {
    double time = iTime*deltaTime;
    double forceGlobal = stiffHigh*u0; 
    for (int iMw = 0; iMw < nMaxwell; ++iMw) {
      forceGlobal -= u0*stiffFacMw[iMw]*stiffHigh*(1-exp(-time/tauMw[iMw]));
    }

    // dump state in quasi-log steps
    if (iTime == 64*increment) increment *= 2;
    if (iTime%increment==0) output << time << "\t" << forceGlobal << "\n";
  }
  output.close();
}

void dispGFMD(double F0, string filename) {
  //cerr << "  dispGFMD()" << endl;//FLOW

  unsigned int nTime = RELAX_TIME_FAC*vmax(tauMw)*(stiffHigh/stiffness)/dTime;
  vector<double> dispNow(nMaxwell);
  double dispGlobal, dispGlobOld, forceGlobal;
  double damping = dampGlobal*dTime;

  // initial conditions: relax from dispGlobal=forceGlobal/stiffness
  if (RELAX_MODE) {
    for (int iMw = 0; iMw < nMaxwell; ++iMw) dispNow[iMw] = F0/stiffness;
    dispGlobal = F0/stiffness;
    forceGlobal = 0;
  }
  // initial conditions: relax toward dispGlobal=forceGlobal/stiffness
  else {
    for (int iMw = 0; iMw < nMaxwell; ++iMw) dispNow[iMw] = 0;
    dispGlobal = dispGlobOld = 0;
    forceGlobal = F0;
  }
  
  // start time loop
  ofstream output;
  output.open(filename);
  unsigned int increment = 1;

  // time step 0 and 1
  double forceMw = 0; 
  output << 0 << "\t" << dispGlobal << "\t" << forceGlobal+forceMw << "\n";
  dispGlobal += F0*dTime*dTime/stiffHigh/2;
  //output << dTime << "\t" << dispGlobal << "\n";//DEBUG

  for (int iTime = 1; iTime < nTime; ++iTime) {
    // dump state in quasi-log steps
    if (iTime == 64*increment) increment *= 2;
    if (iTime%increment==0) output << iTime*dTime << "\t" << dispGlobal << "\t" << forceGlobal+forceMw << "\n";

    // update Maxwell elements via improved Euler's method (a.k.a. Heun's method)
    for (int iMw = 0; iMw < nMaxwell; ++iMw) {
      double slopeNow = invTauMw[iMw]*(dispGlobal - dispNow[iMw]);
      double dispFnew = dispNow[iMw] + slopeNow*dTime;
      double slopeNew = invTauMw[iMw]*(dispGlobal - dispFnew);
      dispNow[iMw] += 0.5*(slopeNow + slopeNew) * dTime;
    }

    // update dispGlobal
    forceMw = 0;
    for (int iMw = 0; iMw < nMaxwell; ++iMw) forceMw += stiffFacMw[iMw]*(1*stiffHigh)*dispNow[iMw];
    double acc = (forceGlobal + forceMw - (1*stiffHigh)*dispGlobal)/(1*massGFMD);
    double dispNew = (2.-damping)*dispGlobal - (1.-damping)*dispGlobOld + acc*dTime*dTime;

    // step forward 
    dispGlobOld = dispGlobal;
    dispGlobal = dispNew;
  }
  output.close();
}


void runRelax() {
  //cerr << "  runRelax()" << endl;//FLOW

  if (fRunRelax==2) {
    cout << "Do you want to calculate analytic relaxations? (y/n): ";
    string cmd; cin >> cmd;
    if (cmd != "y") return;
  }
  else if (fRunRelax==0) return;
  
  string fileName = "MWstressAna"+to_string(ID)+".dat";
  stressAna(1./stiffHigh,fileName); 
  cout << fileName << "\n";
  
  fileName = "MWdispAna"+to_string(ID)+".dat";
  dispAna(stiffness,fileName);
  cout << fileName << "\n";
  
  fileName = "MWdispGFMD"+to_string(ID)+".dat";
  dispGFMD(stiffness,fileName);
  cout << fileName << "\n";
}