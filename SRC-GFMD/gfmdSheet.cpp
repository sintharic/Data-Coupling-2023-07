// MM2CM: http://gnuplot.sourceforge.net/docs_4.2/node103.html
// ... zu bedenken, wenn wieder (2D statt 1D) Movies gemacht werden

// MM2All: Think about nVeloStopStep, nVeloStartStep

#include "header.h"
#include "gfmdSheet.h"
#include "globals.h"

// functions called from main

void gfmdSheet::initParams(int newID) {

  if (newID>1) cerr << "# ID might be too large in gfmdSheet.";
  ID = newID;

  // setting parameters
  initParamsDefault();
  initParamsFile();
  //NEW-MW
  if (fMaxwell) initMaxwell(); 
  else {for (int iLayer=0; iLayer<nLayer; ++iLayer) stiffHigh[iLayer] = stiffness[iLayer];}

  // sanity checks and post-processing of parameters

  if (ny==1 && nx>1) termination("Do not use ny=1, but rather nx=1!");

  nxH = nx/2; nyHP1 = (ny/2) + 1;
  nFour = nx*nyHP1;
  nReal = 2*nFour;
  nxny = nx*ny;
  areaPP = areaXY/nxny;

  dx = lengthX / nx; dqx = TWOPI / lengthX; dqx2 = dqx*dqx;
  dy = lengthY / ny; dqy = TWOPI / lengthY; dqy2 = dqy*dqy;

  // do post processing of roughness parameters first

  fTopography = 0;

  // add 4 to fTopoAdd, if this bit is not set in fTopoAdd
  if ( fAddSWR && !(fTopoAdd&4) ) fTopoAdd += 4;

  if ( fTopoAdd || fTopoRead ) fTopography = 1;


  fSheetType = 0;
  if (fTopography) fSheetType += 1;
  if (nElast) {
    fSheetType += 2;
    if (nElast>=2) fDispY = 1;
    if (nElast==3) fDispX = 1;
  }

  if ( rRoughNorm<=0 ) {
    rRoughNorm = 1;
    cerr << "# Forcing rRoughNorm=1" << endl;
  }

  if (peklenik<=0) {
    cerr << "# Forcing peklenik=1" << endl;
    peklenik = 1;
  }

  if ( fSteppedRamp && (vzConstCOM!=0) ) {
    cerr << "# Setting vzConstCOM=0 due to fSteppedRamp" << endl;
    vzConstCOM = 0;
  }

  //NEW-MW
  if (fDumpMaxwell && !fMaxwell) {
    cerr << "# Setting fDumpMaxwell to 0\n";
    fDumpMaxwell = 0;
  }
  if (fMaxwell && !fMassWeightg) termination("fMaxwell only works with fMassWeightg.");
  if (fMaxwell && fKelvinVoigt) termination("fMaxwell and fKelvinVoigt are exclusive.");

  if (fFire&&fKelvinVoigt) termination("fFire AND fKelvinVoigt not allowed!");

  if (fKelvinVoigt) {
    if (frictRelax != 1) cerr << "# CAUTION: sheet" << ID << " fKelvinVoigt and frictRelax." << endl;
    if (fMassWeightg) termination("fMassWeightg AND fKelvinVoigt not allowed!");
    if (rKV_LPF < 0) rKV_LPF = 0;
  }

  //if (fLateral&&fDispX) termination("fLateral AND fDispX not allowed!");//TEMP: is this now allowed?
  //if (fLateral&&fDispY) termination("fLateral AND fDispY not allowed!");//TEMP: is this now allowed?
  //if (fDispX && !fDispY) termination("fDispX without fDispY not allowed!");
  if (stiffness[1]!=0 && nElast>1) termination("Lateral displacements not implemented for multilayers!");
  if (fThickness[0]==1 && nElast>1) termination("Lateral displacements not implemented for fThickness=1!");

  konfigName = "konfig" + to_string(ID);
  if (fTopography) konfigName = konfigName + "E"; 
  if (nElast) konfigName = konfigName + "D";
  if (nElast>=2) konfigName = konfigName + "zy"; 
  if (nElast==3) konfigName = konfigName + "x";

  konfigName = konfigName + ".";

  // finish function with processing of elasticity-related parametrizations
  if (!nElast) return;

  if ( (fSteppedRamp==1) && (!fConstCOM) ) fConstCOM = 1; 
  if ( (fSteppedRamp==2) && fConstCOM) termination("fSteppedRamp=2 incompatible with fConstCOM!");

  if (zOpposite) zOpposite *= nxny;
  if (fExtrapolateInf && fzOpposite) 
    cerr << "CAUTION: fExtrapolateInf with fzOpposite is prone to instabilities!";

  if ( (fConstCOM&&fFire) && (vzConstCOM!=0) ) {
    cerr << "# Switching off FIRE becaucse of vzConstCOM!=0\n";
    fFire = 0;
  }

  if ( (fSteppedRamp)&&(fFire) ) {
    if ( (fireIncrmt!=1)||(fireDecrmt!=1)||(fireRedrct!=0) )
      cerr << "# Resetting FIRE variables because of ramp, which may be unnecessary.\n";
    fireIncrmt = 1;
    fireDecrmt = 1;
    fireRedrct = 0;
  }

  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    if ( (elastExpnt[iLayer]!=1) && (fThickness[iLayer]!=0) )
    termination("Exponent of layer violates thickness.");
  }

  if (fSteppedRamp==1) {
    ddzRamp = nxny * dzRamp / rampSteps;
    pressInit = pressFinal = 0;
  } 
  else if (fSteppedRamp==2) {
    ddpRamp = dpRamp / rampSteps;
  }

  if (forceInit != 0) {
    cerr << "# using forceInit rather than pressInit.\n";
    pressInit = forceInit/(lengthX*lengthY);
  }
  if (forceFinal != 0) {pressFinal = forceFinal/lengthX/lengthY;}
  pressure = pressInit;
  if (checkTurnPress && (forceTurnaround != 1e+85)) {
    cerr << "# using forceTurnaround rather than pressTurnaround.\n";
    pressTurnaround = forceTurnaround/(lengthX*lengthY);
  }
  if (rVeloTurnStep>0) nVeloTurnStep = rVeloTurnStep*nTime;
  if ((nVeloTurnStep>0) && (nVeloTransition>0)) {
    tStartTransition = nVeloTurnStep - nVeloTransition/2;
    tEndTransition = nVeloTurnStep + 3*nVeloTransition/2;
    t0Transition = nVeloTurnStep + nVeloTransition/2;
  }
  else if (!checkTurnPress) {
    nVeloTransition = -1;
    tStartTransition = -1;
    tEndTransition = -1;
  }
  sgnTurnPress = (pressTurnaround < 0) - (pressTurnaround > 0);// it's a pressure, therefore sign change
  absTurnPress = abs(pressTurnaround);
  vzInit = vzConstCOM;


  if (fOnSitePotential!=0) fRealSpaceStress = 1;
  //if (fOnSitePotential==1) fConstCOM = 1;
  if ((fOnSitePeriod!=0)&&(fOnSitePotential==1)) fOnSiteFreq = 2*PI/fOnSitePeriod;

  if (scalKV<1+1e-5) scalKV = 1+1e-5;

  if (tauKV<=1e-10) cerr << "# Check units of tauKV!\n";

} // initParams

// function called from within object

void gfmdSheet::initParamsDefault() {

  nx = nxGlobal;
  ny = nyGlobal;

  fTopoAdd = fTopoRead = nElast = 0;

  fDispX = fDispY = 0;

  relXXpre = relYYpre = 1.5;
  relXYpre = 0.5;
  relXZpre = relYZpre = Complex(0.,1./8);

  if (ID==0) fTopoAdd = 1;
  if (ID==1) nElast = 1; 
  if (nSheet==1) {fTopoAdd = 0;}
  else if (nSheet!=2) cerr << "# Be aware of default setting for nSheet = " << nSheet << endl;

  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    stiffness[iLayer] = stiffHigh[iLayer] = 0;//NEW-MW
    elastExpnt[iLayer] = 1;
    poisson[iLayer] = 0.25;
    fThickness[iLayer] = 0;
    thickness[iLayer] = lengthY;
  } 
  stiffness[0] = 0.5; // makes contact modulus of regular sheet = 1

  fMassWeightg = 0;
  zeroModeMass = 1;

  // default values for elastomer
  pressure = pressInit = pressFinal = 0.01;
  forceInit = forceFinal = 0;

  fConstCOM = 0;	// 1 (2) : constr velocity of top (bottom) surface
  zConstCOM = 0; 
  vzConstCOM = 0;
  fSteppedRamp = 0;
  dzRamp = 1.e-2;
  rampSteps = 130;
  rampRelax = 170;
  rVeloTurnStep = -1;
  nVeloTurnStep = -1;
  pressTurnaround = 1e+85;
  forceTurnaround = 1e+85;
  checkTurnPress = 0;
  nVeloTransition = 0;

  fzOpposite = fExtrapolateInf = 0;
  zOpposite = zInfOffset = zInfErrorOld = 0;  

  fConstCOMx = 0; xConstCOM = 0; vxConstCOM = 0;
  fConstCOMy = 0; yConstCOM = 0; vyConstCOM = 0;

  fLateral = fFriction = vX = vY = 0;
  frictRelax = 1;
  fSteadySlide = false;

  fOnSitePotential = 0;
  fOnSitePeriod = 0;
  pressPhi = 90;
  pressTheta = 0;

  vXOnSite = 0;
  vYOnSite = 0;
  frictionCoeffOS = 0.3;
  surfEnergOS = 1.0e-05;//relatively strong
  potCurveRelOS = 0.4;

  fKelvinVoigt = 0;
  tauKV = 1.;
  scalKV = 1000.;
  rKV_LPF = 0;

  //NEW-MW
  fMaxwell = 0;
  nMaxwell = 2;
  fDumpMaxwell = 0;
  massGFMD = 1;

  // default roughness values

  // Hertz roughness 	fTopoAdd&1
  rXhertz = rYhertz = 1;
  hertzExpnt = 2;

  // random roughness	fTopoAdd&2
  hurst = 0.8;
  lambdaR = 0.5;
  lambdaS = 0.02;
  fRollOff = 1;
  fRoughNorm = 1;	// 1: normalize to rms gradient; 2: ... to rms height
  rRoughNorm = 1.;
  peklenik = 1.;
  fBoxMuller = 0;

  // SWR: single-wavelength roughness	fTopoAdd&4
  fAddSWR = 0;
  nqxAddSWR = 2; nqyAddSWR = nqxAddSWR;
  heightSWR = 0.1;

  // flat punch:	fTopoAdd&8
  radiusFlatPunch = 0.2;
  heightFlatPunch = 1;

  // hemisphere:	fTopoAdd&16
  rSphere = 0.2;

  // discrete steps
  dzStep = -1./200;

  fRealSpaceStress = 0;

  // movie parameters
  f3dMovie = 0;
  resolMovie = 512; // max number of points per line
}

bool gfmdSheet::initParamsFile() {
  ifstream readParams;
  readParams.open("params.in");
  if ( !readParams.is_open() ) return 0;

  // overwrite following default setting
  nElast = fTopoAdd = fAddSWR = 0;

  int fReadParams = 0;
  while ( !readParams.eof() ) {

    double param;
    std::string ROL; // rest of line
    std::size_t NIS = std::string::npos; // NIS == Not In String
    if (readParams.eof()) break;
    readParams >> param; getline(readParams,ROL);

    if ( (param==ID) && (ROL.find("# sheet start") !=NIS) ) fReadParams = 1;
    if (fReadParams==0) continue;
    if (ROL.find("# sheet end")   !=NIS) break;

    // topography
    if (ROL.find("# fTopoAdd #") !=NIS) fTopoAdd = param;
    if (ROL.find("# fRoughAdd #") !=NIS) fTopoAdd = param;
    if (ROL.find("# fTopoRead #") !=NIS) fTopoRead = param;
    if (ROL.find("# fRoughRead #") !=NIS) fTopoRead = param;
    if (ROL.find("# rXhertz #") !=NIS) {rXhertz = param; rYhertz = param;}
    if (ROL.find("# rYhertz #") !=NIS) rYhertz = param;
    if (ROL.find("# hertzExpnt #") !=NIS) hertzExpnt = param;
    if (ROL.find("# rSphere #") !=NIS) rSphere = param;
    if (ROL.find("# hurst #") !=NIS) hurst = param;
    if (ROL.find("# lambdaR #") !=NIS) lambdaR = param;
    if (ROL.find("# lambdaS #") !=NIS) lambdaS = param;
    if (ROL.find("# rRoughNorm #") !=NIS) rRoughNorm = param;
    if (ROL.find("# peklenik #") !=NIS) peklenik = param;
    if (ROL.find("# fRollOff #") !=NIS) fRollOff = param;
    if (ROL.find("# fRoughNorm #") !=NIS) fRoughNorm = param;
    if (ROL.find("# fBoxMuller #") !=NIS) fBoxMuller = param;
    if (ROL.find("# fAddSWR #") !=NIS) fAddSWR = param;
    if (ROL.find("# nqxAddSWR #") !=NIS) {nqxAddSWR = param; nqyAddSWR = nqxAddSWR;}
    if (ROL.find("# nqyAddSWR #") !=NIS) nqyAddSWR = param;
    if (ROL.find("# heightSWR #") !=NIS) heightSWR = param;
    if (ROL.find("# radiusFlatPunch #") !=NIS) radiusFlatPunch = param;
    if (ROL.find("# heightFlatPunch #") !=NIS) heightFlatPunch = param;
    if (ROL.find("# dzStep #") !=NIS) dzStep = param;

    // (visco-)elastic properties
    if (ROL.find("# nElast #") !=NIS) nElast = param;
    if (ROL.find("# fDispX #") !=NIS) {
      if (param) {
        nElast = 3;
        cerr << "# WARNING: fDispX is deprecated. Use nElast=3 instead.\n";
      }
    }
    if (ROL.find("# fDispY #") !=NIS) {
      if (nElast<2 && param>0) {
        nElast = 2;
        cerr << "# WARNING: fDispY is deprecated. Use nElast>=2 instead.\n";
      }
    }
    for (int iLayer=0; iLayer<nLayer; ++iLayer) {
      std::string text = "# stiffness" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) stiffness[iLayer] = param; 
      text = "# stiffHigh" + to_string(iLayer) + " #";//NEW-MW
      if (ROL.find(text) != NIS) stiffHigh[iLayer] = param;//NEW-MW
      text = "# elastExpnt" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) elastExpnt[iLayer] = param; 
      text = "# fThickness" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) fThickness[iLayer] = param; 
      text = "# poisson" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) poisson[iLayer] = param; 
      text = "# thickness" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) thickness[iLayer] = param; 
    }
    if (ROL.find("# fKelvinVoigt #") !=NIS) fKelvinVoigt = param;
    if (ROL.find("# tauKV #") !=NIS) tauKV = param;
    if (ROL.find("# scalKV #") !=NIS) scalKV = param;
    if (ROL.find("# rKV_LPF #") !=NIS) rKV_LPF = param;
    if (ROL.find("# fMassWeightg #") !=NIS) fMassWeightg = param;
    if (ROL.find("# zeroModeMass #") !=NIS) zeroModeMass = param;
    //NEW-MW
    if (ROL.find("# fMaxwell #") !=NIS) fMaxwell = param;
    if (ROL.find("# nMaxwell #") !=NIS) {
      nMaxwell = param;
      stiffFacMw.resize(nMaxwell,0); invTauMw.resize(nMaxwell,0);
    }
    if (ROL.find("# massGFMD #") !=NIS) massGFMD = param;
    if (ROL.find("# fDumpMaxwell #") !=NIS) fDumpMaxwell = param;

    // pressure/displacement/velocity control
    if (ROL.find("# pressInit #") !=NIS) {pressInit  = param; pressFinal=param;}
    if (ROL.find("# pressFinal #") !=NIS) pressFinal = param;
    if (ROL.find("# forceInit #") !=NIS) {forceInit  = param; forceFinal=param;}
    if (ROL.find("# forceFinal #") !=NIS) forceFinal = param;
    if (ROL.find("# fConstCOM #") !=NIS) fConstCOM = param;
    if (ROL.find("# zConstCOM #") !=NIS) zConstCOM = param;
    if (ROL.find("# vzConstCOM #") !=NIS) vzConstCOM = param;
    if (ROL.find("# fConstCOMx #") !=NIS) fConstCOMx = param;
    if (ROL.find("# xConstCOM #") !=NIS) xConstCOM = param;
    if (ROL.find("# vxConstCOM #") !=NIS) vxConstCOM = param;
    if (ROL.find("# fConstCOMy #") !=NIS) fConstCOMy = param;
    if (ROL.find("# yConstCOM #") !=NIS) yConstCOM = param;
    if (ROL.find("# vyConstCOM #") !=NIS) vyConstCOM = param;
    if (ROL.find("# fSteppedRamp #") !=NIS) fSteppedRamp = param;
    if (ROL.find("# rampSteps #") !=NIS) rampSteps = param;
    if (ROL.find("# rampRelax #") !=NIS) rampRelax = param;
    if (ROL.find("# dzRamp #") !=NIS) dzRamp = param;
    if (ROL.find("# dpRamp #") !=NIS) dpRamp = param;
    if (ROL.find("# rVeloTurnStep #") !=NIS) rVeloTurnStep = param;
    if (ROL.find("# nVeloTurnStep #") !=NIS) nVeloTurnStep = param;
    if (ROL.find("# nVeloTransition #") !=NIS) nVeloTransition = param;
    if (ROL.find("# pressTurnaround #") !=NIS) {
      pressTurnaround = param;
      checkTurnPress = 1;
    }
    if (ROL.find("# veloTurnPress #") !=NIS) {
      cerr << "# WARNING: veloTurnPress is deprecated. Use pressTurnaround instead.\n";
      pressTurnaround = param;
      checkTurnPress = 1;
    }
    if (ROL.find("# forceTurnaround #") !=NIS) {
      forceTurnaround = param;
      checkTurnPress = 1;
    }
    if (ROL.find("# fzOpposite #") !=NIS) fzOpposite = param;
    if (ROL.find("# fExtrapolateInf #") !=NIS) fExtrapolateInf = param;
    if (ROL.find("# zOpposite #") !=NIS)  zOpposite = param;
    if (ROL.find("# fLateral #") !=NIS) fLateral = param; 
    if (ROL.find("# vX #") !=NIS) vX = param;
    if (ROL.find("# vY #") !=NIS) vY = param;
    if (ROL.find("# fSteadySlide #") !=NIS) fSteadySlide = param; 

    // OnSite methods
    if (ROL.find("# fOnSitePotential #") !=NIS) fOnSitePotential = param;
    if (ROL.find("# fOnSitePeriod #") !=NIS) fOnSitePeriod = param;
    if (ROL.find("# frictRelax #") !=NIS) frictRelax = param;
    if (ROL.find("# frictionCoeffOS #") !=NIS) frictionCoeffOS = param;
    if (ROL.find("# surfEnergOS #") !=NIS) surfEnergOS = param;
    if (ROL.find("# potCurveRelOS #") !=NIS) potCurveRelOS = param;
    if (ROL.find("# vXOnSite #") !=NIS) vXOnSite = param;
    if (ROL.find("# vYOnSite #") !=NIS) vYOnSite = param;
    if (ROL.find("# xOnSite #") !=NIS) xOnSite = param;
    if (ROL.find("# yOnSite #") !=NIS) yOnSite = param;
    if (ROL.find("# pressPhi #") !=NIS) pressPhi = param;
    if (ROL.find("# pressTheta #") !=NIS) pressTheta = param;

    // output options
    if (ROL.find("# f3dMovie #") !=NIS) f3dMovie = param;
    if (ROL.find("# resolMovie #") !=NIS) resolMovie = param;
  } 

  readParams.close();
  return 1;
}

//NEW-MW
void gfmdSheet::initMaxwell(){
  string maxwellFile = "maxwell"+to_string(ID)+".in";
  ifstream input(maxwellFile);
  int nMaxwell_test;
  double dTime_test=0, kLow_test=0, kHigh_test=0, damping_test=0, mass_test=0;

  if(!input.is_open()) termination(maxwellFile+" does not exist.");

  cerr << "# Reading "+maxwellFile+".\n";

  size_t pos; // current position in file stream
  std::string ROL; //REST OF LINE
  std::size_t NIS = std::string::npos; //NIS = NOT IN STRING
  
  while (input.peek()!=EOF){

    // skip empty lines and lines starting with '#'
    char firstChar;
    pos = input.tellg();
    input >> firstChar;
    input.seekg(pos,input.beg);
    if (firstChar=='#' || firstChar=='\n') {
      getline(input, ROL);
      //cout << "skipping: _" << ROL << "_\n";//DEBUG
      continue;
    }

    // read stiffness, stiffHigh, nMaxwell
    pos = input.tellg();
    double param;
    input >> param;
    getline(input, ROL);
    if (ROL.find("# stiffness #") !=NIS) {
      kLow_test = param;
      if (kLow_test != stiffness[0]) termination("stiffness0 incompatible with "+maxwellFile);
    }
    if (ROL.find("# stiffHigh #") !=NIS) {
      kHigh_test = param;
      if (kHigh_test != stiffHigh[0]) termination("stiffHigh0 stiffness incompatible with "+maxwellFile);
    }
    if (ROL.find("# nMaxwell #") !=NIS) {
      nMaxwell_test = param;
      if (nMaxwell_test != nMaxwell) termination("nMaxwell incompatible with "+maxwellFile);
    }
    if (ROL.find("# dTime #") !=NIS) dTime_test = param;
    if (ROL.find("# massGFMD #") !=NIS) {
      mass_test = param;
      if (mass_test != massGFMD) termination("massGFMD incompatible with "+maxwellFile);
    }
    if (ROL.find("# dampGlobal #") !=NIS) damping_test = param;
    if (ROL.find("# MWelement #") !=NIS) {
      double iMw=-1, stiffMw=-1, tauMw=-1;

      input.seekg(pos,input.beg);
      input >> iMw >> stiffMw >> tauMw;
      getline(input,ROL);
      invTauMw[iMw] = 1./tauMw;
      stiffFacMw[iMw] = stiffMw/stiffHigh[0];
      //cout << iMw << " " << stiffMw << " " << tauMw << endl;//DEBUG
      //cout << iMw << " " << stiffFacMw[iMw] << " " << invTauMw[iMw] << endl;//DEBUG
    }
  }
  /*//DEBUG
  cout << "dTime = \t" << dTime_test << "\t" << dTime << "\n";
  cout << "mass = \t" << mass_test << "\t" << massGFMD << "\n";
  cout << "damping = \t" << damping_test << "\t" << dampGlobal << "\n";
  cout << "stiffness = \t" << kLow_test << "\t" << stiffness[0] << "\n";
  cout << "stiffHigh = \t" << kHigh_test << "\t" << stiffHigh[0] << "\n";//*/

  /*//DEBUG
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    cout << "#" << iMw << " " << stiffFacMw[iMw]*stiffHigh[0] << " " << 1./(invTauMw[iMw]*stiffHigh[0]) << endl;
  }//*/

  // sanity checks
  if (dTime_test < 0.99*dTime) {
    cerr << "CAUTION! dTime likely too large. Consider using this one or smaller:\n";
    cerr << dTime_test << "\t\t# dTime #\n";
  } 
  if (damping_test != dampGlobal) cerr << "CAUTION! You are using a different dampGlobal than suggested by maxwell.cpp!\n";

  //cerr << "Finished reading maxwell.in.\n";//DEBUG
}

void gfmdSheet::writeParams() {
  if (!ID) writeParamsDefault();

  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet start\n";

  if (nElast) {
    output << nElast << "\t\t# nElast #\n";
    for (int iLayer=0; iLayer<nLayer; ++iLayer) {
      if (stiffness[iLayer] == 0)  continue;
      output << stiffness[iLayer] << "\t\t# stiffness" << iLayer << " #\n";
      if (stiffness[iLayer] != stiffHigh[iLayer])
        output << stiffHigh[iLayer] << "\t\t# stiffHigh" << iLayer << " #\n";//NEW-MW
      if (fThickness[iLayer]) {
        output << fThickness[iLayer] << "\t\t# fThickness" << iLayer << " #\n";
        output << poisson[iLayer]    << "\t\t# poisson"    << iLayer << " #\n";
        output << thickness[iLayer]  << "\t\t# thickness"  << iLayer << " #\n";
      }
      if (elastExpnt[iLayer] != 1)
        output << elastExpnt[iLayer] << "\t\t# elastExpnt" << iLayer << " #\n";
    }
    if (nElast>=2) {
      if (!fThickness[0]) output << poisson[0] << "\t\t# poisson0 #\n";
      //output << fDispY << "\t\t# fDispY #\n";
    }
    if (fMaxwell) {//NEW-MW
      output << fMaxwell << "\t\t# fMaxwell #\n";
      output << nMaxwell << "\t\t# nMaxwell #\n";
      output << fDumpMaxwell << "\t\t# fDumpMaxwell #\n";
    }
    if (massGFMD != 1) output << massGFMD << "\t\t# massGFMD #\n";
    if (fKelvinVoigt) {
      output << fKelvinVoigt <<"\t\t# fKelvinVoigt #\n";
      output << tauKV <<"\t\t# tauKV #\n";
      output << scalKV <<"\t\t# scalKV #\n";
      output << rKV_LPF <<"\t\t# rKV_LPF #\n";
    }
    if (fMassWeightg) {
      output << fMassWeightg << "\t\t# fMassWeightg" << " #\n";
      output << zeroModeMass << "\t\t# zeroModeMass" << " #\n";
    }

    if (forceInit != 0) {
      output << forceInit << "\t\t# forceInit #\n";
      if (forceInit!=forceFinal) output << forceFinal << "\t\t# forceFinal #\n";
    }
    else {
      output << pressInit << "\t\t# pressInit #\n";
      if (pressInit!=pressFinal) output << pressFinal << "\t\t# pressFinal #\n";
    }
    if (fConstCOM) {
      output << fConstCOM << "\t\t# fConstCOM #\n";
      output << zConstCOM << "\t\t# zConstCOM #\n";
      output << vzConstCOM << "\t\t# vzConstCOM #\n";
    }
    if (nElast>=2) {
      if (fConstCOMy) {
        output << fConstCOMy << "\t\t# fConstCOMy #\n";
        output << yConstCOM << "\t\t# yConstCOM #\n";
        output << vyConstCOM << "\t\t# vyConstCOM #\n";
      }
    }
    if (nElast==3) {
      if (fConstCOMx) {
        output << fConstCOMx << "\t\t# fConstCOMx #\n";
        output << xConstCOM << "\t\t# xConstCOM #\n";
        output << vxConstCOM << "\t\t# vxConstCOM #\n";
      }
    }

    if (fSteppedRamp) {
      output << fSteppedRamp << "\t\t# fSteppedRamp #\n";
      output << rampSteps << "\t\t# rampSteps #\n";
      output << rampRelax << "\t\t# rampRelax #\n";
      if (fSteppedRamp==1) output << dzRamp << "\t\t# dzRamp #\n";
      else if (fSteppedRamp==2) output << dpRamp << "\t\t# dpRamp #\n";
    }

    if (rVeloTurnStep>0) output << rVeloTurnStep << "\t\t# rVeloTurnStep #\n";
    if (nVeloTurnStep>0) output << nVeloTurnStep <<"\t\t# nVeloTurnStep #\n";
    if (checkTurnPress) {
      if (forceTurnaround != 1e+85) output << forceTurnaround <<"\t\t# forceTurnaround #\n";
      else output << pressTurnaround <<"\t\t# pressTurnaround #\n";
    }
    if (nVeloTransition>0) output << nVeloTransition <<"\t\t# nVeloTransition #\n";
    if (fzOpposite) {
      output << fzOpposite <<"\t\t# fzOpposite #\n";
      output << zOpposite/nxny <<"\t\t# zOpposite #\n";
      if ( (!fThickness[0])&&nElast ) 
        output << thickness[0] << "\t\t# thickness0 #\n";
    }
    if (fExtrapolateInf) output << fExtrapolateInf <<"\t\t# fExtrapolateInf #\n";

    if (fOnSitePotential) {
      output << fOnSitePotential <<"\t\t# fOnSitePotential #\n";
      output << fOnSitePeriod <<"\t\t# fOnSitePeriod #\n";
      output << xOnSite <<"\t\t# xOnSite #\n";
      output << yOnSite <<"\t\t# yOnSite #\n";
      if (fOnSitePotential&4) {//NEW-OnHertz
        output << rXhertz <<"\t\t# rXhertz #\n";
        output << rYhertz <<"\t\t# rYhertz #\n";
        output << vXOnSite <<"\t\t# vXOnSite #\n";
        output << vYOnSite <<"\t\t# vYOnSite #\n";
        output << frictionCoeffOS <<"\t\t# frictionCoeffOS #\n";
        output << surfEnergOS <<"\t\t# surfEnergOS #\n";
        output << potCurveRelOS <<"\t\t# potCurveRelOS #\n";
      }
    }
    if (pressPhi!=90) output << pressPhi <<"\t\t# pressPhi #\n";
    if (pressTheta!=0) output << pressTheta <<"\t\t# pressTheta #\n";
      
  } // end if (nElast)

  if (fTopoRead) output << fTopoRead << "\t\t# fTopoRead #\n";

  if (fTopoAdd) {
    output << fTopoAdd << "\t\t# fTopoAdd #\n";
    if (fTopoAdd&1) {
      output << rXhertz << "\t\t# rXhertz #\n";
      if (rYhertz!=rXhertz) output << rYhertz << "\t\t# rYhertz #\n";
      if (hertzExpnt!=2) output << hertzExpnt << "\t\t# hertzExpnt #\n";
    }
    if (fTopoAdd&2) {
      output << hurst << "\t\t# hurst #\n";
      output << lambdaR << "\t\t# lambdaR #\n";
      output << lambdaS << "\t\t# lambdaS #\n";
      if (fRoughNorm!=1) output << fRoughNorm << "\t\t# fRoughNorm #\n";
      if (rRoughNorm!=1) output << rRoughNorm << "\t\t# rRoughNorm #\n";
      if (peklenik!=1)output << peklenik << "\t\t# peklenik #\n";
      output << fRollOff << "\t\t# fRollOff #\n";
      if (fBoxMuller) output << fBoxMuller << "\t\t# fBoxMuller #\n";
    }
    if (fTopoAdd&4) {
      output << fAddSWR << "\t\t# fAddSWR #\n";
      output << nqxAddSWR << "\t\t# nqxAddSWR #\n";
      if (nqxAddSWR!=nqyAddSWR) output << nqyAddSWR << "\t\t# nqyAddSWR #\n";
      output << heightSWR << "\t\t# heightSWR #\n";
    }
    if (fTopoAdd&8) {
      output << radiusFlatPunch << "\t\t# radiusFlatPunch #\n";
      output << heightFlatPunch << "\t\t# heightFlatPunch #\n";
    }
    if (fTopoAdd&16) {
      output << rSphere << "\t\t# rSphere #\n";
    }

    if (dzStep>0) output << dzStep << "\t\t# dzStep #\n";
  }

  if (resolMovie != 512) output << resolMovie << "\t\t# resolMovie #\n";
  if (f3dMovie) output << f3dMovie << "\t\t# f3dMovie #\n";

  if (fLateral) {
    output << fLateral << "\t\t# fLateral #\n";
    if (vX!=0) output << vX << "\t\t# vX #\n";
    if (vY!=0) output << vY << "\t\t# vY #\n";
  }
  if (frictRelax != 1) output << frictRelax << "\t\t# frictRelax #\n";

  if (fSteadySlide) output << fSteadySlide << "\t\t# fSteadySlide #\n";


  output << ID << "\t# sheet end" << endl << endl;
  output.close();

} // writeParams

void gfmdSheet::writeParamsDefault() {
  ofstream output("params.def",ofstream::app);
  output << "! ! ! ! ! ! ! ! ! sheet (default) values\n\n";

  output << 0    << "\t\t# nElast #\n";
  output << 0.5  << "\t\t# stiffnessN #\n";
  output << 1    << "\t\t# elastExpntN #\n\n";
  output << 0    << "\t\t# fThicknessN #\n";
  output << 0.25 << "\t\t# poissonN #\n";
  output << 1    << "\t\t# thicknessN #\n\n";
  output << 1    << "\t\t# fKelvinVoigt #\n";
  output << 1    << "\t\t# tauKV #\n";
  output << 1000 << "\t\t# scalKV #\n";
  output << 1    << "\t\t# rKV_LPF #\n\n";
  //NEW-MW
  output << 1    << "\t\t# fMaxwell #\n";
  output << 2    << "\t\t# nMaxwell #\n";
  output << 1    << "\t\t# massGFMD #\n";
  output << 1    << "\t\t# fDumpMaxwell #\n\n";

  output << 1 << "\t\t# fMassWeightg #\n";
  output << 1 << "\t\t# zeroModeMass #\n\n";

  output << 1     << "\t\t# fConstCOM #\n";
  output << 1     << "\t\t# zConstCOM #\n";
  output << 0.001 << "\t\t# vzConstCOM #\n\n";
  output << 1     << "\t\t# fConstCOMx #\n";
  output << 1     << "\t\t# xConstCOM #\n";
  output << 0.001 << "\t\t# vxConstCOM #\n\n";
  output << 1     << "\t\t# fConstCOMy #\n";
  output << 1     << "\t\t# yConstCOM #\n";
  output << 0.001 << "\t\t# vyConstCOM #\n\n";
  output << 1     << "\t\t# fSteppedRamp #\n";
  output << 0.01  << "\t\t# dzRamp #\n";
  output << 0.001 << "\t\t# dpRamp #\n";
  output << 130   << "\t\t# rampSteps #\n";
  output << 170   << "\t\t# rampRelax #\n\n";
  output << 0.5   << "\t\t# rVeloTurnStep #\n";
  output << 1000  << "\t\t# nVeloTurnStep #\n\n";
  output << 1.    << "\t\t# pressInit #\n";
  output << 1.    << "\t\t# pressFinal #\n";
  output << 1.    << "\t\t# forceInit #\n";
  output << 1.    << "\t\t# forceFinal #\n";
  output << 0.1   << "\t\t# pressTurnaround #\n\n";
  output << 0.1   << "\t\t# forceTurnaround #\n\n";
  output << 50    << "\t\t# nVeloTransition #\n";
  output << 1	    << "\t\t# fzOpposite #\n";
  output << 1     << "\t\t# fExtrapolateInf #\n";
  output << 0	    << "\t\t# zOpposite #\n\n";
  output << 1 	  << "\t\t# fOnSitePotential #\n";
  output << 1 	  << "\t\t# fOnSitePeriod #\n";
  output << 0     << "\t\t# xOnSite #\n";
  output << 0     << "\t\t# yOnSite #\n";
  output << 0     << "\t\t# vXOnSite #\n";
  output << 0     << "\t\t# vYOnSite #\n";
  output << 0.3   << "\t\t# frictionCoeffOS #\n";
  output << 1e-4  << "\t\t# surfEnergOS #\n";
  output << 0.4   << "\t\t# potCurveRelOS #\n";
  output << 45    << "\t\t# pressPhi # (deg)\n";
  output << 45    << "\t\t# pressTheta # (deg)\n";

  output << 1 << "\t\t# fLateral #\n";
  output << 1 << "\t\t# fSteadySlide #\n";
  output << 1 << "\t\t# vX #\n";
  output << 1 << "\t\t# vY #\n\n";
  output << 1 << "\t\t# frictRelax #\n";

  output << 1 << "\t\t# fTopoAdd #\n";
  output << 1 << "\t\t# fTopoRead #\n\n";
  
  output << 1 << "\t\t# rXhertz #\n";
  output << 1 << "\t\t# rYhertz #\n";
  output << 2 << "\t\t# hertzExpnt #\n\n";

  output << 0.2 << "\t\t# rSphere #\n\n";

  output << 0.8 << "\t\t# hurst #\n";
  output << 0.5 << "\t\t# lambdaR #\n";
  output << .05 << "\t\t# lambdaS #\n";
  output << 1.0 << "\t\t# rRoughNorm #\n";
  output << 1.0 << "\t\t# peklenik #\n";
  output << 1   << "\t\t# fRollOff #\n";
  output << 1   << "\t\t# fRoughNorm #\n";
  output << 1   << "\t\t# fBoxMuller #\n\n";

  output << 1  << "\t\t# fAddSWR #\n";
  output << 2  << "\t\t# nqxAddSWR #\n";
  output << 2  << "\t\t# nqyAddSWR #\n";
  output << .1 << "\t\t# heightSWR #\n\n";

  output << 0.2 << "\t\t# radiusFlatPunch #\n";
  output << 1.  << "\t\t# heightFlatPunch #\n";

  output << 0.2 << "\t\t# rSphere #\n\n";

  output << .01 << "\t\t# dzStep #\n\n";
  
  output << 0   << "\t\t# f3dMovie #\n";
  output << 512 << "\t\t# resolMovie #\n\n";

  output.close();
}

void gfmdSheet::initSystem() {

  // allocate all fields needed by FFTW, first "commodity fields";

  // fieldRFFT not always needed, but for ease of coding always allocated
  fieldRFFT = (double *)  fftw_malloc( nReal*sizeof(double) );
  fieldFFFT = (Complex *) fftw_malloc( nFour*sizeof(Complex));

  if (fTopography) initConfig();
  if (nElast) prepareDispZ();
  if (nElast>=2) prepareDispY();
  if (nElast==3) prepareDispX();
  if (fOnSitePotential==4) initFriction();
  if (fLateral) initLateral();
  if (nElast) prepareGFMD();

  kFastest = getLinCIndex(nx/2, ny/2);
  rCentral = getLinRIndex(nx/2, ny/2);

  if (fTopography) dumpConfig();
}

void gfmdSheet::initConfig() {

  equilPos = (double *) fftw_malloc( nReal * sizeof(double) );
  for (int k = 0; k < nReal; ++k) equilPos[k] = 0;

  if (fTopoRead&1) addConfigR();
  if (fTopoRead&2) addConfigF();

  if (fTopoAdd&1) addHertz();
  if (fTopoAdd&2) addSelfAffine();
  if (fTopoAdd&4) addSWR();
  if (fTopoAdd&8) addFlatPunch();
  if (fTopoAdd&16) addSphere();

  if (dzStep > 0) makeSteps();

  if (ID) {// (rigid indenter might come from below if ID>0)
    for (int k = 0; k < nReal; ++k) equilPos[k] *= -1;
  }
}

void gfmdSheet::addConfigR() {

  string fileName = konfigName + "real";
  ifstream test(fileName);
  if (!test.is_open()) termination("File "+fileName+ " does not exist.");
  test.close();
  cerr << "# Reading " + fileName + "\n";
  if (nx==1) readReal(fieldRFFT, nx, ny, fileName, 2);
  else readReal(fieldRFFT, nx, ny, fileName, 3);
  for (Lint k = 0; k < nReal; ++k) equilPos[k] += fieldRFFT[k];

  // make surface analysis
  vector<double> props(6);
  realSpaceAnalysis(fieldRFFT, &props, nx, ny);
  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet [in] roughness info start\n";
  output << sqrt(props[0]) << "\t" << props[5] << "\t# rms and max height\n";
  output << sqrt(props[1]) << "\t" << sqrt(props[2]) << "\t# x and y rms grad\n";
  output << sqrt(props[3]) << "\t" << sqrt(props[4]) << "\t# x and y rms curv\n";
  output << sqrt(props[1]+props[2]) << "\t"  << 1. / sqrt( (props[3]+props[4])/2 )
         << "\t# rms grad and radOfCurv\n";
  output << ID << "\t# sheet [in] roughness info end\n\n";
  output.close();
  
}

void gfmdSheet::addConfigF() {
  return;

  // MM2Change: infrastructure for reading Fourier
  string fileName = konfigName + "four";
  ifstream test(fileName);
  if (!test.is_open()) termination("File "+fileName+ " does not exist.");
  test.close();
  cerr << "# Reading " + fileName + "\n";
  // readFour(fieldFFFT, nx, ny, fileName, 3);
  for (Lint k = 0; k < nReal; ++k) equilPos[k] += fieldRFFT[k];

}

void gfmdSheet::addHertz() {
  for (Lint k=0; k<nxny; ++k) {
    get_irxiry(k);
    double deltaX = (irx-nxH)*dx, deltaY = (iry-ny/2)*dy;
    double height = deltaX*deltaX/rXhertz + deltaY*deltaY/rYhertz;
    height = pow(height, hertzExpnt/2);
    height /= hertzExpnt;
    equilPos[k] += height;
  }
}

void gfmdSheet::addSphere() {
  for (Lint k=0; k<nxny; ++k) {
    get_irxiry(k);
    double deltaX2 = (irx-nxH)*dx; deltaX2 *= deltaX2;
    double deltaY2 = (iry-ny/2)*dy; deltaY2 *= deltaY2;
    double rSphere2 = rSphere*rSphere;
    if (deltaX2 + deltaY2 < rSphere2)
      equilPos[k] += rSphere - sqrt(rSphere2 - deltaX2 - deltaY2);
    else equilPos[k] += rSphere;
  }
}

void gfmdSheet::addSelfAffine() {

  if (lambdaS>lengthX) cerr << "# lambdaS>lengthX\n";
  if (lambdaS>lengthY) cerr << "# lambdaS>lengthY\n";

  // NX, NY of a system (just) large enough to accommodate complete spectrum
  Lint NX = 0.1+2*lengthX/lambdaS;
  Lint NY = 0.1+2*lengthY/lambdaS;

  // redefine hurst for one-dimensional interface
  if (nx==1) {NX = 1; hurst -= 0.5;}
  if (ny==1) {NY = 1; hurst -= 0.5;}

  Lint NYHP1 = NY/2 + 1;
  Lint nFOUR = NX * NYHP1;

  // fftw_plan for heightF (fieldFFFT) --> heightR (fieldRFFT)
  fftw_plan heightF2R =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, fieldRFFT, FFTW_ESTIMATE);
  for (Lint k = 0; k < nFour; ++k) fieldFFFT[k] = 0;

  double qRoll = TWOPI / lambdaR;
  double qMax  = TWOPI / lambdaS;

  double varHeight = 0, gradX2 = 0, gradY2 = 0, curvX2 = 0, curvY2 = 0;

  for (int iWarm = 0; iWarm<100; ++iWarm) rrand();

  for (Lint K = 1; K < nFOUR; ++K) {

    // compute (hypothetical) indices
    Lint IQX = K/NYHP1, IQY = K%NYHP1;
    Lint JQX = abs(NX-IQX);
    IQX = (IQX<JQX) ? IQX : JQX;

    // check if wavenumber within cutoff
    double q = sqrt( peklenik*pow(IQX,2)*dqx2 + pow(IQY,2)*dqy2/peklenik );
    if (q>qMax) continue;

    // draw random numbers first
    double phase = rrand()*TWOPI;
    double bm = 1;
    if (fBoxMuller) bm = sqrt(-2.*log(rrand()))*cos(2*PI*rrand());
    
    // compute contribution to spectrum 
    double absHeight = sqrtSpec(q);
    double specLocal = absHeight*absHeight;
    double weight=2;
    if ( (IQY==0) || (IQY==NY/2) ) weight = 1;
    specLocal *= weight;

    // compute height characteristics

    // height variance
    varHeight += specLocal;

    // slope variance
    double qxN = IQX*dqx2*IQX;
    double qyN = IQY*dqx2*IQY;
    gradX2 += qxN * specLocal;
    gradY2 += qyN * specLocal;

    // curvature variance
    qxN *= qxN; qyN *= qyN;
    curvX2 += qxN * specLocal;
    curvY2 += qyN * specLocal;

    // assign to FT(height) if it fits
    if ( (IQX>nx/2) || (IQY>ny/2) ) continue;

    if (IQX!=JQX) iqx = IQX; 
    else iqx = nx-IQX;
    iqy = IQY;

    Lint k = iqx*nyHP1 + iqy;

    fieldFFFT[k] = bm * absHeight * exp(Complex(0,1)*phase);

    // enforce fieldFFFT({\bf q}) = fieldFFFT*(-{\bf q})
    if ( (iqy==0) || (iqy==ny/2) ) {
      if ( (iqx==0) || (iqx==nxH) ) {   // Fourier transform is purely real
        phase = (phase>PI) ? -1 : 1;
        fieldFFFT[k] = absHeight * phase;
      } 
      else if (iqx>nxH) {             // complex conjugate already drawn
        Lint kConj = getLinCIndex(nx-iqx,iqy);
        fieldFFFT[k]  = conj(fieldFFFT[kConj]);
      }
    }

  } // loop over k

  // normalize surface in desired way (height or rms gradient)
  double scaleWeight = rRoughNorm*rRoughNorm, grad2 = gradX2 + gradY2;
  if (fRoughNorm==1) scaleWeight /= grad2;
  else if (fRoughNorm==2) scaleWeight /= varHeight;
  else cerr << "# no rescaling in addSelfAffine of sheet " << ID << endl;

  varHeight *= scaleWeight;
  grad2 *= scaleWeight; gradX2 *= scaleWeight; gradY2 *= scaleWeight;
  curvX2 *= scaleWeight; curvY2 *= scaleWeight;

  scaleWeight = sqrt(scaleWeight);
  for (int k=0; k < nFour; ++k) fieldFFFT[k] *= scaleWeight;
  
  // transform to real space
  fftw_execute(heightF2R);

  // shift bottom layer to zero
  double minHeight = fieldRFFT[0], maxHeight = fieldRFFT[0];
  double height1 = 0, height2 = 0;
  for (Lint k = 1; k < nxny; ++k) {
    minHeight = min(minHeight, fieldRFFT[k]);
    maxHeight = max(maxHeight, fieldRFFT[k]);
    height1 += fieldRFFT[k];
    height2 += fieldRFFT[k]*fieldRFFT[k];
  } 
  maxHeight -= minHeight;
  height1 /= nxny;
  height2 /= nxny; height2 -= height1*height1;

  for (Lint k = 0; k < nxny; ++k) {
    double dummy = fieldRFFT[k] - minHeight;
    equilPos[k] += dummy;
  }

  // output measurements
  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet [self-affine] roughness info start\n";
  output << sqrt(height2)   << "\t\t# rms real-space height\n";
  output << sqrt(varHeight) << " " << maxHeight << "\t# rms and max height (Fourier)\n";
  output << sqrt(gradX2) << "  " << sqrt(gradY2) << "\t# x and y rms slope\n";
  output << sqrt(curvX2) << "  " << sqrt(curvY2) << "\t\t# x and y rms curve\n";
  output << sqrt(grad2) << "\t"  << 1. / sqrt( (curvX2+curvY2)/2 )
         << "\t# rms grad and radOfCurv\n";
  output << ID << "\t# sheet [self-affine] roughness info end\n\n";
  output.close();

  fftw_destroy_plan(heightF2R);

} // addSelfAffine

void gfmdSheet::addSWR() {

// fAddSWR = ...
// 1 single wave parallel to y
// 2 single wave parallel to x
// 3 square roughness, see Eq. (1) in Dapp and Muser, EPL 109, 44001 (2015)
// 4 hexagonal, see Eq. (2) in same paper
// 5 triangular, see Eq. (3)
// 6 hexagonal  again but x <--> y
// 7 triangular again but x <--> y
// 4 ... 7, correct lengthX/lengthY ratio must be imposed to yield correct results

  const double qx = nqxAddSWR*2*PI/lengthX, qy = nqyAddSWR*2*PI/lengthY;
  int int_nqxAddSWR = nqxAddSWR, int_nqyAddSWR = nqyAddSWR;

  if (nqxAddSWR != int_nqxAddSWR) cerr << "### potential commensurability issue with nqx in sheet " << ID << "\n";
  else if ( ((fAddSWR==4)||(fAddSWR==5)) && (int_nqxAddSWR%2) ) {
    cerr << "### potential commensurability issue with nqx in sheet " << ID << "\n";
  }
  if (nqyAddSWR != int_nqyAddSWR) cerr << "### potential commensurability issue with nqy in sheet " << ID << "\n";
  else if ( ((fAddSWR==6)||(fAddSWR==7)) && (int_nqyAddSWR%2) ) {
    cerr << "### potential commensurability issue with nqy in sheet " << ID << "\n";
  }

  for (Lint ix=0; ix < nx; ++ix) { double x = ix*dx;//-lengthX/2 + ix*dx;//TEMP
  for (Lint iy=0; iy < ny; ++iy) { double y = iy*dy;//-lengthY/2 + iy*dy;//TEMP
    double height;
    if (fAddSWR==1) height =  1 + cos(qx*x) ;
    if (fAddSWR==2) height =  1 + cos(qy*y) ;
    if (fAddSWR==3) height =  2 + cos(qx*x) + cos(qy*y);
    if (fAddSWR==4) height = sqrt(2./3) * ( 1.5 + 2*cos(qx*x)*cos(qy*y/2) + cos(qy*y));
    if (fAddSWR==5) height = sqrt(13.5) 
			 - sqrt(2./3) * ( 1.5 + 2*cos(qx*x)*cos(qy*y/2) + cos(qy*y));
    if (fAddSWR==6) height = sqrt(2./3) * ( 1.5 + 2*cos(qy*y)*cos(qx*x/2) + cos(qx*x));
    if (fAddSWR==7) height = sqrt(13.5)
                          - sqrt(2./3) * ( 1.5 + 2*cos(qy*y)*cos(qx*x/2) + cos(qx*x));
    equilPos[ix*ny+iy] += heightSWR * height;
  } }

  if (fAddSWR>3) { // report relative deviations of wavelengths in x and y
    ofstream output("params.out", ofstream::app);
    output << ID << "\t# sheet [addSWR] roughness info start\n";
    double lambdaRatio = (lengthX/nqxAddSWR) / (lengthY/nqyAddSWR);
    if ( (fAddSWR==4)||(fAddSWR==5) ) lambdaRatio /= (2./sqrt(3.));
    if ( (fAddSWR==6)||(fAddSWR==7) ) lambdaRatio *= (2./sqrt(3.));
    double violation = 1 - lambdaRatio;
    output << violation << "\t\t# wavelength mismatch\n";
    if (abs(violation)>0.1) 
    cerr << "\n# " << violation << " wavelength mismatch in ID "
	 << ID << endl;
    output << ID << "\t# sheet [addSWR] roughness info end\n\n";
    output.close();
  }

} // addSWR

void gfmdSheet::addFlatPunch() {

  int fFlatPunch = 1; // MM2change (allow for read in and dump)
  double r2cut = (1.+1e-6)*radiusFlatPunch*radiusFlatPunch;
  for (Lint k=0; k < nReal; ++k) {
    get_irxiry(k);
    double deltaX = (irx-nxH)*dx;
    double deltaY = (iry-ny/2)*dy;
    double r2 = deltaX*deltaX+deltaY*deltaY;
//  if (fIndMode==0) r2 += deltaX*deltaX;
    if (r2>2*r2cut) equilPos[k] += heightFlatPunch;
    else if (r2>r2cut) {
      if (fFlatPunch==1) equilPos[k] += heightFlatPunch;
      else equilPos[k] += 
      heightFlatPunch*(cos(PI*(r2-2*r2cut)/(r2cut))+1)/2;
    }
  }


}


void gfmdSheet::makeSteps() {
  for (Lint k = 0; k < nReal; ++k) {
    double zLocal = equilPos[k];
    zLocal /= dzStep;
    zLocal = dzStep*round(zLocal);
    equilPos[k] = zLocal;
  }
}


double gfmdSheet::sqrtSpec(double q) {
  double qRoll = TWOPI / lambdaR;
  double qRed = q/qRoll;
  if (fRollOff==1) return( pow(1./sqrt(1.+qRed*qRed),1+hurst) ); // smooth roll-off
  if (qRed<1.-1e-9) {
    if (fRollOff==0) return(0.);        // cut-off
    else if (fRollOff==2) return(1.);   // roll-off with cusp
  }
  return(pow(1./qRed,1+hurst));
}


void gfmdSheet::prepareDispZ() {
  cerr << "   gfmdSheet" << ID << "::prepareDispZ()\n";//FLOW

  // initialize fields
  dispZR = (double *)  fftw_malloc(nReal*sizeof(double) );
  dispZF = (Complex *) fftw_malloc(nFour*sizeof(Complex));
  int iDummy = 0;
  if (fConstCOM) iDummy = 1;
  for (Lint k = 0; k < nReal; ++k) dispZR[k] = iDummy*zConstCOM;

  stressZR = (double *)  fftw_malloc( nReal*sizeof(double) );
  for (Lint k = 0; k < nReal; ++k) stressZR[k] = 0;

  dispZF    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  dispZFold = (Complex *) fftw_malloc( nFour*sizeof(Complex) );

  for (Lint k = 1; k < nFour; ++k) dispZFold[k] = dispZF[k] = 0.;

  stressZF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  if (fKelvinVoigt) rhsKV_LPF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );

  if (fSteadySlide) preFacCorrSS = (Complex *) fftw_malloc( nFour*sizeof(Complex) );

  // read displacemenmt

  string fileName = konfigName + "old";
  ifstream test(fileName);
  if (test.is_open()) {
    cerr << "# Reading Z displacements\n";
    test.close();
    if ( (nx!=1)&&(ny!=1) ) {
      if (fTopography) readReal(dispZR, nx, ny, fileName, 4);
      else readReal(dispZR, nx, ny, fileName, 3);
    } 
    else {
      if (fTopography) readReal(dispZR, nx, ny, fileName, 3);
      else readReal(dispZR, nx, ny, fileName, 2);
    }
  }

  // make all plans

  dispZR2F =
  fftw_plan_dft_r2c_2d(nx, ny, dispZR, (fftw_complex*) dispZF, FFTW_ESTIMATE);

  dispZF2R =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, dispZR, FFTW_ESTIMATE);

  stressZR2F =
  fftw_plan_dft_r2c_2d(nx, ny, stressZR, (fftw_complex*) stressZF, FFTW_ESTIMATE);

  // initialize (Fourier) displacements from file

  fftw_execute(dispZR2F);
  for (int k = 0; k < nFour; ++k) dispZFold[k] = dispZF[k];
  if (fConstCOM) {
    dispZFold[0] = dispZF[0] = zConstCOM*nxny;
    dispF2R();
  }

  if (fSteadySlide) {
    fConstCOM = 1;
    zConstCOM = dispZF[0].real()/nxny;
  }

  // change MM2SS: Better way to initialize?
  if (fKelvinVoigt) {
    for (Lint k = 0; k < nFour; ++k) rhsKV_LPF[k] = 0;
  }

  //NEW-MW
  if (fMaxwell){
    dispZMw.resize(nMaxwell);
    for (int iMw = 0; iMw < nMaxwell; ++iMw){
      dispZMw[iMw] = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
      for (Lint k = 0; k < nFour; ++k) dispZMw[iMw][k] = 0.;
      //if (fConstCOM) dispZMw[iMw][0] = zConstCOM*nxny;
    }

    // read Maxwell displacements
    string fileName = "konfig" + to_string(ID) + "Mw.old";
    ifstream test(fileName);
    if (test.is_open()) {
      test.close();
      for (int iMw=0; iMw<nMaxwell; ++iMw) {
        cerr << "# Reading Maxwell element " << iMw << "\n";
        readFour(dispZMw[iMw], nx, ny, fileName, 3+2*iMw);
      }
    }
  }

  // initialize steady state sliding prefactors
  if (fSteadySlide) {
    for (int k = 0; k < nFour; ++k) {
      // compute frequency of mode
      get_iqxiqy(k);
      if (iqx>nxH)   iqx -= nxH;
      if (iqy>nyHP1) iqy -= nyHP1;
      double qx = iqx*dqx, qy = iqy*dqy;
      double omega = qx*vX + qy*vY;
      if (fKelvinVoigt) {
        preFacCorrSS[k]  = (1.+Complex(0,1)*omega*tauKV/scalKV);
        preFacCorrSS[k] /= (1.+Complex(0,1)*omega*tauKV);
      } 
      else {
        termination("Steady-state sliding only implemented for KV.");
      }
    }
  }

  if (fExtrapolateInf) {
    zInfOffset = nxny*dispInf() - dispZF[0].real();
  }
} // prepareDispZ

void gfmdSheet::prepareDispY() {
  cerr << "   gfmdSheet" << ID << "::prepareDispY()\n";//FLOW

  if (ny==1) termination("nElast>1 and ny==1 not meaningful!");

  // initialize fields
  dispYR = (double *)  fftw_malloc(nReal*sizeof(double) );
  dispYF = (Complex *) fftw_malloc(nFour*sizeof(Complex));
  int iDummy = 0;
  if (fConstCOMy) iDummy = 1;
  for (Lint k = 0; k < nReal; ++k) dispYR[k] = iDummy*yConstCOM;

  stressYR = (double *)  fftw_malloc( nReal*sizeof(double) );
  for (Lint k = 0; k < nReal; ++k) stressYR[k] = 0;

  dispYF    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  dispYFold = (Complex *) fftw_malloc( nFour*sizeof(Complex) );

  for (Lint k = 1; k < nFour; ++k) dispYFold[k] = dispYF[k] = 0.;
  //TODO: 0-mode?

  stressYF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  if (fKelvinVoigt) {
    termination("No Kelvin Voigt for tangential displacements implemented.");
    //rhsKV_LPF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  }

  // read displacemenmt

  string fileName = konfigName + "old";
  ifstream test(fileName);
  int iColumn = 4;//3;//CM-dumpReal
  if (test.is_open()) {
    cerr << "# Reading Y positions\n";
    test.close();
    if ( (nx!=1)&&(ny!=1) ) {
      if (fTopography) readReal(dispYR, nx, ny, fileName, iColumn+2);
      else readReal(dispYR, nx, ny, fileName, iColumn+1);
    } 
    else {
      if (fTopography) readReal(dispYR, nx, ny, fileName, iColumn+1);
      else readReal(dispYR, nx, ny, fileName, iColumn);
    }
  }

  // make plans

  dispYR2F =
  fftw_plan_dft_r2c_2d(nx, ny, dispYR, (fftw_complex*) dispYF, FFTW_ESTIMATE);

  dispYF2R =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, dispYR, FFTW_ESTIMATE);

  stressYR2F =
  fftw_plan_dft_r2c_2d(nx, ny, stressYR, (fftw_complex*) stressYF, FFTW_ESTIMATE);

  // initialize (Fourier) displacements from file

  fftw_execute(dispYR2F);
  for (int k = 0; k < nFour; ++k) dispYFold[k] = dispYF[k];
  if (fConstCOMy) {
    dispYFold[0] = dispYF[0] = yConstCOM*nxny;
    dispF2R();
  }

  //TODO: KV in the future

}

void gfmdSheet::prepareDispX() {
  cerr << "   gfmdSheet" << ID << "::prepareDispX()\n";//FLOW

  if (nx==1) termination("nElast==3 and nx==1 not meaningful!");

  // initialize fields
  dispXR = (double *)  fftw_malloc(nReal*sizeof(double) );
  dispXF = (Complex *) fftw_malloc(nFour*sizeof(Complex));
  int iDummy = 0;
  if (fConstCOMx) iDummy = 1;
  for (Lint k = 0; k < nReal; ++k) dispXR[k] = iDummy*xConstCOM;

  stressXR = (double *)  fftw_malloc( nReal*sizeof(double) );
  for (Lint k = 0; k < nReal; ++k) stressXR[k] = 0;

  dispXF    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  dispXFold = (Complex *) fftw_malloc( nFour*sizeof(Complex) );

  for (Lint k = 1; k < nFour; ++k) dispXFold[k] = dispXF[k] = 0.;
  //TODO 0-mode?

  stressXF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  if (fKelvinVoigt) {
    termination("No Kelvin Voigt for tangential displacements implemented.");
    //rhsKV_LPF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  }

  // read displacemenmt

  string fileName = konfigName + "old";
  ifstream test(fileName);
  int iColumn = 4; if (nElast>=2) iColumn += 2;
  if (test.is_open()) {
    cerr << "# Reading X positions\n";//TEMP
    test.close();
    if ( (nx!=1)&&(ny!=1) ) {
      if (fTopography) readReal(dispXR, nx, ny, fileName, iColumn+2);
      else readReal(dispXR, nx, ny, fileName, iColumn+1);
    } 
    else {
      if (fTopography) readReal(dispXR, nx, ny, fileName, iColumn+1);
      else readReal(dispXR, nx, ny, fileName, iColumn);
    }
  }

  dispXR2F =
  fftw_plan_dft_r2c_2d(nx, ny, dispXR, (fftw_complex*) dispXF, FFTW_ESTIMATE);

  dispXF2R =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, dispXR, FFTW_ESTIMATE);

  stressXR2F =
  fftw_plan_dft_r2c_2d(nx, ny, stressXR, (fftw_complex*) stressXF, FFTW_ESTIMATE);

  // initialize (Fourier) displacements from file

  fftw_execute(dispXR2F);
  for (int k = 0; k < nFour; ++k) dispXFold[k] = dispXF[k];
  if (fConstCOMx) {
    dispXFold[0] = dispXF[0] = xConstCOM*nxny;
    dispF2R();
  }

  //TODO:  KV in the future

}

void gfmdSheet::initFriction() {
  fFriction = 1;
  if (fOnSitePotential==4 && vYOnSite==0 && vXOnSite==0) fFriction = 0;

  if (nElast>=2) {
    int iDummy = 0;
    if (fConstCOMy) iDummy = 1;
  }

  if (nElast==3) {
    int iDummy = 0;
    if (fConstCOMx) iDummy = 1;
  }
}


void gfmdSheet::prepareGFMD() {
  cerr << "   gfmdSheet" << ID << "::prepareGFMD()\n";//FLOW

  // memory allocation
  if (!fKelvinVoigt) invMass = (double *) fftw_malloc( nFour*sizeof(double) );
  weightQ = (short int *) fftw_malloc( nFour*sizeof(short int) );
  stressZZPre = (double *) fftw_malloc( nFour*sizeof(double) );
  for (Lint k = 0; k < nFour; ++k) stressZZPre[k] = 0;
  
  // initialize weights
  weightQ[0] = 1;
  for (Lint k = 1; k < nFour; ++k) {
    get_iqxiqy(k);
    weightQ[k] = 2;
    if ( (iqy==0) || (iqy==ny/2) ) weightQ[k] = 1;
  }

  // initialize Green's functions (stiffnesses)
  if (nElast==1) prepareStiffZ(); // GFMD only in Z direction
  else if (nElast==2) prepareStiffZY(); // GFMD in Z and Y directions
  else if (nElast==3) prepareStiffZYX(); // GFMD in all three directions

  // determine min and max stiffness
  double q2Max = 0, q2Min = min(dqx2,dqx2);
  if (nx>1) q2Max += pow(PI/dx,2);
  if (ny>1) q2Max += pow(PI/dy,2);
  if (nx==1) q2Min = dqy2;
  if (ny==1) q2Min = dqx2;
  double qMax = sqrt(q2Max), qMin = sqrt(q2Min);
  stiffMax = contModEff = stiffMin = 0;
  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    stiffMax += stiffness[iLayer]*pow(sqrt(q2Max),elastExpnt[iLayer]);
    stiffMin += stiffness[iLayer]*pow(sqrt(q2Min),elastExpnt[iLayer]);
  }
  contModEff = 2*stiffMax / sqrt(q2Max);

  // initialize masses and damping
  if (!fKelvinVoigt) {
    for (Lint k=0; k < nFour; ++k) invMass[k] = 1;
  }
  if (fMassWeightg) {
    if (fMaxwell) {
      for (int k = 0; k < nFour; ++k)
        invMass[k] = stiffHigh[0]/(stressZZPre[k]*massGFMD);
      if (stressZZPre[0] == 0) invMass[0] = 2*invMass[1];
    }
    else if (zeroModeMass==0) {
      for (int k = 0; k < nFour; ++k) {
        double stiffnessZ = stressZZPre[k];
        //if (nElast > 1) stiffnessZ = stressZZPre[k] + pow(ImStressXZPre[k],2)/stressXXPre[k];//OPTION
        invMass[k] = dampGlobal*dampGlobal/(4*stiffnessZ);
      }
      if (stressZZPre[0] == 0) invMass[0] = 2*invMass[1];
    }
    else {
      for (int k = 0; k < nFour; ++k) {
        double stiffnessZ = stressZZPre[k];
        //if (nElast > 1) stiffnessZ = stressZZPre[k] + pow(ImStressXZPre[k],2)/stressXXPre[k];//OPTION
        invMass[k] = pow(zeroModeMass*stiffMin,2.) + pow(stiffnessZ,2);
        invMass[k] = 1./sqrt(invMass[k]);
      }
    }
    damping = dampGlobal; 
  } 
  else {
    double stiffEff = sqrt(0*pressure*pressure/areaXY + stiffMin*stiffMin); 
    damping = dampGlobal*sqrt(stiffEff); 
  }
  kFastest = getLinCIndex(nx/2, ny/2);


  // dump info to params.out
  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet elasticity info start\n";
  if (!fConstCOM) output << pressure*areaXY << "\t\t# total force\n";
  output << stiffMin << "\t\t# stiffMin\n";
  output << stiffMax << "\t\t# stiffMax\n";
  output << stressZZPre[1] << "\t\t# stressZZPre[1]\n";
  output << stressZZPre[kFastest] << "\t\t# stressZZPre[kFastest]\n";
  if (nElast>=2) {
    output << stressYYPre[1] << "\t\t# stressYYPre[1]\n";
    output << stressYYPre[kFastest] << "\t\t# stressYYPre[kFastest]\n";
    output << ImStressYZPre[1] << "\t\t# ImStressYZPre[1]\n";
    output << ImStressYZPre[kFastest] << "\t\t# ImStressYZPre[kFastest]\n";
  }
  if (nElast==3) {
    output << stressXXPre[1] << "\t\t# stressXXPre[1]\n";
    output << stressXXPre[kFastest] << "\t\t# stressXXPre[kFastest]\n";
    output << stressXYPre[1] << "\t\t# stressXYPre[1]\n";
    output << stressXYPre[kFastest] << "\t\t# stressXYPre[kFastest]\n";
    output << ImStressXZPre[1] << "\t\t# ImStressXZPre[1]\n";
    output << ImStressXZPre[kFastest] << "\t\t# ImStressXZPre[kFastest]\n";
  }
  if (!fKelvinVoigt) {
    output << 1./invMass[0] << "\t\t# massCOM\n";
    output << 1./invMass[kFastest] << "\t\t# massFast\n";
  }
  output << damping  << "\t\t# damping\n";
  if (stressZZPre[0]!=0) output << stressZZPre[0] << "\t\t# stressZZPre[0]\n";
  output << ID << "\t# sheet elasticity info end\n\n";
  output.close();


  /*//DEBUG
  if (!fKelvinVoigt) {
    cerr << "\n# Check for critical damping condition:";
    cerr << "\n# damping^2       = " << damping*damping;
    cerr << "\n# 4 * stiffHigh/m[1]   = " << 4*stressZZPre[1] * invMass[1];
    cerr << "\n# 4 * stiffHigh/m[-1]  = " << 4*stressZZPre[nFour-1] * invMass[nFour-1];
    cerr << "\n";
  }//*/

  /*//DEBUG
  cerr << "# Exporting stress prefactors.\n";
  ofstream debug;
  debug.open("stressPre.dat");
  for (Lint k=0; k < nFour; ++k) {
    get_iqxiqy(k);
    debug << iqx << "\t" << iqy << "\t";
    debug << getQ(k);
    debug << "\t" << stressXXPre[k];
    debug << "\t" << stressXYPre[k];
    debug << "\t" << stressYYPre[k];
    debug << "\t" << stressZZPre[k];
    debug << "\t" << ImStressXZPre[k];
    debug << "\t" << ImStressYZPre[k];
    debug << "\n";
  }
  debug.close();//*/

}


void gfmdSheet::prepareStiffZ() {
  cerr << "   gfmdSheet" << ID << "::prepareStiffZ()\n";//FLOW

  // initialize 0-mode stiffness (if fThickness==2)
  stressZZPre[0] = 0;
  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    if (fThickness[iLayer] == 2) {
      double nu = poisson[iLayer], thick = thickness[iLayer], stiff = stiffHigh[iLayer];
      stressZZPre[0] += (stiff/thick) * (2*(1.-nu)*(1.-nu))/(1.-2*nu); // Carbone
      //stressZZPre[0] += (stiff/thick) * (1.-nu*nu); // E_Young
    }
  }

  // initialize finite-q stiffnesses
  for (Lint k = 1; k < nFour; ++k) {
    double q = getQ(k);
    for (int iLayer=0; iLayer<nLayer; ++iLayer) {
      
      // semi-infinite
      double preFac = 1;
      double qh = thickness[iLayer]*q;
      double width2 = qh*qh;
      const double POISSON = poisson[iLayer];

      // from Carbone et al. J. Mech. Phys. Solids 56 (2008)

      // finite thickness, bottom boundary constant stress
      if (fThickness[iLayer]==1) {
        preFac = cosh(2*qh)-2*width2-1;
        preFac /= sinh(2*qh)+2*qh;
        if ( qh>40 ) preFac = 1;// sinh(qh),cosh(qh) become infinite for large qh
      } 

      // finite thickness, bottom boundary constant displacement
      else if (fThickness[iLayer]==2) {
        preFac  = (3-4*POISSON)*cosh(2*qh)+2*width2-4*POISSON*(3-2*POISSON)+5;
        preFac /= (3-4*POISSON)*sinh(2*qh)-2*qh;
        if (qh>40) preFac = 1;// sinh(qh),cosh(qh) become infinite for large qh
      }

      // Green's function with qh-dependet q-Factor
      stressZZPre[k] += preFac * pow(q,elastExpnt[iLayer]) * stiffHigh[iLayer];
    }
  }

  if (fConstCOM && fzOpposite && stressZZPre[0]==0) stressZZPre[0] = stressZZPre[1];//CM-FixCOM
  if (fKelvinVoigt && stressZZPre[0]==0) stressZZPre[0] = stressZZPre[1];//CM-FixCOM
  //cout << "k[0]=" << stressZZPre[0] << ", k[1]=" << stressZZPre[1] << "\n";//DEBUG
}


void gfmdSheet::prepareStiffZY() {
  cerr << "   gfmdSheet" << ID << "::prepareStiffZY()\n";//FLOW
  if (fThickness[0] == 1) termination("(1+1)D GFs for fThickness=1 not implemented!");
  if (fKelvinVoigt) termination("(1+1)D doesn't do Kelvin Voigt.");

  // memory allocation
  stressYYPre = (double *) fftw_malloc( nFour*sizeof(double) );
  ImStressYZPre = (double *) fftw_malloc( nFour*sizeof(double) );
  stressYYPre[0] = 0;
  ImStressYZPre[0] = 0;

  const double POISSON = poisson[0];
  const double STIFFNESS = stiffHigh[0];
  const double THICKNESS = thickness[0];

  stressYYPre[0] = STIFFNESS/(1-POISSON)/8/THICKNESS;
  if (fThickness[0] == 2) {
    stressZZPre[0] = STIFFNESS/(1-2*POISSON)/4/THICKNESS;
  }

  // Construct Green's tensor
  double GreenXX, GreenZZ, ImGreenXZ;
  for (Lint k = 1; k < nFour; ++k) {
    get_iqxiqy(k);
    double q = dqy * iqy;

    // from Menga. Int. J. Solids Struct. 164, 212-220 (2019)

    // semi-infinite
    double preFacXX = 1;
    double preFacXZ = 0.5/(1.-POISSON);
    double preFacZZ = 1;
    if (fThickness[0]==0) preFacXZ *= (2*POISSON-1);
    //TODO: Add LPF here, too?

    // finite thickness, bottom boundary displacement = 0
    else if (fThickness[0]==2) {
      double qh = THICKNESS*q;
      double width2 = qh*qh;

      if (qh < 40) {// sinh(qh),cosh(qh) become infinite for large qh
        double normG = 5+2*width2-4*POISSON*(3-2*POISSON)+(3-4*POISSON)*cosh(2*qh);

        preFacXX *= +2*qh + (3-4*POISSON)*sinh(2*qh); //Menga2019 (A.9)
        preFacXX /= normG;

        preFacXZ *= +3 + 2*width2 - 2*POISSON*(5-4*POISSON) - (3-2*POISSON*(5-4*POISSON))*cosh(2*qh); //Menga2019 (B.3)
        preFacXZ /= normG;

        preFacZZ *= -2*qh + (3-4*POISSON)*sinh(2*qh); //Carbone2008 (B.5)
        preFacZZ /= normG;
      }
      else preFacXZ *= (2*POISSON-1);

    }

    // Green's function with qh-dependet q-Factor
    GreenXX = preFacXX * pow(q,-elastExpnt[0]) / STIFFNESS;
    ImGreenXZ = preFacXZ * pow(q,-elastExpnt[0]) / STIFFNESS;
    GreenZZ = preFacZZ * pow(q,-elastExpnt[0]) / STIFFNESS;

    // invert Green's tensor to obtain stiffness tensor Phi
    double PhiXX, PhiZZ, ImPhiXZ;
    double detG = GreenXX*GreenZZ - ImGreenXZ*ImGreenXZ;
    PhiXX = GreenZZ/detG;
    ImPhiXZ = -ImGreenXZ/detG;
    PhiZZ = GreenXX/detG;

    // save to array
    stressYYPre[k] = PhiXX;
    ImStressYZPre[k] = ImPhiXZ;
    stressZZPre[k] = PhiZZ;

  }
}


void gfmdSheet::prepareStiffZYX() {
  cerr << "   gfmdSheet" << ID << "::prepareStiffZYX()\n";//FLOW
  if (fThickness[0] == 1) termination("(2+1)D GFs for fThickness=1 not implemented!");
  if (fKelvinVoigt) termination("(2+1)D doesn't do Kelvin Voigt.");

  // memory allocation
  stressXXPre = (double *) fftw_malloc( nFour*sizeof(double) );
  stressYYPre = (double *) fftw_malloc( nFour*sizeof(double) );
  ImStressXZPre = (double *) fftw_malloc( nFour*sizeof(double) );
  ImStressYZPre = (double *) fftw_malloc( nFour*sizeof(double) );
  stressXYPre = (double *) fftw_malloc( nFour*sizeof(double) );

  const double POISSON = poisson[0];
  const double STIFFNESS = stiffHigh[0];
  const double THICKNESS = thickness[0];

  int nyH = nyHP1 - 1;
  for (Lint k = 1; k < nFour; ++k) {
    int iqxLoc = k/nyHP1; 
    int iqyLoc = k%nyHP1;
    if (iqxLoc > nxH) iqxLoc = iqxLoc - nx;

    double q_x = iqxLoc*dqx, q_y = iqyLoc*dqy;
    double qXmax = PI/dx, qYmax = PI/dy;
    double q = sqrt(q_x*q_x + q_y*q_y);

    //CM2MM: should the LPF range (i.e. fix for imaginary components) be adjustable?
    double stiffLPF = 1;

    //OPTION: set only the edges of the Fourier domain to 0
    //int hack_nxH = 1;
    //int hack_nyH = 1;
    //if (iqxLoc == nxH) hack_nxH = 0; // XY, XZ
    //if (iqyLoc == nyH) hack_nyH = 0; // XY, YZ

    //OPTION: filter making stiffnesses globally differentiable, except sign changes
    //q_x = qXmax*sin(PI*q_x/qXmax/2);
    //q_y = qYmax*sin(PI*q_y/qYmax/2);
    
    //OPTION: low-pass filter with four-fold symmetry
    //stiffLPF *= 0.5*(1+cos(PI*q_x/qXmax)); // full range LPF
    //stiffLPF *= 0.5*(1+cos(PI*q_y/qYmax)); // full range LPF
    //if (abs(q_x) > qXmax/2) stiffLPF *= 0.5*(1-cos(2*PI*q_x/qXmax)); // half range LPF
    //if (abs(q_y) > qYmax/2) stiffLPF *= 0.5*(1-cos(2*PI*q_y/qYmax)); // half range LPF
    
    //OPTION: low-pass filter with rotational symmetry (i.e. isotropic)
    /*double qMax = min(qXmax,qYmax);
    stiffLPF *= 0.5*(1+cos(PI*q/qMax)); // full range LPF
    if (q > qMax) stiffLPF = 0;         // full range LPF*/
    //if (q > qMax) stiffLPF = 0;                               // half range LPF
    //else if (q > qMax/2) stiffLPF *= 0.5*(1-cos(2*PI*q/qMax)) // half range LPF

    //cout << iqxLoc << "\t" << iqyLoc << "\t" << stiffLPF << "\n";//DEBUG

    // semi-infinite
    double preFacXX = 4*(1-POISSON)*(1-POISSON)/(3-4*POISSON);
    double preFacXZ = 2*(1-POISSON)*(1-2*POISSON)/(3-4*POISSON);
    double preFacYY = 1-POISSON;
    double preFacZZ = 4*(1-POISSON)*(1-POISSON)/(3-4*POISSON);

    // finite thickness, bottom boundary displacement = 0
    if (fThickness[0]==2) {
      double qh = THICKNESS*q;

      if (qh < 40) {// sinh(qh),cosh(qh) become infinite for large qh
        double denominator = pow((3-4*POISSON)*sinh(qh),2) - pow(qh,2);

        preFacXX = 2*pow(1-POISSON,2)*(-2*qh + (3-4*POISSON)*sinh(2*qh));
        preFacXX /= denominator;

        preFacYY = (1-POISSON)/tanh(qh);

        preFacXZ = 2*(1-POISSON)*(-pow(qh,2) + (1-2*POISSON)*(3-4*POISSON)*pow(sinh(qh),2));
        preFacXZ /= denominator;

        preFacZZ = 2*pow(1-POISSON,2)*( 2*qh + (3-4*POISSON)*sinh(2*qh));
        preFacZZ /= denominator;
      }
    }

    // rotation using Gxx, Gyy, Gxz, Gzz
    double cosFac = q_x/q, sinFac = q_y/q;
    double stiffXX = preFacXX*cosFac*cosFac + preFacYY*sinFac*sinFac;
    double stiffXY = stiffLPF*(preFacXX - preFacYY)*cosFac*sinFac;
    double stiffXZ = stiffLPF*preFacXZ*cosFac;
    double stiffYY = preFacXX*sinFac*sinFac + preFacYY*cosFac*cosFac;
    double stiffYZ = stiffLPF*preFacXZ*sinFac;
    double stiffZZ = preFacZZ;

    double qFac = pow(q,elastExpnt[0]);
    stressXXPre[k] = stiffXX * qFac * STIFFNESS;
    stressXYPre[k] = stiffXY * qFac * STIFFNESS;
    stressYYPre[k] = stiffYY * qFac * STIFFNESS;
    stressZZPre[k] = stiffZZ * qFac * STIFFNESS;

    ImStressXZPre[k] = stiffXZ * qFac * STIFFNESS;
    ImStressYZPre[k] = stiffYZ * qFac * STIFFNESS;
  }

  // stiffness of 0-mode: always assigns finite stiffness in lateral directions for friction simulations!
  stressXXPre[0] = (1-POISSON)*STIFFNESS/THICKNESS;
  stressXYPre[0] = 0;
  stressYYPre[0] = (1-POISSON)*STIFFNESS/THICKNESS;
  ImStressXZPre[0] = 0;
  ImStressYZPre[0] = 0;
  if (fThickness[0]==2) stressZZPre[0] += 2*STIFFNESS*(1-POISSON)*(1-POISSON)/(1-2*POISSON)/THICKNESS;
  else stressZZPre[0] = 0;
}


void gfmdSheet::initLateral() {

  if ( (vX!=0) || (vY!=0) ) {
    fSliding = true;		// reset class variable
    fExternalDriving = true;	// reset global variable
  }

  if (fTopography) {
    equilPos0F = (Complex *) fftw_malloc( nFour * sizeof(Complex) );
    for (int k=0; k < nReal; ++k) fieldRFFT[k] = equilPos[k]; 
    fftw_plan equilPosR2F =
    fftw_plan_dft_r2c_2d(nx, ny, fieldRFFT, (fftw_complex*) equilPos0F, FFTW_ESTIMATE);
    fftw_execute(equilPosR2F);
    fftw_destroy_plan(equilPosR2F);
    equilPosF2R =
    fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, equilPos , FFTW_ESTIMATE);
  }
}


void gfmdSheet::moveLateral() {
  double xShiftRel = -vX*dTime;
  double yShiftRel = -vY*dTime;
  double xShiftAbs;
  double yShiftAbs;
  if (frictRelax != 1) {
    xShiftAbs = xShiftRel * ((int) iTime/frictRelax);
    yShiftAbs = yShiftRel * ((int) iTime/frictRelax);
  }
  else {// Fire in-/decrements are allowed, so mdTime might not be dTime*iTime
    xShiftAbs = vX*mdTime;
    yShiftAbs = vY*mdTime;
  }

  if (fTopography) {
    for (Lint k = 0; k < nFour; ++k) {
      get_iqxiqy(k);
      if (iqx >= nx/2) iqx-=nx;
      double qx = iqx*dqx, qy = iqy*dqy;
      fieldFFFT[k] = exp(Complex(0,1)*(qx*xShiftAbs+qy*yShiftAbs)) * equilPos0F[k];
      fieldFFFT[k] /= (double) nxny;
    }
    fftw_execute(equilPosF2R);
  }

  if (nElast) {
    termination("dispZFold might have to be shifted too!"); 
    for (Lint k = 0; k < nFour; ++k) {
      get_iqxiqy(k);
      if (iqx >= nx/2) iqx-=nx;
      double qx = iqx*dqx, qy = iqy*dqy;
      fieldFFFT[k] = exp(Complex(0,1)*(qx*xShiftRel+qy*yShiftRel)) * dispZF[k];
      fieldFFFT[k] /= (double) nxny;
    }
    fftw_execute(dispZF2R);
  }

  return;

}


double gfmdSheet::getQ(Lint k) {
  get_iqxiqy(k);
  int jqx = abs(iqx-nx);
  jqx = (iqx<jqx) ? iqx : jqx;
  double q2 = jqx*jqx*dqx2 + iqy*iqy*dqy2;
  return (sqrt(q2));
}


void gfmdSheet::dumpFrame(int iFrame) {

  // no need to dump frame of static equilPos
  if ((fSheetType<2)&&(!fLateral)) return; 

  int xStep = (nx>resolMovie) ? (nx-1)/resolMovie + 1: 1;
  int yStep = (ny>resolMovie) ? (ny-1)/resolMovie + 1: 1;

  string fileName = "frames/" + konfigName + to_string(iFrame) + ".dat";

  if (nElast && !fRealSpaceStress) stressF2R();

  for (int k = 0; k < nReal; ++k) {
    fieldRFFT[k] = 0;
    if (nElast) fieldRFFT[k]  = dispZR[k];
    if (fTopography) fieldRFFT[k] += equilPos[k];
  }

  // fields to be dumped
  vector<double*> arrays = {fieldRFFT};
  vector<string> arrayNames = {"Z"};
  if (fSheetType>1) {
    arrays.push_back(stressZR);
    arrayNames.push_back("stressZ");
  }
  if (nElast>=2) {
    arrays.push_back(dispYR);
    arrayNames.push_back("dispY");
    arrays.push_back(stressYR);
    arrayNames.push_back("stressY");
  }
  if (nElast==3) {
    arrays.push_back(dispXR);
    arrayNames.push_back("dispX");
    arrays.push_back(stressXR);
    arrayNames.push_back("stressX");
  }

  //memcpy(fieldRFFT, stressZR,  nReal*sizeof(double));
  //memcpy(fieldFFFT, stressZF,  nFour*sizeof(Complex));
  //updateStress();

  if (f3dMovie) dumpReal(arrays, nx, ny, xStep, yStep, fileName);
  dumpRealCS(arrays, arrayNames, nx, ny, xStep, yStep, fileName+"H");

}


void gfmdSheet::dumpConfig() {

  string fileName = konfigName + "dat";

  if ( (!fTopoAdd) && (!nElast) ) { // dump cross section of input surface
    vector<double*> arrays = {equilPos};
    vector<string> arrayNames = {"Z_eq"};
    dumpRealCS(arrays, arrayNames, nx, ny, 1, 1, fileName+"H");
    return;
  }

  cerr << "   gfmdSheet" << ID << "::dumpConfig()\n";//FLOW

  // fields to be dumped
  vector<double*> arrays;
  vector<string> arrayNames;
  if (fSheetType&1) {
    arrays.push_back(equilPos);
    arrayNames.push_back("Z_eq");
  }
  if (fSheetType&2) {
    arrays.push_back(dispZR);
    arrayNames.push_back("dispZ");
    arrays.push_back(stressZR);
    arrayNames.push_back("stressZ");
  }
  if (nElast>=2) {
    arrays.push_back(dispYR);
    arrayNames.push_back("dispY");

    arrays.push_back(stressYR);
    arrayNames.push_back("stressY");
  }
  if (nElast==3) {
    arrays.push_back(dispXR);
    arrayNames.push_back("dispX");

    arrays.push_back(stressXR);
    arrayNames.push_back("stressX");
  }

  if (fSheetType==1) {
    dumpReal(arrays, nx, ny, 1, 1, fileName);
    dumpRealCS(arrays, arrayNames, nx, ny, 1, 1, fileName+"H");
    return;
  }

  updateStress();

  dumpReal(arrays, nx, ny, 1, 1, fileName);
  dumpRealCS(arrays, arrayNames, nx, ny, 1, 1, fileName+"H");

  dumpMaxwell();//NEW-MW

}

void gfmdSheet::dumpMaxwell() {//NEW-MW
  if (!fMaxwell) return;

  string fileName = "konfig" + to_string(ID) + "Mw.dat";
  dumpFour(dispZMw, nx, ny, fileName);
}


void gfmdSheet::updateStress() {
  zeroStress();

  if (fKelvinVoigt) stressKV();
  else stressGFMD();

  // external stress as part of elastic stress
  if (!fConstCOM) {
    double theta = pressTheta*2*PI/360;
    double phi = pressPhi*2*PI/360;
    double sinTheta = sin(theta);
    double sinPhi = sin(phi);
    double cosTheta = cos(theta);
    double cosPhi = cos(phi);

    stressZF[0] += pressure*cosTheta * nxny;
    vGravit -= pressure * dispZF[0].real() * areaPP;
    if (nElast>=2) {
      stressYF[0] += pressure*sinTheta*sinPhi * nxny;
      vGravit -= pressure * dispYF[0].real() * areaPP;
    }
    if (nElast==3) {
      stressXF[0] += pressure*sinTheta*cosPhi * nxny;
      vGravit -= pressure * dispXF[0].real() * areaPP;
    }
  }

  // Fourier based filtering of stress field would go here

  stressF2R();
}


void gfmdSheet::zeroStress() {
  // memset better? 
  if (!nElast) return;

  for (int k = 0; k < nFour; ++k) stressZF[k] = 0;
  if (fRealSpaceStress) for (int k = 0; k < nReal; ++k) stressZR[k] = 0;

  if (nElast>=2) {
    for (int k = 0; k < nFour; ++k) stressYF[k] = 0;
    if (fRealSpaceStress) for (int k = 0; k < nReal; ++k) stressYR[k] = 0;
  }

  if (nElast==3) {
    for (int k = 0; k < nFour; ++k) stressXF[k] = 0;
    if (fRealSpaceStress) for (int k = 0; k < nReal; ++k) stressXR[k] = 0;
  }
}


double gfmdSheet::addStress() {

  if (!nElast) return(0.);
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::addStress()\n";//FLOW

  vOnSitePot = 0.;

  // complete real-space stress with on-site potential
  if (fOnSitePotential&1) vOnSitePot += deltaIndenter();
  if (fOnSitePotential&2) vOnSitePot += einsteinSolid();
  if (fOnSitePotential==4) vOnSitePot += OnSiteHertzPotential();//NEW-OnHertz
  if (fOnSitePotential==8) OnSiteHertzConstraint();//NEW-OnHertz

  // initialize Fourier stress with FFTW from real-space stress
  for (Lint k=0; k<nFour; ++k) stressZF[k] = 0;
  if (fRealSpaceStress) {
    fftw_execute(stressZR2F);

    if (nElast>=2) {
      fftw_execute(stressYR2F);
      stressInterY = stressYF[0].real()/nxny;
    }

    if (nElast==3) {
      fftw_execute(stressXR2F);
      stressInterX = stressXF[0].real()/nxny;
    }
  }
  stressInterZ = stressZF[0].real()/nxny;

  // add elastic and thermal stresses
  vElastic = stressGFMD(); 
  if (fzOpposite) vElastic += stressOpposite();
  if (fMaxwell) stressMaxwell();//NEW-MW
  if ( fLangevin && (temp!=0) ) thermostat();
  //TEMP-stressRelax: output stresses of individual modes
  //cout << mdTime << "\t" << stressZF[1].real() << "\t" << abs(stressZF[kFastest]) << "\n";

  // add external pressure
  vGravit = 0;
  if (!fConstCOM) {
    double theta = pressTheta*2*PI/360;
    double phi = pressPhi*2*PI/360;
    double sinTheta = sin(theta);
    double sinPhi = sin(phi);
    double cosTheta = cos(theta);
    double cosPhi = cos(phi);

    stressZF[0] += pressure*cosTheta * nxny;
    vGravit -= pressure * dispZF[0].real() * areaPP;
    if (nElast>=2) {
      stressYF[0] += pressure*sinTheta*sinPhi * nxny;
      vGravit -= pressure * dispYF[0].real() * areaPP;
    }
    if (nElast==3) {
      stressXF[0] += pressure*sinTheta*cosPhi * nxny;
      vGravit -= pressure * dispXF[0].real() * areaPP;
    }
  }
  stressCOM = stressZF[0].real()/nxny;

  if (fSteadySlide) chi2SteadySlide = computeChi2SS();
  if (iTime==0) {
    chi2SteadySlideOld = chi2SteadySlide;
    chi2SSRestart = 1e+50*chi2SteadySlide;
    if (chi2SteadySlide<1e-10) chi2SSRestart = 1e+50*lengthX*lengthY;
    weightSS = 1e-5;
    invStiffSS_COMRef = 1./stressZZPre[1];
    invStiffSS_COM = invStiffSS_COMRef;
    oldSSPress = 1;
    prefacSSnxny = pow(nxny,0.9);
  }

  vPotTot = vOnSitePot+vElastic+vGravit;

  return(vPotTot);
}

//NEW-MW
void gfmdSheet::stressMaxwell(){
  if (iTime <= 2) cout << "   gfmdSheet" << ID << ".stressMaxwell()\n";//FLOW

  // add Maxwell stress
  for (Lint k = 0; k < nFour; ++k) {
    Complex stressMw = 0;
    for (int iMw = 0; iMw < nMaxwell; ++iMw) {
      stressMw += stiffFacMw[iMw]*stressZZPre[k]*dispZMw[iMw][k];
    }
    stressZF[k] += stressMw;//TEMP isolate conventional GFMD element by commenting this out
  }

  // update Maxwell elements via improved Euler's method (a.k.a. Heun's method)
  for (Lint k = 1; k < nFour; ++k) {
    for (int iMw = 0; iMw < nMaxwell; ++iMw) {
      Complex slopeNow = invTauMw[iMw]*(dispZF[k] - dispZMw[iMw][k]);
      Complex dispFnew = dispZMw[iMw][k] + slopeNow*dTime;
      Complex slopeNew = invTauMw[iMw]*(dispZF[k] - dispFnew);
      dispZMw[iMw][k] += 0.5*(slopeNow + slopeNew) * dTime;
    }
  }
  // treat 0-mode separately
  if (!fConstCOM || (fConstCOM&&fzOpposite)) {
    for (int iMw = 0; iMw < nMaxwell; ++iMw) {
      Complex slopeNow = invTauMw[iMw]*(dispZF[0]-zOpposite - dispZMw[iMw][0]);
      Complex dispFnew = dispZMw[iMw][0] + slopeNow*dTime;
      Complex slopeNew = invTauMw[iMw]*(dispZF[0]-zOpposite - dispFnew);
      dispZMw[iMw][0] += 0.5*(slopeNow + slopeNew) * dTime;
    }
  }
}

void gfmdSheet::fireHalt() {
  if (!nElast) return;
//memcpy(dispZFcopy, dispZFold, nFour*sizeof(Complex));
  memcpy(dispZF, dispZFold,  nFour*sizeof(Complex));

  //CM-FireXY
  if (nElast>=2) memcpy(dispYF, dispYFold,  nFour*sizeof(Complex));
  if (nElast==3) memcpy(dispXF, dispXFold,  nFour*sizeof(Complex));

  dispF2R();
}


void gfmdSheet::fireRedirect() {

  if ( (!nElast) || (fireRedrct==0) ) return;

  Complex vDotV = 0., vDotF = 0., fDotF = 0;
  for (Lint k = 0; k < nFour; ++ k) {
    Complex velocity = dispZF[k] - dispZFold[k];
    vDotV += weightQ[k] * 1. * velocity   * conj(velocity);
    vDotF += weightQ[k] * 1. * velocity   * conj(stressZF[k]);
    fDotF += weightQ[k] * 1. * stressZF[k] * conj(stressZF[k]);
  }

  double delta_PH = (1. - fireRedrct) * vDotF.real() / fDotF.real();
  double delta_Q  = fireRedrct*(fireRedrct - 2) * vDotV.real() / fDotF.real();

  double scaleV = (1. - fireRedrct);
  double scaleF = -delta_PH + sqrt(delta_PH*delta_PH - delta_Q);

  for (Lint k = 0; k < nFour; ++ k) {
    dispZF[k]  = dispZFold[k] + scaleV * (dispZF[k] - dispZFold[k]);
    dispZF[k] += scaleF * stressZF[k];
  }

  //CM-FireXY
  if (nElast>=2) {
    vDotV = vDotF = fDotF = 0;
    for (Lint k = 0; k < nFour; ++ k) {
      Complex velocity = dispYF[k] - dispYFold[k];
      vDotV += weightQ[k] * 1. * velocity   * conj(velocity);
      vDotF += weightQ[k] * 1. * velocity   * conj(stressYF[k]);
      fDotF += weightQ[k] * 1. * stressYF[k] * conj(stressYF[k]);
    }

    delta_PH = (1. - fireRedrct) * vDotF.real() / fDotF.real();
    delta_Q  = fireRedrct*(fireRedrct - 2) * vDotV.real() / fDotF.real();

    scaleV = (1. - fireRedrct);
    scaleF = -delta_PH + sqrt(delta_PH*delta_PH - delta_Q);

    for (Lint k = 0; k < nFour; ++ k) {
      dispYF[k]  = dispYFold[k] + scaleV * (dispYF[k] - dispYFold[k]);
      dispYF[k] += scaleF * stressYF[k];
    }
  }

  if (nElast==3) {
    vDotV = vDotF = fDotF = 0;
    for (Lint k = 0; k < nFour; ++ k) {
      Complex velocity = dispXF[k] - dispXFold[k];
      vDotV += weightQ[k] * 1. * velocity   * conj(velocity);
      vDotF += weightQ[k] * 1. * velocity   * conj(stressXF[k]);
      fDotF += weightQ[k] * 1. * stressXF[k] * conj(stressXF[k]);
    }

    delta_PH = (1. - fireRedrct) * vDotF.real() / fDotF.real();
    delta_Q  = fireRedrct*(fireRedrct - 2) * vDotV.real() / fDotF.real();

    scaleV = (1. - fireRedrct);
    scaleF = -delta_PH + sqrt(delta_PH*delta_PH - delta_Q);

    for (Lint k = 0; k < nFour; ++ k) {
      dispXF[k]  = dispXFold[k] + scaleV * (dispXF[k] - dispXFold[k]);
      dispXF[k] += scaleF * stressXF[k];
    }
  }

  if (fRealSpaceStress) dispF2R();

}


double gfmdSheet::deltaIndenter() {
  //const int iPoint = 1*rCentral;
  Lint ixIndenter = round( (double) (xOnSite+lengthX/2)/dx );
  ixIndenter = ixIndenter%nx;
  Lint iyIndenter = round( (double) (yOnSite+lengthY/2)/dy );
  iyIndenter = iyIndenter%ny;
  const Lint iPoint = getLinRIndex(ixIndenter,iyIndenter);

  double stress = pressure * nxny;
  if (fOnSitePeriod!=0) stress *= sin(fOnSiteFreq*mdTime);
  double extWork;

  double theta = pressTheta*2*PI/360;
  double phi = pressPhi*2*PI/360;
  double sinTheta = sin(theta);
  double sinPhi = sin(phi);
  double cosTheta = cos(theta);
  double cosPhi = cos(phi);
  if (nElast>=2) {
    stressYR[iPoint] -= stress*sinTheta*sinPhi;
    extWork += stressYR[iPoint]*dispYR[iPoint]*areaPP;
  }
  if (nElast==3) {
    stressXR[iPoint] -= stress*sinTheta*cosPhi;
    extWork += stressXR[iPoint]*dispXR[iPoint]*areaPP;
  }
  stressZR[iPoint] -= stress*cosTheta;
  extWork += stressZR[iPoint]*dispZR[iPoint]*areaPP;

  return (extWork);
}


void gfmdSheet::OnSiteHertzConstraint() {
  const Lint iPoint = 1*rCentral;
  get_irxiry(iPoint);
  double ixCentral = irx + xOnSite/dx;
  double iyCentral = iry + yOnSite/dy;
  
  /*// apply lateral velocity to indenter
  double shiftX = vXOnSite*mdTime;
  double shiftY = vYOnSite*mdTime;
  // make periodic
  while (shiftX > lengthX - ixCentral*dx) shiftX -= lengthX;
  while (shiftX < -lengthX + ixCentral*dx) shiftX += lengthX;
  while (shiftY > lengthY - iyCentral*dy) shiftY -= lengthY;
  while (shiftY < -lengthY + iyCentral*dy) shiftY += lengthY;
  ixCentral -= shiftX/dx;
  iyCentral -= shiftY/dy;//*/

  double lengthXH = lengthX/2;
  double lengthYH = lengthY/2;
  
  // constraint
  for (Lint k = 0; k<nxny; ++k) {
    get_irxiry(k);
    double deltaX = (irx-ixCentral)*dx; //deltaX2 *= deltaX2;
    double deltaY = (iry-iyCentral)*dy; //deltaY2 *= deltaY2;
    
    // make periodic
    if (deltaX > lengthXH) deltaX = (irx-ixCentral-nx)*dx;
    else if (deltaX < -lengthXH) deltaX = (irx-ixCentral+nx)*dx;
    if (deltaY > lengthYH) deltaY = (iry-iyCentral-ny)*dy;
    else if (deltaY < -lengthYH) deltaY = (iry-iyCentral+ny)*dy;

    double surface = deltaX*deltaX/rXhertz + deltaY*deltaY/rYhertz;
    surface /= 2;
    if (dispZR[k] > surface) dispZR[k] = surface;
  }
  dispR2F();
}


double gfmdSheet::OnSiteHertzPotential() {
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::OnSiteHertzPotential()\n";//FLOW

  // interaction parameters
  double potCurve = potCurveRelOS*stiffMax;
  double freq = 0, potRange = 0;
  if (surfEnergOS > 0) {
    freq = sqrt(2*potCurve/surfEnergOS); // potParam[0]
    potRange = PI / freq; // potParam[1]
  }

  // Hertzian parameters
  const Lint iPoint = 1*rCentral;
  get_irxiry(iPoint);
  double ixCentral = irx + xOnSite/dx;
  double iyCentral = iry + yOnSite/dy;

  // apply lateral velocity to indenter
  /*double shiftX = vXOnSite*dTime * ((int) iTime/frictRelax);
  double shiftY = vYOnSite*dTime * ((int) iTime/frictRelax);
  while (shiftX > lengthX - ixCentral*dx) shiftX -= lengthX;
  while (shiftX < -lengthX + ixCentral*dx) shiftX += lengthX;
  while (shiftY > lengthY - iyCentral*dy) shiftY -= lengthY;
  while (shiftY < -lengthY + iyCentral*dy) shiftY += lengthY;
  ixCentral += shiftX/dx;
  iyCentral += shiftY/dy;//*/

  double lengthXH = lengthX/2;
  double lengthYH = lengthY/2;
  double energLocal = 0;
  stressXOnSite = stressYOnSite = 0;

  // to ensure stability, the lateral component is switched on gradually.
  double frictionNow = frictionCoeffOS;
  /*if (iTime < frictRelax/2) frictionNow = 0;
  else if (iTime < frictRelax) frictionNow *= 0.5*(1+cos(iTime*2*PI/frictRelax));//*/

  for (Lint k = 0; k<nxny; ++k) {
    get_irxiry(k);
    double deltaX = (irx-ixCentral)*dx;
    double deltaY = (iry-iyCentral)*dy;
    
    // make periodic
    if (deltaX > lengthXH) deltaX = (irx-ixCentral-nx)*dx;
    else if (deltaX < -lengthXH) deltaX = (irx-ixCentral+nx)*dx;
    if (deltaY > lengthYH) deltaY = (iry-iyCentral-ny)*dy;
    else if (deltaY < -lengthYH) deltaY = (iry-iyCentral+ny)*dy;

    double surface = deltaX*deltaX/rXhertz + deltaY*deltaY/rYhertz;
    surface /= 2;
    double gapZ = surface - dispZR[k];
    
    // Z potential2
    double stressZ;
    if (gapZ<0) { 
      stressZ = potCurve*gapZ; 
      energLocal += -surfEnergOS + stressZ*gapZ/2;
    }
    else if (gapZ>potRange) {
      stressZ = 0; 
      energLocal += 0;
    }
    else { 
      double zEff = freq*gapZ;
      energLocal  += -surfEnergOS * (1+cos(zEff)) / 2;
      stressZ =  surfEnergOS * freq * sin(zEff) / 2;
    }//*/

    stressZR[k] += stressZ;

    // true X and Y coupling
    double velXnow = 0, stressX = 0;
    double velYnow = 0, stressY = 0;
    if (nElast>=2) {
      stressY = -stressZ*deltaY/rYhertz;
      stressYOnSite += stressY;
      
      Lint kR = (irx+1)%nx*ny + iry;    // right neighbor with pbc
      Lint kL = (irx+nx-1)%nx*ny + iry; // left neighbor with pbc
      Lint kT = irx*ny + (iry+1)%ny;    // top neighbor with pbc
      Lint kB = irx*ny + (iry+ny-1)%ny; // bottom neighbor with pbc
      double duy_dx = (dispYR[kR] - dispYR[kL])/(2*dx);
      double duy_dy = (dispYR[kT] - dispYR[kB])/(2*dy);
      velYnow = -duy_dx*vXOnSite - duy_dy*vYOnSite;

      velYnow -= vYOnSite;

      // friction in X and Y direction
      if (nElast==3) {
        stressX = -stressZ*deltaX/rXhertz;
        stressXOnSite += stressX;

        double dux_dx = (dispXR[kR] - dispXR[kL])/(2*dx);
        double dux_dy = (dispXR[kT] - dispXR[kB])/(2*dy);
        velXnow = -dux_dx*vXOnSite - dux_dy*vYOnSite;

        velXnow -= vXOnSite;
        double velLateral = velXnow*velXnow + velYnow*velYnow;
        velLateral = sqrt(velLateral);
        if (velLateral==0) velLateral = 1;// avoid division by 0.
        if (gapZ < 0) { // no anti-friction allowed.
          stressY = frictionNow*stressZ*velYnow/velLateral;
          stressX = frictionNow*stressZ*velXnow/velLateral;
          stressYR[k] += stressY;
          stressXR[k] += stressX;
        }
      }

      // friction in Y direction only
      else {
        int velSign = (velYnow>0) - (velYnow<0);
        if (gapZ < 0) { // no anti-friction allowed.
          stressY = velSign*frictionNow*stressZ;
          stressYR[k] += stressY;
        }
      }
    }
    // measure geometric X and Y stress components
    else {
      stressY = -stressZ*deltaY/rYhertz;
      stressYOnSite += stressY;
      stressX = -stressZ*deltaX/rXhertz;
      stressXOnSite += stressX;
    }
  }
  stressYOnSite /= nxny;
  stressXOnSite /= nxny;

  energLocal *= areaPP;
  return(energLocal);

} 


double gfmdSheet::einsteinSolid() {
  double vPot = 0;
  double zEq = 1, kEinstein = sqrt(stiffMax*stiffMin);
  if (iTime==0) cerr << "# kEinstein = " << kEinstein
		                 << " (stress)\t" << kEinstein*areaPP << " (per CGA)" << endl;
  for (Lint k = 0; k < nxny; ++k) {
    double dz = dispZR[k]-zEq;
    stressZR[k] -= kEinstein*dz;
    vPot += dz*dz;
  } 
  vPot *= kEinstein*areaPP/2;
  return(vPot); 
}

double gfmdSheet::stressOpposite() {
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::stressOpposite()\n";//FLOW
  double prefac = (fKelvinVoigt)? scalKV : 1.;
  Complex stressLoc = prefac*stressZZPre[0] * (dispZF[0] - zOpposite);
  stressZF[0] -= stressLoc;
  Complex dummy = stressLoc * conj(dispZF[0] - zOpposite);
  return weightQ[0]*dummy.real();
}

double gfmdSheet::stressGFMD() {
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::stressGFMD()\n";//FLOW
  
  if (fKelvinVoigt) return(0);//CM-FixCOM
  
  double vElaLoc = 0;
  if (nElast>=2) {
    if (nElast==3) vElaLoc += stressGFMDzyx(); 
    else vElaLoc += stressGFMDzy();
  }
  else {
    for (Lint k = 1; k < nFour; ++k) {
      Complex stressLoc = stressZZPre[k] * dispZF[k];
      stressZF[k] -= stressLoc;
      Complex dummy = stressLoc * conj(dispZF[k]);
      vElaLoc += weightQ[k]*dummy.real();
    }
  }

  vElaLoc *= areaPP / 2;
  vElaLoc /= nxny;

  return(vElaLoc);
}


double gfmdSheet::stressGFMDzy() {
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::stressGFMDzy()\n";//FLOW
  double vElaLoc = 0.;

  for (Lint k = 1; k < nFour; ++k) { 
    Complex dummy = 0.;
    
    // stiffness components
    double stiffYY = stressYYPre[k];
    Complex stiffYZ = Complex(0, ImStressYZPre[k]);
    double stiffZZ = stressZZPre[k];

    // YZ
    Complex stressLoc = stiffYZ * dispZF[k];
    stressYF[k] -= stressLoc;
    dummy += stressLoc * conj(dispYF[k]);

    // ZY
    stressLoc = conj(stiffYZ) * dispYF[k];
    stressZF[k] -= stressLoc;
    dummy += stressLoc * conj(dispZF[k]);

    // YY
    stressLoc = stiffYY * dispYF[k];
    stressYF[k] -= stressLoc;
    dummy += stressLoc * conj(dispYF[k]);

    // ZZ
    stressLoc = stiffZZ * dispZF[k];
    stressZF[k] -= stressLoc;
    dummy += stressLoc * conj(dispZF[k]);

    vElaLoc += weightQ[k]*dummy.real();
  }

  return(vElaLoc);
}


double gfmdSheet::stressGFMDzyx() {
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::stressGFMDzyx()\n";//FLOW
  double vElaLoc = 0.;

  // treat 0-mode in x- and y-direction separately
  if (!fConstCOMx) {
    Complex stressLoc = stressXXPre[0] * dispXF[0];
    stressXF[0] -= stressLoc;
    vElaLoc += (stressLoc * conj(dispXF[0])).real();
  }
  if (!fConstCOMy) {
    Complex stressLoc = stressYYPre[0] * dispYF[0];
    stressYF[0] -= stressLoc;
    vElaLoc += (stressLoc * conj(dispYF[0])).real();
  }

  for (Lint k = 1; k < nFour; ++k) { 
    Complex dummy = 0., stressLoc = 0.;

    double q = getQ(k);
    get_iqxiqy(k);
    Lint jqx = abs(iqx-nx);
    iqx = (iqx<jqx) ? iqx : -jqx;
    double q_x = iqx*dqx, q_y = iqy*dqy;
    double cosFac = q_x/q, sinFac = q_y/q;
    
    /*// rotate stiffness tensor Phi to local q-space system
    double stiffXX = stressXXPre[k] - 2*stressXYPre[k]*sinFac*cosFac;
    double stiffXY = stressXYPre[k]*(2*cosFac*cosFac - 1);
    double stiffYY = stressXXPre[k] + 2*stressXYPre[k]*sinFac*cosFac;
    Complex stiffXZ = Complex(0, ImStressXZPre[k]*(cosFac - sinFac));
    Complex stiffYZ = Complex(0, ImStressXZPre[k]*(cosFac + sinFac));
    double stiffZZ = stressZZPre[k];//*/
    
    // no rotation ( rotation already happens in prepareGFMD() )
    double stiffXX = stressXXPre[k];
    double stiffXY = stressXYPre[k];
    double stiffYY = stressYYPre[k];
    Complex stiffXZ = Complex(0, ImStressXZPre[k]);
    Complex stiffYZ = Complex(0, ImStressYZPre[k]);
    double stiffZZ = stressZZPre[k];//*/

    // XX
    stressLoc = stiffXX * dispXF[k];
    stressXF[k] -= stressLoc;
    dummy += stressLoc * conj(dispXF[k]);

    // XY
    stressLoc = stiffXY * dispYF[k];
    stressXF[k] -= stressLoc;
    dummy += stressLoc * conj(dispXF[k]);

    // XZ
    stressLoc = stiffXZ * dispZF[k];
    stressXF[k] -= stressLoc;
    dummy += stressLoc * conj(dispXF[k]);

    // YX
    stressLoc = stiffXY * dispXF[k];
    stressYF[k] -= stressLoc;
    dummy += stressLoc * conj(dispYF[k]);

    // YY
    stressLoc = stiffYY * dispYF[k];
    stressYF[k] -= stressLoc;
    dummy += stressLoc * conj(dispYF[k]);

    // YZ
    stressLoc = stiffYZ * dispZF[k];
    stressYF[k] -= stressLoc;
    dummy += stressLoc * conj(dispYF[k]);

    // ZX
    stressLoc = conj(stiffXZ) * dispXF[k];
    stressZF[k] -= stressLoc;
    dummy += stressLoc * conj(dispZF[k]);

    // ZY
    stressLoc = conj(stiffYZ) * dispYF[k];
    stressZF[k] -= stressLoc;
    dummy += stressLoc * conj(dispZF[k]);

    // ZZ
    stressLoc = stiffZZ * dispZF[k];
    stressZF[k] -= stressLoc;
    dummy += stressLoc * conj(dispZF[k]);

    // dummy adds up to be real since stiffness tensor is Hermitian.
    vElaLoc += weightQ[k]*dummy.real();
  }

  return(vElaLoc);
}


double gfmdSheet::stressKV() {

  double vElaLoc = 0;

  for (Lint k = 1; k != nFour; ++k) {
    stressZF[k] = -stressZZPre[k]*dispZF[k];
    Complex dummy = -stressZF[k]*conj(dispZF[k]);
    vElaLoc += weightQ[k]*dummy.real();
  }
  vElaLoc *= areaPP/2;
  vElaLoc /= nxny;
  return(vElaLoc);
    
}


void gfmdSheet::thermostat() {

  double preFac = sqrt(1.5*nxny*damping*temp*0.5/areaPP) / dTime; 
  for (Lint k = 0; k != nFour; ++k) {
    double randZR = 2*mix64() - 1;
    double randZI = 2*mix64() - 1;
    Complex randZC = Complex(randZR, randZI);
    stressZF[k] += preFac*randZC / sqrt(2*invMass[k]);
    get_iqxiqy(k);
  }
}


void gfmdSheet::initMeasure() {
  if (!nElast) return;

  tKinetic1 = tKinetic2 = 0;

  string num = to_string(ny);
  while (num.length() < 4) num = "0" + num;
  string moniName = "moni" + to_string(ID) + "-" + num + ".dat";
  string rampName = "ramp" + to_string(ID) + "-" + num + ".dat";
  string frictName = "frict" + to_string(ID) + "-" + num + ".dat";
  string maxwellName = "maxwell" + to_string(ID) + "-" + num + ".dat";

  moni.open(moniName, ofstream::out);
  if (fSteppedRamp) {
    ramp.open(rampName, ofstream::out);
    ramp << "# displacement  stress  time\n";
  }
  if (fFriction) {
    frictout.open(frictName, ofstream::out);
    frictout << "# time  stressZ";
    if (nElast>=2) frictout << "  stressY";
    if (nElast==3) frictout << "  stressX";
    if (fOnSitePotential==4) {
      frictout << "  geomStressY  geomStressX";
    }
    frictout << "\n";
  }
  if (fMaxwell && fDumpMaxwell) {
    maxwellout.open(maxwellName);
    maxwellout << "# time";
    for (int iMw = 0; iMw < nMaxwell; ++iMw) {
      maxwellout << "  dispZMw[" << iMw << "][0]";
      maxwellout << "  dispZMw[" << iMw << "][1]";
      maxwellout << "  dispZMw[" << iMw << "][kFastest]";
    }
    maxwellout << "\n";
  }
  moni << "# mdTime"
  << "\tdispZF[0].real()"
  << "  abs(dispZF[1])"
  << "  abs(dispZF[kFastest])"
  << "  dispZR[0]"
  << "  dizpZR[central]";
  if (fConstCOM) moni << "  stressCOM"; 
  if ( (pressFinal!=pressInit) || (fSteppedRamp==2) ) moni << "  pressure";
  if (fzOpposite) moni << "  zOpposite";
  if (fExtrapolateInf) moni << "  zInf";
  if (nElast>=2) {
    moni << "\tdispYF[0].real()"
         << "  abs(dispYF[1])"
         << "  abs(dispYF[kFastest])"
         << "  dispYR[0]"
         << "  dizpYR[central]"
         << "  stressYF[0]";
  }
  if (nElast==3) {
    moni << "\tdispXF[0].real()"
         << "  abs(dispXF[1])"
         << "  abs(dispXF[kFastest])"
         << "  dispXR[0]"
         << "  dizpXR[central]"
         << "  stressXF[0]";
  }
  moni << endl;

  //DEBUG
  //cout << "# stressInterZ\tstressCOM\tdispZF[0]\tdispZFold[0]\tzOpposite\tstressZZPre[0]\n";
}


void gfmdSheet::measure() {
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::measure()\n";//FLOW
  if (!nElast) return;

  tKinetic1 += tKinetic;
  tKinetic2 += tKinetic*tKinetic;

  if (iTime) {
    moni << mdTime << "\t";
    moni << dispZF[0].real()/nxny
         << "\t" << dispZF[1].real()/nxny
         << "\t" << dispZF[kFastest].real()/nxny
         << "\t" << dispZR[0]
         << "\t" << dispZR[rCentral];
    if (fConstCOM) {
      if (fzOpposite) moni << "\t" << stressInterZ;
      else moni << "\t" << stressCOM;
    }
    if ( (pressFinal!=pressInit) || (fSteppedRamp==2) ) moni << "\t" << pressure;
    if (fzOpposite) moni << "\t" << zOpposite/nxny;
    if (fExtrapolateInf) {
      //if (fzOpposite) moni << "\t" << (zOpposite - dispZF[0].real())/nxny + dispInf();
      moni << "\t" << dispInf();
    }
    if (nElast>=2) {
      moni << "\t" << dispYF[0].real()/nxny
           << "\t" << abs(dispYF[1])/nxny
           << "\t" << abs(dispYF[kFastest])/nxny
           << "\t" << dispYR[0]
           << "\t" << dispYR[rCentral]
           << "\t" << stressYF[0].real()/nxny;
      if (nElast==3) {
        moni << "\t" << dispXF[0].real()/nxny
             << "\t" << abs(dispXF[1])/nxny
             << "\t" << abs(dispXF[kFastest])/nxny
             << "\t" << dispXR[0]
             << "\t" << dispXR[rCentral]
             << "\t" << stressXF[0].real()/nxny;
      }
    }
    moni << endl;

    if (fFriction) {
      frictout << mdTime << "\t" << stressInterZ;
      if (nElast>=2) frictout << "\t" << stressYF[0].real()/nxny - stressInterY;
      if (nElast==3) frictout << "\t" << stressXF[0].real()/nxny - stressInterX;
      if (fOnSitePotential==4)
        frictout << "\t" << stressYOnSite << "\t" << stressXOnSite;
      frictout << endl;
    }

    if (fDumpMaxwell) {
      maxwellout << mdTime;
      for (int iMw = 0; iMw < nMaxwell; ++iMw) {
        maxwellout << "\t" << abs(dispZMw[iMw][0])/nxny;
        maxwellout << "\t" << abs(dispZMw[iMw][1])/nxny;
        maxwellout << "\t" << abs(dispZMw[iMw][kFastest])/nxny;
      }
      maxwellout << "\n";
    }
  }
  //cout << iTime << "\t" << stressZF[0].real() << "\n";//--DEBUG
  //cout << abs(stressZF[1298])/nxny << "\t" << abs(stressYF[1298])/nxny << "\t" << abs(stressXF[1298])/nxny << "\n";//--DEBUG

  //DEBUG
  //cout << stressInterZ << "\t" << stressCOM << "\t" << dispZF[0].real() << "\t"
  //     << dispZFold[0].real() << "\t" << zOpposite << "\t" << stressZZPre[0] << "\n";
}


double gfmdSheet::tKinHalfStep() {

  if (!nElast) return(0);

  // initialize prefactors
  double dt2 = dTime2;
  double weightNow = 2, weightOld = 1;

  tKinetic = 0;
  int kStart = 0;
  if (fConstCOM) kStart = 1;
  for (Lint k = kStart; k < nFour; ++k) {
    Complex vDiffN = dispZF[k] - dispZFold[k];
    vDiffN = vDiffN * conj(vDiffN);
    tKinetic += weightQ[k] * vDiffN.real() / invMass[k];
  }
  //CM-FireXY
  kStart = 0;
  if (fConstCOMy) kStart = 1;
  if (nElast>=2) {
    for (Lint k = kStart; k < nFour; ++k) {
      Complex vDiffN = dispYF[k] - dispYFold[k];
      vDiffN = vDiffN * conj(vDiffN);
      tKinetic += weightQ[k] * vDiffN.real() / invMass[k];
    }
  }
  kStart = 0;
  if (fConstCOMx) kStart = 1;
  if (nElast==3) {
    for (Lint k = kStart; k < nFour; ++k) {
      Complex vDiffN = dispXF[k] - dispXFold[k];
      vDiffN = vDiffN * conj(vDiffN);
      tKinetic += weightQ[k] * vDiffN.real() / invMass[k];
    }
  }
  tKinetic *= areaPP / (2*dt2*nxny);
  return(tKinetic);
}


double gfmdSheet::propagate() {
  if (iTime < 3) cerr << "   gfmdSheet" << ID << "::propagate()\n";//FLOW

  // MM2MM: need to add moveLateral on displacements (for nElast==1 sheets)
  if (fSliding&&!fSteadySlide) moveLateral(); 

  if (!nElast) return(0.);

  //switch off FIRE while being ramped
  if (fSteppedRamp) {
    rampTime = iTime % (rampSteps+rampRelax);
    if (fFire && (rampTime<rampSteps)) fFireOn = 0;
  }

  // update pressure
  if (fSteppedRamp != 2) {
    pressure = pressInit;
    if (iTime>nRelax) pressure += (iTime-nRelax)*(pressFinal-pressInit)/(nTime-nRelax);
  } 
  else if (iTime>nRelax) {
    if ( (rampTime<rampSteps) && (iTime>rampSteps) ) {
      double prefacRamp = sin(rampTime*PI/rampSteps);
      prefacRamp *= 2*prefacRamp;
      pressure += ddpRamp*prefacRamp;
    } 
    else if (rampTime+1==rampSteps+rampRelax) { // dump info at last relax step
      ramp << dispZF[0].real()/nxny << "\t";
      ramp << stressInterZ << "\t" << mdTime << endl;
    }
  }

  // viscoelasticity is treated by propagateKV instead
  if (fKelvinVoigt) {
    if (fSteadySlide) {
      propagateSteadySlide();
      return(0); 
    } 
    propagateKV();
    return(0); 
  }

  // initialize prefactors
  double dt2 = dTime2;
  double weightNow = 2, weightOld = 1;

  if (!fFireOn) {
    double dampDL = damping*dTime; // damping in units of dTime
    dt2 /= (1.0 + dampDL);
    double weightDiff = 1-exp(-dampDL);
    weightNow -= weightDiff;
    weightOld -= weightDiff;
  }
  
  // external propagation of COM/Opposite
  int kStart = 0;
  if (fConstCOM) {
    if (!fzOpposite) {
      kStart = 1;
      if (iTime > nRelax) moveCOM();
    }
    else if (iTime > nRelax) moveOpposite();
  }//TODO: is this equivalent and more readable?
  //if (kStart==0) cerr << invMass[0]*stressZF[0]*dt2 << "\n";//DEBUG

  // propagate internal q modes with Verlet
  tKinetic = 0;
  for (Lint k = kStart; k < nFour; ++k) {
    Complex dispFnew = weightNow*dispZF[k] - weightOld*dispZFold[k]
		     + invMass[k]*stressZF[k]*dt2;
    Complex vDiffN = dispFnew - dispZFold[k];
    vDiffN = vDiffN * conj(vDiffN);
    tKinetic += weightQ[k] * vDiffN.real() / invMass[k]; 
    dispZFold[k] = dispZF[k];
    dispZF[k] = dispFnew;//TEMP-stressRelax: comment this out to look at constant displacement
  }

  if (nElast>=2) {
    kStart = 0;
    if (fConstCOMy) kStart = 1;
    for (Lint k = kStart; k < nFour; ++k) {
      Complex dispFnew = weightNow*dispYF[k] - weightOld*dispYFold[k]
    		       + invMass[k]*stressYF[k]*dt2;
      Complex vDiffN = dispFnew - dispYFold[k];
      vDiffN = vDiffN * conj(vDiffN);
      tKinetic += weightQ[k] * vDiffN.real() / invMass[k]; 
      dispYFold[k] = dispYF[k];
      dispYF[k] = dispFnew;
    }
  }

  if (nElast==3) {
    kStart = 0;
    if (fConstCOMx) kStart = 1;
    for (Lint k = kStart; k < nFour; ++k) {
      Complex dispFnew = weightNow*dispXF[k] - weightOld*dispXFold[k]
               + invMass[k]*stressXF[k]*dt2;
      Complex vDiffN = dispFnew - dispXFold[k];
      vDiffN = vDiffN * conj(vDiffN);
      tKinetic += weightQ[k] * vDiffN.real() / invMass[k]; 
      dispXFold[k] = dispXF[k];
      dispXF[k] = dispFnew;
    }
  }

  tKinetic *= areaPP / (8*dt2*nxny);

  if (fFire&&(fireRedrct!=0)) fireRedirect();

  // turn around movement if necessary
  turnaroundCOM();

  dispF2R();
  return(tKinetic);
}


double gfmdSheet::dispInf() {
  return 6*dispZR[0] - 2.5*(dispZR[getLinRIndex(nx/2,0)]+ dispZR[getLinRIndex(0,ny/2)]);
}


void gfmdSheet::moveCOM() {
  if (fExtrapolateInf && iTime >= 1) {
    double zInfTarget, zInfOld = dispInf();
    double deltaOld = dispZF[0].real() - dispZFold[0].real();
    if (nVeloTurnStep > 0) {
      if (nVeloTransition > 0) {
        if (iTime > tEndTransition)
          zInfTarget = zConstCOM + 2*vzInit*t0Transition*dTime - vzInit*(iTime-1)*dTime;
        else if (iTime > tStartTransition) {
          zInfTarget = zConstCOM + vzInit*dTime*(tStartTransition+t0Transition)/2;
          zInfTarget -= vzInit*dTime/(2*nVeloTransition)*pow(iTime-1-t0Transition, 2);
        }
        else zInfTarget = zConstCOM + vzInit*(iTime-1)*dTime;
      }
      else if (iTime > nVeloTurnStep) 
        zInfTarget = zConstCOM + 2*vzInit*nVeloTurnStep*dTime - vzInit*(iTime-1)*dTime;
    }
    else zInfTarget = zConstCOM + vzInit*(iTime-1)*dTime;
    
    // add zInfError to dispZF[0] with LPF:
    const double weightLPFold = 0.3, weightLPFnow = 0.7;
    double zInfErrorNow = zInfTarget - zInfOld;
    dispZF[0] += nxny*(vzInit*dTime + weightLPFnow*zInfErrorNow + weightLPFold*zInfErrorOld);
    zInfErrorOld = zInfErrorNow;
  }
  else dispZF[0] += vzConstCOM*dTime*nxny;

  if (fSteppedRamp) {
    if ( (rampTime<rampSteps) && (iTime>rampSteps) ) {
      double prefacRamp = sin(rampTime*PI/rampSteps);
      prefacRamp *= 2*prefacRamp;
      dispZF[0] += ddzRamp*prefacRamp;
    } 
    else if (rampTime+1==rampSteps+rampRelax) {
      ramp << dispZF[0].real()/nxny << "\t";
      ramp << stressInterZ << "\t" << mdTime << endl;
    }
    dispZFold[0] = dispZF[0];
  }
}


void gfmdSheet::moveOpposite() {
  if (fExtrapolateInf && iTime >= 1) {
    double zInfTarget, zInfOld = dispInf();
    double deltaOld = dispZF[0].real() - dispZFold[0].real();
    double weightLPFold = 0.3, weightLPFnow = 0.7;
    if (nVeloTurnStep > 0) {
      if (nVeloTransition > 0) {
        if (iTime > tEndTransition) 
          zInfTarget = zConstCOM + 2*vzInit*t0Transition*dTime - vzInit*(iTime-1)*dTime;
        else if (iTime > tStartTransition) {
          zInfTarget = zConstCOM + vzInit*dTime*(tStartTransition+t0Transition)/2;
          zInfTarget -= vzInit*dTime/(2*nVeloTransition)*pow(iTime-1-t0Transition, 2);
        }
        else zInfTarget = zConstCOM + vzInit*(iTime-1)*dTime;
      }
      else if (iTime > nVeloTurnStep) 
        zInfTarget = zConstCOM + 2*vzInit*nVeloTurnStep*dTime - vzInit*(iTime-1)*dTime;
    }
    else zInfTarget = zConstCOM + vzInit*(iTime-1)*dTime;

    // add zInfError to dispZF[0] with LPF:
    double zInfErrorNow = zInfTarget - zInfOld;
    zOpposite += nxny*(vzInit*dTime + weightLPFnow*zInfErrorNow + weightLPFold*zInfErrorOld);
    zInfErrorOld = zInfErrorNow;
  }
  else zOpposite += vzConstCOM*dTime*nxny;

  if (fSteppedRamp) {
    if ( (rampTime<rampSteps) && (iTime>rampSteps) ) {
      double prefacRamp = sin(rampTime*PI/rampSteps);
      prefacRamp *= 2*prefacRamp;
      zOpposite += ddzRamp*prefacRamp;
    } 
    else if (rampTime+1==rampSteps+rampRelax) {
      ramp << zOpposite/nxny << "\t";
      ramp << stressInterZ << "\t" << mdTime << endl;
    }
  }
}

void gfmdSheet::turnaroundCOM() {
  if (checkTurnPress && (sgnTurnPress*stressInterZ > 0.995*absTurnPress)) {
    //cerr << "# preload has been reached.\n";//DEBUG
    absTurnPress = 1e+85;
    checkTurnPress = 0;
    if (fSteppedRamp) {
      ddzRamp *= -1;
      ddpRamp *= -1;
    }
    else {
      tStartTransition = iTime + 1;
      nVeloTurnStep = tStartTransition + nVeloTransition/2;
      tEndTransition = nVeloTurnStep + 3*nVeloTransition/2;
      t0Transition = nVeloTurnStep + nVeloTransition/2;
      /*//DEBUG
      cout << "# nVeloTransition = " << nVeloTransition << "\n";
      cout << "# tStartTransition = " << tStartTransition << "\n";
      cout << "# t0Transition = " << t0Transition << "\n";
      cout << "# tEndTransition = " << tEndTransition << "\n";//*/
    }
  }
  else if ((nVeloTransition>0) && (iTime>=tStartTransition) && (iTime<tEndTransition)) { 
    //if (iTime == t0Transition) t0Transition -= 1; // skip time step where vz would be 0
    vzConstCOM = vzInit*(t0Transition - iTime)/nVeloTransition;
  }
  else if (iTime==nVeloTurnStep) {
    vzConstCOM *= -1;
    ddzRamp *= -1;
    ddpRamp *= -1;
  }
}

double gfmdSheet::propagateKV() {

  // initialize prefactors
  double weightNew = 1./(1+rKV_LPF);
  double weightOld = 1-weightNew;
  double prefacDot = weightNew*tauKV/scalKV;

  // external propagation of COM/Opposite
  int kStart = 0;
  if (fConstCOM) {
    if (!fzOpposite) {
      kStart = 1;
      if (iTime > nRelax) moveCOM();
    }
    else if (iTime > nRelax) moveOpposite();
  }
  //if (iTime == 1) cout << "kStart=" << kStart << "\n";//DEBUG

  // propagate internal q modes with KV-SLS
  double dissipKV = 0.; // MM2fix: dissipKV not yet computed
  for (Lint k = kStart; k < nFour; ++k) {

    // NOTE: dispZFold is actually stressZFold here.
    Complex extStressDot = (stressZF[k]-dispZFold[k])/dTime; 
    Complex dispLocOld = dispZF[k];

    rhsKV_LPF[k] *= weightOld;
    rhsKV_LPF[k] += weightNew*stressZF[k];
    if (iTime>1) rhsKV_LPF[k] += prefacDot*extStressDot;

    //0 mode fix
    Complex effectiveStress;
    if (k==0) effectiveStress = rhsKV_LPF[k];//CM-TODO: is this necessary?
    else effectiveStress = rhsKV_LPF[k] - stressZZPre[k] * dispZF[k];
    Complex  velocity = effectiveStress / (stressZZPre[k]*tauKV);
    dispZF[k] += velocity * dTime;

    dispZFold[k] = stressZF[k];

  }

  // turn around movement if necessary
  turnaroundCOM();

  dispF2R();
  
  return(dissipKV);

}


double gfmdSheet::computeChi2SS() {
  Complex chi2 = 0;
  for (int k = 1; k < nFour; ++k) {
    stressZF[k] *= preFacCorrSS[k]/stressZZPre[k];
    Complex diff = stressZF[k] - dispZF[k];
    chi2 += diff*conj(diff);
  }
  return(chi2.real()); 
}


void gfmdSheet::propagateSteadySlide() {

  if (chi2SteadySlide >= chi2SteadySlideOld) {
    chi2SteadySlide = chi2SteadySlideOld;
    weightSS /= 4;
    if (weightSS<1e-10) {
      stressCOM = stressZF[0].real()/nxny;
      double pressDiff = pressInit + stressCOM;
      if (abs(pressDiff) < pressInit*1.e-5) {
        dispF2R();
        cerr << "# relative pressure error<1e-5 induces finish";
        ofstream finish("finish");
        finish.close();
      } 
      else {
        // adjusting the effective mass
        if (oldSSPress*pressDiff > 0) invStiffSS_COM*=pow(2.,1./4);
        else {
	  if (invStiffSS_COM > invStiffSS_COMRef) invStiffSS_COM = invStiffSS_COMRef;
	  invStiffSS_COM /= 2;
          oldSSPress *= -1;
        }
        dispZF[0] += invStiffSS_COM * pressDiff * prefacSSnxny;
        
        chi2SteadySlideOld = chi2SSRestart;
        weightSS = 1.e-5; // restart the weight
        for (int k = 1; k < nFour; ++k) dispZFold[k] = dispZF[k];
      }
    }
  } 
  else {
    chi2SteadySlideOld = chi2SteadySlide;
    for (int k = 1; k < nFour; ++k) {
      Complex newDispZF;
      newDispZF = weightSS*stressZF[k] + (1.-weightSS)*dispZFold[k];
      dispZFold[k] = dispZF[k];
      dispZF[k] = newDispZF;
    }
    weightSS *= pow(2.,1./4);
    if (weightSS>4) weightSS = 4;
  }
 
  dispF2R();
}


void gfmdSheet::stressF2R() {
  memcpy(fieldFFFT, stressZF, nFour*sizeof(Complex));
  fftw_plan stressZF2R  =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, stressZR, FFTW_ESTIMATE);
  fftw_execute(stressZF2R);
  fftw_destroy_plan(stressZF2R);
  for (Lint k = 0; k < nxny; ++k) stressZR[k] /= nxny;

  // same for Y
  if (nElast>=2) {
    if (!fConstCOMy) stressYF[0] += pressY * nxny;
    memcpy(fieldFFFT, stressYF, nFour*sizeof(Complex));
    fftw_plan stressYF2R  =
    fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, stressYR, FFTW_ESTIMATE);
    fftw_execute(stressYF2R);
    fftw_destroy_plan(stressYF2R);
    for (Lint k = 0; k < nxny; ++k) stressYR[k] /= nxny;
  }

  // same for X
  if (nElast==3) {
    if (!fConstCOMx) stressXF[0] += pressX * nxny;
    memcpy(fieldFFFT, stressXF, nFour*sizeof(Complex));
    fftw_plan stressXF2R  =
    fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, stressXR, FFTW_ESTIMATE);
    fftw_execute(stressXF2R);
    fftw_destroy_plan(stressXF2R);
    for (Lint k = 0; k < nxny; ++k) stressXR[k] /= nxny;
  }
}


void gfmdSheet::dispF2R() {
  // transform from Fourier to real space with back-up copy
  memcpy(fieldFFFT, dispZF, nFour*sizeof(Complex));
  fftw_execute(dispZF2R);
  for (Lint k = 0; k < nxny; ++k) dispZR[k] /= nxny;

  // same for Y direction
  if (nElast>=2) {
    memcpy(fieldFFFT, dispYF, nFour*sizeof(Complex));
    fftw_execute(dispYF2R);
    for (Lint k = 0; k < nxny; ++k) dispYR[k] /= nxny;
  }

  // same for X direction
  if (nElast==3) {
    memcpy(fieldFFFT, dispXF, nFour*sizeof(Complex));
    fftw_execute(dispXF2R);
    for (Lint k = 0; k < nxny; ++k) dispXR[k] /= nxny;
  }
}


void gfmdSheet::dispR2F() {
  fftw_execute(dispZR2F);
  if (nElast>=2) fftw_execute(dispYR2F);
  if (nElast==3) fftw_execute(dispXR2F);
  // MM2fix: if we ever start with finite velocities, then we need
  //		a real copy of dispZF at iTime==0. 
  if (iTime==0) {
    memcpy(dispZFold, dispZF, nFour*sizeof(Complex));
    if (nElast>=2) memcpy(dispYFold, dispYF, nFour*sizeof(Complex));
    if (nElast==3) memcpy(dispXFold, dispXF, nFour*sizeof(Complex));
  }
}


void gfmdSheet::outMeasure() {

  if (!nElast) return;

  ofstream output("params.out",ofstream::app);
  output << ID << "\t\t# sheet measure start\n";

  if (vElastic!=0) output << vElastic << "\t\t# vElastic\n"; 
  if (vOnSitePot!=0) output << vOnSitePot << "\t\t# vOnSitePot\n"; 
  if (vGravit!=0) output << vGravit << "\t\t# vGravit\n"; 
  if (areaXY!=1) {
    if (vElastic!=0) output << vElastic/areaXY << "\t\t# vElastic/areaXY\n"; 
    if (vOnSitePot!=0) output << vOnSitePot/areaXY << "\t\t# vOnSitePot/areaXY\n"; 
    if (vGravit!=0) output << vGravit/areaXY << "\t\t# vGravit/areaXY\n"; 
  }

  tKinetic1 /= (iTime-nRelax);
  tKinetic2 /= (iTime-nRelax);
  tKinetic2 -= tKinetic1*tKinetic1;
  if (fLangevin || fExternalDriving) {
    output << tKinetic1 << "\t" << tKinetic2 << "\t# tKinetic1, tKinetic2\n";
  }

  output << dispZF[0].real()/nxny << "\t\t# zConstCOM #\n";
  if (fzOpposite) output << zOpposite/nxny << "\t\t# zOpposite #\n";

  output << ID << "\t\t# sheet measure end\n\n";
  output.close();
}


void gfmdSheet::outSystem() {
  if (fOnSitePotential&8) OnSiteHertzConstraint();//NEW-OnHertz
  if (nElast||fLateral) dumpConfig();
}


gfmdSheet::~gfmdSheet() {

  if (nElast) {
    moni.close();
    if (fSteppedRamp) ramp.close();
    if (fFriction) frictout.close();
    if (fDumpMaxwell) maxwellout.close();
  }

  //cout << "ID=" << ID << ": Destroying topo/commodity fields: " << endl;//DEBUG
  fftw_free(fieldRFFT);
  fftw_free(fieldFFFT);
  fftw_free(equilPos);
  if (nElast) {
    //cout << "ID=" << ID << ": Destroying fields in Z direction: " << endl;//DEBUG
    fftw_free(dispZR);
    fftw_free(dispZF);
    fftw_free(dispZFold);
    fftw_free(stressZR);
    fftw_free(stressZF);
    fftw_free(weightQ);
    if (!fKelvinVoigt) fftw_free(invMass);
    fftw_free(stressZZPre);
    fftw_destroy_plan(dispZR2F);
    fftw_destroy_plan(dispZF2R);
    fftw_destroy_plan(stressZR2F);
  }
  if (nElast > 1) {
    //cout << "ID=" << ID << ": Destroying fields in Y direction: " << endl;//DEBUG
    fftw_destroy_plan(dispYR2F);
    fftw_destroy_plan(dispYF2R);
    fftw_destroy_plan(stressYR2F);
  }
  if (nElast==3) {
    //cout << "ID=" << ID << ": Destroying fields in X direction: " << endl;//DEBUG
    fftw_destroy_plan(dispXR2F);
    fftw_destroy_plan(dispXF2R);
    fftw_destroy_plan(stressXR2F);
  }
  if (fTopography && fSliding) {
    fftw_free(equilPos0F);
    fftw_destroy_plan(equilPosF2R);
  }

  //cout << "Destroyed sheet " << ID << endl;//DEBUG
}
