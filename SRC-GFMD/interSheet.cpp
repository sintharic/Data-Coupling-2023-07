#include "header.h"
#include "gfmdSheet.h"
#include "interSheet.h"
#include "globals.h"

// additional global variables
extern vector <gfmdSheet> sheet;
extern bool fFireConstraint;


void interSheet::initParams(int newID){

  ID = newID;

  initParamsDefault();
  initParamsFile();

  // ensure correct order of sheets
  if (sheetID0==sheetID1) termination("Do not couple sheet to itself!");
  else if (sheetID0>sheetID1) {
    // order sheets so that upper sheet has lower number (can be reversed later)
    int dummy = sheetID0;
    sheetID0 = sheetID1;
    sheetID1 = dummy;
    // MM2XX: What is fPotentialTest. Permanent or temporary?
    if (fPotential||fPotentialTest) {
      sheet[sheetID0].fRealSpaceStress = true;
      sheet[sheetID1].fRealSpaceStress = true;
    }
  }

  // sanity checks 
  if ( fElastic && (fConstraint||fPotential||fReynolds||fShifted) )
  termination("Do not use fElastic with other inter flags!");

  if ( fShifted && fReynolds )
  termination("Do not use fShifted with fReynolds."); 

  if (fConstraint) fFireConstraint = true;

  if (frictionCoeff != 0) fDumpLateral = 1;

} // initParams


void interSheet::initParamsFile() {
  ifstream readParams("params.in");
  if (readParams.is_open()) fConstraint = 0;    // reverse original default
  int fReadParams = 0;
  fFirstCurve = -1; // defines what parameter is used for curvature
  while (readParams.is_open()) {
    double param;
    std::string ROL; // rest of line
    std::size_t NIS = std::string::npos; // NIS == Not In String
    if (readParams.eof()) break;
    readParams >> param; getline(readParams,ROL);
    if ( (param==ID) && (ROL.find("# inter start") !=NIS) ) fReadParams = 1;
    if ( fReadParams==0 ) continue;
    if (ROL.find("# inter end")  !=NIS) break;
    else if (ROL.find("# sheetID0 #") !=NIS) sheetID0 = param;
    else if (ROL.find("# sheetID1 #") !=NIS) sheetID1 = param;
    else if (ROL.find("# fConstraint #")!=NIS) fConstraint = param;
    else if (ROL.find("# fPotential #") !=NIS) fPotential = param;
    else if (ROL.find("# fReynolds #")  !=NIS) fReynolds = param;
    else if (ROL.find("# fElastic #") !=NIS) fElastic = param;
    else if (ROL.find("# surfEnerg #")  !=NIS) surfEnerg = param;
    else if ( (ROL.find("# potCurveRel #")!=NIS) && (fFirstCurve==-1) )
      {potCurveRel = param; fFirstCurve = 1;}
    else if ( (ROL.find("# potRange #") !=NIS) && (fFirstCurve==-1) ) 
      {potRange = param; fFirstCurve = 3;}
    else if ( (ROL.find("# potCurve #") !=NIS) && (fFirstCurve==-1) )
      {potCurve = param; fFirstCurve = 2;}
    else if (ROL.find("# frictionCoeff #") !=NIS) frictionCoeff = param;
    else if (ROL.find("# relFrictPeak #") !=NIS) relFrictPeak = param;
    else if (ROL.find("# velFrictThresh #") !=NIS) velFrictThresh = param; 
    else if (ROL.find("# frictRelax #") !=NIS) frictRelax = param; 
    else if (ROL.find("# fPotentialTest #") !=NIS) fPotentialTest = param; 
    else if (ROL.find("# fDumpGap #") !=NIS) fDumpGap = param;
    else if (ROL.find("# resolMovie #") !=NIS) resolMovie = param;
    else if (ROL.find("# fDumpFrame #") !=NIS) fDumpFrame = param;
    else if (ROL.find("# fDumpLateral #") !=NIS) fDumpLateral = param;
    else if (ROL.find("# nTimeOn #") !=NIS) nTimeOn = param;
    else if (ROL.find("# nTimeOff #") !=NIS) nTimeOff = param;
  }
  readParams.close();
}


void interSheet::initSystem(){

  // further sanity checks needing information of sheets

  if (fConstraint) {
    if ( sheet[sheetID0].nElast && sheet[sheetID0].fMassWeightg ) 
    termination("fMassWeightg in sheetID0 needs to be 0 when using fConstraint.");
    if ( sheet[sheetID1].nElast && sheet[sheetID1].fMassWeightg ) 
    termination("fMassWeightg in sheetID1 needs to be 0 when using fConstraint.");
  }

  // make sure elastic sheets fft stressZR2F
  if (fPotential) {
    if (sheet[sheetID0].nElast) sheet[sheetID0].fRealSpaceStress = 1;
    if (sheet[sheetID1].nElast) sheet[sheetID1].fRealSpaceStress = 1;
  }

  // determine size and allocate fields 
  nx = sheet[sheetID0].nx;
  ny = sheet[sheetID0].ny;
  nxH = nx/2; nyHP1 = (ny/2) + 1;
  nFour = nx*nyHP1;
  nReal = 2*nFour;
  nxny = nx*ny;

  if (fElastic) { // elastic coupling between two sheets via Green's functions
    nx = min(nx, sheet[sheetID1].nx);
    ny = min(ny, sheet[sheetID1].ny);
    nGap = 1;
  } 
  else if (fShifted) {
    termination("fShifted not yet implemented");
  } 
  else { // regular contact mechanics
    if ( (nx!=sheet[sheetID1].nx) || (ny!=sheet[sheetID1].ny) ) {
      termination("Different resolutions for real-space interaction not implemented");
    }
    nGap = nx * ny;

    // MM2MM: change:  should only be called in sliding simulations 
    gap	= (double *) fftw_malloc( nGap*sizeof(double) );
    gapStress = (double *) fftw_malloc( nGap*sizeof(double) );

    // define fInterType
    if ((sheet[sheetID0].fSheetType==1) && (sheet[sheetID1].fSheetType==2)) {
      // top: rough, bottom: elastic
      fInterType = 1;
    } 
    else
    if ((sheet[sheetID0].fSheetType==2) && (sheet[sheetID1].fSheetType==1)) {
      // top: elastic, bottom: rough
      fInterType = 2;
    } 
    else
    if ((sheet[sheetID0].fSheetType==1) && (sheet[sheetID1].fSheetType==3)) {
      // both: rough, bottom: elastic
      fInterType = 3;
    } 
    else
    if ((sheet[sheetID0].fSheetType==3) && (sheet[sheetID1].fSheetType==1)) {
      // both: rough, top: elastic
      fInterType = 4;
    } 
    else
    if ((sheet[sheetID0].fSheetType==3) && (sheet[sheetID1].fSheetType==3)) {
      // both: rough and elastic
      fInterType = 5;
    } 
    else fInterType = 1;//termination("fInterType undefined in inter ",ID); //TEMP
    //change: for some reason, the code terminates 30 lines BEFORE this, if this line is changed to terminate ???

    // initialization of gap and potentials depends on whether one or two sheets are elastic

    if (fInterType==1) {
      (*this).makeGapX = &interSheet::makeGap1;
      (*this).constraintX = &interSheet::constraint1;
    } 
    else if (fInterType==2) {
      (*this).makeGapX = &interSheet::makeGap2;
      (*this).constraintX = &interSheet::constraint2;
    } 
    else if (fInterType==3) {
      (*this).makeGapX = &interSheet::makeGap3;
      (*this).constraintX = &interSheet::constraint3;
    } 
    else if (fInterType==4) {
      (*this).makeGapX = &interSheet::makeGap4;
      (*this).constraintX = &interSheet::constraint4;
    } 
    else if (fInterType==5) {
      (*this).makeGapX = &interSheet::makeGap5;
      (*this).constraintX = &interSheet::constraint5;
    } 
    else termination ("fInterType could not be assigned");

    if (fPotential) potentialInit();

  } // end if (fElastic vs. fShifted vs. regular cont mech)

  areaPP = areaXY/nGap;

  double fullCEE = fullContactElastEnergy();

  // fftw fields and plans to compute lateral forces and gradient corrections
  if (fDumpLateral||(fPotential==3)) 
  fieldRFFT = (double *)  fftw_malloc( nReal*sizeof(double) );

  ofstream output("params.out",ofstream::app);
  output << ID << "\t# inter info start\n";
  if (fPotential>1) {
    output << surfEnerg << "\t\t# surfEnerg #\n"; // always redundant
    output << potRange << "\t\t# potRange #\n"; // redundant if fFirstCurve==3
    output << potCurve << "\t\t# potCurve #\n"; // redundant if fFirstCurve==2
    output << potCurveRel << "\t\t# potCurveRel #\n"; // redundant if fFirstCurve==1
  }
  if (frictionCoeff) {
    output << frictionCoeff << "\t\t# frictionCoeff #\n";
    output << frictRelax << "\t\t# frictRelax #\n";
  }
  output << 0 << "\t# rmsFullContactStress\n"; // MM2MM to be done
  output << fullCEE << "\t# fullContactElastEnergy\n";
  output << ID << "\t# inter info end\n\n";
  output.close();

  if (frictionCoeff != 0) {
    // sanity checks
    if (fFire) {
      if ( (fireIncrmt!=1)||(fireDecrmt!=1)||(fireRedrct!=0) )
        cerr << "# Resetting FIRE variables because of friction.\n";
      fireIncrmt = 1;
      fireDecrmt = 1;
      fireRedrct = 0;
    }

    double vXeff = sheet[sheetID0].vX - sheet[sheetID1].vX;
    double vYeff = sheet[sheetID0].vY - sheet[sheetID1].vY;
    if (vXeff==0 && vYeff==0) termination("frictionCoeff != 0 but no sliding velocities");
    
    // update sheet's frictRelax
    cerr << "# inter" << ID << " overwrites frictRelax in sheets." << endl;//DEBUG
    sheet[sheetID0].frictRelax = frictRelax;
    sheet[sheetID1].frictRelax = frictRelax;
    
    // initialize old lateral displacements in real space necessary for friction
    if (sheet[sheetID0].nElast) sheet[sheetID0].initFriction();
    if (sheet[sheetID1].nElast) sheet[sheetID1].initFriction();
  }      

} // initSystem


double interSheet::fullContactElastEnergy(){

  int nFour0 = sheet[sheetID0].nFour;
  if (nFour0!= sheet[sheetID1].nFour){
    cerr << "Different resolutions in fullContactElastEnergy() not implemented";
    return(0.);
  }

  // assign gap in real space first 
  double *fullDispZR = (double *) fftw_malloc( nGap*sizeof(double) );
  if (fInterType==1) {
    for (int iGap=0; iGap<nGap; ++iGap) 
    fullDispZR[iGap] = sheet[sheetID0].equilPos[iGap];
  } 
  else if (fInterType==2) {
    for (int iGap=0; iGap<nGap; ++iGap)
    fullDispZR[iGap] = sheet[sheetID1].equilPos[iGap];
  } 
  else if ( (fInterType>=3) && (fInterType<=5) ) {
    for (int iGap=0; iGap<nGap; ++iGap)
    fullDispZR[iGap] = sheet[sheetID1].equilPos[iGap] 
		     - sheet[sheetID0].equilPos[iGap];
  } 

  // transform gap to Fourier
  Complex *fullDispZF = (Complex *) fftw_malloc( nFour0*sizeof(Complex) );
  for (Lint k = 1; k < nFour0; ++k) fullDispZF[k] = 0;
  fftw_plan fullDispZR2F =
  fftw_plan_dft_r2c_2d(nx, ny, fullDispZR, (fftw_complex*) fullDispZF, FFTW_ESTIMATE);
  fftw_execute(fullDispZR2F);

  double fullCEE = 0;
  if ( (fInterType==1) || (fInterType==3) ) {
    for (Lint k = 1; k < nFour0; ++k) {
      Complex stressLoc = sheet[sheetID1].stressZZPre[k] * fullDispZF[k];
      Complex dummy = stressLoc * conj(fullDispZF[k]);
      fullCEE += sheet[sheetID1].weightQ[k]*dummy.real();
    }
  } 
  else if ( (fInterType==2) || (fInterType==4) ) {
    for (Lint k = 1; k < nFour0; ++k) {
      Complex stressLoc = sheet[sheetID0].stressZZPre[k] * fullDispZF[k];
      Complex dummy = stressLoc * conj(fullDispZF[k]);
      fullCEE += sheet[sheetID0].weightQ[k]*dummy.real();
    }
  } 
  else if ( fInterType==5 ) {
    for (Lint k = 1; k < nFour0; ++k) {
      double stressPre =  1./sheet[sheetID0].stressZZPre[k] 
		 	+ 1./sheet[sheetID1].stressZZPre[k];
      stressPre = 1./stressPre;
      Complex stressLoc = stressPre * fullDispZF[k];
      Complex dummy = stressLoc * conj(fullDispZF[k]);
      fullCEE += sheet[sheetID0].weightQ[k]*dummy.real();
    }
  } 
  fullCEE *= areaPP / (2*nGap);

  // free memeory
  fftw_free(fullDispZR);
  fftw_free(fullDispZF);
  fftw_destroy_plan(fullDispZR2F);

  return(fullCEE);
}


void interSheet::dumpGap(){
  // ((*this).*(makeGapX))();
  vector<double*> arrays = {gap};
  vector<string> arrayNames = {"gapZ"};
  if (fPotential) {
    arrays.push_back(gapStress);
    arrayNames.push_back("stressZ");
  }
  string fileName = "gap" + to_string(ID) + ".dat";
  dumpReal(arrays, nx, ny, 1, 1, fileName);
  dumpRealCS(arrays, arrayNames, nx, ny, 1, 1, fileName+"H");
}


void interSheet::initParamsDefault(){

  sheetID0 = ID;
  sheetID1 = ID+1;

  fPotential = fReynolds = fElastic = fPotentialTest = 0;
  fConstraint = 1;

  // potential
  surfEnerg = 1.e-4;
  potRange = 0;
  potCurveRel = 0.4;
  frictionCoeff = 0.0;
  frictRelax = 1;
  relFrictPeak = 1;
  velFrictThresh = 0;

  nGap = 1;

  // output handling
  fDumpGap = 1;
  fDumpFrame = 0;
  resolMovie = 512;
  fDumpLateral = 0;
  
  nTimeOn = 0;
  nTimeOff = 0;

}


void interSheet::writeParams(){

  ofstream output("params.out",ofstream::app);
  output << ID << "\t# inter start\n";

  output << sheetID0 << "\t\t# sheetID0 #\n";
  output << sheetID1 << "\t\t# sheetID1 #\n";

  if (fConstraint) {
    output << fConstraint << "\t\t# fConstraint #\n";
  }

  if (fPotential) {
    output << fPotential << "\t\t# fPotential #\n";
    if (surfEnerg>0) output << surfEnerg  << "\t\t# surfEnerg #\n";
    if (potRange>0)  output << potRange << "\t\t# potRange #\n";
    if (potCurveRel>0) output << potCurveRel << "\t\t# potCurveRel #\n";
    if (potCurve>0) output << potCurve << "\t\t# potCurve #\n";
    if (fPotentialTest) output << fPotentialTest << "\t\t# fPotentialTest #\n";
  }
  if (frictionCoeff>0) {
    output << frictionCoeff  << "\t\t# frictionCoeff #\n";
    output << frictRelax  << "\t\t# frictRelax #\n";
    output << relFrictPeak << "\t\t# relFrictPeak #\n";
    output << velFrictThresh << "\t\t# velFrictThresh #\n";
  }

  if (fReynolds) {
    output << fReynolds << "\t\t# fReynolds #\n";
  }

  if (fElastic) {
    output << fElastic << "\t\t# fElastic #\n";
  }

  output << fDumpGap << "\t\t# fDumpGap #\n";

  if (fDumpFrame) {
    output << fDumpFrame << "\t\t# fDumpFrame #\n";
    output << resolMovie << "\t\t# resolMovie #\n";
  }
   
  if (fDumpLateral) {
    output << fDumpLateral << "\t\t# fDumpLateral #\n";
  }

  if (nTimeOn > 0) {
    output << nTimeOn << "\t\t# nTimeOn #\n";
  }
  if (nTimeOff) {
    output << nTimeOff << "\t\t# nTimeOff #\n";
  }
    

  output << ID << "\t# inter end" << endl << endl;

  if (!ID) writeParamsDefault();

}


void interSheet::writeParamsDefault(){
  ofstream output("params.def",ofstream::app);
  output << "! ! ! ! ! ! ! ! ! selected interSheet (default) values\n\n";
  output << 1 << "\t\t# fConstraint #\n";
  output << 1 << "\t\t# fPotential #\n";
  output << 0.001<< "\t\t# surfEnerg #\n";
  output << 0.4 << "\t\t# potCurveRel #\n";
  output << 0.01 << "\t\t# potRange #\n";
  output << 0.3 << "\t\t# frictionCoeff #\n";
  output << 1.5 << "\t\t# relFrictPeak #\n";
  output << 0.001 << "\t\t# velFrictThresh #\n";
  output << 250 << "\t\t# frictRelax #\n";
  output << fPotentialTest << "\t\t# fPotentialTest #\n";
  output << 1 << "\t\t# fDumpGap #\n";
  output << 1 << "\t\t# fDumpFrame #\n";
  output << 1 << "\t\t# fDumpLateral #\n";
  output << 512 << "\t\t# resolMovie #\n";
  output << 1000 << "\t\t# nTimeOff #\n\n";
  output << 250 << "\t\t# nTimeOn #\n\n";
  output.close();
}


// make gap options

void interSheet::makeGap1(){
  for (int iGap=0; iGap<nGap; ++iGap) 
  gap[iGap] = sheet[sheetID0].equilPos[iGap] - sheet[sheetID1].dispZR[iGap];
}

void interSheet::makeGap2(){
  for (int iGap=0; iGap<nGap; ++iGap)
  gap[iGap] = sheet[sheetID0].dispZR[iGap] - sheet[sheetID1].equilPos[iGap];
}

void interSheet::makeGap3(){
  for (int iGap=0; iGap<nGap; ++iGap)
  gap[iGap] = sheet[sheetID0].equilPos[iGap] 
	    - sheet[sheetID1].equilPos[iGap] - sheet[sheetID1].dispZR[iGap];
}

void interSheet::makeGap4(){
  for (int iGap=0; iGap<nGap; ++iGap)
  gap[iGap] = sheet[sheetID0].equilPos[iGap] + sheet[sheetID0].dispZR[iGap]
	    - sheet[sheetID1].equilPos[iGap];
}

void interSheet::makeGap5(){
  for (int iGap=0; iGap<nGap; ++iGap)
  gap[iGap] = sheet[sheetID0].equilPos[iGap] + sheet[sheetID0].dispZR[iGap]
	    - sheet[sheetID1].equilPos[iGap] - sheet[sheetID1].dispZR[iGap];
}

// make constraint options

void interSheet::constraint(){

  ((*this).*(makeGapX))();
  ((*this).*(constraintX))();

  // compute kinetic energies
//if (sheet[sheetID0].nElast!=0) sheet[sheetID0].tKinHalfStep();
//if (sheet[sheetID1].nElast!=0) sheet[sheetID1].tKinHalfStep();
}

void interSheet::constraint1(){
  for (Lint k = 0; k < nGap; ++k){
    if (gap[k]<0) {
      sheet[sheetID1].dispZR[k] = sheet[sheetID0].equilPos[k] ;
      gap[k] = 0;
    }
  }
  // up-date Fourier coefficients after imposing constraint
  sheet[sheetID1].dispR2F();
}

void interSheet::constraint2(){
  for (Lint k = 0; k < nGap; ++k){
    if (gap[k]<0) {
      sheet[sheetID0].dispZR[k] = sheet[sheetID1].equilPos[k] ;
      gap[k] = 0;
    }
  }
  // up-date Fourier coefficients after imposing constraint
  sheet[sheetID1].dispR2F();
}

void interSheet::constraint3(){
  for (Lint k = 0; k < nGap; ++k){
    if (gap[k]<0) {
      sheet[sheetID1].dispZR[k] += gap[k];
      gap[k] = 0;
    }
  }
  // up-date Fourier coefficients after imposing constraint
  sheet[sheetID1].dispR2F();
}

void interSheet::constraint4(){
  for (Lint k = 0; k < nGap; ++k){
    if (gap[k]<0) {
      sheet[sheetID0].dispZR[k] -= gap[k];
      gap[k] = 0;
    }
  }
  // up-date Fourier coefficients after imposing constraint
  sheet[sheetID0].dispR2F();
}

void interSheet::constraint5(){
  double stiffInt = 1./sheet[sheetID0].stiffMax + 1./sheet[sheetID1].stiffMax;
  stiffInt = 1. / stiffInt;
  double shift0 = stiffInt / sheet[sheetID0].stiffMax;
  double shift1 = stiffInt / sheet[sheetID1].stiffMax;
  for (Lint k = 0; k < nGap; ++k){
    if (gap[k]<0) {
      sheet[sheetID0].dispZR[k] -= shift0 * gap[k];
      sheet[sheetID1].dispZR[k] += shift1 * gap[k];
      gap[k] = 0;
    }
  }
  // up-date Fourier coefficients after imposing constraint
  sheet[sheetID0].dispR2F();
  sheet[sheetID1].dispR2F();
}


void interSheet::makeGradientPotential(){

  double *topo;
  if (fInterType==1) topo = sheet[sheetID0].equilPos;
  else if (fInterType==2) topo = sheet[sheetID1].equilPos;
  else termination("makeGradientPotential() not implemented for fInterType>2!");
  
  // grad calculated in real space 
  for (Lint k = 0; k < nGap; ++k) fieldRFFT[k] = 0;
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      Lint k = ix*ny + iy;
      Lint kR = (ix+1)%nx*ny + iy;
      Lint kL = (ix+nx-1)%nx*ny + iy;
      Lint kT = ix*ny + (iy+1)%ny;
      Lint kB = ix*ny + (iy+ny-1)%ny;
      double gradX = (topo[kR] - topo[kL])*nx/(2*lengthX);
      double gradY = (topo[kT] - topo[kB])*ny/(2*lengthY);
      fieldRFFT[k] = (gradX*gradX + gradY*gradY);
    }
  }

  vPotTot = 0;
  for (int k = 0; k < nGap; ++k) {
    double f = 1;
    double energLocal = surfEnerg*sqrt(1 + fieldRFFT[k]);
    energLocal = min(energLocal, 2*surfEnerg);
    // harmonic repulsion, cosine adhesion
    if (gap[k]<0) {
      gapStress[k] = potCurve*gap[k];
      vPotTot += -energLocal + (gapStress[k])*gap[k]/2;
    } 
    else if (gap[k]>potRange) {
      gapStress[k] = 0;
    } 
    else {
      double zEff = qCosinePot*gap[k];
      vPotTot -= energLocal * (1+cos(zEff)) / 2;
      gapStress[k] =  energLocal * qCosinePot * sin(zEff) / 2;
    }
  } vPotTot *= areaPP;

}


double interSheet::addStress(){
  if (!fPotential) return(0.);

  if (nTimeOff) {
    if (iTime>=nTimeOff) return(0.);
  }

  ((*this).*(makeGapX))();

  double totStress = 0;
  if (fPotential==3) {
     makeGradientPotential();
  } 
  else {
    vPotTot = 0;
    for (int k = 0; k < nGap; ++k) {
      double f = 1;
      vPotTot += ((*this).*(potentialX))(gap[k], &f);
      gapStress[k] = f;
      totStress += f;
    } vPotTot *= areaPP; 
  }
  totStress *= lengthY/ny;

  if (sheet[sheetID1].nElast){
    for (int k = 0; k < nGap; ++k) sheet[sheetID1].stressZR[k] += gapStress[k];
  }

  if (sheet[sheetID0].nElast){
    for (int k = 0; k < nGap; ++k) sheet[sheetID0].stressZR[k] -= gapStress[k];
  }

  if (nTimeOn) {if (iTime<=nTimeOn) return(vPotTot);}

  if (frictionCoeff != 0) vPotTot += addFriction();

  return(vPotTot);

} // addStress


double interSheet::addFriction() {
  //if (iTime < 3) cerr << "  interSheet" << ID << "::addFriction()\n";//FLOW

  double energLocal = 0;

  // to ensure stability, the lateral component is switched on gradually.
  double frictionNow = frictionCoeff;
  if (iTime < nTimeOn) return(0.);
  else if (iTime < 2*nTimeOn) frictionNow *= 0.5*(1+cos(iTime*PI/nTimeOn));

  //change: use Newton's action/reaction axiom and true relative velocities if both surfaces are elastic!
  vector<int> elastIDs = {};
  if (sheet[sheetID0].nElast) elastIDs.push_back(sheetID0);
  if (sheet[sheetID1].nElast) elastIDs.push_back(sheetID1);

  for (int thisID:elastIDs) {
    
    // determine relative sheet velocities
    int oppositeID = sheetID0;
    if (thisID==sheetID0) oppositeID = sheetID1;
    double vXmeanRel = sheet[thisID].vX, vYmeanRel = sheet[thisID].vY;
    vXmeanRel -= sheet[oppositeID].vX;
    vYmeanRel -= sheet[oppositeID].vY;
    
    /*//DEBUG
    if (iTime==nTimeOn+1) {
      cout << "elastic sheet" << thisID << ":\n";
      cout << " frictRelax=" << sheet[thisID].frictRelax << "\n";
      cout << " vXmeanRel=" << vXmeanRel << "\n vYmeanRel=" << vYmeanRel << "\n";
    }//*/

    // apply friction
    double dx = lengthX/nx;
    double dy = lengthY/ny;
    for (Lint k=0; k<nGap; ++k) {
      if (gap[k] >= 0) continue; // no anti-friction allowed.

      double velXnow = vXmeanRel, velYnow = vYmeanRel;
      if (sheet[thisID].nElast>=2) {
        Lint irx = k/ny; 
        Lint iry = k%ny;
        Lint kR = (irx+1)%nx*ny + iry;    // right neighbor with pbc
        Lint kL = (irx+nx-1)%nx*ny + iry; // left neighbor with pbc
        Lint kT = irx*ny + (iry+1)%ny;    // top neighbor with pbc
        Lint kB = irx*ny + (iry+ny-1)%ny; // bottom neighbor with pbc
        double duy_dx = (sheet[thisID].dispYR[kR] - sheet[thisID].dispYR[kL])/(2*dx);
        double duy_dy = (sheet[thisID].dispYR[kT] - sheet[thisID].dispYR[kB])/(2*dy);
        velYnow += duy_dx*vXmeanRel + duy_dy*vYmeanRel;
        
        // friction in X and Y direction
        if (sheet[thisID].nElast==3) {
          double dux_dx = (sheet[thisID].dispXR[kR] - sheet[thisID].dispXR[kL])/(2*dx);
          double dux_dy = (sheet[thisID].dispXR[kT] - sheet[thisID].dispXR[kB])/(2*dy);
          velXnow += dux_dx*vXmeanRel + dux_dy*vYmeanRel;
          double velLateral = velXnow*velXnow + velYnow*velYnow;
          velLateral = sqrt(velLateral);
          double xi = 1.;
          // rate dependent friction model, see eq. 7 in Menga2022 (doi.org/10.1007/s11071-022-07997-0)
          if (velFrictThresh > 0) {
            double dummy = 2*(1 - asin(1./relFrictPeak)/M_PI);
            xi = relFrictPeak*sin(dummy*atan(velLateral/velFrictThresh));
          }
          if (velLateral==0) velLateral = 1;// avoid division by 0.
          sheet[thisID].stressYR[k] += xi*frictionNow*sheet[thisID].stressZR[k]*velYnow/velLateral;
          sheet[thisID].stressXR[k] += xi*frictionNow*sheet[thisID].stressZR[k]*velXnow/velLateral;
        }
        
        // friction in Y direction only
        else {
          int velSign = (velYnow>0) - (velYnow<0);
          sheet[thisID].stressYR[k] += frictionNow*velSign*sheet[thisID].stressZR[k];
        }
      }
    }
  }
  energLocal *= areaPP; //change: if both are elastic and Newton is applied, make sure not to count twice.

  return(0.);
}


void interSheet::potentialInit(){

  if (fFirstCurve==-1) termination("define potCurve in interSheet ", ID);

  // determine maximum interacial stiffness first
  interStiffMax = 0; 
  if ( (fInterType==1) || (fInterType==3) ) {
    interStiffMax = sheet[sheetID1].stiffMax;
  } 
  else if ( (fInterType==2) || (fInterType==4) ) {
    interStiffMax = sheet[sheetID0].stiffMax;
  } 
  else if (fInterType==5) {
    interStiffMax = 1./sheet[sheetID0].stiffMax + 1./sheet[sheetID1].stiffMax;
    interStiffMax = 1./interStiffMax;
  } 
  else termination("fInterType not defined in potentialInit");

  if (fFirstCurve==1) potCurve = potCurveRel * interStiffMax;

  if (fFirstCurve==3) {
    if (fPotential != 2) termination("fPotential not compatible with potRange");
    potCurve = PI*PI*surfEnerg/(2*potRange*potRange);
    potCurveRel = potCurve / interStiffMax;
    if (potCurve > interStiffMax) {
      cerr << "# CAUTION: Potential stiffness " << potCurve << " vs "
      << interStiffMax << " elastic stiffness in inter " << ID << "\n";
    }
  }

  // define pointer
  if (fPotential==1) {
    // harmonic hard-wall: only needs potCurve
    (*this).potentialX = &interSheet::potential1;
  } 
  else if (fPotential==2 || fPotential==3) {
    // harmonic repulsion, cosine adhesion: 
    (*this).potentialX = &interSheet::potential2;
    qCosinePot = sqrt(2*potCurve/surfEnerg);	// q
    potRange = PI / qCosinePot;		// cutOffRadius
  } 
  else termination ("fPotential not implemented");

  if (fPotentialTest) potentialTest();

}


double interSheet::potential1(double z, double *s){
  // harmonic hard wall
  *s = 0;
  if (z>0) return(0.);
  *s = potCurve*z;
  return( (*s)*z/2 );
}

double interSheet::potential2(double z, double *s){
  // harmonic repuslsion, cosine adhesion
  double e;
  if (z<0) { *s = potCurve*z; e = -surfEnerg + (*s)*z/2 ;}
  else if (z>potRange) {*s = 0; e = 0;}
  else { 
    double zEff = qCosinePot*z;
    e  = -surfEnerg * (1+cos(zEff)) / 2;
    *s =  surfEnerg * qCosinePot * sin(zEff) / 2;
  } 
  return(e);
}

double interSheet::potential3(double z, double *s){
  // Morse
  return(0.);
}


void interSheet::potentialTest(){

  if (fPotential==3) return;

  ofstream potTest("potTest.dat");
  potTest << "# potential:= " << fPotential << "\n\n";
  potTest << "potRange: " << potRange << "\n";
  potTest << "qCosinePot: " << qCosinePot << "\n";
  potTest << "\n";

  double zScale = 0;
  if (nx==1) zScale = lengthY/ny;
  else if (ny==1) zScale = lengthX/nx;
  else zScale = min(lengthX/nx, lengthY/ny);

  for (double z = -4*zScale; z <= 12*zScale; z += zScale/20) {
    double f, e = ((*this).*(potentialX))(z, &f);
    potTest << z << "\t" << e << "\t" << f << "\n";
  }

  potTest.close();
  if (fPotentialTest!=1) termination("potentialTest calls termination");
}


void interSheet::initMeasure(){
  string moniName = "iMoni" + to_string(ID) + "-";
  if (ny<1000) {moniName = moniName + "0";}
  if (ny<100)  {moniName = moniName + "0";}
  if (ny<10)   {moniName = moniName + "0";}
  moniName = moniName + to_string(ny) + ".dat";
  moni.open(moniName, ofstream::out);
  moni << "# mdTime";
  moni << "\tstressZ";
  if (fDumpLateral) moni << "\tstressY\tstressX";
  if (fDumpGap >= 2 && (fPotential||fConstraint)) {
    moni << "\tmeanGap\tmeanPosGap\trelContArea";
    if (fPotential > 1) moni << "\trelAttractArea";
  }

#ifdef _ROTSYM_
  // quick fix for Bo
  moni << "\tcontRadius(maxStress)\tcontRadius(gap=0)";
#endif

  moni << endl;
}


void interSheet::get_iqxiqy(Lint k){ iqx = k/nyHP1; iqy = k%nyHP1; };


void interSheet::measure(){

  if (iTime%10) return;
  moni << mdTime;

  ((*this).*(makeGapX))(); // MM2MM: Why is this needed?

 // determine contact topography from a rough sheet
  int roughID = (fInterType==2)? sheetID1 : sheetID0;


  if (fDumpLateral||fPotential==3) {
    for (int k = 0; k < nxny; ++k) fieldRFFT[k] = sheet[roughID].equilPos[k];
    if (sheet[roughID].nElast) {
      for (int k = 0; k < nxny; ++k) fieldRFFT[k] += sheet[roughID].dispZR[k];
    }
  }

  double stressZ = 0;
  for (int k=0; k<nxny; ++k) stressZ += gapStress[k];

  stressZ /= nxny;
  moni << "\t" << stressZ;

  if (fDumpLateral) {
    // calculate stressX and stressY from gapStress and topography gradient
    double normX = nx/(2*lengthX), normY = ny/(2*lengthY);
    double stressX = 0, stressY = 0;
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        Lint k = ix*ny + iy;
        Lint kR = (ix+1)%nx*ny + iy;    // right neighbor with pbc
        Lint kL = (ix+nx-1)%nx*ny + iy; // left neighbor with pbc
        Lint kT = ix*ny + (iy+1)%ny;    // top neighbor with pbc
        Lint kB = ix*ny + (iy+ny-1)%ny; // bottom neighbor with pbc
        double slopeX = (fieldRFFT[kR] - fieldRFFT[kL])*normX;
        double slopeY = (fieldRFFT[kT] - fieldRFFT[kB])*normY;
        stressY -= slopeY*gapStress[k];
        stressX -= slopeX*gapStress[k];
      }
    }
    stressY /= nxny; stressX /= nxny;
    moni << "\t" << stressY << "\t" << stressX;
  }

  // compute mean gap and contact area
  if (fDumpGap >= 2 && (fPotential||fConstraint)) {
    double threshDist = lengthY*1.e-10;
    contactArea = 0, meanGap = 0, meanPosGap = 0;
    for (int iGap=0; iGap<nGap; ++iGap) {
      meanGap += gap[iGap];
      if (gap[iGap] > threshDist) meanPosGap += gap[iGap];
      else contactArea += 1;
    }
    meanGap /= nGap; attractArea /= nGap; meanPosGap /= nGap; contactArea /= nGap;
    moni << "\t" << meanGap;
    moni << "\t" << meanPosGap;
    moni << "\t" << contactArea;

    // adhesive contact area (if applicable)
    if (fPotential > 1) {
      attractArea = 0;
      threshDist = 0.75*potRange; //CM-Adjust
      for (int iGap=0; iGap<nGap; ++iGap) {
        if ((gap[iGap] < threshDist) && (gap[iGap] > 0)) attractArea += 1;
      }
      attractArea /= nGap;
      moni << "\t" << attractArea;
    }
  }

#ifdef _ROTSYM_// formerly _QUICKFIXBO_
  // quick fix for Bo
  Lint nxH = nx/2;
  Lint iyMaxTension = ny/2, iyGapZero = ny/2, nxnyH = nxH*ny;
  Lint k = nxnyH + iyMaxTension;
  double stressMax = gapStress[k];
  double gapNow = gap[k];
  for (int iy=ny/2+1; iy<ny; ++iy) {
    k = nxnyH + iy;
    if ( gapStress[k] > stressMax) {
      stressMax = -gapStress[k];
      iyMaxTension = iy;
    } 
    if (gap[k]>=0 && gapNow<0) iyGapZero = iy;
    gapNow = gap[k];
  } 
  moni << "\t" << (iyMaxTension-ny/2)*lengthY/ny << "\t" << (iyGapZero-ny/2)*lengthY/ny;
#endif

  moni << endl;
}


void interSheet::dumpFrame(int iFrame){
  if (!fDumpFrame) return;
  string fileName = "frames/inter"+to_string(ID)+"contact."+to_string(iFrame)+".dat";
  ofstream output;
  output.open(fileName);

  int xStep = (nx>resolMovie) ? (nx-1)/resolMovie + 1 : 1;
  int yStep = (ny>resolMovie) ? (ny-1)/resolMovie + 1 : 1;

  output << "#" << nx/xStep << "\t" << ny/yStep << "\n";
  //output << "#ix\tiy (in contact)\n";
  output << "-1\t" << contactArea << "\n"; // make sure file contains points for gnuplot

  // dump positions of points in repulsive contact
  double thresh = lengthY*1.e-10;
  for (int ix = 0; ix < nx; ix += xStep) {
    for (int iy = 0; iy < ny; iy += yStep) {
      Lint k = ix*ny + iy;
      if (gap[k] < thresh) output << ix/xStep << "\t" << iy/yStep << "\n";
    }
  }
  output.close();

  // dump positions of points in attractive contact
  if (fPotential > 1) {
    fileName = "frames/inter"+to_string(ID)+"attract."+to_string(iFrame)+".dat";
    output.open(fileName);
    output << "#" << nx/xStep << "\t" << ny/yStep << "\n";
    //output << "#ix\tiy (in contact)\n";
    output << "-1\t" << attractArea << "\n"; // make sure file contains points for gnuplot

    thresh = 1.*potRange;
    for (int ix = 0; ix < nx; ix += xStep) {
      for (int iy = 0; iy < ny; iy += yStep) {
        Lint k = ix*ny + iy;
        if ((gap[k] < thresh) && (gap[k] > 0)) output << ix/xStep << "\t" << iy/yStep << "\n";
      }
    }
    output.close();
  }
}


void interSheet::outMeasure(){

  ofstream output("params.out",ofstream::app);
  ((*this).*(makeGapX))();
  output << ID << "\t# inter measure start\n";

  if (fPotential) {
    output << vPotTot << "\t\t# vPotTot\n";
    if (areaXY!=1) output << vPotTot/areaXY << "\t\t# vPotTot/areaX\n";
  }

  if (fPotential||fConstraint) {
    double thresh = lengthY*1.e-10;
    contactArea = 0, meanGap = 0, meanPosGap = 0;
    for (int iGap=0; iGap<nGap; ++iGap) {
      meanGap += gap[iGap];
      if (gap[iGap]>thresh) meanPosGap += gap[iGap];
      else contactArea += 1;
    }
    meanGap /= nGap; meanPosGap /= nGap; contactArea /= nGap;
    output << meanGap << "\t\t# mean gap\n";
    if (!fConstraint) { 
      output << meanPosGap << "\t\t# mean positive gap\n";
      output << 2*meanPosGap - meanGap << "\t\t# estimated converged gap\n";
    }
    output << contactArea << "\t\t# relative contact area\n";
    if (areaXY!=1) output << contactArea*areaXY << "\t\t# absolute contact area\n";
  }

  output << ID << "\t# inter measure end\n\n";
  output.close();
}


void interSheet::outSystem(){
  if (!fDumpGap) return; //CM2CM: check if makeGapX() is needed regardless.
  ((*this).*(makeGapX))();
  dumpGap();
}


interSheet::~interSheet(){
  moni.close();
  fftw_free(gap);
  fftw_free(gapStress);
  fftw_free(fieldRFFT);
}
