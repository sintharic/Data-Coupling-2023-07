#ifndef GFMDSHEET_H
#define GFMDSHEET_H
//endif at EOF

#include "header.h"

class gfmdSheet{

private:

  int ID;

  double areaPP;
  double zeroModeMass;

  double pressure, pressInit, pressFinal;
  double forceInit, forceFinal;
  double pressX, pressY;

  int fConstCOMx, fConstCOMy;
  double xConstCOM, yConstCOM, vxConstCOM, vyConstCOM;

  int fConstCOM; // = 1 controlled COM
  int fSteppedRamp, rampSteps, rampRelax;
  double zConstCOM, vzConstCOM, stressCOM, dzRamp, ddzRamp, dpRamp, ddpRamp;
  double stressInterZ, stressInterY, stressInterX;
  double rVeloTurnStep; // relative number of steps before inversion of ramp
  int nVeloTurnStep;	// absolute number of steps ...
  int nVeloTransition, tStartTransition, tEndTransition, t0Transition;
  double forceTurnaround, pressTurnaround, absTurnPress;
  bool checkTurnPress;
  short int sgnTurnPress;
  double vzInit;

  int fzOpposite, fExtrapolateInf;
  double zOpposite, zInfOffset, zInfErrorOld;

  double dxOffset, dyOffset;

  int fOnSitePotential;
  double fOnSitePeriod, fOnSiteFreq;
  double vXOnSite, vYOnSite, stressXOnSite, stressYOnSite;
  double xOnSite, yOnSite;
  double surfEnergOS, potCurveRelOS;
  double frictionCoeffOS;
  double pressPhi, pressTheta;
  int fKelvinVoigt;
  double tauKV, scalKV;

  //NEW-MW
  int fMaxwell, nMaxwell;
  double massGFMD;
  vector<double> invTauMw, stiffFacMw;
  vector<Complex*> dispZMw; 
  void stressMaxwell();
  ofstream maxwellout;
  int fDumpMaxwell;

  double rXhertz, rYhertz, hertzExpnt;
  double rSphere;

  double hurst, lambdaR, lambdaS, rRoughNorm, peklenik;
  int fRollOff, fRoughNorm, fBoxMuller;

  int fTopoRead, fTopoAdd;

  double fAddSWR, nqxAddSWR, nqyAddSWR; // SWR = single-wavelength roughness
  double heightSWR;

  int fFlatPunch;
  double radiusFlatPunch, heightFlatPunch;

  double dzStep;

  // space-holder fields to be used in FFTW
  double *fieldRFFT;
  Complex *fieldFFFT;

  Lint kFastest, rCentral;

  int rampTime;
  double damping; 

  string konfigName;
  int f3dMovie, resolMovie;

  ofstream moni, ramp, frictout;
  int fFriction = 0;

  void initParamsDefault();
  bool initParamsFile();
  void initMaxwell();
  void writeParamsDefault();

  void addConfigR();
  void addConfigF();

  void initConfig();
  void addHertz();
  void addSphere();
  void addSelfAffine(), selfAffineAnalysis(vector<double> *derivs);
  void addSWR();
  void addFlatPunch();

  void makeSteps();//CM-dzSteps
  double sqrtSpec(double);

  void prepareDispZ();
  void prepareDispY();
  void prepareDispX();
  void prepareStiffZ();
  void prepareStiffZY();
  void prepareStiffZYX();
  void prepareGFMD();
  
  void initLateral();
  void moveLateral();

  double deltaIndenter(), einsteinSolid();
  void OnSiteHertzConstraint();
  double OnSiteHertzPotential();
  double stressOpposite();
  double stressGFMD();
  double stressGFMDzy();
  double stressGFMDzyx();
  double stressKV();
  void updateStress();
  void thermostat();

  void moveCOM(), moveOpposite(), turnaroundCOM();
  double dispInf();

  double computeChi2SS();
  void propagateSteadySlide(); 

  // divisions and modulos could be optimized with bit-shift operations!
  Lint irx, iry, iqx, iqy;
  inline Lint getLinRIndex(Lint ix, Lint iy) {return(ix*ny+iy);};
  inline void get_irxiry(Lint k){ irx = k/ny; iry = k%ny; };
  inline Lint getLinCIndex(Lint ix, Lint iy) {return(ix*nyHP1+iy);};
  inline void get_iqxiqy(Lint k){ iqx = k/nyHP1; iqy = k%nyHP1; };

  double getQ(Lint);

public:

  Lint nx, ny;
  Lint nxH, nyHP1, nxny, nReal, nFour;

  double dx, dqx, dqx2, dy, dqy, dqy2;

  int fTopography, nElast;
  int fMassWeightg;
  int fLateral, fSteadySlide;
  bool fSliding;
  double vX, vY;
  int frictRelax;
  double chi2SteadySlide, chi2SteadySlideOld, weightSS, invStiffSS_COM, invStiffSS_COMRef, chi2SSRestart, prefacSSnxny;
  int oldSSPress;

  int fSheetType; //	if rough: += 1; if nElast: += 2; makes 3 different types

  // energy: v = stiffness * q^elastExpnt * u^2 / 2
  // contactMod = 2 * stiffness; ==> default: stiffness = 0.5
  static const int nLayer=4;
  double stiffness[nLayer], stiffHigh[nLayer], elastExpnt[nLayer];
  double poisson[nLayer], thickness[nLayer];
  int fThickness[nLayer];
  double stiffMax, stiffMin, contModEff, massScal;

  double tKinHalfStep();
  double tKinetic, vPotTot, tKinetic1, tKinetic2;
  double vOnSitePot, vElastic, vGravit;

  // config fields
  double *equilPos;
  Complex *equilPosF, *equilPos0F;

  double  *dispZR, *stressZR, rKV_LPF;
  Complex *dispZF, *dispZFold, *stressZF, *rhsKV_LPF; 

  int fDispX, fDispY;
  Complex relXXpre, relXYpre, relXZpre, relYYpre, relYZpre;
  double  *dispXR, *stressXR, *dispYR, *stressYR;
  Complex *dispXFold, *dispYFold;
  Complex  *dispXF, *stressXF, *dispYF, *stressYF;

  short int *weightQ;
  double *invMass, *stressZZPre;
  double *stressXXPre, *stressXYPre, *stressYYPre;
  //double *stressOffsetPre;//TODO-remove
  double *ImStressXZPre, *ImStressYZPre;
  //Complex *stressXZPre, *stressYZPre;//TODO-remove
  Complex *preFacCorrSS;

  // fftw_plans
  fftw_plan dispZR2F, dispZF2R, stressZR2F;
  fftw_plan dispYR2F, dispYF2R, stressYR2F;
  fftw_plan dispXR2F, dispXF2R, stressXR2F;
  fftw_plan equilPosF2R;

  // called from contMech
  void initParams(int);
  void initSystem();
  void initFriction();
  void writeParams();
  void initMeasure();
  double propagate();
  double propagateKV();
  void zeroStress();
  double addStress();
  void measure();
  void fireHalt();
  void fireRedirect();
  void outMeasure();
  void outSystem();
  void dumpFrame(int);
  void dumpConfig();
  void dumpMaxwell();
  ~gfmdSheet();

  // called / set from inter
  bool fRealSpaceStress;
  void dispF2R();
  void dispR2F();
  void stressF2R();

};

#endif