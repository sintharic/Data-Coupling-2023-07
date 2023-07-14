# ROUGH
#inpaths = ["/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press1p0e-03/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press1p6e-03/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press2p5e-03/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press4p0e-03/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press6p3e-03/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press1p0e-02/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press1p6e-02/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press2p5e-02/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press4p0e-02/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press6p3e-02/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press1p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press2p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press3p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press4p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press5p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press6p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press7p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press8p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press9p0e-01/frict1p0",
#"/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/press1p0e+00/frict1p0"]
#rigidpath = "/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/konfig0E.dat"
#outpath = "/Users/christian/sim/coupled/H0p8-n0p49-h0p10-final/resultsV3/res-frict1p0.dat"
#OnSiteHertz = False

# HERTZ
#inpaths = ["/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-04/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force1p6e-04/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force2p5e-04/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force4p0e-04/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force6p3e-04/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-03/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force1p6e-03/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force2p5e-03/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force4p0e-03/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force6p3e-03/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-02/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force2p0e-02/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force3p0e-02/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force4p0e-02/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force5p0e-02/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force6p0e-02/geometric",
#"/Users/christian/sim/coupled/Hertz-final/Nvar/force7p0e-02/geometric"]
##"/Users/christian/sim/coupled/Hertz-final/Nvar/force8p0e-02/geometric",
##"/Users/christian/sim/coupled/Hertz-final/Nvar/force9p0e-02/geometric",
##"/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-01/geometric"]
#outpath = "/Users/christian/sim/coupled/Hertz-final/Nvar/resultsV3/res-geometric.dat"
inpaths = ["/Users/christian/sim/coupled/Hertz-final/N1024-test/force1p0e-02/stat-geom",
"/Users/christian/sim/coupled/Hertz-final/N1024-test/force1p5e-02/stat-geom",
"/Users/christian/sim/coupled/Hertz-final/N1024-test/force2p3e-02/stat-geom",
"/Users/christian/sim/coupled/Hertz-final/N1024-test/force3p5e-02/stat-geom",
"/Users/christian/sim/coupled/Hertz-final/N1024-test/force5p3e-02/stat-geom",
"/Users/christian/sim/coupled/Hertz-final/N1024-test/force8p0e-02/stat-geom"]
outpath = "/Users/christian/sim/coupled/Hertz-final/N1024-test/resultsV3/res-stat-geom.dat"
rigidpath = ""
OnSiteHertz = True



import numpy as np
import cmutils as cm 
import os, sys
cm.WARN = False
import matplotlib.pyplot as plt 
plt.close('all')


def process(paths):
  result = np.zeros((len(paths),15))  


  # sanity checks
  for path in paths:
    if not os.path.isdir(path): sys.exit(path + " does not exist")
  outfolder = os.path.split(outpath)[0]
  if not os.path.isdir(outfolder): 
    print(f"Making new folder: {outfolder}")
    os.makedirs(outfolder)



  # import rigid surface
  if not OnSiteHertz:
    dat0 = cm.readConfig(rigidpath,2)
    meanRough = dat0.mean()



  for iPath,path in enumerate(paths):

    # read params
    press = force = 0
    mu = 0
    frictRelax = 0
    nx = ny = 0
    lengthX = lengthY = 0
    relCont = -1
    poisson = 1
    thick = 0
    rXhertz = rYhertz = 0
    fullContElaEnergy = 0
    vX = vY = 0

    if path[-1] != "/": path += "/"
    print(path)
    with open(path+"params.out","r") as params:
      for line in params:
        if "# nxGlobal #" in line: nx = int(line.split()[0])
        elif "# nyGlobal #" in line: ny = int(line.split()[0])
        elif "# lengthX #" in line: lengthX = float(line.split()[0])
        elif "# lengthY #" in line: lengthY = float(line.split()[0])
        elif "# rXhertz #" in line: rXhertz = float(line.split()[0])
        elif "# rYhertz #" in line: rYhertz = float(line.split()[0])
        elif "# pressInit #" in line: press = float(line.split()[0])
        elif "# forceInit #" in line: force = float(line.split()[0])
        elif "# vX #" in line: vX = float(line.split()[0])
        elif "# vY #" in line: vY = float(line.split()[0])
        elif "# vXOnSite #" in line: vX = float(line.split()[0])
        elif "# vYOnSite #" in line: vY = float(line.split()[0])
        elif "# frictionCoeff #" in line: mu = float(line.split()[0])
        elif "# frictionCoeffOS #" in line: mu = float(line.split()[0])
        elif "# frictRelax #" in line: frictRelax = int(line.split()[0])
        elif "# relative contact area" in line: relCont = float(line.split()[0])
        elif "# poisson0" in line: poisson = float(line.split()[0])
        elif "# thickness0" in line: thick = float(line.split()[0])
        elif "# fullContactElastEnergy" in line: fullContElaEnergy = float(line.split()[0])
    if ny == 0: ny = nx
    if lengthY == 0: lengthY = lengthX
    if rYhertz == 0: rYhertz = rXhertz
    if (force !=0) and (press==0): press = force/lengthX/lengthY
    if (press !=0) and (force==0): force = press*lengthX*lengthY
    dx = lengthX/nx
    dy = lengthY/ny


    # import config(s)
    if OnSiteHertz:
      #X,Y = cm.XY((nx,ny))
      X,Y = np.ogrid[:nx,:ny]
      X = (X-nx/2)*dx
      Y = (Y-ny/2)*dy
      dat0 = np.square(X)/(2*rXhertz) + np.square(Y)/(2*rYhertz)
      dat1 = cm.readConfig(path+"konfig0Dzyx.dat",2)
    else:
      dat1 = cm.readConfig(path+"konfig1Dzyx.dat",2)
    
    

    # system descriptors
    result[iPath,0] = poisson
    result[iPath,1] = thick
    result[iPath,2] = mu
    if OnSiteHertz: result[iPath,3] = press*lengthX*lengthY
    else: result[iPath,3] = press


    # macroscopic observables
    mask = dat1 >= dat0
    if OnSiteHertz: relCont = 1.*mask.sum()/nx/ny
    dh_dx = ( np.diff(dat0,append=dat0[:1,:],axis=0) + np.diff(dat0,prepend=dat0[-1:,:],axis=0) ) / (2*dx)
    dh_dy = ( np.diff(dat0,append=dat0[:,:1],axis=1) + np.diff(dat0,prepend=dat0[:,-1:],axis=1) ) / (2*dy)
    sq_grad = np.square(dh_dx) + np.square(dh_dy)
    rmsgrad_global = np.sqrt(sq_grad.mean())
    rmsgrad_cont = np.sqrt(sq_grad[mask].mean())
    norm = np.copy(dat0)
    norm[norm==0] = 1e-36
    gap = (dat0 - dat1)#/norm
    if OnSiteHertz: result[iPath,6] = relCont*lengthX*lengthY
    else: result[iPath,6] = relCont
    result[iPath,7] = gap.mean()
    gap[gap<0] = 0
    result[iPath,8] = gap.mean()
    result[iPath,9] = rmsgrad_cont
    mGradX = dh_dx[mask].mean() # gradX in contact
    mGradY = dh_dy[mask].mean() # gradY in contact
    result[iPath,10] = mGradY
    result[iPath,11] = mGradX
    mGradX = dh_dx[np.logical_not(mask)].mean() # gradX in non-contact
    mGradY = dh_dy[np.logical_not(mask)].mean() # gradY in non-contact
    frac_pos_gradX = np.logical_and(mask, dh_dx>0).sum()/mask.sum()
    frac_pos_gradY = np.logical_and(mask, dh_dy>0).sum()/mask.sum()
    result[iPath,12] = frac_pos_gradY
    result[iPath,13] = frac_pos_gradX
    result[iPath,-1] = fullContElaEnergy


    # geometric stress in lateral directions
    if mu==0:
      result[iPath,4] = 0 #geomStressY
      result[iPath,5] = 0 #geomStressX
      continue

    if OnSiteHertz:
      geomStressLateral = np.loadtxt(path+"frict0-"+str(ny).zfill(4)+".dat", usecols=(4,5))[-1,:] # geomStressY, geomStressX
      stressCOMLateral = np.loadtxt(path+"frict0-"+str(ny).zfill(4)+".dat", usecols=(2,3))[-1,:] # stressYF[0], stressXF[0]
      result[iPath,4] = geomStressLateral[0]*lengthX*lengthY
      result[iPath,5] = geomStressLateral[1]*lengthX*lengthY
    else:
      geomStressLateral = np.loadtxt(path+"iMoni0-"+str(ny).zfill(4)+".dat", usecols=(2,3))[-1,:] # geomStressY, geomStressX
      stressCOMLateral = np.loadtxt(path+"frict1-"+str(ny).zfill(4)+".dat", usecols=(2,3))[-1,:] # stressYF[0], stressXF[0]
      result[iPath,4] = geomStressLateral[0]
      result[iPath,5] = geomStressLateral[1]
      
    

    # deviation from expectation value
    target = mu*press
    frictStress = np.sqrt(np.square(stressCOMLateral).sum())
    err = np.abs(frictStress - target)/target
    print(f"relative error in friction force: {err:.4e}")

  return result

def save(filepath, data):
  if OnSiteHertz: header = "poisson\tthickness\tfrictCoeff\tforce\tgeomForceY\tgeomForceX\tabsContArea\tmeanGap\tmeanPosGap\trmsGrad_cont\tmGradY_cont\tmGradX_cont\tfrac_pos_gradY\tfrac_pos_gradX\tfullContElaEnergy"
  else: header = "poisson\tthickness\tfrictCoeff\tpressure\tgeomStressY\tgeomStressX\trelContArea\tmeanGap\tmeanPosGap\trmsGrad_cont\tmGradY_cont\tmGradX_cont\tmGradY_noncont\tmGradX_noncont\tfrac_pos_gradY\tfrac_pos_gradX\tfullContElaEnergy"

  np.savetxt(filepath, data, fmt="%g", header=header)

def main():
  data = process(inpaths)
  save(outpath, data)

if __name__ == "__main__": main()