inpath = "/Users/christian/sim/coupled/Hertz-final/Nvar/force5p0e-02/none/"

do_check = True
do_export = True
zoom_factor = 3
relative_pressure_tolerance = 1e-06



import numpy as np 
import cmutils as cm 
import matplotlib.pyplot as plt 
plt.close("all")
import mpl 
mpl.RevTex()



def zoom(array, factor):
  assert factor>1
  start0 = int(array.shape[0]*(1 - 1./factor)/2)
  end0 = array.shape[0] - start0
  start1 = int(array.shape[1]*(1 - 1./factor)/2)
  end1 = array.shape[1] - start1
  return array[start0:end0,start1:end1]



def process(path):
  if path[-1] != "/": path += "/"
  print(path)

  vX = vY = 0
  lengthX = lengthY = 0
  nx = ny = 0
  mu = press = force = 0 
  rXhertz = rYhertz = 1
  fid = open(path+"params.out","r")
  for line in fid.readlines():
    if "# nxGlobal #" in line: nx = int(line.split()[0])
    elif "# nyGlobal #" in line: ny = int(line.split()[0])
    elif "# lengthX #" in line: lengthX = float(line.split()[0])
    elif "# lengthY #" in line: lengthY = float(line.split()[0])
    elif "# vXOnSite #" in line: vX = float(line.split()[0])
    elif "# vYOnSite #" in line: vY = float(line.split()[0])
    elif "# rXhertz #" in line: rXhertz = float(line.split()[0])
    elif "# rYhertz #" in line: rYhertz = float(line.split()[0])
    elif "# vX #" in line: vX = float(line.split()[0])
    elif "# vY #" in line: vY = float(line.split()[0])
    elif "# frictionCoeff" in line: mu = float(line.split()[0])
    elif "# pressInit #" in line: press = float(line.split()[0])
    elif "# forceInit #" in line: force = float(line.split()[0])
  fid.close()
  if ny==0: ny = nx
  if lengthY==0: lengthY = lengthX
  dx = lengthX/nx
  dy = lengthY/ny 
  vExt = np.hypot(vX,vY)
  if (press==0) and (force!=0): press = force/(lengthX*lengthY)

  uz = cm.readConfig(path+"konfig0Dzyx.dat",2)
  uy = cm.readConfig(path+"konfig0Dzyx.dat",4)
  ux = cm.readConfig(path+"konfig0Dzyx.dat",6)
  pz = cm.readConfig(path+"konfig0Dzyx.dat",3)
  py = cm.readConfig(path+"konfig0Dzyx.dat",5)
  px = cm.readConfig(path+"konfig0Dzyx.dat",7)
  pthresh = pz.mean()*relative_pressure_tolerance
  mask = pz > pthresh

  # compute all spatial derivatives with periodic boundary conditions
  dux_dx = ( np.diff(ux,append=ux[:1,:],axis=0) + np.diff(ux,prepend=ux[-1:,:],axis=0) ) / (2*dx)
  dux_dy = ( np.diff(ux,append=ux[:,:1],axis=1) + np.diff(ux,prepend=ux[:,-1:],axis=1) ) / (2*dy)
  duy_dx = ( np.diff(uy,append=uy[:1,:],axis=0) + np.diff(uy,prepend=uy[-1:,:],axis=0) ) / (2*dx)
  duy_dy = ( np.diff(uy,append=uy[:,:1],axis=1) + np.diff(uy,prepend=uy[:,-1:],axis=1) ) / (2*dy)
  duz_dx = ( np.diff(uz,append=uz[:1,:],axis=0) + np.diff(uz,prepend=uz[-1:,:],axis=0) ) / (2*dx)
  duz_dy = ( np.diff(uz,append=uz[:,:1],axis=1) + np.diff(uz,prepend=uz[:,-1:],axis=1) ) / (2*dy)

  # compute local lateral velocities
  vx_loc = - dux_dx*vX - dux_dy*vY
  vy_loc = - duy_dx*vX - duy_dy*vY
  v_rel = np.hypot(vX-vx_loc, vY-vy_loc)

  # compute "geometric" stresses
  stress_geom_x = pz*duz_dx
  stress_geom_y = pz*duz_dy


  if mu != 0:

    # check if stress in lateral directions is calculated correctly
    sigma_x = mu*pz*(vX - vx_loc)/v_rel
    sigma_y = mu*pz*(vY - vy_loc)/v_rel
    print(f"max error in sigma_x: {np.abs(sigma_x+px).max():.4e} (mean: {sigma_x.mean():.4e}, max: {sigma_x.max():.4e})")
    print(f"max error in sigma_y: {np.abs(sigma_y+py).max():.4e} (mean: {sigma_y.mean():.4e}, max: {sigma_y.max():.4e})")

    F_x = sigma_x.sum()*dx*dy
    F_y = sigma_y.sum()*dx*dy
    wdiss1 = mu*pz*v_rel
    
    F_l = mu*pz.sum()*dx*dy

    F_plus_Fgeom_x = (sigma_x + stress_geom_x).sum()*dx*dy
    F_plus_Fgeom_y = (sigma_y + stress_geom_y).sum()*dx*dy
    

    if abs(vX) > abs(vY):
      wdiss2 = mu*pz*vX*(1 + dux_dx)
      geom_l = stress_geom_x.sum()*dx*dy
    elif abs(vY) > abs(vX):
      wdiss2 = mu*pz*vY*(1 + duy_dy)
      geom_l = stress_geom_y.sum()*dx*dy

    print("\nDissipated power:")
    print(f" P_true    : {wdiss1.sum()*dx*dy :.6e}")
    print(f" P_long    : {wdiss2.sum()*dx*dy :.6e}")
    print(f" F_true*v0 : {np.hypot(F_x,F_y)*vExt :.6e}")
    print(f" F_long*v0 : {F_l*vExt :.6e}")
    print(f" (F+F_geom)_true*v0 : {np.hypot(F_plus_Fgeom_x,F_plus_Fgeom_y)*vExt :.6e}")
    print(f" (F+F_geom)_long*v0 : {(F_l+geom_l)*vExt :.6e}")

    # dissipated power
    pressEff = press*nx*ny/mask.sum()
    Pdiss = wdiss1/(mu*pressEff*vExt) # undimensionalized

  else:
    Pdiss = np.zeros_like(pz)

  # make gap
  X,Y = cm.XY((nx,ny))
  X = (X-nx/2)*dx
  Y = (Y-ny/2)*dy
  Hertz = np.square(X)/(2*rXhertz) + np.square(Y)/(2*rYhertz)
  gap = Hertz - uz

  # export zoomed versions of fields
  if do_export and (zoom_factor > 1):
    Pdiss = zoom(Pdiss, zoom_factor)
    mask = zoom(mask, zoom_factor)
    pz = zoom(pz, zoom_factor)
    gap = zoom(gap, zoom_factor)
    ux = zoom(ux, zoom_factor)
    uy = zoom(uy, zoom_factor)
    vx_loc = zoom(vx_loc, zoom_factor)
    vy_loc = zoom(vy_loc, zoom_factor)
    gap[mask] = 0

    cm.dumpConfig(Pdiss, path+"Pdiss_norm.dat", Lx=Pdiss.shape[0]*dx, Ly=Pdiss.shape[1]*dy)
    #cm.dumpConfig(pz/pressEff, path+"pressZ_norm.dat", Lx=pz.shape[0]*dx, Ly=pz.shape[1]*dy)
    cm.dumpConfig(pz, path+"pressZ.dat", Lx=pz.shape[0]*dx, Ly=pz.shape[1]*dy)
    cm.dumpConfig([vx_loc/vExt,vy_loc/vExt], path+"vel_loc_norm.dat", Lx=pz.shape[0]*dx, Ly=pz.shape[1]*dy)
    cm.dumpConfig([ux,uy], path+"ux_uy.dat", Lx=pz.shape[0]*dx, Ly=pz.shape[1]*dy)
    cm.dumpConfig(gap, path+"gap_zoom.dat", Lx=Pdiss.shape[0]*dx, Ly=Pdiss.shape[1]*dy)
    #cm.dumpConfig([vx_loc+vX,vy_loc+vY], path+"vel_tot.dat", Lx=pz.shape[0]*dx, Ly=pz.shape[1]*dy)
    #plt.savefig(path+"Pdiss.pdf")

  # plot for checks
  if do_check:
    cm.plotImg(Pdiss, rot=True)
    cm.plotImg(mask, rot=True)



if __name__ == '__main__':
  process(inpath)