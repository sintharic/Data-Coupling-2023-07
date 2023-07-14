#inpath = "/Users/christian/sim/coupled/H0p8-n0p30-fix/press0p02/vX/frict1p0/"
#inpath = "/Users/christian/sim/coupled/H0p8-n0p30-fix/press0p02/frict0p0/"
#inpath = "/Users/christian/sim/coupled/H0p8-n0p49-fix/press0p05/vX/frict1p0/"
#inpath = "/Users/christian/sim/coupled/H0p8-n0p49-fix/press0p05/frict0p0/"

inpath = "../friction/rough/material/press2p0e-02/"

do_rot = False
do_flip = False
fake_mu = 0.0



import os, sys
import numpy as np
import cmutils as cm
import cmparams as cp



def diffX(uz):
  return ( np.diff(uz, append=uz[:1,:], axis=0) + np.diff(uz, prepend=uz[-1:,:], axis=0) ) / 2

def diffY(uz):
  return ( np.diff(uz, append=uz[:,:1], axis=1) + np.diff(uz, prepend=uz[:,-1:], axis=1) ) / 2



def process(inpath):
  if inpath[-1] != os.sep: inpath += os.sep
  
  params = None
  lengthX = 0; lengthY = 0
  nu = 0.5; thick = 0
  F = 0; E = 0
  R = 0
  mu = 0

  konfig = ""
  desc = ""

  x = np.zeros(1)
  y = np.zeros(1)
  sxx = np.zeros(1)
  sxy = np.zeros(1)
  sxz = np.zeros(1)
  syy = np.zeros(1)
  syz = np.zeros(1)
  szz = np.zeros(1)
  eig0 = np.zeros(1)
  eig1 = np.zeros(1)
  eig2 = np.zeros(1)
  VonMises = np.zeros(1)

  # read params from params file
  params = cp.read(inpath + "params.out")
  nx = params.nxGlobal; ny = params.nyGlobal
  lengthX = params.lengthX; lengthY = params.lengthY
  area = params.lengthX * params.lengthY
  dx, dy = lengthX/nx, lengthY/ny
  R = params.SHEET[0].rXhertz
  if params.SHEET[0].fOnSitePotential:
    konfig = "konfig0Dzyx.dat"
    nu = params.SHEET[0].poisson0; nu = nu
    F = params.SHEET[0].pressInit * area
    E = 2*params.SHEET[0].stiffness0; E = E
    thick = params.SHEET[0].thickness0; thick = thick
    mu = params.SHEET[0].frictionCoeffOS; mu = mu
  else:
    konfig = "konfig1Dzyx.dat"; konfig = konfig
    nu = params.SHEET[1].poisson0; nu = nu
    F = params.SHEET[1].pressInit * area
    E = 2*params.SHEET[1].stiffness0; E = E
    thick = params.SHEET[1].thickness0; thick = thick
    mu = params.INTER[0].frictionCoeff; mu = mu

  # read konfig file
  data = cm.readMulti(inpath+konfig)
  if do_rot: data = cm.rot270Multi(data)
  if do_flip: data = cm.flipMulti(data, axis=1)
  [uz, szz, uy, syz, ux, sxz] = data

  # apply engineering convention
  szz = -szz
  syz = -syz
  sxz = -sxz

  # add fake friction
  if fake_mu:
    print("Applying simple fake friction.")
    konfig = "FAKE_" + konfig
    sxz = - fake_mu * szz


  # calculate slopes for strains
  dux_dx = diffX(ux)/dx
  dux_dy = diffY(ux)/dy
  duy_dx = diffX(uy)/dx
  duy_dy = diffY(uy)/dy
  if do_flip:
    dux_dx = -dux_dx
    dux_dy = -dux_dy
    duy_dx = -duy_dx
    duy_dy = -duy_dy

  # calculate missing components of the stress tensor
  sxx = dux_dx*E*(1 + nu*(1-2*nu))/2 + duy_dy*E*nu*(1-nu) + 2*nu*szz
  sxy = E*(1-nu)*(dux_dy + duy_dx)/4
  syy = duy_dy*E*(1 + nu*(1-2*nu))/2 + dux_dx*E*nu*(1-nu) + 2*nu*szz

  # calculate stress eigenvalues
  eig0 = np.zeros_like(sxx)
  eig1 = np.zeros_like(sxx)
  eig2 = np.zeros_like(sxx)
  for ix in range(nx):
    for iy in range(ny):
      M = np.array([[sxx[ix,iy], sxy[ix,iy], sxz[ix,iy]], [sxy[ix,iy], syy[ix,iy], syz[ix,iy]], [sxz[ix,iy], syz[ix,iy], szz[ix,iy]]])
      eig0[ix,iy], eig1[ix,iy], eig2[ix,iy] = sorted(np.linalg.eigvals(M))
  #VonMises = (np.square(sxx-syy) + np.square(sxx-szz) + np.square(syy-szz))/2
  #VonMises += 3*(np.square(sxy) + np.square(sxz) + np.square(syz))
  #VonMises = np.sqrt(VonMises)
  VonMises = np.square(eig0-eig1) + np.square(eig0-eig2) + np.square(eig1-eig2)
  VonMises = np.sqrt(VonMises/2)


  filename = inpath+konfig[:-4]
  if do_flip: filename = inpath+"FLIPPED"+konfig[:-4]
  cm.dumpConfig(VonMises, filename+"_VonMises.dat", Lx=lengthX, Ly=lengthY)
  cm.dumpConfig(eig2, filename+"_eig2.dat", Lx=lengthX, Ly=lengthY)
  cm.dumpConfig(eig1, filename+"_eig1.dat", Lx=lengthX, Ly=lengthY)
  cm.dumpConfig(eig0, filename+"_eig0.dat", Lx=lengthX, Ly=lengthY)



if __name__ == '__main__': 
  if len(sys.argv) > 1:
    newpath = sys.argv[1]
    if os.path.isdir(newpath): inpath = newpath
    else: print("ERROR: invalid path. Using "+inpath)

  process(inpath)