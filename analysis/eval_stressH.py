#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-none/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-mat/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-geom/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-geomV5/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/none/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/material/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/geometric/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/geomV5/"

#inpath = "/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-02/geometric/"
#inpath = "/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-02/material/"
inpath = "/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-02/none/"
#inpath = "/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-02/stat-mat/"
#inpath = "/Users/christian/sim/coupled/Hertz-final/Nvar/force1p0e-02/stat-geom/"

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
  if inpath[-1]!=os.sep: inpath += os.sep
  if len(sys.argv) > 1:
    newpath = sys.argv[1]
    if os.path.isdir(newpath): inpath = newpath
    else: print("ERROR: invalid path. Using "+inpath)
  if inpath[-1] != os.sep: inpath += os.sep


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
  szz = -szz
  syz = -syz
  sxz = -sxz

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

  # calculate eigenvalues along profile
  x = np.linspace(-lengthX/2, lengthX/2, nx, endpoint=False)
  y = np.linspace(-lengthY/2, lengthY/2, ny, endpoint=False)
  sxx_x = sxx[:,ny//2]; sxx_x = sxx_x
  sxy_x = sxy[:,ny//2]; sxy_x = sxy_x
  sxz_x = sxz[:,ny//2]; sxz_x = sxz_x
  syy_x = syy[:,ny//2]; syy_x = syy_x
  syz_x = syz[:,ny//2]; syz_x = syz_x
  szz_x = szz[:,ny//2]; szz_x = szz_x

  # calculate stress eigenvalues
  evecs0 = np.zeros((len(x),3))
  evecs1 = np.zeros((len(x),3))
  evecs2 = np.zeros((len(x),3))
  eig0_x = np.zeros_like(x)
  eig1_x = np.zeros_like(x)
  eig2_x = np.zeros_like(x)
  for i in range(len(x)):
    M = np.array([[sxx_x[i], sxy_x[i], sxz_x[i]], [sxy_x[i], syy_x[i], syz_x[i]], [sxz_x[i], syz_x[i], szz_x[i]]])
    vals, vecs = np.linalg.eigh(M)
    f = np.argsort(vals)
    vals = vals[f]
    vecs = vecs[:,f]
    eig0_x[i], eig1_x[i], eig2_x[i] = vals
    evecs0[i,:] = vecs[:,0]
    evecs1[i,:] = vecs[:,1]
    evecs2[i,:] = vecs[:,2]

  Mises_x = np.square(eig0_x-eig1_x) + np.square(eig0_x-eig2_x) + np.square(eig1_x-eig2_x)
  Mises_x = np.sqrt(Mises_x/2)
  #Mises_x = (np.square(sxx_x-syy_x) + np.square(sxx_x-szz_x) + np.square(syy_x-szz_x))/2
  #Mises_x += 3*(np.square(sxy_x) + np.square(sxz_x) + np.square(syz_x))
  #Mises_x = np.sqrt(Mises_x)


  # export results
  header = "x\tsigmaI(x)\tsigmaII(x)\tsigmaIII(x)"
  filename = "stressEigenvalues.dat"
  if do_flip: filename = "FLIPPED" + filename
  filename = os.path.split(inpath)[0] + os.sep + filename
  np.savetxt(filename, 
             np.vstack([x, eig0_x, eig1_x, eig2_x]).transpose(), 
             fmt="%g", header=header)

  header = "x\tvecI_x(x)\tvecI_y(x)\tvecI_z(x)\tvecII_x(x)\tvecII_y(x)\tvecII_z(x)\tvecIII_x(x)\tvecIII_y(x)\tvecIII_z(x)"
  filename = "stressEigenvects.dat"
  if do_flip: filename = "FLIPPED" + filename
  filename = os.path.split(inpath)[0] + os.sep + filename
  np.savetxt(filename, 
             np.vstack([x, evecs0[:,0], evecs0[:,1], evecs0[:,2], evecs1[:,0], evecs1[:,1], evecs1[:,2], evecs2[:,0], evecs2[:,1], evecs2[:,2]]).transpose(), 
             fmt="%g", header=header)

  filename = "stressComponents.dat"
  if do_flip: filename = "FLIPPED" + filename
  filename = os.path.split(inpath)[0] + os.sep + filename
  header = "x\tsigma_XX(x)\tsigma_XY(x)\tsigma_XZ(x)\tsigma_YY(x)\tsigma_YZ(x)\tsigma_ZZ(x)"
  np.savetxt(filename, 
             np.vstack([x, sxx_x, sxy_x, sxz_x, syy_x, syz_x, szz_x]).transpose(), 
             fmt="%g", header=header)


if __name__ == '__main__': process(inpath)