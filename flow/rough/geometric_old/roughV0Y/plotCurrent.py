params_mu0 = "/Users/christian/sim/coupled/H0p8-n0p49-fix/press0p05/frict0p0/params.out"
params_mu1 = "/Users/christian/sim/coupled/H0p8-n0p49-fix/press0p05/frict1p0/params.out"

path_mu0 = "/Users/christian/sim/flow/roughV0Y/mu0/currentCalc.x_01024.current.000000.dat"
path_mu0R = "/Users/christian/sim/flow/roughV0Y/mu0_rot90/currentCalc.x_01024.current.000000.dat"
path_mu1 = "/Users/christian/sim/flow/roughV0Y/mu1/currentCalc.x_01024.current.000000.dat"
path_mu1R = "/Users/christian/sim/flow/roughV0Y/mu1_rot90/currentCalc.x_01024.current.000000.dat"
path_p0 = "/Users/christian/sim/flow/roughV0Y/p0/currentCalc.x_01024.current.000000.dat"
path_p0R = "/Users/christian/sim/flow/roughV0Y/p0_rot90/currentCalc.x_01024.current.000000.dat"

CLIM = (0, 1.52e-9)

import numpy as np
import cmutils as cm
import cmparams as cp
import os

def readFlow(path, shape):
  with open(path) as file:
    lines = [line for line in file.readlines() if line[0].isalnum()]
  array = np.array([float(line.split()[2]) for line in lines])
  return array.reshape(shape)

sim = cp.read(params_mu1)
dist = (sim.nTime//sim.SHEET[0].frictRelax) * sim.dTime * sim.SHEET[0].vY
rel_shift = dist/sim.lengthY - int(dist/sim.lengthY)
idxShift = round(rel_shift*sim.nyGlobal)//2 # reduced from 2048 to 1024!

def shift(array, xshift=0, yshift=0):
  data = np.copy(array)
  if yshift>0:
    data[:,:-yshift] = array[:,yshift:]
    data[:,-yshift:] = array[:,:yshift]
  if xshift>0:
    data_copy = np.copy(data)
    data[:-xshift,:] = data_copy[xshift:,:]
    data[-xshift:,:] = data_copy[:xshift,:]
  return data

flow_mu0 = readFlow(path_mu0, (1024,1024))
flow_mu0 = shift(flow_mu0,0,0)
flow_mu1 = readFlow(path_mu1, (1024,1024))
#flow_mu1[flow_mu1==0] = np.nan
#flow_mu1 = cm.smoothPBC(flow_mu1, 42)
flow_mu1 = shift(flow_mu1,0,415)
#flow_p0 = readFlow(path_p0, (1024,1024))
#flow_p0 = cm.rot270(flow_p0)

flow_mu0R = readFlow(path_mu0R, (1024,1024))
flow_mu0R = cm.rot270(flow_mu0R)
flow_mu1R = readFlow(path_mu1R, (1024,1024))
flow_mu1R = cm.smoothPBC(flow_mu1R, 42)
flow_mu1R = cm.rot270(shift(flow_mu1R,610,0))
#flow_mu1R[flow_mu1R==0] = np.nan
#flow_p0R = readFlow(path_p0R, (1024,1024))


cm.plotImg(flow_mu0, rot=True, clim=CLIM, title="flow in x")
cm.plotImg(flow_mu0R, rot=True, clim=CLIM, title="flow in y")

total_mu0 = flow_mu0.mean()
total_mu1 = flow_mu1.mean()
title = "flow in x, frict. in y (%.2f)" % (total_mu1/total_mu0)
cm.plotImg(flow_mu1, rot=True, clim=CLIM, title=title)
print(title, ":", total_mu1/total_mu0)

total_mu0R = flow_mu0R.mean()
total_mu1R = flow_mu1R.mean()
title = "flow in y, frict. in y (%.2f)" % (total_mu1R/total_mu0R)
cm.plotImg(flow_mu1R, rot=True, clim=CLIM, title=title)
print(title, ":", total_mu1R/total_mu0R)
#cm.plotImg(flow_mu0R, rot=True)
#cm.plotImg(flow_mu1R, rot=True)
#cm.plotImg(flow_p0, rot=True)
#cm.plotImg(flow_p0R, rot=True)

#DEBUG
#print(os.path.isfile(params_mu0))
#print(os.path.isfile(params_mu1))
#print(os.path.isfile(path_mu0))
#print(os.path.isfile(path_mu0R))
#print(os.path.isfile(path_mu1))
#print(os.path.isfile(path_mu1R))
#print(os.path.isfile(path_p0))
#print(os.path.isfile(path_p0R))
