path = "../../friction/wave/geometric/frict0p0/"
#path = "/Users/christian/sim/coupled/squareSWR/press0p50/frict1p0"


import numpy as np
import cmparams as cp
import cmutils as cm

def expand(array, n=8):
  result = np.zeros((array.shape[0]+2*n, array.shape[1]+2*n))
  
  # center
  result[n:-n, n:-n] = array
  
  # edges
  result[:n, n:-n] = array[-n:, :]
  result[n:-n, :n] = array[:, -n:]
  result[-n:, n:-n] = array[:n, :]
  result[n:-n, -n:] = array[:, :n]
  
  # corners
  result[:n, :n] = array[-n:, -n:]
  result[:n, -n:] = array[-n:, :n]
  result[-n:, :n] = array[:n, -n:]
  result[-n:, -n:] = array[:n, :n]

  return result



if path[-1]!="/": path += "/"
sim = cp.read(path+"params.in")
dat0 = cm.readConfig(path+sim.sheet[0].konfigName)
dat1 = cm.readConfig(path+sim.sheet[1].konfigName)

dat0M = np.copy(dat0)
refpos = np.maximum(dat0[dat0.shape[0]//2,0], dat1)
dat0M[dat0M<refpos] = refpos[dat0M<refpos]

dat0M = expand(dat0M, 768)[:1536,:1536]
dat1M = expand(dat1, 768)[:1536,:1536]

# plot both surfaces
#ax = cm.plotSurf(dat1M, kind="color", stride=16)
ax = cm.plotSurf(dat1M, kind="color", cmap="inferno", stride=8)
#ax = cm.plotSurf(expand(dat0, 768)[:1536,:1536]+0.001, stride=16, alpha=0.3, kind="color", cmap="gray", axis=ax)
ax.view_init(elev=30, azim=-90)
ax.set_zlim(0,0.1)

# plot only the (flipped) rough surface
ax = cm.plotSurf(expand(dat0,512), alpha=0.8, stride=16, cmap="turbo")
ax.view_init(elev=30, azim=-90)
ax.set_zlim(-0.3, 0.5)
