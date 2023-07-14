import numpy as np
import cmutils as cm

def plot_constrict(path="."):
  if path[-1]!="/": path+="/"
  press = cm.readConfig(path+"konfig1D.dat",3)
  nx, ny = press.shape
  mask = press > 1e-8
  img = np.zeros_like(mask)
  img[:nx//2,:] = mask[nx//2:,:]
  img[nx//2:,:] = mask[:nx//2,:]
  return cm.plotImg(img[3*nx//8:5*nx//8, 3*ny//8:5*ny//8])
