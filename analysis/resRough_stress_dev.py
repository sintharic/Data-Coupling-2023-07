paths = ["../friction/rough/geometric1/press4p0e-02/frict1p0/",
         "../friction/rough/geometric1/press4p0e-02/seed1_vX/",
         "../friction/rough/geometric1/press4p0e-02/seed2_vX/",
         "../friction/rough/geometric1/press4p0e-02/seed3_vX/",
         "../friction/rough/geometric1/press4p0e-02/frict1p0_vY/",
         "../friction/rough/geometric1/press4p0e-02/seed1_vY/",
         "../friction/rough/geometric1/press4p0e-02/seed2_vY/",
         "../friction/rough/geometric1/press4p0e-02/seed3_vY/"]


import os
import numpy as np
import matplotlib.pyplot as plt
import cmutils as cm
import eval_stress
import mpl
mpl.RevTex()

do_split = False
plot_hist = True
do_log = True


# generate non-existent data first
for path in paths:
  if path[-1] != os.sep: path += os.sep
  if os.path.isfile(path+"konfig1Dzyx_eig2.dat"): continue
  print("processing %s"%path)
  if "_vY" in path: eval_stress.do_rot = True
  else: eval_stress.do_rot = False
  eval_stress.process(path)


eig = []
for path in paths:
  eig.append(cm.readConfig(path+"konfig1Dzyx_eig2.dat"))

vM = []
for path in paths:
  vM.append(cm.readConfig(path+"konfig1Dzyx_VonMises.dat"))


def split4(data):
  nxH = data.shape[0]//2
  nyH = data.shape[1]//2
  return data[:nxH,:nyH], data[:nxH,nyH:], data[nxH:,:nyH], data[nxH:,nyH:]



# ----- Eigenstress ----- #

if do_split:
  ms_eig = np.zeros(4*len(eig))
  for idata,eigdata in enumerate(eig):
    ms_eig[(4*idata):4*(idata+1)] = np.array([np.square(data).mean() for data in split4(eigdata)])
else:
  ms_eig = np.array([np.square(data).mean() for data in eig])
print("%.2g +/- %.2g" % (ms_eig.mean(), ms_eig.std(ddof=1)))

if plot_hist:
  bin_edges = np.linspace(-0.2, 0.5, 51)
  bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

  plt.figure()
  if do_log: 
    plt.semilogy()
    plt.ylim(5e-7,5e-1)
  plt.ylabel("probability")
  plt.xlabel("eigenstress")
  for idata,data in enumerate(eig):
    hist,_ = np.histogram(data.flatten(), bins=bin_edges)
    if do_log: plt.bar(bin_centers, hist/hist.sum(), label="seed%i"%idata, alpha=0.3, width=0.9*(bin_edges[1]-bin_edges[0]))
    else: plt.plot(bin_centers, hist/hist.sum(), label="seed%i"%idata, alpha=0.3)
  plt.legend()
  plt.savefig("eig_deviation.png")



# ----- Von Mises ----- #

if do_split:
  ms_vM = np.zeros(4*len(vM))
  for idata,vMdata in enumerate(vM):
    ms_vM[(4*idata):4*(idata+1)] = np.array([np.square(data).mean() for data in split4(vMdata)])
else:
  ms_vM = np.array([np.square(data).mean() for data in vM])
print("%.2g +/- %.2g" % (ms_vM.mean(), ms_vM.std(ddof=1)))

if plot_hist:
  bin_edges = np.linspace(0, 0.5, 51)
  bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

  plt.figure()
  if do_log: 
    plt.semilogy()
    plt.ylim(1e-3,5e-1)
  plt.ylabel("probability")
  plt.xlabel("VonMises stress")
  for idata,data in enumerate(vM):
    hist,_ = np.histogram(data.flatten(), bins=bin_edges)
    if do_log: plt.bar(bin_centers, hist/hist.sum(), label="seed%i"%idata, alpha=0.3, width=0.9*(bin_edges[1]-bin_edges[0]))
    else: plt.plot(bin_centers, hist/hist.sum(), label="seed%i"%idata, alpha=0.3)
  plt.legend()
  plt.savefig("VonMises_deviation.png")