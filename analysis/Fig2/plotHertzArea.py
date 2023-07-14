path_mat0p0 = "../results/hertz/res-material-frict0p0_red.dat"
path_mat1p0 = "../results/hertz/res-material-frict1p0_red.dat"
path_geom0p0 = "../results/hertz/res-geometric2-frict0p0_red.dat"
path_geom1p0 = "../results/hertz/res-geometric2-frict1p0_red.dat"

#figsize = (1.2, 1.2)
figsize = (1.8,1.2)
outpath = "HertzContArea.pdf"


import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import ticker
import mpl 
mpl.RevTex(2)
mpl.grid(False)


from scipy.optimize import curve_fit as fit
def func(x,a): return a*np.power(x,1/3)

class result():
  valid = np.zeros(1)
  thickness = np.zeros(1)
  force = np.zeros(1)
  frictCoeff = 0
  absContArea = 0
  geomForce = np.zeros(1)
  meanGap = np.zeros(1)
  meanPosGap = np.zeros(1)

  force = np.zeros(1)
  delta_mu = np.zeros(1)
  popt, pcov = (np.zeros(1), np.zeros(1))

  def __init__(self, path, do_filter=False):
    data = np.loadtxt(path)
    
    self.thickness = data[:,1]
    self.frictCoeff = data[:,2]
    self.force = data[:,3]
    #self.geomForce = data[:,4] #geomForceY
    self.geomForce = data[:,5] #geomForceX
    self.absContArea = data[:,6]

    # quality of data can be determined by the overlap of the surfaces
    if do_filter:
      self.meanGap = data[:,7]
      self.meanPosGap = data[:,8]
      self.valid = self.meanGap > 0.8*self.meanPosGap

      self.thickness = self.thickness[valid]
      self.frictCoeff = self.frictCoeff[valid]
      self.force = self.force[valid]
      self.geomForce = self.geomForce[valid]
      self.absContArea = self.absContArea[valid]
    else:
      self.valid = np.ones(data.shape[0], dtype=bool)

    self.delta_mu = self.geomForce/self.force/self.frictCoeff
    if (self.frictCoeff > 0).sum() == len(self.frictCoeff):
      popt, pcov = fit(func, self.force[-4:], self.delta_mu[-4:])
      (self.popt, self.pcov) = (popt, pcov)
    
    


geom0p0 = result(path_geom0p0)
geom1p0 = result(path_geom1p0)
mat0p0 = result(path_mat0p0)
mat1p0 = result(path_mat1p0)


# set up broken axis plot
fig, ax = plt.subplots()
ax.set_xlim(0.0,0.07)
ax.set_ylim(-0.00,0.021)
ax.set_xlabel(r"$F_\mathrm{z} / E^\mathrm{*} R^2$")
ax.set_ylabel(r"$\Delta A_\mathrm{c} / A_\mathrm{c}(\mu\!=\!0)$")
# pmat, = ax.plot(mat1p0.force, np.sqrt(mat1p0.absContArea/mat0p0.absContArea)-1, "--", label=r"$\nu=0.25$", color=mpl.color(0))
# pgeom, = ax.plot(geom1p0.force, np.sqrt(geom1p0.absContArea/geom0p0.absContArea)-1, label=r"$h=0.2 R$", color=mpl.color(0))
pmat, = ax.plot(mat1p0.force, mat1p0.absContArea/mat0p0.absContArea-1, "--", label=r"$\nu=0.25$", color=mpl.color(0))
pgeom, = ax.plot(geom1p0.force, geom1p0.absContArea/geom0p0.absContArea-1, label=r"$h=0.2 R$", color=mpl.color(0))
plt.legend(loc="upper left")
#lmat = plt.legend([pmat], [pmat.get_label()], loc="upper left")
#lgeom = plt.legend([pgeom], [pgeom.get_label()], loc="lower right")
#ax.add_artist(lmat)

mpl.format(fig, figsize, right=0.05, yticl=1.0)
mpl.label(fig, "b)")
mpl.digitalize(fig, "HertzContArea")
print("figsize (in): ", fig.get_size_inches())

if outpath: plt.savefig(outpath)