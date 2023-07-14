# Hertzian
path_geom1p0 = "../results/hertz/res-geometric2-frict1p0.dat"
path_mat1p0 = "../results/hertz/res-material-frict1p0.dat"

do_percentages = False # on y axis
#figsize = (1.2, 1.2)
figsize = (1.8,1.2)

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import ticker
import mpl 
mpl.RevTex(2)
#mpl.set_fontsize(9)
mpl.grid(False)

color = [mpl.color(1), mpl.color(2), mpl.color(0)]


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
    
    


#geom0p0 = result(path_geom0p0)
geom1p0 = result(path_geom1p0)
#mat0p0 = result(path_mat0p0)
mat1p0 = result(path_mat1p0)

# set up broken axis plot (from https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/broken_axis.html)
fig, (ax_geom, ax_mat) = plt.subplots(2, 1, sharex=True)
ax_geom.semilogx(); ax_geom.set_yscale("symlog", linthresh=1e-4); #mpl.setnonsci(ax_geom, "y")
ax_mat.semilogx(); ax_mat.set_yscale("symlog", linthresh=1e-4); #mpl.setnonsci(ax_mat, "y")
ax_geom.spines['bottom'].set_visible(False)
ax_mat.spines['top'].set_visible(False)
ax_geom.xaxis.tick_top()
ax_geom.tick_params(labeltop=False)  # don't put tick labels at the top
ax_mat.xaxis.tick_bottom()
locmaj = ticker.LogLocator(base=10, numticks=12)
ax_geom.xaxis.set_major_locator(locmaj)
locmin = ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8), numticks=12)
ax_geom.xaxis.set_minor_locator(locmin)
ax_geom.xaxis.set_minor_formatter(ticker.NullFormatter())
d = .020  # how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax_geom.transAxes, color='k', clip_on=False, ls="-", marker="None")
ax_geom.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax_geom.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax_mat.transAxes)  # switch to the bottom axes
ax_mat.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax_mat.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


# format x axis
ax_mat.set_xlabel(r"$F_\mathrm{z} / E^\mathrm{*} R^2$")
ax_mat.set_xlim(1e-4,5e-2)
ax_geom.set_xlim(1e-4,5e-2)

# format y axis with percentages
if do_percentages:
  ax_geom.set_yticks([1, 2, 4, 6, 8, 10])
  ax_geom.set_ylim(0.2, 20)
  ax_mat.set_yticks([-1, -2, -4, -6, -8, -10])
  ax_mat.set_ylim(-10, -0.4)
  ax_geom.text(-0.20, -0.05, r"$\Delta \mu/\mu_\mathrm{c}$ (\%)", rotation="vertical", transform=ax_geom.transAxes, va="center", ha="right")
# format y axis with absolute
else:
  ax_geom.set_yticks([0.0004, 0.0006, 0.0008, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10], minor=True)
  ax_geom.set_ylim(0.00015, 0.20)
  ax_mat.set_yticks([-0.01, -0.02, -0.04, -0.06, -0.08, -0.10], minor=True)
  ax_mat.set_ylim(-0.11, -0.005)
  ax_geom.text(-0.24, 0.00, r"$\Delta \mu/\mu_\mathrm{c}$", rotation="vertical", transform=ax_geom.transAxes, va="center", ha="right")


# plot data
ax_geom.plot(geom1p0.force, geom1p0.delta_mu, ls=mpl.line(0), color=color[2])#, marker=mpl.marker(5))
ax_geom.text(3e-4, 0.015, "$h = 0.2 R$", va = "bottom", ha="left")
ax_mat.plot(mat1p0.force[1:], mat1p0.delta_mu[1:], ls=mpl.line(1), color=color[2])#, marker=mpl.marker(5))
ax_mat.text(4e-4, -0.03, "$\\nu = 0.25$", va = "top", ha="left")

# plot fit functions
#ax_geom.text(1e-2, 3.4, r"$\propto + F_\mathrm{z}^{1/3}$", va="top", ha="left", fontsize=9)
ax_mat.text(4.5e-2, -2.5e-2, r"$\propto F_\mathrm{z}^{1/3}$", va="bottom", ha="right", fontsize=9)

# format figure
mpl.format(fig, figsize, right=0.05, yticl=1.0)
mpl.label(fig, "d)")
plt.subplots_adjust(hspace=0.05)
mpl.digitalize(fig, "HertzGeomFrict")
print("figsize (in): ", fig.get_size_inches())

plt.savefig("HertzGeomFrict.pdf")