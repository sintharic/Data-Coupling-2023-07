path_mat1p0 = "../../friction/rough/material/press2p0e-02/frict1p0/"
path_mat0p0 = "../../friction/rough/material/press2p0e-02/frict0p0/"
path_geom1p0 = "../../friction/rough/geometric1/press4p0e-02/frict1p0/"
path_geom0p0 = "../../friction/rough/geometric1/press4p0e-02/frict0p0/"
fake_mu = 1.0

do_flip = False
box_props = dict(facecolor='white', edgecolor='none')
cb_props = dict(aspect=10, shrink=0.85, pad=0.02, fraction=0.2)

CLIM = [(0,0.25), (0,1.7)]#[(0,1.3), (0,0.26)]
ylabel = [r"$\sigma_{\alpha\beta} / p_\mathrm{\scriptsize H}$", r"$\sigma / E^\mathrm{*}$"]


import numpy as np
import cmutils as cm
import cmparams as cp
import matplotlib.pyplot as plt
import mpl
#import reader
mpl.RevTex(2)
mpl.set_fontsize(9)
mpl.grid(False)
mpl.rcp["text.latex.preamble"] = "\\usepackage{bm}\n\\usepackage{amsmath}"



class simResult:
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

  def __init__(self, path, add_mu=0):

    # read params from params file
    self.params = cp.read(path + "params.out")
    self.nx = self.params.nxGlobal; self.ny = self.params.nyGlobal
    self.lengthX = self.params.lengthX; self.lengthY = self.params.lengthY
    area = self.params.lengthX * self.params.lengthY
    dx, dy = self.lengthX/self.nx, self.lengthY/self.ny
    self.R = self.params.SHEET[0].rXhertz
    if self.params.SHEET[0].fOnSitePotential:
      self.konfig = "konfig0Dzyx.dat"
      nu = self.params.SHEET[0].poisson0; self.nu = nu
      self.F = self.params.SHEET[0].pressInit * area
      E = 2*self.params.SHEET[0].stiffness0; self.E = E
      thick = self.params.SHEET[0].thickness0; self.thick = thick
      mu = self.params.SHEET[0].frictionCoeffOS; self.mu = mu
    else:
      konfig = "konfig1Dzyx.dat"; self.konfig = konfig
      nu = self.params.SHEET[1].poisson0; self.nu = nu
      self.F = self.params.SHEET[1].pressInit * area
      E = 2*self.params.SHEET[1].stiffness0; self.E = E
      thick = self.params.SHEET[1].thickness0; self.thick = thick
      mu = self.params.INTER[0].frictionCoeff; self.mu = mu
    
    self.konfig = konfig
    if add_mu: konfig = "FAKE_" + konfig
    if do_flip: konfig = "FLIPPED" + konfig

    # read konfig file
    #self.eig0 = cm.readConfig(path+konfig[:-4]+"_eig0.dat") # not needed
    #self.eig1 = cm.readConfig(path+konfig[:-4]+"_eig1.dat") # not needed
    self.eig2 = cm.readConfig(path+konfig[:-4]+"_eig2.dat")
    self.VonMises = cm.readConfig(path+konfig[:-4]+"_VonMises.dat")


    # update description
    self.desc = "$\\mu = %.1f$\\phantom{1.}\n" % mu
    if thick >= 1: self.desc += "$h \\gg R$\\phantom{0.1}\n"
    else: self.desc += "$h = %.1f R$\n" % thick
    self.desc += "$\\nu = %.2f$\\phantom{.}" % nu

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



geom0p0 = simResult(path_geom0p0, fake_mu)
geom1p0 = simResult(path_geom1p0, 0)
mat0p0 = simResult(path_mat0p0, fake_mu)
mat1p0 = simResult(path_mat1p0, 0)


# plot with friction
fig = plt.figure()
mpl.format(fig, boxsize=(6.3,1.5), left=-0.08, right=-0.02, bottom=0.0, top=0.0, xlabel=0, yticl=0.0, xticl=0)
spec = plt.GridSpec(1,4, width_ratios=[1.0,1.0,1.0,1.28])


ax = fig.add_subplot(spec[0])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(4*mat0p0.eig2.transpose(), cmap="turbo", clim=CLIM[0], origin="lower")
ax.text(0.98, 0.02, "(x4)", color="w", va="bottom", ha="right", transform=ax.transAxes)
mpl.label(ax, "\\Large{a)}", xoffset=0.0, yoffset=0.0, color="w")#, bbox=box_props)

ax = fig.add_subplot(spec[1])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(mat1p0.eig2.transpose(), cmap="turbo", clim=CLIM[0], origin="lower")
mpl.label(ax, "\\Large{b)}", xoffset=0.0, yoffset=0.0, color="w")#, bbox=box_props)
#cb = fig.colorbar(im, ax=ax, **cb_props)
#cb.ax.minorticks_on()


ax = fig.add_subplot(spec[2])
ax.set_xticks(()); ax.set_yticks(())
#im = ax.imshow(shift(geom0p0.VonMises.transpose(), 0, 0), cmap="turbo", clim=CLIM[1], origin="lower")
im = ax.imshow(3*geom0p0.eig2.transpose(), cmap="turbo", clim=CLIM[0], origin="lower")
ax.text(0.98, 0.02, "(x3)", color="w", va="bottom", ha="right", transform=ax.transAxes)
mpl.label(ax, "\\Large{c)}", xoffset=0.0, yoffset=0.0, color="w")#, bbox=box_props)

ax = fig.add_subplot(spec[3])
ax.set_xticks(()); ax.set_yticks(())
#im = ax.imshow(shift(geom1p0.VonMises.transpose(), 0, 835), cmap="turbo", clim=CLIM[1], origin="lower")
im = ax.imshow(geom1p0.eig2.transpose(), cmap="turbo", clim=CLIM[0], origin="lower")
mpl.label(ax, "\\Large{d)}", xoffset=0.0, yoffset=0.0, color="w")#, bbox=box_props)
cb = fig.colorbar(im, ax=ax, **cb_props)
cb.ax.minorticks_on()

mpl.digitalize(fig, "RoughStresses")
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.savefig("RoughStresses.pdf")