plotDots = True
ROWS = [0,1,2,3]
shrink = 2
cols = [2,3,2,2]
do_circles = True
do_rot = True
do_flip = False
do_relative = True
do_none = True
outpath = "Hertz_dissip.pdf"

forceZ = 0.01
rHertz = 1.
E = 1.

pathsNone = ["../../friction/hertz/none/frict1p0/pressZ.dat",
"../../friction/hertz/none/frict1p0/vel_loc_norm.dat",
"../../friction/hertz/none/frict1p0/vel_loc_norm.dat",
"../../friction/hertz/none/frict1p0/konfig0Dzyx_eig2_zoom.dat"]


pathsMat = ["../../friction/hertz/material/frict1p0/pressZ.dat",
"../../friction/hertz/material/frict1p0/vel_loc_norm.dat",
"../../friction/hertz/material/frict1p0/vel_loc_norm.dat",
"../../friction/hertz/material/frict1p0/konfig0Dzyx_eig2_zoom.dat"]


pathsGeom = ["../../friction/hertz/geometric2/frict1p0/pressZ.dat",
"../../friction/hertz/geometric2/frict1p0/vel_loc_norm.dat",
"../../friction/hertz/geometric2/frict1p0/vel_loc_norm.dat",
"../../friction/hertz/geometric2/frict1p0/konfig0Dzyx_eig2_zoom.dat"]

paramsGeom = "../../friction/hertz/geometric2/frict1p0/params.out"


CLIM = [(0,1.26), (-0.45,0.45), (-0.055,0.055), (-0.14,0.14)]
if do_relative: CLIM[1] = (0.7, 1.3)#(-1.26, -0.74)
ylabel = [r"$p_\mathrm{z}/p_\mathrm{H}$",
          r"$\dot{u}_\mathrm{l} / v_0$",
          r"$\dot{u}_\mathrm{t} / v_0$",
          r"$\sigma_\mathrm{I}^\mathrm{max} / E^\mathrm{*}$"]
if do_relative:
  ylabel[1] = r"$v_\mathrm{l}^\mathrm{rel} / v_0$"



import cmutils as cm
cm.WARN = False
import cmparams as cp
import numpy as np
import matplotlib.pyplot as plt
import mpl 
mpl.RevTex(2)
from matplotlib.colors import TwoSlopeNorm
#mpl.set_fontsize(10)


aHertz = (0.75 * forceZ * rHertz / E)**(1./3)
simGeom = cp.read(paramsGeom)
thickness = simGeom.SHEET[0].thickness0

def better_ticks(axlim, n=5, delta=40):
  fac = delta//n
  if delta%n != 0:
    fac = delta//n + 1
    delta = fac*n
  dv = float("%.0e" % (axlim[1]-axlim[0]))/delta
  if 0 in axlim: ticks = dv*np.arange(n)*fac
  else: ticks = float("%.1e" % (0.5*(axlim[0]+axlim[1]))) + dv*np.linspace(-fac//2,fac//2,n)

def read(path, iPath):
  data = cm.readConfig(path,cols[iPath])
  if iPath == 1:
    if do_rot: data = cm.rot270(data)
    if do_flip: data = cm.flip(data, axis=1, direction=0)
  elif iPath == 2:
    if do_rot: data = -cm.rot270(data)
    if do_flip: data = cm.flip(data, axis=1, direction=1)
  else: 
    if do_rot: data = cm.rot270(data)
    if do_flip: data = cm.flip(data, axis=1)

  return data



# Read data
datNone = []
lengthY = []
for iPath,path in enumerate(pathsNone):
  data = read(path, iPath)
  datNone.append(data)
  lengthY.append(cm.YLIM[1]-cm.YLIM[0])
mask = datNone[0] > datNone[0].mean()/2 #adjust
noneCenter = cm.center(mask)
noneRadius = np.sqrt(mask.sum()/np.pi)
#noneLine = plt.Circle((noneCenter[1], noneCenter[0]), noneRadius, color='w', fill=False, linewidth=0.5)


datGeom = []
for iPath,path in enumerate(pathsGeom):
  data = read(path, iPath)
  datGeom.append(data)

datMat = []
for iPath,path in enumerate(pathsMat):
  data = read(path, iPath)
  datMat.append(data)

# normalize pressure
datGeom[0] = datGeom[0] / datNone[0].max()
datMat[0] = datMat[0] / datNone[0].max()
datNone[0] = datNone[0] / datNone[0].max()

# normalize gap
#datGeom[1] = (datGeom[1]-datNone[1]) / datNone[1]
#datMat[1] = (datMat[1]-datNone[1]) / datNone[1]
#datNone[1] = datNone[1]/datNone[1].max()

# calculate relative velocities
if do_relative:
  datNone[1] = 1 - datNone[1]
  datGeom[1] = 1 - datGeom[1]
  datMat[1] = 1 - datMat[1]

# plot
fig = plt.figure()
if do_none: 
  mpl.format(fig, boxsize=(3.3*0.743, 0.743*len(ROWS)), left=0.02, right=0.35, bottom=0.0, top=0.20, xlabel=0, xticl=0, ylabel=1.1, yticl=0)
  spec = plt.GridSpec(len(ROWS), 3, width_ratios=[1.03,1,1.25])
else: 
  mpl.format(fig, boxsize=(2,1.0*len(ROWS)), left=0.02, right=0.32, bottom=-0.1, top=0.1, xlabel=0, yticl=0.1, xticl=0)
  spec = plt.GridSpec(len(ROWS), 2, width_ratios=[1,1.25])
print("figure size (inches): %.2f, %.2f" % tuple(fig.get_size_inches()))
for iRow in ROWS:
  cmap_norm = None
  if iRow==3: cmap_norm = TwoSlopeNorm(0)
  
  # no coupling
  if do_none:
    axNone = fig.add_subplot(spec[iRow,0])
    axNone.set_xticks(())
    axNone.set_yticks(())
    axNone.label_outer()
    if iRow==0: axNone.set_title("reference")
    data = datNone[iRow].transpose()
    shape = data.shape
    # zoom in in some cases
    if iRow != 666:
      xstart = shape[0]//8
      xend = int(7*shape[0]/8)
      ystart = shape[1]//8
      yend = int(7*shape[1]/8)
      posFac = 0.75
    else:
      xstart = ystart = 0
      xend = yend = -1
      posFac = 1
    im = axNone.imshow(data[xstart:xend, ystart:yend], cmap="turbo", clim=CLIM[iRow], origin="lower", norm=cmap_norm)
    if do_circles: 
      noneLine = plt.Circle((noneCenter[1]*posFac, noneCenter[0]*posFac), noneRadius, color='gray', fill=False, linewidth=0.5)
      axNone.add_patch(noneLine)

    axNone.set_ylabel(ylabel[iRow])


  # material coupling
  axMat = fig.add_subplot(spec[iRow,0+do_none])
  axMat.set_xticks(())
  axMat.set_yticks(())
  axMat.label_outer()
  if iRow==0: axMat.set_title(r"$\nu = 0.25$")
  data = datMat[iRow].transpose()
  shape = data.shape
  # zoom in in some cases
  if iRow != 666:
    xstart = shape[0]//8
    xend = int(7*shape[0]/8)
    ystart = shape[1]//8
    yend = int(7*shape[1]/8)
    posFac = 0.75
  else:
    xstart = ystart = 0
    xend = yend = -1
    posFac = 1
  im = axMat.imshow(data[xstart:xend, ystart:yend], cmap="turbo", clim=CLIM[iRow], origin="lower", norm=cmap_norm)
  if do_circles: 
    noneLine = plt.Circle((noneCenter[1]*posFac, noneCenter[0]*posFac), noneRadius, color='gray', fill=False, linewidth=0.5)
    axMat.add_patch(noneLine)
  
  if not do_none: axMat.set_ylabel(ylabel[iRow])


  # geometric coupling
  axGeom = fig.add_subplot(spec[iRow,1+do_none])
  axGeom.set_xticks(())
  axGeom.set_yticks(())
  axGeom.label_outer()
  if iRow==0: axGeom.set_title("$h \\approx r_\\mathrm{H}$")
  data = datGeom[iRow].transpose()
  shape = data.shape
  # zoom in in some cases
  if iRow < 90:
    xstart = shape[0]//8
    xend = int(7*shape[0]/8)
    ystart = shape[1]//8
    yend = int(7*shape[1]/8)
    posFac = 0.75
  else:
    xstart = ystart = 0
    xend = yend = -1
    posFac = 1
  im = axGeom.imshow(data[xstart:xend, ystart:yend], cmap="turbo", clim=CLIM[iRow], origin="lower", norm=cmap_norm)
  if do_circles:
    noneLine = plt.Circle((noneCenter[1]*posFac, noneCenter[0]*posFac), noneRadius, color='gray', fill=False, linewidth=0.5)
    axGeom.add_patch(noneLine)
  
  ticks = better_ticks(CLIM[iRow], 3, 20)
  if do_none: cb = fig.colorbar(im, ax=axGeom, aspect=5, shrink=0.8)#, ticks=ticks)
  else: cb = fig.colorbar(im, ax=axGeom, aspect=8, shrink=0.65)#, ticks=ticks)
  cb.ax.minorticks_on()
  yl = cb.ax.get_yticklabels()
  #print(yl)#DEBUG
  if do_none:
    if iRow in (1,2):
      for t in yl: t.set_horizontalalignment('right'); t.set_x(3.5+0.2)
    elif iRow==3:
      for t in yl: t.set_horizontalalignment('right'); t.set_x(3.0+0.2)
    else:
      for t in yl: t.set_horizontalalignment('right'); t.set_x(2.1+0.2)
  else:
    if iRow in (0,1):
      for t in yl: t.set_horizontalalignment('right'); t.set_x(3.6+0.2)
    elif iRow==3:
      for t in yl: t.set_horizontalalignment('right'); t.set_x(2.5+0.2)
    else:
      for t in yl: t.set_horizontalalignment('right'); t.set_x(4.2+0.2)

if do_none: plt.subplots_adjust(wspace=0.01, hspace=0.04)
else: plt.subplots_adjust(wspace=0.05, hspace=-0.2)

mpl.label(fig, "a)")
mpl.digitalize(fig, "HertzDissip")

if outpath: plt.savefig(outpath)
