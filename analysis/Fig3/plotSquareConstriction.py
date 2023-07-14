import numpy as np
import matplotlib.pyplot as plt
import mpl
import cmutils as cm

mat0p0 = cm.readConfig("../../friction/wave/material/frict0p0/konfig1D.dat", usecols=3)
mat1p0 = cm.readConfig("../../friction/wave/material/frict1p0/konfig1Dzyx.dat", usecols=3)
geom0p0 = cm.readConfig("../../friction/wave/geometric/frict0p0/konfig1D.dat", usecols=3)
geom1p0 = cm.readConfig("../../friction/wave/geometric/frict1p0/konfig1Dzyx.dat", usecols=3)

scale_bar = dict(arrowstyle="|-|, widthA=0.25, widthB=0.25")

outpath = "Hertz"
cb_props = dict(aspect=10, shrink=0.85, pad=0.02, fraction=0.2)
color_intensity = 0.65


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

mat0p0[mat0p0 >= 1e-6] = 1
mat1p0[mat1p0 >= 1e-6] = 1
geom0p0[geom0p0 >= 1e-6] = 1
geom1p0[geom1p0 >= 1e-6] = 1
mat0p0[mat0p0 < 1e-6] = np.nan
mat1p0[mat1p0 < 1e-6] = np.nan
geom0p0[geom0p0 < 1e-6] = np.nan
geom1p0[geom1p0 < 1e-6] = np.nan


fig = plt.figure()
mpl.format(fig, boxsize=(1.6,1.6), left=0.0, right=0.0, bottom=0.0, top=0.0, xlabel=0, ylabel=0, yticl=0.0, xticl=0)
spec = plt.GridSpec(2,2, width_ratios=[1.0,1.0])

ax = fig.add_subplot(spec[0,0])
ax.set_xticks(()); ax.set_yticks(())
#im = ax.imshow(-color_intensity*shift(mat0p0,0,512)[384:640,384:640].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.9)
#im = ax.imshow(+color_intensity*shift(mat1p0,0,512)[384:640,384:640].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.6)
#ax.text(128, 32, r"$\lambda/8$", va="bottom", ha="center", fontsize=9)
#ax.annotate("", (64, 24), (192, 24), arrowprops=scale_bar, color="k")
im = ax.imshow(-color_intensity*shift(mat0p0,0,512)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.9)
im = ax.imshow(+color_intensity*shift(mat1p0,0,512)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.6)
mpl.label(ax,"j)", color="k")

ax = fig.add_subplot(spec[1,0])
ax.set_xticks(()); ax.set_yticks(())
#im = ax.imshow(-color_intensity*shift(mat0p0,512,0)[384:640,384:640].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.9)
#im = ax.imshow(+color_intensity*shift(mat1p0,512,0)[384:640,384:640].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.6)
#ax.text(128, 32, r"$\lambda/8$", va="bottom", ha="center", fontsize=9)
#ax.annotate("", (64, 24), (192, 24), arrowprops=scale_bar, color="k")
im = ax.imshow(-color_intensity*shift(mat0p0,512,0)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.9)
im = ax.imshow(+color_intensity*shift(mat1p0,512,0)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.6)
ax.text(256, 64, r"$\lambda/4$", va="bottom", ha="center", fontsize=9)
ax.annotate("", (128, 48), (384, 48), arrowprops=scale_bar, color="k")
mpl.label(ax,"l)", color="k")


ax = fig.add_subplot(spec[0,1])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(-color_intensity*shift(geom0p0,0,512)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.9)
im = ax.imshow(+color_intensity*shift(geom1p0,0,512)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.6)
#ax.text(256, 64, r"$\lambda/4$", va="bottom", ha="center", fontsize=9)
#ax.annotate("", (128, 48), (384, 48), arrowprops=scale_bar, color="k")
mpl.label(ax,"k)", color="k")

ax = fig.add_subplot(spec[1,1])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(-color_intensity*shift(geom0p0,512,0)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.9)
im = ax.imshow(+color_intensity*shift(geom1p0,512,0)[256:768,256:768].transpose(), cmap="turbo", clim=(-1,+1), origin="lower", alpha=0.6)
#ax.text(256, 64, r"$\lambda/4$", va="bottom", ha="center", fontsize=9)
#ax.annotate("", (128, 48), (384, 48), arrowprops=scale_bar, color="k")
mpl.label(ax,"m)", color="k")


#cb = fig.colorbar(im, ax=ax, **cb_props)
#cb.ax.minorticks_on()

plt.subplots_adjust(wspace=0.05, hspace=0.05)

mpl.digitalize(fig, "HertzConstriction")
if outpath: plt.savefig(outpath+"Constriction.pdf")