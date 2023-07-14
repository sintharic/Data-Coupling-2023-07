file = "currentCalc.x_01024.current.000000.dat"
path_flowX = f"../../flow/rough/geometric_old/roughV0Y/mu0/{file}"
path_flowY = f"../../flow/rough/geometric_old/roughV0Y/mu0_rot90/{file}"
path_flowX_muY = f"../../flow/rough/geometric_old/roughV0Y/mu1/{file}"
path_flowY_muY = f"../../flow/rough/geometric_old/roughV0Y/mu1_rot90/{file}"
path_flowX_muX = f"../../flow/rough/geometric_old/roughV0X/mu1/{file}"
path_flowY_muX = f"../../flow/rough/geometric_old/roughV0X/mu1_rot90/{file}"
shape = (1024,1024)
lab_flow = (500,250)
lab_frict = (500,450)
lab_dist = 1.

# file = "currentCalc.x_02048.current.000000.dat"
# path_flowX = f"../../flow/rough/geometric/seed0_v0_flowX/{file}"
# path_flowY = f"../../flow/rough/geometric/seed0_v0_flowY/{file}"
# path_flowX_muX = f"../../flow/rough/geometric/seed0_vX_flowX/{file}"
# path_flowY_muX = f"../../flow/rough/geometric/seed0_vX_flowY/{file}"
# path_flowX_muY = f"../../flow/rough/geometric/seed0_vY_flowX/{file}"
# path_flowY_muY = f"../../flow/rough/geometric/seed0_vY_flowY/{file}"
# shape = (2048,2048)
# lab_flow = (1000,500)
# lab_frict = (1000,900)
# lab_dist = 2.

CLIM = (0, 0.9)#(4.2,1)


boxsize = (5.9, 1.3)

cb_props = dict(aspect=10, shrink=0.85, pad=0.02, fraction=0.2)

import numpy as np
import matplotlib.pyplot as plt
import mpl
mpl.RevTex(2)
mpl.grid(False)
mpl.set_fontsize(9)
import cmutils as cm
import cmparams as cp
import os

def readFlow(path, shape):
  with open(path) as file:
    lines = [line for line in file.readlines() if line[0].isalnum()]
  array = np.array([float(line.split()[2]) for line in lines])
  return array.reshape(shape)

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

# import and rotate data files properly
img_flowX = readFlow(path_flowX, shape)
img_flowX = shift(img_flowX,0,0) #####
#img_flowX = cm.rot270(img_flowX)

img_flowY = readFlow(path_flowY, shape)
img_flowY = cm.rot270(img_flowY) #####
img_flowY = cm.rot90(img_flowY)
#img_flowY = cm.rot180(img_flowY)
max_flow_mu0 = max(img_flowX.max(), img_flowY.max())

img_flowX_muY = readFlow(path_flowX_muY, shape)
img_flowX_muY = shift(img_flowX_muY,0,415) #####
#img_flowX_muY = cm.rot270(img_flowX_muY)

img_flowY_muY = readFlow(path_flowY_muY, shape)
img_flowY_muY = cm.smoothPBC(img_flowY_muY, 42)
img_flowY_muY = cm.rot270(shift(img_flowY_muY,610,0)) #####
#img_flowY_muY = cm.rot180(shift(img_flowY_muY,610,0)) #####
### img_flowY_muY = cm.rot90(img_flowY_muY)

img_flowX_muX = readFlow(path_flowX_muX, shape)
img_flowX_muX = cm.smoothPBC(img_flowX_muX, 42)
img_flowX_muX = shift(img_flowX_muX,415,0) #####
#img_flowX_muX = cm.rot270(img_flowX_muX)

img_flowY_muX = readFlow(path_flowY_muX, shape)
img_flowY_muX = cm.smoothPBC(img_flowY_muX, 42)
img_flowY_muX = shift(img_flowY_muX,0,415) #####
img_flowY_muX = cm.rot270(img_flowY_muX) #####
### img_flowY_muX = cm.rot90(img_flowY_muX)
#img_flowY_muX = cm.rot180(img_flowY_muX)

# plot frictionless cases
#cm.plotImg(img_flowX, rot=True, clim=CLIM, title="flow in x")
#cm.plotImg(img_flowY, rot=True, clim=CLIM, title="flow in y")

# plot with friction
fig = plt.figure()
mpl.format(fig, boxsize=(6.3,1.5), left=-0.08, right=-0.02, bottom=0.0, top=0.0, xlabel=0, yticl=0.0, xticl=0)
spec = plt.GridSpec(1,4, width_ratios=[1.0,1.0,1.0,1.28])

mean_flowX = img_flowX.mean()
mean_flowX_muX = img_flowX_muX.mean()
title = "flow in x, frict. in x (%.2f)" % (mean_flowX_muX/mean_flowX)
ax = fig.add_subplot(spec[0])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(img_flowX_muX.transpose()/max_flow_mu0, cmap="turbo", clim=CLIM, origin="lower")
#ax = cm.plotImg(img_flowX_muX, rot=True, clim=CLIM)#, title=title)
ax.annotate("flow", (lab_flow[0]+(200*lab_dist),lab_flow[1]), (lab_flow[0]-(30*lab_dist),lab_flow[1]), color="w", ha="right", va="center", arrowprops=dict(arrowstyle='-|>',color="w"))
ax.annotate("sliding", (lab_frict[0]+(200*lab_dist),lab_frict[1]), (lab_frict[0]-(30*lab_dist),lab_frict[1]), color="orange", ha="right", va="center", arrowprops=dict(arrowstyle='-|>',color="orange"))
mpl.label(ax,"\\Large{e)}", color="w")
print(title)

mean_flowX = img_flowX.mean()
mean_flowX_muY = img_flowX_muY.mean()
title = "flow in x, frict. in y (%.2f)" % (mean_flowX_muY/mean_flowX)
ax = fig.add_subplot(spec[1])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(img_flowX_muY.transpose()/max_flow_mu0, cmap="turbo", clim=CLIM, origin="lower")
#ax = cm.plotImg(img_flowX_muY, rot=True, clim=CLIM)#, title=title)
ax.annotate("flow", (lab_flow[0]+(200*lab_dist),lab_flow[1]), (lab_flow[0]-(30*lab_dist),lab_flow[1]), color="w", ha="right", va="center", arrowprops=dict(arrowstyle='-|>',color="w"))
ax.annotate("", (lab_frict[0]+(50*lab_dist),lab_frict[1]+(100*lab_dist)), (lab_frict[0]+(50*lab_dist),lab_frict[1]-(100*lab_dist)), color="orange", arrowprops=dict(arrowstyle='-|>',color="orange"))
ax.text(lab_frict[0], lab_frict[1], "sliding", va="center", ha="right", color="orange")
mpl.label(ax,"\\Large{f)}", color="w")
#cb = fig.colorbar(im, ax=ax, **cb_props)
#cb.ax.minorticks_on()
print(title)

mean_flowY = img_flowY.mean()
mean_flowY_muX = img_flowY_muX.mean()
title = "flow in y, frict. in x (%.2f)" % (mean_flowY_muX/mean_flowY)
ax = fig.add_subplot(spec[2])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(img_flowY_muX.transpose()/max_flow_mu0, cmap="turbo", clim=CLIM, origin="lower")
#ax = cm.plotImg(img_flowY_muX, rot=True, clim=CLIM)#, title=title)
ax.annotate("", (lab_flow[0]+(50*lab_dist),lab_flow[1]+(100*lab_dist)), (lab_flow[0]+(50*lab_dist),lab_flow[1]-(100*lab_dist)), color="w", arrowprops=dict(arrowstyle='-|>',color="w"))
ax.text(lab_flow[0], lab_flow[1], "flow", va="center", ha="right", color="w")
ax.annotate("sliding", (lab_frict[0]+(200*lab_dist),lab_frict[1]), (lab_frict[0]-(30*lab_dist),lab_frict[1]), color="orange", ha="right", va="center", arrowprops=dict(arrowstyle='-|>',color="orange"))
mpl.label(ax,"\\Large{g)}", color="w")
print(title)

mean_flowY = img_flowY.mean()
mean_flowY_muY = img_flowY_muY.mean()
title = "flow in y, frict. in y (%.2f)" % (mean_flowY_muY/mean_flowY)
ax = fig.add_subplot(spec[3])
ax.set_xticks(()); ax.set_yticks(())
im = ax.imshow(img_flowY_muY.transpose()/max_flow_mu0, cmap="turbo", clim=CLIM, origin="lower")
#ax = cm.plotImg(img_flowY_muY, rot=True, clim=CLIM)#, title=title)
#ax.annotate("", (lab_flow[0]-70,lab_flow[1]+(100*lab_dist)), (lab_flow[0]-70,lab_flow[1]-(100*lab_dist)), color="w", arrowprops=dict(arrowstyle='-|>',color="w"))
#ax.text(lab_flow[0]-(100*lab_dist), lab_flow[1], "flow", va="center", ha="right", color="w")
#ax.annotate("", (lab_flow[0]+(50*lab_dist),lab_flow[1]+(100*lab_dist)), (lab_flow[0]+(50*lab_dist),lab_flow[1]-(100*lab_dist)), color="orange", arrowprops=dict(arrowstyle='-|>',color="orange"))
#ax.text(lab_flow[0]+60, lab_flow[1], "sliding", va="center", ha="left", color="orange")
ax.annotate("", (lab_flow[0]+(50*lab_dist),lab_flow[1]+(100*lab_dist)), (lab_flow[0]+(50*lab_dist),lab_flow[1]-(100*lab_dist)), color="w", arrowprops=dict(arrowstyle='-|>',color="w"))
ax.text(lab_flow[0], lab_flow[1], "flow", va="center", ha="right", color="w")
ax.annotate("", (lab_frict[0]+(50*lab_dist),lab_frict[1]+(100*lab_dist)), (lab_frict[0]+(50*lab_dist),lab_frict[1]-(100*lab_dist)), color="orange", arrowprops=dict(arrowstyle='-|>',color="orange"))
ax.text(lab_frict[0], lab_frict[1], "sliding", va="center", ha="right", color="orange")
mpl.label(ax,"\\Large{h)}", color="w")
print(title)

cb = fig.colorbar(im, ax=ax, **cb_props)
cb.ax.minorticks_on()
plt.subplots_adjust(wspace=0.05, hspace=0.05)
#mpl.label(fig, "f)")

mpl.digitalize(fig, "RoughFlow")
plt.savefig("RoughFlow.pdf")