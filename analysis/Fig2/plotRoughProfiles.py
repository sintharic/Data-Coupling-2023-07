inpaths = ["../../friction/rough/geometric1/press1p0e-01/frict0p0/",
           "../../friction/rough/geometric1/press1p0e-01/frict0p3/",
           "../../friction/rough/geometric1/press1p0e-01/frict1p0/"]

outpath = "RoughConfigsX.pdf"
mu = [0.0, 0.3, 1.0]
refPath = 0
length = 2

grayLevel = 0.85
top = 1.20
shift = 0

#figsize = (1.2, 1.2)
figsize = (1.8,1.2)



import numpy as np
import mpl 
mpl.RevTex(2, simple=False)
#mpl.set_fontsize(9)
mpl.simple_legend(False)
mpl.grid(False)
import matplotlib.pyplot as plt 
plt.close("all")

color = [mpl.color(2), mpl.color(1), mpl.color(0)]
line  = [mpl.line(0), mpl.line(1), mpl.line(2)]



dat0 = []
dat1 = []
for path in inpaths:
  dat0.append(np.loadtxt(path+"konfig0E.datH",usecols=1))
  dat1.append(np.loadtxt(path+"konfig1Dzyx.datH",usecols=1))

if shift:
  for iPath in range(len(inpaths)):
    data = np.zeros(dat0[iPath].shape)
    data[:-shift] = dat0[iPath][shift:]
    data[-shift:] = dat0[iPath][:shift]
    dat0[iPath] = data

    data = np.zeros_like(dat1[iPath])
    data[:-shift] = dat1[iPath][shift:]
    data[-shift:] = dat1[iPath][:shift]
    dat1[iPath] = data
    mask = dat1[iPath]>dat0[iPath]
    dat1[iPath][mask] = dat0[iPath][mask]

coord = np.linspace(0,length,dat0[0].shape[0])

fig = plt.figure()
plt.xlim(-1.2,0.6)
plt.ylim(0.022, 0.048)
plt.plot(coord, dat0[refPath], "k-")
plt.fill_between(coord,top*dat0[refPath].max(),dat0[refPath],facecolor=(grayLevel,grayLevel,grayLevel))
plt.plot(coord-length, dat0[refPath], "k-")
plt.fill_between(coord-length,top*dat0[refPath].max(),dat0[refPath],facecolor=(grayLevel,grayLevel,grayLevel))

#plt.annotate(text="sliding", xy=(0.6, 0.045), xytext=(0.2, 0.045), arrowprops=dict(arrowstyle='-|>',color="k"), ha="right", va="center")
plt.annotate(text="", xy=(0.4, 0.045), xytext=(-0.2, 0.045), arrowprops=dict(arrowstyle='-|>',color="k"))
plt.text(0.1, 0.044, "sliding", va="top", ha="center")

for i in (2,1,0):
  plt.plot(coord, dat1[i], "-", label="$\\mu_\\mathrm{c}="+str(mu[i])+"$", color=color[i])#, ls=line[i])
  plt.plot(coord-length, dat1[i], "-", color=color[i])#, ls=line[i])

plt.xlabel(r"$x/\lambda_\mathrm{l}$")
plt.ylabel(r"$z/\lambda_\mathrm{l}$")
plt.legend(loc="upper left")
mpl.format(fig, figsize, right=0.05, yticl=1.0)
mpl.label(fig, "c)")
mpl.digitalize(fig, "RoughConfigsX")
print("figsize (in): ", fig.get_size_inches())
plt.savefig(outpath)