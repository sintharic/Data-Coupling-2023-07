import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import ticker
import mpl 
mpl.RevTex(2)
#mpl.set_fontsize(9)

#figsize = (1.2, 1.2)
figsize = (1.8,1.2)


# GEOMETRIC: compare different friction coeffs and thicknesses
inpaths = ["../results/rough/res-geom1-frict0p3.dat",
  "../results/rough/res-geom1-frict1p0.dat",
  "../results/rough/res-geom2-frict0p3.dat",
  "../results/rough/res-geom2-frict1p0.dat"
]
refpaths = ["../results/rough/res-geom1-frict0p3.dat",
  "../results/rough/res-geom2-frict0p3.dat",
]
#inpaths = ["/Users/christian/sim/coupled/H0p8-n0p49-fix/resultsV2/res-frict0p3.dat",
##"/Users/christian/sim/coupled/H0p8-n0p49-fix/resultsV2/res-frict0p6.dat",
#"/Users/christian/sim/coupled/H0p8-n0p49-fix/resultsV2/res-frict1p0.dat",
##"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h065/resultsV2/res-frict0p3.dat",
##"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h065/resultsV2/res-frict0p6.dat",
##"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h065/resultsV2/res-frict1p0.dat",
#"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h080/resultsV2/res-frict0p3.dat",
##"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h080/resultsV2/res-frict0p6.dat",
#"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h080/resultsV2/res-frict1p0.dat"]
#refpaths = ["/Users/christian/sim/coupled/H0p8-n0p49-fix/resultsV2/res-frict0p0.dat",
##"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h065/resultsV2/res-frict0p0.dat",
#"/Volumes/Elements4TB/sim/coupled/H0p8-n0p49-fix-h080/resultsV2/res-frict0p0.dat",]
refnums = [2,2]#[3,3,3,0]
reflabel = [r"$h/\lambda_\mathrm{l} = 0.05$", r"$h/\lambda_\mathrm{l} = 0.1$"]#["$h=0.050$", "$h=0.065$", "$h=0.080$"]
mu = [0.3, 1.0]*2#[0.3, 0.6, 1.0]*3
RMSgrad_global = [0.276452]*4#9
labels = ["$\\mu=%g$" % val for val in mu[:2]]#3]]
labels += [""]*2#6
color = [mpl.color(1), mpl.color(0)]*2#range(3)]*3
line = ["-"]*2 + [":"]*2#["-"]*3 + ["--"]*3 + ["dashdot"]*3
lambdaR = [1]*4#9

outpath = "Rough"

add_inset = False 
filter_valid = False





dataSet = []
dataRef = []
for path in inpaths:
  dataSet.append(np.loadtxt(path))
for path,num in zip(refpaths,refnums):
  dataRef = dataRef + [np.loadtxt(path)]*num

# plot geom stress
fig1,axGeomFrict = plt.subplots()
axGeomFrict.set_xlabel(r"$p h / \lambda_\mathrm{l} E^\mathrm{*}$")
axGeomFrict.set_ylabel(r"$\Delta \mu / \mu_\mathrm{c} \bar{g}$")

if add_inset:
  inset = axGeomFrict.inset_axes([0.2,0.12,0.43,0.43])
  inset.semilogx()
  inset.set_xlim(1e-4,2e-2)
  #inset.set_ylim(1e-8,1e-1)


# format figures
mpl.format(fig1, figsize, right=0.05, yticl=1.0)


for i in range(len(refpaths)):
  if i==0: label = r"$\mu=0$"
  else: label = ""
  j = 0
  for k in range(0,i): j+= refnums[k]
  data = dataRef[j]

  # data validity criterion
  if filter_valid:
    meanGap = data[:,7]
    meanPosGap = data[:,8]
    valid = meanGap > meanPosGap/2
  else:
    valid = range(data.shape[0])

  poisson = data[valid,0]
  thickness = data[valid,1]
  frictCoeff = data[valid,2]
  pressure = data[valid,3]
  #geomStress = data[valid,4] # Y direction
  geomStress = data[valid,5] # X direction
  relCont = data[valid,6]
  meanGap = data[valid,7]
  meanPosGap = data[valid,8]
  rmsgrad = data[valid,9]
  #rmsgrad = RMSgrad_global[i]
  
  pressN = pressure*thickness*2*np.pi/lambdaR[i]/RMSgrad_global[i]




for i,data in enumerate(dataSet):
  # data file header:
  # poisson thickness pressure  frictCoeff  error relContArea geomStressY meanNormGap meanNormPosGap  rmsGrad meanGap fullContElaEnergy
  
  # data validity criterion
  if filter_valid:
    meanGap = data[:,7]
    meanPosGap = data[:,8]
    valid = meanGap > meanPosGap/2
  else:
    valid = range(data.shape[0])

  poisson = data[valid,0]
  thickness = data[valid,1]
  frictCoeff = data[valid,2]
  pressure = data[valid,3]
  #geomStress = data[valid,4] # Y direction
  geomStress = data[valid,5] # X direction
  relCont = data[valid,6]
  meanGap = data[valid,7]
  meanPosGap = data[valid,8]
  rmsgrad = data[valid,9]
  #rmsgrad = RMSgrad_global[i]

  if mu[i] != frictCoeff[0]: 
    print("WARNING: frictCoeff in "+inpaths[i]+" does not match "+__name__)
  
  #pressN = pressure*thickness*2*np.pi/lambdaR[i]/RMSgrad_global[i]
  pressN = pressure*thickness/lambdaR[i]/RMSgrad_global[i]
  geomFrict = geomStress/pressure/frictCoeff/RMSgrad_global[i]

  if mu[i] > 0: 
    axGeomFrict.plot(pressN, geomFrict, label=labels[i],marker="",color=color[i],ls=line[i])
    if add_inset: inset.plot(pressN, geomFrict, marker="", color=color[i], ls=line[i])

#LinFormat = ticker.ScalarFormatter()
#LinFormat = ticker.FormatStrFormatter("%g")
#axGeomFrict.legend(loc='lower right')
axGeomFrict.semilogx()
axGeomFrict.set_xlim(5e-4, 0.13)
axGeomFrict.set_ylim(0.022, 0.22)
#mpl.setsci(axGeomFrict,"y")
axGeomFrict.text(7e-4, 0.12, reflabel[0], va="bottom", ha="left", rotation=30)
axGeomFrict.text(7e-4, 0.035, reflabel[1], va="bottom", ha="left", rotation=15)
axGeomFrict.text(2.0e-2, 0.172, r"$\mu_\mathrm{c}=0.3$", va="bottom", ha="left", color=axGeomFrict.lines[0].get_color())
axGeomFrict.text(2.0e-2, 0.147, r"$\mu_\mathrm{c}=1$", va="bottom", ha="left", color=axGeomFrict.lines[1].get_color())



mpl.label(fig1,"e)")
mpl.digitalize(fig1, "RoughGeom")
print("figsize (in): ", fig1.get_size_inches())
plt.figure(fig1.number); 
plt.savefig(outpath+"GeomStress.pdf")