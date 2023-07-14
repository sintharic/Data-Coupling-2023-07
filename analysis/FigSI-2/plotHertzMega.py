path_none_mu0p0 = "../../friction/hertz/none/frict0p0"
path_none_mu1p0 = "../../friction/hertz/none/frict1p0"
path_mat_mu0p0 = "../../friction/hertz/material/frict0p0"
path_mat_mu1p0 = "../../friction/hertz/material/frict1p0"
path_geomV1_mu0p0 = "../../friction/hertz/geometric1/frict0p0"
path_geomV1_mu1p0 = "../../friction/hertz/geometric1/frict1p0"
path_geomV2_mu0p0 = "../../friction/hertz/geometric2/frict0p0"
path_geomV2_mu1p0 = "../../friction/hertz/geometric2/frict1p0"

paths = [path_none_mu0p0, path_none_mu1p0, path_mat_mu1p0, path_geomV2_mu0p0, path_geomV2_mu1p0, path_geomV1_mu1p0] 
do_rot = True

F = 0.01
R = 1
E = 1

ylabel = [r"$u_\mathrm{z}$", r"$p_\mathrm{z}/p_\mathrm{H}$", r"$w_\mathrm{diss}/p_\mathrm{H} v_0$",
          r"$v^\mathrm{rel}_\mathrm{l} / v_0$", r"$\dot{u}_\mathrm{t} / v_0$", 
          r"$\sigma/E^\mathrm{*}$", r"$\mathbf{v}(\sigma_\mathrm{I}^\mathrm{max})$"]
title  = [r"$h \to \infty, \nu \approx 0.5$", r"$v_0>0, \mu=1$", 
          r"$\nu=0.25$", r"$h \approx a_\mathrm{H}, \nu \approx 0.5$", 
          r"$v_0>0, \mu=1$", r"$h \approx a_\mathrm{H}/4$, $\mu=1$"]
ylims = [(-0.005,0.035), (-0.05,1.35), (-0.05, 1.35), (0.45, 1.35), (-0.065,0.065), (-0.135, 0.155), (-0.05, 1.05)]



import os
import numpy as np
import cmutils as cm
import matplotlib.pyplot as plt
import cmparams as cp
import mpl
mpl.RevTex(2)
mpl.set_fontsize(9)
mpl.simple_legend(False)

def diffX(array):
  return ( np.diff(array, append=array[:1,:], axis=0) + np.diff(array, prepend=array[-1:,:], axis=0) ) / 2

def diffY(array):
  return ( np.diff(array, append=array[:,:1], axis=1) + np.diff(array, prepend=array[:,-1:], axis=1) ) / 2

def d(vect):
  return ( np.diff(vect, append=vect[0]) + np.diff(vect, prepend=vect[-1]) ) / 2

aHertz = (0.75 * F * R / E)**(1./3)
pHertz = 3*F/(2*np.pi*aHertz**2)



x = []
y = []
ux_x = []
ux_y = []
uy_x = []
uy_y = []
uz_x = []
uz_y = []
pz_x = []
pz_y = []
indX = []
indY = []
eig2_x = []
vec2_x_x = []
vec2_x_y = []
vec2_x_z = []
VonMises = []
vxLoc_x = []
vxLoc_y = []
vyLoc_x = []
vyLoc_y = []
Lx = []
Ly = []
nx = []
ny = []

for path in paths:
  print(path)
  sim = cp.read(path + os.sep + "params.out")

  data = cm.readMulti(path+os.sep+sim.sheet[0].konfigName)
  if do_rot: 
    data = cm.rot270Multi(data)
    nx.append( sim.nyGlobal )
    ny.append( sim.nxGlobal )
    Lx.append( sim.lengthY )
    Ly.append( sim.lengthX )
    vx = sim.sheet[0].vYOnSite
    vy = -sim.sheet[0].vXOnSite
    rXhertz = sim.sheet[0].rYhertz
    rYhertz = sim.sheet[0].rXhertz
  else:
    nx.append( sim.nxGlobal )
    ny.append( sim.nyGlobal )
    Lx.append( sim.lengthX )
    Ly.append( sim.lengthY )
    vx = sim.sheet[0].vXOnSite
    vy = sim.sheet[0].vYOnSite
    rXhertz = sim.sheet[0].rXhertz
    rYhertz = sim.sheet[0].rYhertz

  dx = Lx[-1]/nx[-1]
  dy = Ly[-1]/ny[-1]

  x.append( np.linspace(-sim.lengthX/2, sim.lengthX/2, sim.nxGlobal+1)[:-1] )
  y.append( sim.lengthY/sim.lengthX*x[-1] )

  uz = data[0]
  pz = data[1]
  uy = data[2]
  py = data[3]
  ux = data[4]
  px = data[5]

  dux_dx = diffX(ux) / dx
  dux_dy = diffX(ux) / dy
  duy_dx = diffX(uy) / dx
  duy_dy = diffX(uy) / dy
  duz_dx = diffX(uz) / dx
  duz_dy = diffX(uz) / dy

  if vx==0: 
    vxLoc_x.append( -np.ones_like(x[-1]) )
    vxLoc_y.append( -np.ones_like(y[-1]) )
    vyLoc_x.append( -np.ones_like(x[-1]) )
    vyLoc_y.append( -np.ones_like(y[-1]) )
  else:
    vxLoc_x.append( -dux_dx[:,ny[-1]//2] )
    vxLoc_y.append( -dux_dx[nx[-1]//2,:] )
    vyLoc_x.append( -duy_dx[:,ny[-1]//2] )
    vyLoc_y.append( -duy_dx[nx[-1]//2,:] )

  ux_x.append( ux[:,ny[-1]//2] )
  ux_y.append( ux[nx[-1]//2,:] )
  uy_x.append( uy[:,ny[-1]//2] )
  uy_y.append( uy[nx[-1]//2,:] )
  uz_x.append( uz[:,ny[-1]//2] )
  uz_y.append( uz[nx[-1]//2,:] )
  pz_x.append( pz[:,ny[-1]//2] )
  pz_y.append( pz[nx[-1]//2,:] )

  indX.append( np.square(x[-1])/sim.sheet[0].rXhertz/2 )
  indY.append( np.square(y[-1])/sim.sheet[0].rYhertz/2 )
  eig0_x = np.loadtxt(sim.path+os.sep+"stressEigenvalues.dat", usecols=1)
  eig1_x = np.loadtxt(sim.path+os.sep+"stressEigenvalues.dat", usecols=2)
  eig2_x.append( np.loadtxt(sim.path+os.sep+"stressEigenvalues.dat", usecols=3) )
  vec2_x_x.append( np.loadtxt(sim.path+os.sep+"stressEigenvects.dat", usecols=7) )
  vec2_x_y.append( np.loadtxt(sim.path+os.sep+"stressEigenvects.dat", usecols=8) )
  vec2_x_z.append( np.loadtxt(sim.path+os.sep+"stressEigenvects.dat", usecols=9) )
  VonMises.append( np.square(eig0_x-eig1_x) + np.square(eig0_x-eig2_x[-1]) + np.square(eig1_x-eig2_x[-1]) )



fig = plt.figure()
mpl.format(fig, boxsize=(5.6, 4.8), left=0.0, right=-0.25, bottom=0.0, top=0.20, xlabel=1, xticl=1, ylabel=1, yticl=1)
mpl.format(fig, boxsize=(5.6, 5.5), left=0.0, right=-0.25, bottom=0.0, top=0.20, xlabel=1, xticl=1, ylabel=1, yticl=1)
spec = plt.GridSpec(7, 6, width_ratios=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
plt.subplots_adjust(wspace=0.05, hspace=0.05)

for iSim in range(len(paths)):
  dx = Lx[iSim]; dy = Ly[iSim]
  for iaxis in range(len(ylabel)):
    ax = fig.add_subplot(spec[iaxis,iSim])
    if iaxis==0: ax.set_title(title[iSim])
    if iaxis==(len(ylabel)-1): ax.set_xlabel(r"$r/a_\mathrm{H}$")
    if iSim==0: ax.set_ylabel(ylabel[iaxis])
    ax.label_outer()
    ax.set_xlim(-1.8, 1.8)
    
    if iaxis==0:
      ax.set_ylim(ylims[0])
      ax.fill_between(x[iSim]/aHertz, indX[iSim], indX[iSim].max(), color=(0.7,0.7,0.7))
      ax.plot(x[iSim]/aHertz, uz_x[iSim], label="long.")
      ax.plot(y[iSim]/aHertz, uz_y[iSim], "--", label="trans.")
      if iSim==0: ax.legend(loc="lower center")
    if iaxis==1: 
      ax.set_ylim(ylims[1])
      if iSim==(len(paths)-1):
        ax.plot(x[iSim]/aHertz, pz_x[iSim]/pHertz/2.5)
        ax.plot(y[iSim]/aHertz, pz_y[iSim]/pHertz/2.5, "--")
        mpl.label(ax, "(*)", xoffset=0.05, yoffset=0.05)
      else:
        ax.plot(x[iSim]/aHertz, pz_x[iSim]/pHertz)
        ax.plot(y[iSim]/aHertz, pz_y[iSim]/pHertz, "--")
    if iaxis==2:
      ax.set_ylim(ylims[2])
      if iSim in (0,3):
        ax.text(0.5, 0.5, "n.a.", transform=ax.transAxes, va="center", ha="center")
      elif iSim==(len(paths)-1):
        ax.plot(x[iSim]/aHertz, 1.*pz_x[iSim]*(1-vxLoc_x[iSim])/pHertz/2.5)
        ax.plot(y[iSim]/aHertz, 1.*pz_y[iSim]*(1-vxLoc_y[iSim])/pHertz/2.5, "--")
        mpl.label(ax, "(*)", xoffset=0.05, yoffset=0.05)
      else: 
        ax.plot(x[iSim]/aHertz, 1.*pz_x[iSim]*(1-vxLoc_x[iSim])/pHertz)
        ax.plot(y[iSim]/aHertz, 1.*pz_y[iSim]*(1-vxLoc_y[iSim])/pHertz, "--")
    if iaxis==3:
      ax.set_ylim(ylims[3])
      if iSim in (0,3):
        ax.text(0.5, 0.5, "n.a.", transform=ax.transAxes, va="center", ha="center")
      else: 
        ax.plot(x[iSim]/aHertz, 1-vxLoc_x[iSim])
        ax.plot(y[iSim]/aHertz, 1-vxLoc_y[iSim], "--")
    if iaxis==4:
      ax.set_ylim(ylims[4])
      if iSim in (0,3):
        ax.text(0.5, 0.5, "n.a.", transform=ax.transAxes, va="center", ha="center")
      else: 
        ax.plot(x[iSim]/aHertz, vyLoc_x[iSim])
        ax.plot(y[iSim]/aHertz, vyLoc_y[iSim], "--")
    if iaxis==5:
      ax.set_ylim(ylims[5])
      if iSim==(len(paths)-1):
        ax.plot(x[iSim]/aHertz, eig2_x[iSim], "k-")
        ax.plot(x[iSim]/aHertz, VonMises[iSim]/5, "k--", label=r"$\sigma_\mathrm{vM}(l)/5$")
        ax.legend(fancybox=False, frameon=False, facecolor="none", edgecolor="none", loc=(0.02,0.02))
      else:
        ax.plot(x[iSim]/aHertz, eig2_x[iSim], "k-", label=r"$\sigma^\mathrm{max}_\mathrm{I}(l)$")
        ax.plot(x[iSim]/aHertz, VonMises[iSim], "k--", label=r"$\sigma_\mathrm{vM}(l)$")
      if iSim==0: ax.legend(fancybox=False, frameon=False, facecolor="none", edgecolor="none")
    if iaxis==6:
      ax.set_ylim(ylims[6])
      ax.plot(x[iSim]/aHertz, np.abs(vec2_x_x[iSim]), color=mpl.color(0), label=r"$v_\mathrm{l}$")
      ax.plot(x[iSim]/aHertz, np.abs(vec2_x_y[iSim]), color=mpl.color(1), label=r"$v_\mathrm{t}$")
      ax.plot(x[iSim]/aHertz, np.abs(vec2_x_z[iSim]), color=mpl.color(2), label=r"$v_\mathrm{z}$")
      if iSim==0: ax.legend(fancybox=False, frameon=False, facecolor="none", edgecolor="none")

mpl.digitalize(fig, "HertzMega")
plt.savefig("HertzMega.pdf")