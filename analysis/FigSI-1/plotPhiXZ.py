import numpy as np
import matplotlib.pyplot as plt
import mpl
mpl.RevTex(2)
mpl.set_fontsize(9)
#mpl.cycle()

def ImPhiXZ_qE(x, nu=0.5):
  return (1-nu) * ((3-4*nu)*(1-2*nu)*np.sinh(x)**2 - (x)**2) / ((3-4*nu)**2*np.sinh(x)**2 - (x)**2)

x = np.linspace(1e-2, 5)
fig = plt.figure()
mpl.format(fig, (2.6,1.8), right=0.05, yticl=0.8)
plt.xlabel("$q h$")
plt.ylabel(r"$\Phi_\mathrm{13}(qh)/\mathrm{i}qE^\mathrm{*}$")
plt.ylim(-0.5,0.35)
plt.xlim(0,4.5)
plt.plot(x, np.zeros_like(x), "-", color="gray")
plt.plot(x, ImPhiXZ_qE(x,0.0), label=r"$\nu = 0$", ls=mpl.line(0))
plt.plot(x, ImPhiXZ_qE(x,0.25), label=r"$\nu = 0.25$", ls=mpl.line(1))
#plt.plot(x, ImPhiXZ_qE(x,0.33), label=r"$\nu = 0.33$", ls=mpl.line(2))
plt.plot(x, ImPhiXZ_qE(x,0.4), label=r"$\nu = 0.4$", ls=mpl.line(2))
plt.plot(x, ImPhiXZ_qE(x,0.5), label=r"$\nu = 0.5$", ls=mpl.line(3))

plt.legend()
mpl.digitalize(fig, "ImPhiXZ")
plt.savefig("ImPhiXZ.pdf", dpi=300)