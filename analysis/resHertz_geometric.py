import numpy as np
import eval_frict

eval_frict.rigidpath = ""
eval_frict.OnSiteHertz = True



# without friction
inpaths0p0 = ["../friction/hertz/geometric2/frict0p0/force1p0e-04",
"../friction/hertz/geometric2/frict0p0/force1p6e-04",
"../friction/hertz/geometric2/frict0p0/force2p5e-04",
"../friction/hertz/geometric2/frict0p0/force4p0e-04",
"../friction/hertz/geometric2/frict0p0/force6p3e-04",
"../friction/hertz/geometric2/frict0p0/force1p0e-03",
"../friction/hertz/geometric2/frict0p0/force1p6e-03",
"../friction/hertz/geometric2/frict0p0/force2p5e-03",
"../friction/hertz/geometric2/frict0p0/force4p0e-03",
"../friction/hertz/geometric2/frict0p0/force6p3e-03",
"../friction/hertz/geometric2/frict0p0/force1p0e-02",
"../friction/hertz/geometric2/frict0p0/force2p0e-02",
"../friction/hertz/geometric2/frict0p0/force3p0e-02",
"../friction/hertz/geometric2/frict0p0/force4p0e-02",
"../friction/hertz/geometric2/frict0p0/force5p0e-02",
"../friction/hertz/geometric2/frict0p0/force6p0e-02",
"../friction/hertz/geometric2/frict0p0/force7p0e-02",
"../friction/hertz/geometric2/frict0p0/force8p0e-02",
"../friction/hertz/geometric2/frict0p0/force9p0e-02",
"../friction/hertz/geometric2/frict0p0/force1p0e-01"]
output0p0 = "results/hertz/res-geometric2-frict0p0.dat"
frict0p0 = eval_frict.process(inpaths0p0)
# manipulate first data point
# frict0p0 = frict0p0[1:,:]
# frict0p0 = np.vstack((frict0p0[:1,:], frict0p0))
# frict0p0[0,3] = 1e-4
# save
eval_frict.save(output0p0, frict0p0)


# with friction
inpaths1p0 = ["../friction/hertz/geometric2/frict1p0/force1p0e-04",
"../friction/hertz/geometric2/frict1p0/force1p6e-04",
"../friction/hertz/geometric2/frict1p0/force2p5e-04",
"../friction/hertz/geometric2/frict1p0/force4p0e-04",
"../friction/hertz/geometric2/frict1p0/force6p3e-04",
"../friction/hertz/geometric2/frict1p0/force1p0e-03",
"../friction/hertz/geometric2/frict1p0/force1p6e-03",
"../friction/hertz/geometric2/frict1p0/force2p5e-03",
"../friction/hertz/geometric2/frict1p0/force4p0e-03",
"../friction/hertz/geometric2/frict1p0/force6p3e-03",
"../friction/hertz/geometric2/frict1p0/force1p0e-02",
"../friction/hertz/geometric2/frict1p0/force2p0e-02",
"../friction/hertz/geometric2/frict1p0/force3p0e-02",
"../friction/hertz/geometric2/frict1p0/force4p0e-02",
"../friction/hertz/geometric2/frict1p0/force5p0e-02",
"../friction/hertz/geometric2/frict1p0/force6p0e-02",
"../friction/hertz/geometric2/frict1p0/force7p0e-02",
"../friction/hertz/geometric2/frict1p0/force8p0e-02",
"../friction/hertz/geometric2/frict1p0/force9p0e-02",
"../friction/hertz/geometric2/frict1p0/force1p0e-01"]
output1p0 = "results/hertz/res-geometric2-frict1p0.dat"
frict1p0 = eval_frict.process(inpaths1p0)
# manipulate first data point
# frict1p0 = frict1p0[1:,:]
# frict1p0 = np.vstack((frict0p0[:1,:], frict1p0))
# save
eval_frict.save(output1p0, frict1p0)
