import numpy as np
import eval_frict

eval_frict.rigidpath = ""
eval_frict.OnSiteHertz = True



# without friction
inpaths0p0 = ["../friction/hertz/material2/frict0p0/force1p0e-02",
"../friction/hertz/material2/frict0p0/force1p5e-02",
"../friction/hertz/material2/frict0p0/force2p3e-02",
"../friction/hertz/material2/frict0p0/force3p5e-02",
"../friction/hertz/material2/frict0p0/force5p3e-02",
"../friction/hertz/material2/frict0p0/force8p0e-02"]
output0p0 = "results/hertz/res-mat-frict0p0_red.dat"
frict0p0 = eval_frict.process(inpaths0p0)
# manipulate first data point
frict0p0 = np.vstack((frict0p0[:1,:], frict0p0))
frict0p0[0,3] = 1e-4
# save
eval_frict.save(output0p0, frict0p0)



# with friction
inpaths1p0 = ["../friction/hertz/material2/frict1p0/force1p0e-02",
"../friction/hertz/material2/frict1p0/force1p5e-02",
"../friction/hertz/material2/frict1p0/force2p3e-02",
"../friction/hertz/material2/frict1p0/force3p5e-02",
"../friction/hertz/material2/frict1p0/force5p3e-02",
"../friction/hertz/material2/frict1p0/force8p0e-02"]
output1p0 = "results/hertz/res-mat-frict1p0_red.dat"
frict1p0 = eval_frict.process(inpaths1p0)
# manipulate first data point
frict1p0 = np.vstack((frict1p0[:1,:], frict1p0))
frict1p0[0,3] = 1e-4
# save
eval_frict.save(output1p0, frict1p0)