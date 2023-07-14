import numpy as np
import eval_frict

eval_frict.rigidpath = "../friction/rough/konfig0E.dat"
eval_frict.OnSiteHertz = False


# without friction
inpaths0p3 = ["../friction/rough/geometric1/press1p0e-03/frict0p3",
"../friction/rough/geometric1/press1p6e-03/frict0p3",
"../friction/rough/geometric1/press2p5e-03/frict0p3",
"../friction/rough/geometric1/press4p0e-03/frict0p3",
"../friction/rough/geometric1/press6p3e-03/frict0p3",
"../friction/rough/geometric1/press1p0e-02/frict0p3",
"../friction/rough/geometric1/press1p6e-02/frict0p3",
"../friction/rough/geometric1/press2p5e-02/frict0p3",
"../friction/rough/geometric1/press4p0e-02/frict0p3",
"../friction/rough/geometric1/press6p3e-02/frict0p3",
"../friction/rough/geometric1/press1p0e-01/frict0p3",
"../friction/rough/geometric1/press2p0e-01/frict0p3",
"../friction/rough/geometric1/press3p0e-01/frict0p3",
"../friction/rough/geometric1/press4p0e-01/frict0p3",
"../friction/rough/geometric1/press5p0e-01/frict0p3",
"../friction/rough/geometric1/press6p0e-01/frict0p3",
"../friction/rough/geometric1/press7p0e-01/frict0p3",
"../friction/rough/geometric1/press8p0e-01/frict0p3",
"../friction/rough/geometric1/press9p0e-01/frict0p3",
"../friction/rough/geometric1/press1p0e+00/frict0p3"]
output0p3 = "results/rough/res-geom1-frict0p3.dat"
frict0p3 = eval_frict.process(inpaths0p3)
# manipulate first data point
# frict0p3 = frict0p3[1:,:]
# frict0p3 = np.vstack((frict0p3[:1,:], frict0p3))
# frict0p3[0,3] = 1e-4
# save
eval_frict.save(output0p3, frict0p3)


# with friction
inpaths1p0 = ["../friction/rough/geometric1/press1p0e-03/frict1p0",
"../friction/rough/geometric1/press1p6e-03/frict1p0",
"../friction/rough/geometric1/press2p5e-03/frict1p0",
"../friction/rough/geometric1/press4p0e-03/frict1p0",
"../friction/rough/geometric1/press6p3e-03/frict1p0",
"../friction/rough/geometric1/press1p0e-02/frict1p0",
"../friction/rough/geometric1/press1p6e-02/frict1p0",
"../friction/rough/geometric1/press2p5e-02/frict1p0",
"../friction/rough/geometric1/press4p0e-02/frict1p0",
"../friction/rough/geometric1/press6p3e-02/frict1p0",
"../friction/rough/geometric1/press1p0e-01/frict1p0",
"../friction/rough/geometric1/press2p0e-01/frict1p0",
"../friction/rough/geometric1/press3p0e-01/frict1p0",
"../friction/rough/geometric1/press4p0e-01/frict1p0",
"../friction/rough/geometric1/press5p0e-01/frict1p0",
"../friction/rough/geometric1/press6p0e-01/frict1p0",
"../friction/rough/geometric1/press7p0e-01/frict1p0",
"../friction/rough/geometric1/press8p0e-01/frict1p0",
"../friction/rough/geometric1/press9p0e-01/frict1p0",
"../friction/rough/geometric1/press1p0e+00/frict1p0"]
output1p0 = "results/rough/res-geom1-frict1p0.dat"
frict1p0 = eval_frict.process(inpaths1p0)
# manipulate first data point
# frict1p0 = frict1p0[1:,:]
# frict1p0 = np.vstack((frict0p3[:1,:], frict1p0))
# save
eval_frict.save(output1p0, frict1p0)
