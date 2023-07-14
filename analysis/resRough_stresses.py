import eval_stress

eval_stress.fake_mu = 0.0
eval_stress.process("../friction/rough/material/press2p0e-02/frict1p0/")

eval_stress.fake_mu = 1.0
eval_stress.process("../friction/rough/material/press2p0e-02/frict0p0/")

eval_stress.fake_mu = 0.0
eval_stress.process("../friction/rough/geometric1/press4p0e-02/frict1p0/")

eval_stress.fake_mu = 1.0
eval_stress.process("../friction/rough/geometric1/press4p0e-02/frict0p0/")