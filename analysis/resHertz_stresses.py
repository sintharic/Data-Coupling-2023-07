#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-none/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-mat/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-geom/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict0p0/press100/stat-geomV5/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/none/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/material/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/geometric/"
#inpath = "/Users/christian/sim/coupled/Hertz-fix/N1024/gamma0/frict1p0/press100/geomV5/"

import eval_stressH

eval_stressH.do_flip = False
eval_stressH.do_rot = True


eval_stressH.process("../friction/hertz/none/frict0p0/")
eval_stressH.process("../friction/hertz/none/frict1p0/")
eval_stressH.process("../friction/hertz/material/frict0p0/")
eval_stressH.process("../friction/hertz/material/frict1p0/")
eval_stressH.process("../friction/hertz/geometric1/frict0p0/")
eval_stressH.process("../friction/hertz/geometric1/frict1p0/")
eval_stressH.process("../friction/hertz/geometric2/frict0p0/")
eval_stressH.process("../friction/hertz/geometric2/frict1p0/")