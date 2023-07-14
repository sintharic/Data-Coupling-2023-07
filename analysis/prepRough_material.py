import prepare
prepare.src = "../friction/rough/material/press2p0e-02"
prepare.dst = "../flow/rough/material"

import prepare_flowX
prepare_flowX.src = "../friction/rough/material/press2p0e-02"
prepare_flowX.dst = "../flow/rough/material"

prepare_flowX.main()



import prepare_flowY
prepare_flowY.src = "../friction/rough/material/press2p0e-02"
prepare_flowY.dst = "../flow/rough/material"

prepare_flowY.main()