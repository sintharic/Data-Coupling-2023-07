import eval_dissip

eval_dissip.do_check = True
eval_dissip.do_export = True

path = "../friction/hertz/none/frict1p0"
eval_dissip.zoom_factor = 2.5
eval_dissip.relative_pressure_tolerance = 1e-06
eval_dissip.process(path)

path = "../friction/hertz/material/frict1p0"
eval_dissip.zoom_factor = 2.5
eval_dissip.relative_pressure_tolerance = 1e-06
eval_dissip.process(path)

path = "../friction/hertz/geometric2/frict1p0"
eval_dissip.zoom_factor = 2.5
eval_dissip.relative_pressure_tolerance = 1e-04
eval_dissip.process(path)