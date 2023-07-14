seeds = ["seed%i"%i for i in range(4)]
velocities = ["_v0_", "_vX_", "_vY_"]

#filter
#velocities = ["_v0_"]


import os
import cmutils as cm
from prepare import *



def main():
    seeds = [dst+"/"+s for s in seeds]
    
    for seed in seeds:
        for v in velocities:
            src = seed + v + "flowX"
            dst = seed + v + "flowY"
            print(dst)

            if not os.path.isfile(src+"/"+gapfile):
                print(f"ERROR: {src}/{gapfile} does not exist\n")
                continue

            os.makedirs(dst, exist_ok=True)

            print(f"reading gap file from {src}")
            gap = cm.readConfig(src+"/"+gapfile)
            gap = cm.rot270(gap)
            
            print("writing rotated gap file")
            cm.dumpConfig(gap, dst+"/"+gapfile)
        
            # convert gap file for flow calculation
            wd = os.getcwd()
            os.chdir(dst)
            os.system(wd+"/"+convert_gap)
            os.chdir(wd)
            print("")

if __name__ == '__main__': main()