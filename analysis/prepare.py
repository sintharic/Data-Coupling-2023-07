convert_gap = "readGap.exe -i gap0.dat -nx 2048 -ny 2048"
reynolds = "rey -n 2048" 

dst = "../flow/rough/geometric"
src = "../friction/rough/geometric1/press4p0e-02"
konfig1 = "konfig1Dzyx.dat"
konfig0 = "konfig0E.dat"
gapfile = "gap0.dat"

folders_v0 = ["frict0p0", "seed1_v0", "seed2_v0", "seed3_v0"]
folders_vX = ["frict1p0", "seed1_vX", "seed2_vX", "seed3_vX"]
folders_vY = ["frict1p0_vY", "seed1_vY", "seed2_vY", "seed3_vY"]

newname = {
    "frict0p0" : "seed0_v0",
    "frict1p0" : "seed0_vX",
    "frict1p0_vY" : "seed0_vY"
}

import os 
import shutil
import cmutils as cm


def prepare_flowX(folder):
    if folder in newname.keys(): 
        newfolder = newname[folder]
    else: 
        newfolder = folder
    
    # create folder
    newfolder += "_flowX"
    print(newfolder)
    os.makedirs(dst+"/"+newfolder, exist_ok=True)
    
    # import gap data
    print(f"reading data from {src}/{folder}")
    shutil.copy2(f"{src}/{folder}/params.out", dst+"/"+newfolder)
    #return#TEMP
    if os.path.isfile(src+"/"+folder+"/"+"/"+gapfile):
        gap = cm.readCOnfig(src+"/"+folder+"/"+"/"+gapfile)
    else:
        pos0 = cm.readConfig(src+"/"+folder+"/"+konfig0)
        pos1 = cm.readConfig(src+"/"+folder+"/"+konfig1)
        gap = pos0-pos1
    
    # make gap non-zero
    gap[gap<0] = 0
    print("writing (strictly non-zero) gap file")
    cm.dumpConfig(gap, dst+"/"+newfolder+"/gap0.dat")

    # convert gap file for flow calculation
    wd = os.getcwd()
    os.chdir(dst+"/"+newfolder)
    os.system(wd+"/"+convert_gap)
    os.chdir(wd)
    print("")


