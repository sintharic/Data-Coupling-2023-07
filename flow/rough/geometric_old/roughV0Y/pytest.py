path = "/Users/christian/sim/flow/mu0/currentCalc.x_01024.current.000000.dat" # 3 cols
path = "/Users/christian/sim/coupled/H0p8-n0p49-fix/press0p10/frict1p0/konfig1Dzyx.dat" # 8 cols


import numpy as np


def version0(path=path, usecols=0, delimiter=None):
  if isinstance(usecols,int): usecols = (usecols,)
  output = ([],)
  for _ in range(1,len(usecols)): output = output + ([],) 
  with open(path) as file:
    for line in file.readlines():
      if not line[0].isalnum(): continue
      split = line.split(delimiter)
      for i in range(len(usecols)): 
        output[i].append(float(split[usecols[i]]))
    
  return output

def version1(path=path):
  
  with open(path) as file:
    lines = [line for line in file.readlines() if line[0].isalnum()]
  x = [float(line.split()[0]) for line in lines]
  y = [float(line.split()[1]) for line in lines]
  
  return x,y

def version1np(path=path):
  
  with open(path) as file:
    lines = [line for line in file.readlines() if line[0].isalnum()]
  x = np.array([float(line.split()[0]) for line in lines])
  y = np.array([float(line.split()[1]) for line in lines])
  
  return x,y

def version_np(path=path):
  data = np.loadtxt(path, usecols=(0,1))
  x = data[:,0]
  y = data[:,1]
  return x,y

def version2(path=path):
  
  with open(path) as file:
    lines = [line.split()[:2] for line in file.readlines() if line[0].isalnum()]
  x = [float(line[0]) for line in lines]
  y = [float(line[1]) for line in lines]
  
  return x,y