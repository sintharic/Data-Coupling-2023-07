from prepare import *

def main():
  for folder in folders_v0: prepare_flowX(folder)
  for folder in folders_vX: prepare_flowX(folder)
  for folder in folders_vY: prepare_flowX(folder)

if __name__ == '__main__': main()