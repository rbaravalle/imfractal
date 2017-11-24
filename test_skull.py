import numpy as np
import os
import argparse
from imfractal import *


def process(args):
    directory = args.dir[0]
    cant_dims = int(args.cant_dims[0])

    classes = os.listdir(directory)

    mfs = MFS()
    mfs.setDef(1,cant_dims,3)

    total = [[] for c in range(len(classes))]
    j = 0
    for c in classes:

        imgs = os.listdir(directory+c)
        print "Computing MFS in dir", directory+c
        for i in imgs:
            total[j].append(mfs.getFDs(directory+c+'/'+i))
            
        j+=1
        
        np.save("suturas_"+str(cant_dims)+".npy", np.asarray(total))

    print "Finished computing. Saved suturas_"+str(cant_dims)+".npy"

def main():
    
    parser = argparse.ArgumentParser(description='Compute MFS of folder (subfolders are classes)')
    parser.add_argument("-d", dest="dir", type=str, required=True, nargs=1, help="Directory with each folder represents a class")
    parser.add_argument("-cant_dims", dest="cant_dims", type=str, required=True, nargs=1, help="Amount of dimensions in the MFS")
    
    args = parser.parse_args()

    process(args)

if __name__ == '__main__':
    main()
