import numpy as np
import os
import argparse
from imfractal import *


def process(args):
    directory = args.dir[0]

    classes = os.listdir(directory)

    mfs = MFS()
    mfs.setDef(1,20,3)

    total = [[] for c in range(len(classes))]
    j = 0
    for c in classes:

        imgs = os.listdir(directory+c)
        print imgs
        for i in imgs:
            total[j].append(mfs.getFDs(directory+c+'/'+i))
            
        j+=1
        
        np.save("suturas.npy", np.asarray(total))

    print mfss

def main():
    
    parser = argparse.ArgumentParser(description='Compute SOM of classes in folders representing labels')
    parser.add_argument("-d", dest="dir", type=str, required=True, nargs=1, help="Directory with each folder represents a class")
    
    args = parser.parse_args()

    process(args)

if __name__ == '__main__':
    main()
