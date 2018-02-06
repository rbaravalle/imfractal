import os
import numpy as np
import Image
import scipy.misc
import argparse


parser = argparse.ArgumentParser(description='Compare binarized files with GT')

parser.add_argument("-i", dest="images", type=str, required=True, nargs='+', help="binarized images")
parser.add_argument("-g", dest="gt", type=str, required=True, nargs='+', help="ground truth images")

args = parser.parse_args()

gt_dir = args.gt[0]+'/'
classified_dir = args.images[0]+'/'

gt_files = os.listdir(gt_dir)
classified_files = os.listdir(classified_dir)

amount_d = 0
corrects_d = 0

amount_sea = 0
corrects_sea = 0

for c in classified_files:    
 
    print "Processing:", c
    values = c.split('_')
    c_filename = values[0]+'_'+values[1]

    data_c = scipy.misc.imread(classified_dir+c)
    data_gt = []

    is_in_gt = False
    # Search in GT dir
    for gt_f in gt_files:
        if c_filename in gt_f:
            is_in_gt = True
            data_gt = scipy.misc.imread(gt_dir+gt_f)
            break

    if not(is_in_gt): 
        print "Not in gt"
        continue

    for i in range(data_c.shape[0]):
        for j in range(data_c.shape[1]):

            v = data_c[i, j]
            r_gt,g_gt,b_gt = data_gt[i,j]

            # if gt is dolphin
            if r_gt == 255 and g_gt == 255 and b_gt == 255:
                amount_d+=1
                # if the classifier says it is dolphin, accuracy increases
                if v == 255:
                    corrects_d+=1
            else:
                # gt sea
                amount_sea += 1
                if v != 255:
                    corrects_sea+=1

    print "Actual Accuracy Dolphin: ", corrects_d/float(amount_d)
    print "Actual Accuracy Sea: ", corrects_sea/float(amount_sea)
            

print "FINAL RESULT:"
print "Corrects: ", corrects_d
print "Amount: ", amount_d
print "Accuracy: ", corrects_d/float(amount_d)
            
print "SEA"
print "Corrects: ", corrects_sea
print "Amount: ", amount_sea
print "Accuracy: ", corrects_sea/float(amount_sea)
