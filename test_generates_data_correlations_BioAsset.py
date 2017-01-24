from os import listdir
from os.path import isfile, join
import numpy as np
import sys, getopt
import csv
from imfractal import *

test_usage_str = sys.argv[0] + " -p <path_mats>"

meta_pos_filename = 2
meta_pos_start_data = 2
meta_pos_end_data = 20

argv = sys.argv[1:]

path_mats = ''

TRANSFORMED_INPUT_STR = ''

if MFS_HOLDER :
    MFS_STR = '_holder'

if APPLY_LAPLACIAN:
    TRANSFORMED_INPUT_STR = '_laplacian'
else:
    if APPLY_GRADIENT:
        TRANSFORMED_INPUT_STR = '_gradient'

try:
    opts, args = getopt.getopt(argv, "h:p:", ["path_mats="])
except getopt.GetoptError:
    print test_usage_str
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print test_usage_str
        sys.exit()
    elif opt in ("-p", "--path_mats"):
        path_mats = arg

print path_mats

if path_mats == '':
    print "Please specify path to matlab matrices with option -p"
    print test_usage_str
    exit()

slice_files = [f for f in listdir(path_mats) if isfile(join(path_mats, f)) and "Slices" in f]

slice_files.sort()  # = sort(slice_files)

path = 'exps/data/'
# one-to-one with slice_files
#meta = np.load(path + 'bioAsset_meta.npy')

# subset of slice_files
print path + 'mfs'+MFS_STR+TRANSFORMED_INPUT_STR +'_BioAsset.npy'
mfs_data = np.load(path + 'mfs'+MFS_STR+TRANSFORMED_INPUT_STR +'_BioAsset.npy')


result = np.array([])

i = 0

f = open(path + 'mfs'+MFS_STR+TRANSFORMED_INPUT_STR +'_BioAsset.csv', 'wt')

writer = csv.writer(f)

line = ('Filename',)

line += tuple(map(lambda x: "MFS " + str(x), range(mfs_data.shape[1])))

writer.writerow(line)

for slice_filename in slice_files:

    [patient_scan, _] = slice_filename.split('.')
    [meta_patient_scan, _] = meta[i][meta_pos_filename].split('.')

    if patient_scan == meta_patient_scan:
        data_i = np.array(mfs_data[i])

        import math
        if not(math.isnan(float(meta[i][3]))):
            for d in range(meta_pos_start_data, meta_pos_end_data + 1):
                data_i = np.hstack((data_i, np.array(meta[i][d])))

            if len(result) == 0:
                result = data_i
            else:
                result = np.vstack((result, data_i))

            line = (patient_scan,)
            line += tuple(mfs_data[i])
            writer.writerow(line)

        i += 1

print "Shape: ", result.shape
np.save(path+'mfs'+MFS_STR+ TRANSFORMED_INPUT_STR +'_and_standard_params.npy', result)