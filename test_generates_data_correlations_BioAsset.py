from os import listdir
from os.path import isfile, join
import numpy as np
import sys, getopt

test_name = "test_correlations_BioAsset.py"

test_usage_str = test_name + " -p <path_mats>"

meta_pos_filename = 2
meta_pos_start_data = 3
meta_pos_end_data = 20

argv = sys.argv[1:]

path_mats = ''

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

# one-to-one with slice_files
meta = np.load('bioAsset_meta.npy')

# subset of slice_files
mfs_data = np.load('mfs_BioAsset.npy')


result = np.array([])

i = 0
for slice_filename in slice_files:

    [patient_scan, _] = slice_filename.split('.')
    [meta_patient_scan, _] = meta[i][meta_pos_filename].split('.')

    if patient_scan == meta_patient_scan:
        data_i = np.array(mfs_data[i])

        for d in range(meta_pos_start_data, meta_pos_end_data + 1):
            data_i = np.hstack((data_i, np.array(meta[i][d])))

        if len(result) == 0:
            result = data_i
        else:
            result = np.vstack((result, data_i))

        i += 1

print "Shape: ", result.shape
np.save('mfs_and_standard_params.npy', result)