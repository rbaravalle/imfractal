"""
Copyright (c) 2016 Rodrigo Baravalle
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from imfractal import *

from pylab import *

import sys
import os
import scipy.stats

sys.path.append(os.path.join(os.path.dirname(sys.path[0]), 'imfractal', 'imfractal',  "Algorithm"))

import qs3D
from numpy import recfromcsv

def stats_test(_output_filename, input_filename):
    _,  extension = os.path.splitext(input_filename)

    if extension == 'csv':
        data = recfromcsv(input_filename, delimiter=',')
        first_f = 1
    else:
        data = np.load(input_filename)
        first_f = 0

    dims = 15
    mfss = np.zeros([len(data), dims])

    for i in range(len(data)):

        print i

        if extension == 'csv':
            mfs = np.array(tuple(data[i])[first_f:])#.astype(np.float32)
        else:
            mfs = data[i]

        if len(mfs) < 100:
            max_fa = np.max(mfs)
            min_fa = np.min(mfs)
            std_fa = np.std(mfs)
            mean_fa = np.mean(mfs)
            median_fa = np.median(mfs)
            sum_fa = np.sum(mfs)
            variation = scipy.stats.variation(mfs)
            var = scipy.stats.tvar(mfs)
            skew = scipy.stats.skew(mfs)
            kurtosis = scipy.stats.kurtosis(mfs)
            arg_max = np.argmax(mfs)
            arg_min = np.argmin(mfs)
            diff = np.max(mfs) - np.min(mfs)
            first_f = mfs[0]
            last_f = mfs[-1]

            mfss[i] = np.array([max_fa, min_fa,
                                mean_fa, std_fa,
                                median_fa, sum_fa,
                                skew, kurtosis,
                                variation, var,
                                arg_max, arg_min, diff,
                                last_f, first_f])


        else:
            tmp = np.array([])

            for l in range(5):
                mmfs = mfs[l*20 : (l+1)*20 - 1]
                max_fa = np.max(mmfs)
                min_fa = np.min(mmfs)
                std_fa = np.std(mmfs)
                mean_fa = np.mean(mmfs)
                median_fa = np.median(mmfs)
                sum_fa = np.sum(mmfs)
                variation = scipy.stats.variation(mmfs)
                var = scipy.stats.tvar(mmfs)
                skew = scipy.stats.skew(mmfs)
                kurtosis = scipy.stats.kurtosis(mmfs)
                arg_max = np.argmax(mmfs)
                arg_min = np.argmin(mmfs)

                tmp = np.hstack((tmp,
                            np.array([max_fa, min_fa,
                                    mean_fa, std_fa,
                                    median_fa, sum_fa,
                                    skew, kurtosis,
                                    variation, var, arg_max, arg_min])))


           # for l in range(5):
            #    skews[l] = scipy.stats.skew(mfs[l*20 : (l+1)*20 - 1])
            #    kurtosiss[l] = scipy.stats.kurtosis(mfs[l*20 : (l+1)*20 - 1])

            #skew = np.mean(skews)
            #kurtosis = np.mean(kurtosiss)

            mfss[i] = tmp


        print mfss[i], i

        np.save(data_path + _output_filename, mfss)

    print "Saved ", data_path + _output_filename

def do_test(_path, _output_filename):

    dims = 21 # should be odd number! to include q = -x , ... q = 0, ..., q = x
    if MFS_HOLDER:
        dims = 20 # Holder 3D MFS
        if LOCAL:
            dims = 20

    if Stats_MFS:
        dims = 10

    # BioAsset bone's multifractal spectra database


    #slices_str = "slice"
    #masks_str = "mask"

    params = {
        "zero": 1,
        "one": 0.75,
        "two": 3.7,
        "three": 1,
        "four": 15,
        "five": 0,
        "mask_filename": '',
        "seven": "no",
        "eight": 'S',
        "nine": 'M',
        "threshold": 200,
        "total_pixels":6000,
        "adaptive" : False, # adaptive threshold (only for not holder)
        "laplacian": APPLY_LAPLACIAN,
        "gradient" : APPLY_GRADIENT
    }

    from os import listdir
    from os.path import isfile, join
    mask_files = [f for f in listdir(_path) if isfile(join(_path, f)) and "Mask" in f]
    slice_files = [f for f in listdir(_path) if isfile(join(_path, f)) and "Slices" in f]

    mfss = np.zeros([len(mask_files), dims])

    mask_files = sort(mask_files)
    slice_files = sort(slice_files)

    if len(mask_files) != len(slice_files):
        print "The directory should contain the same amount of slices and masks"
        exit()

    i = 0
    for mask_filename in mask_files:
        [patient_scan_str, _] = mask_filename.split("Mask")
        [first_str, scan_str] = patient_scan_str.split("_120_")
        [_, patient_str] = first_str.split("BA")

        mask_filename = _path + mask_filename

        params["five"] = 1 #fix me
        params["mask_filename"] = mask_filename

        # obviously we can directly use slice_files[i], but this adds robustness
        slice_filename = _path + "BA" + patient_str + "_120_" + scan_str + "Slices.mat"
        if slice_filename == _path + slice_files[i]:
            print "MASK: ", mask_filename
            print "SLICE: ", slice_filename
        else:
            print "Cannot process test: filename ", _path + slice_files[i], " should be ", slice_filename
            exit()

        if not(MFS_HOLDER):
            aux = CSandbox3D(dims)
            aux.setDef(40, 1.02, True, params)
            mfss[i] = aux.getFDs(slice_filename)
        else:
            aux = MFS_3D()
            if LOCAL:
                #aux = Local_MFS_3D()
                aux = Local_MFS_Pyramid_3D()
                aux.setDef(1, 20, 3, slice_filename, mask_filename, params)
                mfss[i] = aux.getFDs()
            #if Stats_MFS:
            #    aux = Stats_MFS_3D()
            #    aux.setDef(1, 20, 3, slice_filename, mask_filename, params)
            #    mfss[i] = aux.getFDs()

            #if SLICES_MFS:
            #    aux = MFS_3D_Slices()
            #    aux.setDef(1, dims, 3, slice_filename, mask_filename, params)
            #    ax = 2 # X axis
            #    mfss[i] = aux.getFDs(ax)
            #else:
            #    aux.setDef(1, dims, 3, slice_filename, mask_filename, params)
            #    mfss[i] = aux.getFDs()

        # in case something goes wrong, save computed mfs up to here
        print "Data partially saved to ", data_path + _output_filename + ".npy"

        np.save(data_path + _output_filename, mfss)

        i += 1


    print "Data saved to ", data_path + _output_filename + ".npy"
    np.save(data_path + _output_filename, mfss)


    
