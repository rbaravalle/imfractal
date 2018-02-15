"""
Copyright (c) 2017 Rodrigo Baravalle
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

from imfractal import MFS
from PIL import Image
import time
import csv
import sys
import os
from subprocess import *

from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib
from skimage import exposure
from skimage import io, color

import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import misc
from scipy.stats import mode
import argparse

model_directory = "model/"

# resolution, transform, amount of dfs
filename_model_RF_tt=model_directory+"RF_tt_{}_{}_{}.pkl"
filename_model_SVC_tt=model_directory+"SVC_tt_{}_{}_{}.pkl"

sea_label=1
dolphin_label=2

transformation_values = ["YIQ", "R", "G", "B", "lab"]
true_values = ['t', 'T', '1', 1, 'true', 'True', 'TRUE']



def transform_yiq(filename, ins, equalize):
    im = Image.open(filename)

    data = np.array(im.getdata()).reshape((im.size[0], im.size[1], 3))

    data = data / 255.0

    I = 0.595716*data[:,:,0] - 0.274453*data[:,:,1] - 0.321263*data[:,:,2]

    if equalize in true_values:
      I = exposure.equalize_hist(I)

    return ins.getFDs('', I)


def transform_r(filename, ins, equalize):
    im = Image.open(filename)

    data = np.array(im.getdata()).reshape((im.size[0], im.size[1], 3))

    data = data / 255.0

    R = data[:,:,0]

    if equalize in true_values:
      R = exposure.equalize_hist(R)

    return ins.getFDs('', R)

def transform_g(filename, ins, equalize):
    im = Image.open(filename)

    data = np.array(im.getdata()).reshape((im.size[0], im.size[1], 3))

    data = data / 255.0

    G = data[:,:,1]

    if equalize in true_values:
      G = exposure.equalize_hist(G)

    return ins.getFDs('', G)

def transform_b(filename, ins, equalize):
    im = Image.open(filename)

    data = np.array(im.getdata()).reshape((im.size[0], im.size[1], 3))

    data = data / 255.0

    B = data[:,:,2]

    if equalize in true_values:
      B = exposure.equalize_hist(B)

    return ins.getFDs('', B)

# Lab color space
def transform_lab(filename, ins, equalize):
    rgb = io.imread(filename)

    lab = color.rgb2lab(rgb)

    l = lab[:,:,0]
    a = lab[:,:,1]
    b = lab[:,:,2]

    if equalize in true_values:
      l = exposure.equalize_hist(l)
      a = exposure.equalize_hist(a)
      b = exposure.equalize_hist(b)
      
    return ins.getFDs('', l) + ins.getFDs('', a) + ins.getFDs('', b)

# transform image to other color space given filename
def transform_f(filename, transformation, ins, equalize):

    if transformation == "YIQ":
        return transform_yiq(filename, ins, equalize)

    if transformation == "R":
        return transform_r(filename, ins, equalize)

    if transformation == "G":
        return transform_g(filename, ins, equalize)

    if transformation == "B":
        return transform_b(filename, ins, equalize)

    if transformation == "lab":
        return transform_lab(filename, ins, equalize)

    print "ERROR: inexistent transformation"
    exit()

def compute_MFS(path_sea, path_dolphin, args):
    dfs = args.dfs[0]
    transform = args.transform[0]

    dir_sea  = os.listdir(path_sea)
    dir_dolphin = os.listdir(path_dolphin)

    cant_sea = len(dir_sea)
    cant_dolphin = len(dir_dolphin)

    seatrain_i = [[] for i in range(cant_sea)] #np.zeros((cant_sea, dfs)).astype(np.float32)
    dolphintrain_i = [[] for i in range(cant_dolphin)] #np.zeros((cant_dolphin, dfs)).astype(np.float32)

    ins = MFS()
    ins.setDef(1,dfs,3)

    for i in range(cant_sea):
        filename = path_sea + dir_sea[i]
        if(transform in transformation_values):
            seatrain_i[i] = transform_f(filename, transform, ins, args.equalize[0])
        else:
            seatrain_i[i] = ins.getFDs(filename)

    for i in range(cant_dolphin):
        filename = path_dolphin + dir_dolphin[i]
        if(transform in transformation_values):
            dolphintrain_i[i] = transform_f(filename, transform, ins, args.equalize[0])
        else:
            dolphintrain_i[i] = ins.getFDs(filename)

    return np.array(seatrain_i).reshape(cant_sea,len(seatrain_i[0])), np.array(dolphintrain_i).reshape(cant_dolphin,len(dolphintrain_i[0]))




def prepare_dataset(seatrain, dolphintrain, resolutions, args):
    data_path = args.data_path[0]
    path_images = args.path_images[0]

    print "Using", args.dfs[0], "MFS features"
    print "Using Data Path:", data_path
    print "Using Images Path:", path_images



    if args.transform[0] in transformation_values:
        print "Computing", args.transform[0], "transformation"
    else:
        print "No transformation computed"

    for i in range(len(resolutions)):
        r = resolutions[i]
        #if os.path.isfile(data_path+"seatrain"+str(r)+".npy")  \
        #and os.path.isfile(data_path+"dolphintrain"+str(r)+".npy"):

        #    seatrain[i] = np.load(data_path+"seatrain"+str(r)+".npy")
        #    dolphintrain[i] = np.load(data_path+"dolphintrain"+str(r)+".npy")

        #else:
        print "Computing MFSs at Resolution "+str(r)+"..."
        
        path_sea = path_images+'sample-point-sea'+str(r)+'x'+str(r)+'/'
        path_dolphin = path_images+'sample-point-dolphin'+str(r)+'x'+str(r)+'/'
          
        seatrain[i], dolphintrain[i] = compute_MFS(path_sea, path_dolphin, args)

        #np.save(data_path+"seatrain"+r+".npy", seatrain[i])
        #np.save(data_path+"dolphintrain"+r+".npy", dolphintrain[i])

    
    return seatrain, dolphintrain
     


def test_model_amount(train_percentage, seatrain_i, dolphintrain_i, resol, args):

    if train_percentage == 0:
        print "ZERO!"
        return

    cfr = RandomForestClassifier(n_estimators=100)
    cfr2 = svm.LinearSVC(C=1.0)

    data = np.vstack((seatrain_i, dolphintrain_i))
    labels_sea = np.zeros(len(seatrain_i)) + sea_label
    labels_dolphin = np.zeros(len(dolphintrain_i)) + dolphin_label
    labels = np.hstack((labels_sea, labels_dolphin))
    
    X_train_sea, X_test_sea, y_train_sea, y_test_sea = train_test_split(seatrain_i, labels_sea, train_size=train_percentage, random_state=0)
    X_train_do, X_test_do, y_train_do, y_test_do = train_test_split(dolphintrain_i, labels_dolphin, train_size=train_percentage, random_state=0)

    X_train = np.vstack((X_train_sea, X_train_do))
    X_test = np.vstack((X_test_sea, X_test_do))
    y_train = np.hstack((y_train_sea, y_train_do))
    y_test = np.hstack((y_test_sea, y_test_do))

    tr = args.transform[0]
    if tr in transformation_values:
        transform_str = tr
        if args.equalize[0] in true_values:
            transform_str += '_eq'
    else:
        transform_str = "no_transform"
    dfs_str = str(args.dfs[0])

    clf = cfr.fit(X_train, y_train)
    outdir_clf=filename_model_RF_tt.format(resol, transform_str, dfs_str)
    if not os.path.exists(model_directory):
        os.makedirs(model_directory)
    print "Saving classifier...",  outdir_clf
    joblib.dump(clf, outdir_clf) 
        
    clf2 = cfr2.fit(X_train, y_train)
    outdir_clf2=filename_model_SVC_tt.format(resol, transform_str, dfs_str)
    print "Saving classifier...", outdir_clf2
    joblib.dump(clf2, outdir_clf2)
    
    print "RF, SVM: " + str( round(clf.score(X_test, y_test), 4) ) + " " + str( round(clf2.score(X_test, y_test), 4) )

def test_path(amount, seatrain_i, dolphintrain_i, resol, args):

    cant_sea = len(seatrain_i)
    cant_dolphin = len(dolphintrain_i)

    if amount > 0:
        cant_sea = min(amount, cant_sea)
        cant_dolphin = min(amount, cant_dolphin)

    seatrain_i = seatrain_i[:cant_sea]
    dolphintrain_i = dolphintrain_i[:cant_dolphin]


    cfr = RandomForestClassifier(n_estimators=100)
    cfr2 = svm.LinearSVC(C=1.0)

    data = np.vstack((seatrain_i, dolphintrain_i))

    labels_sea = np.zeros(len(seatrain_i)) + 1
    labels_dolphin = np.zeros(len(dolphintrain_i)) + 2
    labels = np.hstack((labels_sea, labels_dolphin))

    scores_rf = cross_validation.cross_val_score(cfr, data, labels, cv=10)
    scores_svm = cross_validation.cross_val_score(cfr2, data, labels, cv=10)
    
    print "RF, SVM: " + str( round(np.array(scores_rf).mean(),4) ) + " " + str( round(np.array(scores_svm).mean(),4) )

def test_all_resolutions(amount, func, seatrain, dolphintrain, resolutions, args):
    for i in range(len(resolutions)):
        r = resolutions[i]
        print r, 'x', r
        func(amount, seatrain[i], dolphintrain[i], r, args)

def test_train_predict(seatrain, dolphintrain, resolutions, args):
    print " "
    print "########### Training-predicting test"

    percentage_train = args.percentage_train[0]
    test_all_resolutions(percentage_train, test_model_amount, seatrain, dolphintrain, resolutions, args)


def test_cross_val(seatrain, dolphintrain, resolutions, args):
    print " "
    print "########### Results with 10-cross validation (90% train)"
    print ""
    print "All samples"

    test_all_resolutions(0, test_path, seatrain, dolphintrain, resolutions, args)
  
    
def do_test():
    parser = argparse.ArgumentParser(description='Binarize an image using classifiers models')
    parser.add_argument("-imgs", dest="path_images", type=str, required=True, nargs=1, help="Path with images to be tested")
    parser.add_argument("-data", dest="data_path", type=str, required=True, nargs=1, help="Path where data will be saved")
    parser.add_argument("-dfs", dest="dfs", type=int, required=True, nargs=1, help="Amount of MFS dimensions per MFS")
    parser.add_argument("-tr", dest="transform", type=str, required=True, nargs=1, help="Convert data. Specify transformation " + ', '.join(transformation_values) + " or no transformation")
    parser.add_argument("-eq", dest="equalize", type=str, required=True, nargs=1, help="Equalize image or not for training")
    parser.add_argument("-ptrain", dest="percentage_train", type=float, required=True, nargs=1, help="Percentage of dataset to be used for training (float)")
    
    args = parser.parse_args()

    # detect resolutions
    dirs = os.listdir(args.path_images[0])
    resolutions = []

    for d in dirs:
        if 'sea' in d:
            _, resx = d.split('sea')
            res, _ = resx.split('x')
            resolutions.append(res)

    resolutions.sort(reverse = True)


    seatrain = [[] for i in range(len(resolutions))]
    dolphintrain = [[] for i in range(len(resolutions))]

    seatrain, dolphintrain = prepare_dataset(seatrain, dolphintrain, resolutions, args)

    test_cross_val(seatrain, dolphintrain, resolutions, args)
    test_train_predict(seatrain, dolphintrain, resolutions, args)
    
