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

import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import misc
from scipy.stats import mode

model_directory = "model/"
filename_model_RF_tt=model_directory+"RF_tt_{}.pkl"
filename_model_SVC_tt=model_directory+"SVC_tt_{}.pkl"

sea_label=1
dolphin_label=2

dDFs  = 5

def compute_MFS(path_sea, path_dolphin):
    dir_sea  = os.listdir(path_sea)
    dir_dolphin = os.listdir(path_dolphin)

    cant_sea = len(dir_sea)
    cant_dolphin = len(dir_dolphin)

    seatrain = np.zeros((cant_sea, dDFs)).astype(np.float32)
    dolphintrain = np.zeros((cant_dolphin, dDFs)).astype(np.float32)

    ins = MFS()
    ins.setDef(1,dDFs,3)

    for i in range(cant_sea):
        filename = path_sea + dir_sea[i]
        seatrain[i] = ins.getFDs(filename)

    for i in range(cant_dolphin):
        filename = path_dolphin + dir_dolphin[i]
        dolphintrain[i] = ins.getFDs(filename)

    return seatrain, dolphintrain

computeMFS = False
seatrain86 = ""
seatrain64 = ""
seatrain32 = ""
seatrain16 = ""

dolphintrain86 = ""
dolphintrain64 = ""
dolphintrain32 = ""
dolphintrain16 = ""

data_path = "exps/data/"

if os.path.isfile(data_path+"seatrain86.npy")  \
and os.path.isfile(data_path+"dolphintrain86.npy"):
    seatrain86 = np.load(data_path+"seatrain86.npy")
    dolphintrain86 = np.load(data_path+"dolphintrain86.npy")

else:
    print "Computing MFSs at Resolution 86..."
    
    path_sea = '/home/deby/spyderws/imfractal-master/images/sea/sample-point-sea86x86/'
    path_dolphin = '/home/deby/spyderws/imfractal-master/images/sea/sample-point-dolphin86x86/'
      
    seatrain86, dolphintrain86 = compute_MFS(path_sea, path_dolphin)

    np.save(data_path+"seatrain86.npy", seatrain86)
    np.save(data_path+"dolphintrain86.npy", dolphintrain86)

if os.path.isfile(data_path+"seatrain64.npy") \
and os.path.isfile(data_path+"dolphintrain64.npy"):
    seatrain64 = np.load(data_path+"seatrain64.npy")
    dolphintrain64 = np.load(data_path+"dolphintrain64.npy")
else:
    print "Computing MFSs at Resolution 64..."
    
    path_sea = 'images/sea/sample-point-sea64x64/'
    path_dolphin = 'images/sea/sample-point-dolphin64x64/'
    
    seatrain64, dolphintrain64 = compute_MFS(path_sea, path_dolphin)
    
    np.save(data_path+"seatrain64.npy", seatrain64)
    np.save(data_path+"dolphintrain64.npy", dolphintrain64)

if os.path.isfile(data_path+"seatrain32.npy")  \
and os.path.isfile(data_path+"dolphintrain32.npy"):
    seatrain32 = np.load(data_path+"seatrain32.npy")
    dolphintrain32 = np.load(data_path+"dolphintrain32.npy")
else:
    print "Computing MFSs at Resolution 32..."
    
    path_sea = 'images/sea/sample-point-sea32x32/'
    path_dolphin = 'images/sea/sample-point-dolphin32x32/'
    
    seatrain32, dolphintrain32 = compute_MFS(path_sea, path_dolphin)
    
    np.save(data_path+"seatrain32.npy", seatrain32)
    np.save(data_path+"dolphintrain32.npy", dolphintrain32)

if os.path.isfile(data_path+"seatrain16.npy") \
and os.path.isfile(data_path+"dolphintrain16.npy"):
    seatrain16 = np.load(data_path+"seatrain16.npy")
    dolphintrain16 = np.load(data_path+"dolphintrain16.npy")
else:
    print "Computing MFSs at Resolution 16..."
    
    path_sea = 'images/sea/sample-point-sea16x16/'
    path_dolphin = 'images/sea/sample-point-dolphin16x16/'
    
    seatrain16, dolphintrain16 = compute_MFS(path_sea, path_dolphin)
    
    np.save(data_path+"seatrain16.npy", seatrain16)
    np.save(data_path+"dolphintrain16.npy", dolphintrain16)


def test_model_amount(train_size, seatrain, dolphintrain, resol):

    if train_size == 0:
        return

    cfr = RandomForestClassifier(n_estimators=100)
    cfr2 = svm.LinearSVC(C=1.0)

    data = np.vstack((seatrain, dolphintrain))
    labels_sea = np.zeros(len(seatrain)) + sea_label
    labels_dolphin = np.zeros(len(dolphintrain)) + dolphin_label
    labels = np.hstack((labels_sea, labels_dolphin))
    
    X_train_sea, X_test_sea, y_train_sea, y_test_sea = train_test_split(seatrain, labels_sea, train_size=min(train_size,len(seatrain)), random_state=0)
    X_train_do, X_test_do, y_train_do, y_test_do = train_test_split(dolphintrain, labels_dolphin, train_size=min(train_size,len(dolphintrain)), random_state=0)

    X_train = np.vstack((X_train_sea, X_train_do))
    X_test = np.vstack((X_test_sea, X_test_do))
    y_train = np.hstack((y_train_sea, y_train_do))
    y_test = np.hstack((y_test_sea, y_test_do))

    clf = cfr.fit(X_train, y_train)
    outdir_clf=filename_model_RF_tt.format(resol)
    if not os.path.exists(model_directory):
        os.makedirs(model_directory)
    #print outdir_clf
    joblib.dump(clf, outdir_clf) 
        
    clf2 = cfr2.fit(X_train, y_train)
    outdir_clf2=filename_model_SVC_tt.format(resol)
    #print outdir_clf2
    joblib.dump(clf2, outdir_clf2)
    
    print "RF, SVM: " + str( round(clf.score(X_test, y_test), 4) ) + " " + str( round(clf2.score(X_test, y_test), 4) )
    



def test_path(amount, seatrain, dolphintrain, resol):

    cant_sea = len(seatrain)
    cant_dolphin = len(dolphintrain)

    if amount > 0:
        cant_sea = min(amount, cant_sea)
        cant_dolphin = min(amount, cant_dolphin)

    seatrain = seatrain[:cant_sea]
    dolphintrain = dolphintrain[:cant_dolphin]


    cfr = RandomForestClassifier(n_estimators=100)
    cfr2 = svm.LinearSVC(C=1.0)

    data = np.vstack((seatrain, dolphintrain))

    labels_sea = np.zeros(len(seatrain)) + 1
    labels_dolphin = np.zeros(len(dolphintrain)) + 2
    labels = np.hstack((labels_sea, labels_dolphin))

    scores_rf = cross_validation.cross_val_score(cfr, data, labels, cv=4)
    scores_svm = cross_validation.cross_val_score(cfr2, data, labels, cv=4)
    
    print "RF, SVM: " + str( round(np.array(scores_rf).mean(),4) ) + " " + str( round(np.array(scores_svm).mean(),4) )

def test_all_resolutions(amount, func):
    print "86x86"
    func(amount, seatrain86, dolphintrain86, 86)
    
    print "64x64"
    func(amount, seatrain64, dolphintrain64, 64)
 
    print "32x32"
    func(amount, seatrain32, dolphintrain32, 32)
     
    print "16x16"
    func(amount, seatrain16, dolphintrain16, 16)

def test_train_predict():
    print " "
    print "########### Training-predicting test"

    for a in range(5,60,5):
        print ""
        print a, " samples"
        test_all_resolutions(a, test_model_amount)


def test_cross_val():
    print " "
    print "########### Results with 4-cross validation"

#    for a in range(0,60,5):
    for a in range(0,40,5):
        print ""
        if a == 0:
            print "All samples"
        else:
            print a, " samples"
        test_all_resolutions(a, test_path)
  
    
def do_test():

    test_cross_val()
    test_train_predict()
    
do_test()
