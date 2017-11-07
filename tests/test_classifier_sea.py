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
path_images = 'images/sea/first/'
data_path = "exps/data/"

filename_model_RF_tt=model_directory+"RF_tt_{}.pkl"
filename_model_SVC_tt=model_directory+"SVC_tt_{}.pkl"

sea_label=1
dolphin_label=2

dDFs  = 5

YIQ = False

# detect resolutions
dirs = os.listdir(path_images)
resolutions = []

for d in dirs:
    if 'sea' in d:
        _, resx = d.split('sea')
        res, _ = resx.split('x')
        resolutions.append(res)

resolutions.sort(reverse = True)

def transform_and_eq_hist(im):

    data = np.array(im.getdata()).reshape((im.size[0], im.size[1], 3))

    data = data / 255.0

    I = 0.595716*data[:,:,0] - 0.274453*data[:,:,1] - 0.321263*data[:,:,2]

    return exposure.equalize_hist(I)

def compute_MFS(path_sea, path_dolphin):
    dir_sea  = os.listdir(path_sea)
    dir_dolphin = os.listdir(path_dolphin)

    cant_sea = len(dir_sea)
    cant_dolphin = len(dir_dolphin)

    seatrain_i = np.zeros((cant_sea, dDFs)).astype(np.float32)
    dolphintrain_i = np.zeros((cant_dolphin, dDFs)).astype(np.float32)

    ins = MFS()
    ins.setDef(1,dDFs,3)

    for i in range(cant_sea):
        filename = path_sea + dir_sea[i]
        if(YIQ):
            # RGB -> YIQ
            im = Image.open(filename)
            im_eq = transform_and_eq_hist(im)
            
            seatrain_i[i] = ins.getFDs('', im_eq)
        else:
            seatrain_i[i] = ins.getFDs(filename)

    for i in range(cant_dolphin):
        filename = path_dolphin + dir_dolphin[i]
        if(YIQ):
            # RGB -> YIQ
            im = Image.open(filename)
            im_eq = transform_and_eq_hist(im)
            
            seatrain_i[i] = ins.getFDs('', im_eq)
        else:
            dolphintrain_i[i] = ins.getFDs(filename)

    return seatrain_i, dolphintrain_i




def prepare_dataset(seatrain, dolphintrain):
    i = 0 # FIX ME
    for r in resolutions:
        if os.path.isfile(data_path+"seatrain"+str(r)+".npy")  \
        and os.path.isfile(data_path+"dolphintrain"+str(r)+".npy"):
            print i, len(seatrain)
            seatrain[i] = np.load(data_path+"seatrain"+str(r)+".npy")
            dolphintrain[i] = np.load(data_path+"dolphintrain"+str(r)+".npy")

        else:
            print "Computing MFSs at Resolution "+str(r)+"..."
            
            path_sea = path_images+'sample-point-sea'+str(r)+'x'+str(r)+'/'
            path_dolphin = path_images+'sample-point-dolphin'+str(r)+'x'+str(r)+'/'
              
            seatrain[i], dolphintrain[i] = compute_MFS(path_sea, path_dolphin)

            np.save(data_path+"seatrain"+r+".npy", seatrain[i])
            np.save(data_path+"dolphintrain"+r+".npy", dolphintrain[i])

        i+=1

    
    return seatrain, dolphintrain
     


def test_model_amount(train_size, seatrain_i, dolphintrain_i, resol):

    if train_size == 0:
        return

    cfr = RandomForestClassifier(n_estimators=100)
    cfr2 = svm.LinearSVC(C=1.0)

    data = np.vstack((seatrain_i, dolphintrain_i))
    labels_sea = np.zeros(len(seatrain_i)) + sea_label
    labels_dolphin = np.zeros(len(dolphintrain_i)) + dolphin_label
    labels = np.hstack((labels_sea, labels_dolphin))
    
    X_train_sea, X_test_sea, y_train_sea, y_test_sea = train_test_split(seatrain_i, labels_sea, train_size=min(train_size,len(seatrain_i)), random_state=0)
    X_train_do, X_test_do, y_train_do, y_test_do = train_test_split(dolphintrain_i, labels_dolphin, train_size=min(train_size,len(dolphintrain_i)), random_state=0)

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
        



def test_path(amount, seatrain_i, dolphintrain_i, resol):

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

    scores_rf = cross_validation.cross_val_score(cfr, data, labels, cv=4)
    scores_svm = cross_validation.cross_val_score(cfr2, data, labels, cv=4)
    
    print "RF, SVM: " + str( round(np.array(scores_rf).mean(),4) ) + " " + str( round(np.array(scores_svm).mean(),4) )

def test_all_resolutions(amount, func, seatrain, dolphintrain):
    i = 0
    for r in resolutions:
        print r, 'x', r
        func(amount, seatrain[i], dolphintrain[i], r)
        i+=1

def test_train_predict(seatrain, dolphintrain):
    print " "
    print "########### Training-predicting test"

    for a in range(5,60,5):
        print ""
        print a, " samples"
        test_all_resolutions(a, test_model_amount, seatrain, dolphintrain)


def test_cross_val(seatrain, dolphintrain):
    print " "
    print "########### Results with 4-cross validation"

    for a in range(0,40,5):
        print ""
        if a == 0:
            print "All samples"
        else:
            print a, " samples"
        test_all_resolutions(a, test_path, seatrain, dolphintrain)
  
    
def do_test():
    seatrain = [[] for i in range(len(resolutions))]
    dolphintrain = [[] for i in range(len(resolutions))]

    seatrain, dolphintrain = prepare_dataset(seatrain, dolphintrain)
    test_cross_val(seatrain, dolphintrain)
    test_train_predict(seatrain, dolphintrain)
    
