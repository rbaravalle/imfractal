"""
Copyright (c) 2013 Rodrigo Baravalle
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
import Image
import time
import csv
import sys
import os
from subprocess import *
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
from sklearn import svm
import matplotlib.pyplot as plt
from pylab import *

def do_test():
    #cant = 10
    dDFs  = 20
    cantClasses = 2

    pathbtr = 'images/train2/bread/'
    dirListbtr=os.listdir(pathbtr)
    cantTrainB = len(dirListbtr)
    pathnbtr = 'images/train2/nonbread/'
    dirListnbtr=os.listdir(pathnbtr)
    cantTrainNB = len(dirListnbtr)
    pathbte = 'images/test2/bread/'
    dirListbte=os.listdir(pathbte)
    cantTestB = len(dirListbte)
    pathnbte = 'images/test2/nonbread/'
    dirListnbte=os.listdir(pathnbte)
    cantTestNB = len(dirListnbte)

    breadtrain = np.zeros((cantTrainB, dDFs)).astype(np.float32)
    breadtest = np.zeros((cantTestB, dDFs)).astype(np.float32)

    nonbreadtrain = np.zeros((cantTrainNB, dDFs)).astype(np.float32)
    nonbreadtest = np.zeros((cantTestNB, dDFs)).astype(np.float32)
    
    ins = MFS()

    print 'Training: computing MFS for the bread database...'
    ins.setDef(1,20,3,True)
    for i in range(cantTrainB):
        filename = pathbtr+dirListbtr[i]
        breadtrain[i] = ins.getFDs(filename)
    for i in range(cantTestB):
        filename = pathbte+dirListbte[i]
        breadtest[i] = ins.getFDs(filename)
    for i in range(cantTrainNB):
        filename = pathnbtr+dirListnbtr[i]
        nonbreadtrain[i] = ins.getFDs(filename)
    for i in range(cantTestNB):
        filename = pathnbte+dirListnbte[i]
        nonbreadtest[i] = ins.getFDs(filename)


    data = np.vstack((breadtrain,breadtest,nonbreadtrain))
    labelsbtr = np.array([1 for i in range(len(breadtrain))])
    labelsnbtr = np.array([2 for i in range(len(nonbreadtrain))])

    labelsbte = np.array([1 for i in range(len(breadtest))])
    labelsnbte = np.array([2 for i in range(len(nonbreadtest))])

    labels = np.hstack((labelsbtr,labelsbte,labelsnbtr,labelsnbte))

    print "Testing..."
    print "1 = Bread"
    print "2 = Nonbread"

    train = np.vstack((breadtrain,nonbreadtrain))
    labels = np.hstack((labelsbtr,labelsnbtr))
    test = np.vstack((breadtest,nonbreadtest))

    lin_svc = svm.LinearSVC(C=1.0).fit(train, labels)
    predictionsSVM = lin_svc.predict(test)

    cfr = RandomForestClassifier(n_estimators=120)
    cfr.fit(train,labels) # train

    gtruth = np.hstack((labelsbte,labelsnbte))
    predictionsRF = cfr.predict(test) # test

    print dirListbte
    print "Random Forest Prediction:"
    print predictionsRF
    print "SVM Prediction:"
    print predictionsSVM
    print "REAL: "
    print gtruth

    x = np.arange(dDFs)
    fsize = 14

    plt.ylabel(r'$R$',fontsize=fsize)
    plt.xlabel('FD',fontsize=fsize)
    plt.plot(x, breadtrain[3:6,:].T, 'k+--', label='bread train',linewidth=2.0)
    plt.plot(x, breadtest.T, 'r*--',  label='bread test',linewidth=2.0)
    plt.legend(loc = 0)
    plt.show()

    scoreRF = (len(gtruth)-sum(abs(gtruth-predictionsRF)))/float(len(gtruth))
    scoreSVM = (len(gtruth)-sum(abs(gtruth-predictionsSVM)))/float(len(gtruth))

    #scores = cross_validation.cross_val_score(cfr, data, labels, cv=4)
    print "Classification performance (Random Forest classifier): " + str( scoreRF*100 ) + "%"
    print "Classification performance (Support Vector Machine classifier): " + str( scoreSVM*100 ) + "%"



