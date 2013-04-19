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

def do_test():
    cant = 10
    dDFs  = 20
    cantClasses = 2
    breadtrain = np.zeros((cant, dDFs)).astype(np.float32)
    breadtest = np.zeros((cant, dDFs)).astype(np.float32)

    nonbreadtrain = np.zeros((cant, dDFs)).astype(np.float32)
    nonbreadtest = np.zeros((cant, dDFs)).astype(np.float32)

    pathbtr = 'images/train/bread/'
    dirListbtr=os.listdir(pathbtr)
    pathnbtr = 'images/train/nonbread/'
    dirListnbtr=os.listdir(pathnbtr)
    pathbte = 'images/test/bread/'
    dirListbte=os.listdir(pathbte)
    pathnbte = 'images/test/nonbread/'
    dirListnbte=os.listdir(pathnbte)
    #print len(dirListbtr), dirListbtr
    
    ins = MFS()

    print 'Training: calculating MFS for the bread database...'
    for i in range(cant):
        ins.setDef(1,20,3,True)
        filename = pathbtr+dirListbtr[i]
        breadtrain[i] = ins.getFDs(filename)
        filename = pathbte+dirListbte[i]
        breadtest[i] = ins.getFDs(filename)

        filename = pathnbtr+dirListnbtr[i]
        nonbreadtrain[i] = ins.getFDs(filename)
        filename = pathnbte+dirListnbte[i]
        nonbreadtest[i] = ins.getFDs(filename)

    cfr = RandomForestClassifier(n_estimators=100)
    data = np.vstack((breadtrain,breadtest,nonbreadtrain,nonbreadtest))
    labels = np.zeros((len(data),1)) # FIX ME
    for i in range(len(data)):
        labels[i] = i
    labels = map(lambda i: np.floor(i/(2*(cant)))+1, labels)
    labels = np.array(labels)
    labels = np.transpose(labels)[0]   # FIX ME
    print "Testing..."
    scores = cross_validation.cross_val_score(cfr, data, labels, cv=4)
    print "Classification performance (Random Forest classifier): " + str( np.array(scores).mean() )



