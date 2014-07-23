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
import csv

def opencsv(filename):
    array = []
    with open(filename, 'rb') as csvfile:
        spamreader = csv.reader(csvfile)
        i = 0
        for row in spamreader:
            array[i] = row
            i = i+1

    return np.array(array)

# CSV WRITE
def writecsv(filename,array):
    with open(filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(array)

def do_test(minn,vals):
    #cant = 10
    dfs = 10
    dDFs  = 2*dfs+1
    cantClasses = 2
    computeMFS = False

    pathbtr = 'images/train2/bread/'
    dirListbtr=os.listdir(pathbtr)
    cantTrainB = 10#len(dirListbtr)
    #pathnbtr = 'images/train2/nonbread/'
    #dirListnbtr=os.listdir(pathnbtr)
    #cantTrainNB = 1#len(dirListnbtr)
    pathbte = 'images/test4/bread/'
    dirListbte=os.listdir(pathbte)
    cantTestB = 2#len(dirListbte)
    #pathnbte = 'images/test2/nonbread/'
    #dirListnbte=os.listdir(pathnbte)
    #cantTestNB = 1#len(dirListnbte)

    breadtrain = np.zeros((cantTrainB, dDFs)).astype(np.float32)
    print breadtrain[0].shape
    breadtest = np.zeros((cantTestB, dDFs)).astype(np.float32)

    #nonbreadtrain = np.zeros((cantTrainNB, dDFs)).astype(np.float32)
    #nonbreadtest = np.zeros((cantTestNB, dDFs)).astype(np.float32)
    
    ins = Sandbox(dfs)

    if(computeMFS):
        print 'Training: computing sandbox MFS for the bread database...'
        ins.setDef(40,1.02,True)
        #print "Computing " + str(cantTrainB) +" bread train..."
        for i in range(cantTrainB):
            filename = pathbtr+dirListbtr[i]
            breadtrain[i] = ins.getFDs(filename)
    else:
        print "Loading CSV"
        with open('breadtrainS.csv', 'rb') as csvfile:
            spamreader = csv.reader(csvfile)
            i = 0
            for row in spamreader:
                breadtrain[i] = row[1:]
                i = i+1

    #print "Computing "+ str(cantTestB) +" bread test..."
    ins.setDef(40,1.15,False)
    for i in range(cantTestB):
        filename = pathbte+dirListbte[i]
        breadtest[i] = ins.getFDs(filename)
    #ins.setDef(40,1.02,True)
    #print "Computing "+str(cantTrainNB)+" non bread train..."
    #for i in range(cantTrainNB):
    #    filename = pathnbtr+dirListnbtr[i]
    #    nonbreadtrain[i] = ins.getFDs(filename)
    #ins.setDef(40,1.02,True)
    #print "Computing "+ str(cantTestNB) +" non bread test..."
    #for i in range(cantTestNB):
    #    filename = pathnbte+dirListnbte[i]
    #    nonbreadtest[i] = ins.getFDs(filename)

    fsize = 14

    if(computeMFS):
        #data = np.vstack((breadtrain,breadtest))#,nonbreadtrain))
        labelsbtr = np.zeros((len(breadtrain),1)) + 1
        #labelsnbtr = np.zeros((len(nonbreadtrain),1)) + 2

        #labelsbte = np.zeros((len(breadtest),1)) + 1
        #labelsnbte = np.zeros((len(nonbreadtest),1)) + 2

        #labels = np.hstack((labelsbtr[:,0],labelsbte[:,0],labelsnbtr[:,0],labelsnbte[:,0]))

        print "Saving CSVs for SOM"

        print "Shapes: labelsbtr: ", labelsbtr.shape, "breadtrain: ",breadtrain.shape

        writecsv('breadtrainS.csv',np.hstack((labelsbtr,breadtrain)) )
        #writecsv('breadtestS.csv',np.hstack((labelsbte,breadtest)) )
        #writecsv('nonbreadtrainS.csv',np.hstack((labelsnbtr,nonbreadtrain)) )
        #writecsv('nonbreadtestS.csv',np.hstack((labelsnbte,nonbreadtest)) )
        #exit()

    if(False):
        #print "Testing..."
        #print "1 = Bread"
        #print "2 = Nonbread"

        train = np.vstack((breadtrain,nonbreadtrain))
        labels = np.hstack((labelsbtr[:,0],labelsnbtr[:,0]))
        test = np.vstack((breadtest,nonbreadtest))

        lin_svc = svm.LinearSVC(C=1.0).fit(train, labels)
        predictionsSVM = lin_svc.predict(test)

        cfr = RandomForestClassifier(n_estimators=120)
        cfr.fit(train,labels) # train

        gtruth = np.hstack((labelsbte[:,0],labelsnbte[:,0]))
        predictionsRF = cfr.predict(test) # test

        print dirListbte
        print dirListbtr
        print "Random Forest Prediction:"
        print predictionsRF[:cantTestB]
        print predictionsRF[cantTestB:]
        print "SVM Prediction:"
        print predictionsSVM[:cantTestB]
        print predictionsSVM[cantTestB:]
        print "REAL: "
        print gtruth

        x = np.arange(dDFs)

        #plt.ylabel(r'$f(alpha)$',fontsize=fsize)
        #plt.xlabel('alpha',fontsize=fsize)
        #plt.plot(x, breadtrain.T, 'k+--', label='bread train',linewidth=2.0)
        #plt.plot(x, breadtest.T, 'g+--', label='bread test',linewidth=2.0)
        #plt.plot(x, nonbreadtrain.T, 'r*--',  label='non bread train',linewidth=2.0)
        #plt.plot(x, nonbreadtest.T, 'b*--',  label='non bread test',linewidth=2.0)
        #plt.legend(loc = 3)
        #plt.show()
    
    y0 = 1.2
    y1 = 3
    plt.figure()
    #plt.subplot(122)
    plt.ylim((y0,y1))
    #plt.xlabel('Real Breads',fontsize=fsize)
    b = plt.boxplot(np.vstack((breadtrain)),sym="")
    mediansReal = map(lambda i: i.get_data()[1][0],b['medians'])
    x = np.arange(len(mediansReal))



    #plt.subplot(121)
    #plt.ylim((y0, y1))
    plt.xlabel('$q$',fontsize=fsize)
    b = plt.boxplot(np.vstack((breadtest)),sym="")
    mediansSynth = map(lambda i: i.get_data()[1][0],b['medians'])
    
    #plt.show()


    err = sum(abs(np.array(mediansReal)-np.array(mediansSynth)))
    err1 = sum(abs(np.array(mediansReal[:dfs])-np.array(mediansSynth[:dfs])))
    err2 = sum(abs(np.array(mediansReal[dfs:])-np.array(mediansSynth[dfs:])))

    qs = range(-dfs,dfs+1)

    xticks(x+1,qs) # translate

    print "ERROR: ", err

    if(err < minn):
        plt.plot(x+1, mediansReal, 'k+--', label='real',linewidth=2.0)   
        plt.plot(x+1, mediansSynth, 'b+--', label='synthetic',linewidth=2.0)
        plt.legend(loc=4)
        plt.ylabel(r'$D_{q}$',fontsize=fsize)
        plt.xlabel(r'$q$',fontsize=fsize)
        plt.title(vals)
        savefig('best_boxplot'+str(err)+'_'+'_.png')

    return err, err1, err2


    if(False):
        scoreRF = (len(gtruth)-sum(abs(gtruth-predictionsRF)))/float(len(gtruth))
        scoreSVM = (len(gtruth)-sum(abs(gtruth-predictionsSVM)))/float(len(gtruth))

        #scores = cross_validation.cross_val_score(cfr, data, labels, cv=4)
        print "Classification performance (Random Forest classifier): " + str( scoreRF*100 ) + "%"
        print "Classification performance (Support Vector Machine classifier): " + str( scoreSVM*100 ) + "%"





