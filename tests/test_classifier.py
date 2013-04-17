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


import Image
import time
import src.Algorithm.Singularity as Singularity
import src.Algorithm.Sandbox as Sandbox
import csv
import sys
import os
from subprocess import *
import conf # confusion matrix
import numpy

def main():
    cant = 10
    dDFs  = 2*7
    cantClasses = 2
    breadtrain = [['Df' for j in range(dDFs)] for i in range(cant)]
    breadtest = [['Df' for j in range(dDFs)] for i in range(cant)]

    nonbreadtrain = [['Df' for j in range(dDFs)] for i in range(cant)]
    nonbreadtest = [['Df' for j in range(dDFs)] for i in range(cant)]

    pathbtr = 'images/train/bread/'
    dirListbtr=os.listdir(pathbtr)
    pathnbtr = 'images/train/nonbread/'
    dirListnbtr=os.listdir(pathnbtr)
    pathbte = 'images/test/bread/'
    dirListbte=os.listdir(pathbte)
    pathnbte = 'images/test/nonbread/'
    dirListnbte=os.listdir(pathnbte)
    print len(dirListbtr), dirListbtr
    
    ins = Sandbox.Sandbox(dDFs)
    ins.setDef(40,1.15)

    nruns = 10
    votings = [[0 for j in range((cant-1)*cantClasses)] for i in range(nruns)]
    crossVs = [0 for j in range(nruns)]

    for h in range(nruns):
        print "Run number: ", h+1
        for i in range(cant):
            ins.setDef(40,1.15)
            filename = pathbtr+dirListbtr[i]
            print filename
            breadtrain[i] = ins.getFDs(filename)
            filename = pathbte+dirListbte[i]
            print filename
            breadtest[i] = ins.getFDs(filename)

            ins.setDef(50,1.05)
            filename = pathnbtr+dirListnbtr[i]
            print filename
            nonbreadtrain[i] = ins.getFDs(filename)
            filename = pathnbte+dirListnbte[i]
            print filename
            nonbreadtest[i] = ins.getFDs(filename)

        with open('exps/sandboxTrain.csv', 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(breadtrain[0:]+nonbreadtrain[0:])

        with open('exps/sandboxTest.csv', 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(breadtest[0:]+nonbreadtest[0:])

        prog = './scripts/a.out' # convert.c
        csvTrain = 'exps/sandboxTrain.csv'
        csvTest = 'exps/sandboxTest.csv'
        txtTrain = 'exps/sandboxTrain.txt'
        txtTest = 'exps/sandboxTest.txt'
        cmd = '{0} "{1}" > "{2}"'.format(prog, csvTrain, txtTrain)
        Popen(cmd, shell = True, stdout = PIPE).communicate()	
        cmd = '{0} "{1}" > "{2}"'.format(prog, csvTest, txtTest)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        # confusion matrix
        testL = conf.test() # test results and show confusion matrix
        votings[h] = testL
        cr =  conf.getCross()
        print "Actual Cross Validation: "
        crossVs[h] = cr

    finalRes = map(lambda i: most_common(i),zip(*votings))
    print 'FinalRes: '
    print finalRes
    print 'Votings: '
    print votings
    print "New Matrix Conf"
    a = [i for i in range(len(finalRes)-1)]
    a = map(lambda i: i/cant+1, a)
    b = conf.conf_mat(finalRes,a)

    suma = 0
    c = 0
    for row in b:
        suma = suma + row[c]
        print row
        c = c+1

    finalRobustness = float(suma)/((cant-1)*cantClasses)
    finalCrossV = numpy.mean(crossVs);

    print "Different Cross Validations: ", crossVs
    print "Final Robustness: ", finalRobustness*100, "%"
    print "Final Cross Validation: ", finalCrossV


def most_common(lst):
    return max(set(lst), key=lst.count)


main()
