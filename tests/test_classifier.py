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
