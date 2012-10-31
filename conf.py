import sys
from subprocess import *

def conf_mat(test, classes):
        m = [ [0 for i in range(max(classes))] for j in range(max(classes))]
        for i in range(len(classes)):
            m[test[i]-1][classes[i]-1] = m[test[i]-1][classes[i]-1] + 1
        return m

# get cross validation of executing ../exps/easy.py ../exps/sandboxTrain.txt
def getCross():

    easy = 'exps/easy.py'
    scanner = 'exps/sandboxTrain.txt'

    cmd = '{0} {1}'.format(easy, scanner)
    print cmd
    f = Popen(cmd, shell = True, stdout = PIPE).stdout

    line = 1
    while True:
        last_line = line
        line = f.readline()
        if str(line).find("rate") != -1:        
            cross = float(line.split()[-1][5:8])
            break
        if not line: break
    return cross

def test():
    arch = 'sandboxTest.txt.predict'

    easy = 'exps/easy.py'
    scanner = "exps/sandboxTrain.txt"
    camera = 'exps/sandboxTest.txt'

    cmd = '{0} {1} {2}'.format(easy, scanner, camera)
    print cmd
    f = Popen(cmd, shell = True).communicate()
   
    cmd = 'cat "{0}"'.format(arch)
    print cmd
    f = Popen(cmd, shell = True, stdout = PIPE).stdout
    g = f

    c = 0
    line = 1
    while True:
        last_line = line
        line = f.readline()
        c = c+1
        if not line: break
    testL = [0 for i in range(c)]

    cmd = 'cat "{0}"'.format(arch)
    print cmd
    f = Popen(cmd, shell = True, stdout = PIPE).stdout
    c = 0
    line = 1
    while True:
        last_line = line
        line = f.readline()
        testL[c] = int(last_line)
        c = c+1
        if not line: break

    a = [i for i in range(len(testL)-1)]
    a = map(lambda i: i/10+1, a)

    print "Test: ", testL

    b = conf_mat(testL,a)
    for row in b:
        print row

    return testL
