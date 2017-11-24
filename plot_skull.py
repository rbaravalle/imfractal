import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
from sklearn import svm
import argparse

# Nc-cross validation type
# FIXME if you want a parameter
Nc = 4

def cross_v(mfss):

    labels = [0 for j in range(len(mfss[0]))]
    data = mfss[0]
    
    for i in range(1,mfss.shape[0]):
        data = np.vstack((data, mfss[i]))
        labels_i = [i for j in range(len(mfss[i]))]
        labels.extend(labels_i)

    cfr = RandomForestClassifier(n_estimators=100)
    cfr2 = svm.LinearSVC(C=1.0)

    labels = np.array(labels)

    scores_rf = cross_validation.cross_val_score(cfr, data, labels, cv=Nc)
    scores_svm = cross_validation.cross_val_score(cfr2, data, labels, cv=Nc)
    
    print Nc, "- Cross Validation scores"
    print "RForest, SVM: " + str( round(np.array(scores_rf).mean(),Nc) ) + " " + str( round(np.array(scores_svm).mean(),Nc) )

def main():
    parser = argparse.ArgumentParser(description='Test '+str(Nc)+'-cross validation of mfss classes')
    parser.add_argument("-npy", dest="npy", type=str, required=True, nargs=1, help="npy you want to test")
    parser.add_argument("--plot", nargs='?', default='F', dest="plot", help="Plot if specified")

    args = parser.parse_args()

    mfss = np.load(args.npy[0])
    cross_v(mfss)

    #FIXME
    if type(args.plot) != str:
        plot(args.npy[0])

def plot(npy):
    mfss = np.load(npy)
    for i in range(mfss.shape[0]):
        for j in range(len(mfss[i])):
            x = np.array(mfss[i][j])
            if i > 10:
                plt.plot(x, '--o', color='C'+str(i%10))
            else:
                plt.plot(x, '-x', color='C'+str(i%10))
    plt.title('All samples')
    plt.show()

if __name__ == '__main__':
    main()
