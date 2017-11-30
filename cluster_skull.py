import numpy as np
import matplotlib.pyplot as plt
import argparse
from sklearn.cluster import KMeans
import pandas as pd


def cluster(mfss, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(mfss)
    print "Fitted", n_clusters, "clusters"
    print ""
    
    #DF2=pd.DataFrame({'real':Y,'pred':cluster_predicted})

    #gg=ggplot(DF2,aes(x='factor(real)',fill='factor(pred)'))+geom_bar()
    #print(gg)

def flat_array(mfss):
    summ = np.array(mfss[0])
    for i in range(1,len(mfss)):
        summ = np.vstack((summ, np.array(mfss[i])))

    return summ


def main():
    parser = argparse.ArgumentParser(description='Clustering algorithms for the skull database')
    parser.add_argument("-npy", dest="npy", type=str, required=True, nargs=1, help="npy you want to test")
    parser.add_argument("-k", dest="nclusters", type=int, required=True, nargs=1, help="amount of clusters for clustering")
    #parser.add_argument("--plot", nargs='?', default='F', dest="plot", help="Plot if specified")

    args = parser.parse_args()

    mfss = np.load(args.npy[0])
    mfss = flat_array(mfss)

    cluster(mfss, args.nclusters[0])

    #FIXME
    #if type(args.plot) != str:
    #    plot(args.npy[0])

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
