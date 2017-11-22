import numpy as np
import matplotlib.pyplot as plt

def main():
    mfss = np.load('suturas.npy')
    for i in range(mfss.shape[0]):
        for j in range(len(mfss[i])):
            x = mfss[i][j]
            plt.plot(x, color='C'+str(i%10))
        plt.show()

if __name__ == '__main__':
    main()
