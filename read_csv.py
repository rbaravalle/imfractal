import sys
import csv
import numpy as np

f = open(sys.argv[1], 'rb') # opens the csv file

from numpy import recfromcsv
my_data = recfromcsv(sys.argv[1], delimiter=',')


print my_data

np.save("exps/data/bioAsset_meta_adaptive.npy", my_data)

exit()

try:
    reader = csv.reader(f)  # creates the reader object
    for row in reader:   # iterates the rows of the file in orders
        print row    # prints each row
        np.vstack((data, row))

finally:
    f.close()      # closing


print data