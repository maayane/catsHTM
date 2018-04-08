
"""*******************************************************
draw 1000 random positions in the sky, make a 1000 cone search wih a loop on these positions, and calculate the average time for one cone search
******************************************************"""
print(__doc__)

import class_HDF5
import numpy as np
import pdb
#from class_HDF5 import HDF5
import catsHTM_work_in_progress
#import catsHTM
import time
import math


DEC=np.random.choice(np.linspace(-math.pi/2.,math.pi/2.,10000),size=1000)# north only
#print(np.linspace(-math.pi/2.,math.pi/2.,10000))
#pdb.set_trace()
RA=np.random.choice(np.linspace(0,math.pi*2.,10000),size=1000)

print RA
print DEC
coord=zip(RA,DEC)
np.savetxt('./random_ra_dec.txt',coord)

