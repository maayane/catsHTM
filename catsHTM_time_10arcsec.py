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

#import catsHTM
#import time
#import numpy as np
#import math

path='/Users/maayanesoumagnac/PostDoc/projects/catsHTM/data'#'/Users/maayanesoumagnac/PostDoc/projects/catsHTM/data'
#catsHTM.cone_search('FIRST',0,0,500,catalogs_dir=path)

'''
start_time = time.time()
for i in range(1000):
        cat,cocell,colu=catsHTM_work_in_progress.cone_search('FIRST',6.283,-7.2431e-5,10,catalogs_dir=path)
end_time=time.time()
print("time elapsed: {:.5f}s".format(end_time - start_time))
print("average time per cone search: {:.5f}s".format((end_time - start_time)/1000.))
'''

coord=np.genfromtxt('random_ra_dec.txt')

radius=10#np.random.choice(np.linspace(10.,10000,10000),size=1000)

'''
start_time = time.time()
for i,j in enumerate(RA):
    catsHTM_work_in_progress.cone_search('FIRST',coord[i,0],coord[i,1],radius,catalogs_dir=path)
end_time=time.time()
print("time elapsed: {:.5f}s".format(end_time - start_time))
print("average time per cone search: {:.5f}s".format((end_time - start_time)/1000.))

pdb.set_trace()
'''

print('FIRST catalog')
start_time = time.time()
for i in range(1000):
#for i,j in enumerate(RA):
    #catsHTM_work_in_progress.cone_search('FIRST',coord[i,0],coord[i,1],radius,catalogs_dir=path)
    #cat,cocell,colu=catsHTM_work_in_progress.cone_search('FIRST',6.283,-7.2431e-5,10,catalogs_dir=path)
    catsHTM_work_in_progress.cone_search('FIRST',coord[i,0],coord[i,1],radius,catalogs_dir=path)
final_time=time.time()
print("Elapsed time is {0} seconds".format(final_time - start_time))

#print("time elapsed for 1000 cone searches: {:.5f}s".format(final_time - start_time))
#print("average time per cone search: {:.5f}s".format((final_time - start_time)/1000.))
#pdb.set_trace()

