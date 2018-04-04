"""*******************************************************
draw 1000 random positions in the sky, make a 1000 cone search wih a loop on these positions, and calculate the average time for one cone search
******************************************************"""
print(__doc__)

import catsHTM
import time
import numpy as np
import math
import pdb


path='/Users/maayanesoumagnac/PostDoc/projects/catsHTM/data'
#cProfile.run(catsHTM.cone_search('FIRST',0,0,500,catalogs_dir=path))
'''
start_time = time.time()
cat,cocell,colu=catsHTM.cone_search('FIRST',0,0,500,catalogs_dir=path)
end_time=time.time()
print("time elapsed on check: {:.5f}s".format(end_time - start_time))

print(cat)
print(cocell)
print(colu)
'''
#pdb.set_trace()
DEC=np.random.choice(np.linspace(0,math.pi/2.,10000),size=1000)# north only
#print(np.linspace(-math.pi/2.,math.pi/2.,10000))
#pdb.set_trace()
RA=np.random.choice(np.linspace(0,math.pi*2.,10000),size=1000)
radius=1000#np.random.choice(np.linspace(10.,10000,10000),size=1000)
#print(RA)
#print(DEC)
#print(radius)
#pdb.set_trace()

start_time = time.time()
for i,j in enumerate(RA):
    catsHTM.cone_search('FIRST',RA[i],DEC[i],radius,catalogs_dir=path)
end_time=time.time()
print("time elapsed: {:.5f}s".format(end_time - start_time))
print("average time per cone search: {:.5f}s".format((end_time - start_time)/1000.))

