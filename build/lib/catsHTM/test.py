import catsHTM
path='../../data'
import numpy as np
#catsHTM.xmatch_2cats('FIRST','FIRST',catalogs_dir=path)

#catsHTM.read_ztf_HDF_matched(815,[10,25],ColCell=None,path=path)

Cat,ColCel=catsHTM.read_ztf_HDF_matched(815,[10,25],ColCell=None,path=path)
print(np.shape(Cat))
print(Cat)
print(ColCel)