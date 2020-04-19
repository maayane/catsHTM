import catsHTM
path='../../data'
#catsHTM.xmatch_2cats('FIRST','FIRST',catalogs_dir=path)

Cat,ColCel=catsHTM.read_ztf_HDF_matched(815,[10,25],path=path)
print(Cat)
print(ColCel)