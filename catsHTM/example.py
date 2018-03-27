
"""*******************************************************
Example cone search
******************************************************"""
print __doc__


import catsHTM
import time

start_time = time.time()
cat,colcell,colunits=catsHTM.cone_search('FIRST',0,0,500)

print "time elapsed: {:.5f}s".format(time.time() - start_time)

print '**** columns names ****'
print colcell
print '**** columns units ****'
print colunits
print '**** Result of cone search ****'
print(cat)
