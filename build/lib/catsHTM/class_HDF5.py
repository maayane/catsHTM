"""*******************************************************
A class that allows to access the variables in HDF5 files datasets by index of the elements in the array.
The class also include some basic tools to save and read HDF5 files.

Input  : null
Output : null
Tested :
By : Maayane T. Soumagnac                    February 2018
URL :
Reliable:


called this way from HDF5 import HDF5
******************************************************"""

import h5py
#import pdb
import numpy as np
import pdb
import os

class HDF5(object):
    def __init__(self,FileName):#,VarName):
        self.FileName = FileName
        #self.VarName=VarName
    def info(self):
        """Description: goves info on the HDF5 file
            Input:

            Output : """
        filename = self.FileName
        f = h5py.File(filename, 'r')
        print('Filename:{0}'.format(filename))
        #print(f.keys())
        print('Number of datasets in the HDF5 file {0}: {1}'.format(self.FileName,np.shape(list(f.keys()))))
        print('{0} pairs of [data file, indexfile]'.format(np.shape(list(f.keys()))[0]/2.))
        print('****')
        print("The keys are:")
        print(list(f.keys()))
        print('****')

    def load(self,dataset_name,numpy_array=False,Offset=None,Block=None,Verbose=True):
        """Description: load an array from a HDF5
        Input:
            -HDF5 file name
            -Dataset name
            -numpy_array: if true, the output is given as a numpy array. Otherwise, it is given as in a HDF5 dataset format
            -[I,J] offset in the array from which to start uploading (?) If not given than get the entire array.
            -[I,J] block size to upload from array.
        Output : - a numpy array
        Example: class_HDF5.HDF5('FIRST_htm_010900.hdf5').load('htm_010921',numpy_array=True,Offset=[0,0],Block=[2,4])
        gives the transpose of Cat = HDF5.load('FIRST_htm_010900.hdf5','htm_010921',[1,1],[2,4]) in matlab"""
        filename = self.FileName
        if os.path.exists(filename)==False:
            if Verbose == True:
                print("Cannot open the file <%s>, it does not exist, the trixel must be empty of sources" % (filename))
            data = np.array([])
        else:
            f = h5py.File(filename, 'r')
            if dataset_name not in f.keys():
                if Verbose==True:
                    print("Cannot read Dataset <%s> from hdf5 file <%s>, the trixel must be empty of sources" % (dataset_name, f))
                f.close()
                data=np.array([])
            else:
                if numpy_array==True:
                    datax=np.array(f[dataset_name])
                    if Block is not None:
                        datay=datax.T
                        #print('Offset is',Offset)
                        #print('Block is',Block)
                        if (Offset[0]-int(Offset[0])>0):
                            print('*********** Warning!!! Offset[0] is not an integer! ***********')
                        if (Offset[1]-int(Offset[1])>0):
                            print('*********** Warning!!! Offset[1] is not an integer! ***********')
                        if (Block[0]-int(Block[0])>0):
                            print('*********** Warning!!! Block[0] is not an integer! ***********')
                        if (Block[1]-int(Block[1])>0):
                            print('*********** Warning!!! Block[1] is not an integer! ***********')
                        dataz=datay[int(Offset[0]):int(Offset[0])+int(Block[0]),int(Offset[1]):int(Offset[1])+int(Block[1])]
                        data=dataz.T
                    else:
                        data=datax
                else:
                    if Block is not None:
                        print('stop, Block not none unsupported')
                        #pdb.set_trace()
                        data=f[dataset_name]
                    #MemSpaceId=h5py.h5s.create_simple(TUPLE dims_tpl, TUPLE max_dims_tpl)
        return data

    def load_colnames(self,Filename,):
        f = h5py.File(Filename, 'r')
        #data = f.get('data/variable1')
        #data = np.array(data)  # For converting to numpy array

