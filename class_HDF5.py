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
        print('Filename:'.format(self.FileName))
        print('Number of datasets in the HDF5 file {0}: {1}'.format(self.FileName,np.shape(f.keys())))
        print('{0} pairs of [data file, indexfile]'.format(np.shape(f.keys())[0]/2.))
        print('****')
        print("The keys are:")
        print(f.keys())
        print('****')
    def load(self,dataset_name,numpy_array=False,Offset=None,Block=None):
        """Description: load an array from a HDF5
        Input:
            -HDF5 file name
            -Dataset name
            -numpy_array: if true, the output is given as a numpy array. Otherwise, it is given as in a HDF5 dataset format
            -[I,J] offset in the array from which to start uploading (?) If not given than get the entire array.
            -[I,J] block size to upload from array.
        Output : - either a numpy array """
        filename = self.FileName
        f = h5py.File(filename, 'r')
        if numpy_array==True:
            data=np.array(f[dataset_name])
        else:
            data=f[dataset_name]
        return data

    def load_colnames(self,Filename,):
        f = h5py.File(Filename, 'r')
        #data = f.get('data/variable1')
        #data = np.array(data)  # For converting to numpy array

