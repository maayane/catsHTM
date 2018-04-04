
"""*******************************************************
A python implementation of catsHTM.m
******************************************************"""
#print __doc__

import math
import numpy as np
import celestial
#import class_HDF5
import scipy.io as sio
import params
import os.path
import sys
import h5py
import pdb
#import time
#class params(object):
#    def __init__(self,path_catalogs,IndexFileTemplate,CatFileTemplate,htmTemplate,NcatinFile,ColCelFile):
#	self.path_catalogs='/Users/maayanesoumagnac/PostDoc/projects/catsHTM/data/'
#	self.IndexFileTemplate='%s_htm.hdf5'
#	self.CatFileTemplate='%s_htm_%06d.hdf5'
#	self.htmTemplate='htm_%06d'
#	self.NcatinFile=100.
#	self.ColCelFile = '%s_htmColCell.mat'
#root_to_data=params.path_catalogs

'''
cata_dic = dict()
cata_dic['TMASS'] = '2MASS'
cata_dic['TMASSxsc'] = 'TMASSxsc'
cata_dic['DECaLS'] = 'DECaLS/DR5'
cata_dic['GAIADR1'] = 'GAIA/DR1'
cata_dic['GALEX'] = 'GALEX/DR6Plus7'
cata_dic['PS1'] = '/PS1'
cata_dic['SDSSDR10'] = 'SDSS/DR10'
cata_dic['SDSSoffset'] = 'SDSSoffset'
cata_dic['UKIDSS'] = 'UKIDSS/DR10'
cata_dic['VISTAviking'] = 'VISTA/Viking/DR2'
cata_dic['VSTatlas'] = 'VST/ATLAS/DR3'
cata_dic['VSTkids'] = 'VST/KiDS/DR3'
cata_dic['FIRST']='FIRST'
'''

def cone_search(CatName,RA,Dec,Radius,catalogs_dir='./data',RadiusUnits='arcsec',IndexFileTemplate=params.IndexFileTemplate,CatFileTemplate=params.CatFileTemplate
                ,htmTemplate=params.htmTemplate,NcatinFile=params.NcatinFile,IndexVarname=None,ColRa = 0,ColDec=1,OnlyCone=True,
                ColCelFile = params.ColCelFile,OutType= 'np_array',verbose=False):
    """Description: Perform a cone search around RA/Dec on a local catalog in HDF5 format sorted into HTM.
    Input  : - Catalog name (e.g., 'GAIADR1'). 
            - J2000.0 R.A. [radians, [H M S], or sexagesimal string].
            - J2000.0 Dec. [radians, [sign D M S], or sexagesimal string].
            - Search radius [arcsec].
            - Optionnal:RadiusUnits - Radius units. Default is 'arcsec'. DO NOT CHANGE THIS DEFAULT
                        IndexFileTemplate - Index Catalog name template. Default is '%s_htm.hdf5'.
                        CatFileTemplate - Catalog name template. Default is '%s_htm_%06d.hdf5'.
                        htmTemplate - HTM dataset template name. Default is 'htm_%06d'.
                        NcatInFile  - Maximum number of Datasets in file.Default is 100.
                        IndexVarName - Default is None.
                        ColRA       - Default is 1.
                        ColDec      - Default is2.
                        OnlyCone    - Return only sources within cone. If false will return also some objects outside cone. Default is true.
                        ColCellFile - Default is '%s_htmColCell.mat'.
    By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018
    Output  : a numpy array where each line is the catalog line for the sources inside the cone """
    #start_time=time.time()
    if verbose==True:
        print('*************')
        print('Catalog: {0}; cone radius: {1} arcsec; cone center: (RA,DEC)=({2},{3})'.format(CatName,Radius,RA,Dec))
        print('*************')

    root_to_data=catalogs_dir+'/'
    if CatName=='TMASS':
        CatDir='2MASS'
    elif CatName=='TMASSxsc':
        CatDir='2MASSxsc'
    elif CatName=='DECaLS':
        CatDir='DECaLS/DR5'
    elif CatName=='GAIADR1':
        CatDir='GAIA/DR1'
    elif CatName=='GALEX':
        CatDir='GALEX/DR6Plus7'
    elif CatName=='PS1':
        CatDir='/PS1'
    elif CatName=='SDSSDR10':
        CatDir='SDSS/DR10'
    elif CatName=='SDSSoffset':
        CatDir='SDSS/DR14offset'
    elif CatName=='UKIDSS':
        CatDir='UKIDSS/DR10'
    elif CatName=='VISTAviking':
        CatDir='VISTA/Viking/DR2'
    elif CatName=='VSTatlas':
        CatDir='VST/ATLAS/DR3'
    elif CatName=='VSTkids':
        CatDir='VST/KiDS/DR3'
    else:
        CatDir=CatName

    Rad = 180. / math.pi

    #if RadiusUnits=='arcsec':
    Radius=Radius/(Rad*3600) #converts arcsec radius into radians radius

    ColCelFile=ColCelFile % CatName
    IndexFilename=IndexFileTemplate % CatName
    if os.path.isfile(root_to_data+CatDir+'/'+ColCelFile)==True:
        test = sio.loadmat(root_to_data+CatDir+'/'+ColCelFile)
        Ncol=np.shape(test['ColCell'])[1]
    else:
        print('ERROR: you need to specify a valid path for the HDF5 catalogs location')
        sys.exit()
    ### computes the list of index of the trixels which intercept the cone
    ID=search_htm_ind(IndexFilename,RA,Dec,Radius,catalogs_dir,VarName=IndexVarname,CatDir=CatDir) #list of IDs of winners leaf
    #print(ID)
    #pdb.set_trace()
    ### computes the catalog with the sources located in those trixels
    ID_matlab=ID+1
    FileID=np.floor(ID_matlab/NcatinFile)*NcatinFile
    Nid=len(ID_matlab) #number of leaf intercepting the circle
    #print(Nid)
    #pdb.set_trace()
    #cat=[]
    if Nid==0:#if none of the catalog's trixel intercept the cone
        if verbose==True:
            print('INFO: the cone does not intercept the catalog')
        cat_onlycone=np.array([])
    else:
        for Iid in range(Nid):
            FileName=CatFileTemplate % (CatName, FileID[Iid])
            DataName=htmTemplate % ID_matlab[Iid]
            filename = root_to_data + CatDir + '/' + FileName
            f = h5py.File(filename, 'r')
            if Iid==0:#optimize
                cat=np.array(f[DataName]).T
            else:
                cat = np.vstack((cat, np.array(f[DataName]).T))
                #cat=np.concatenate((cat,np.array(f[DataName]).T),axis=0) #takes a bit of time

        #if OnlyCone==True:
        D=celestial.sphere_distance_fast(RA,Dec,cat[:,ColRa],cat[:,ColDec])
        cat_onlycone=cat[D<Radius,:]

    ### a colomne with the cell names:
    if cat_onlycone.ndim>1:
        ColCell=np.empty((np.shape(cat_onlycone)[1]),dtype=object)
        ColUnits=np.empty((np.shape(cat_onlycone)[1]),dtype=object)
    else:
        ColCell=np.empty((Ncol),dtype=object)
        ColUnits=np.empty((Ncol),dtype=object)

    for i,j in enumerate(test['ColCell'][0,:]):
        ColCell[i]=str(test['ColCell'][0,i][0])

    for i,j in enumerate(test['ColUnits'][0,:]):
        if len(test['ColUnits'][0,i])>0:
            ColUnits[i]=(test['ColUnits'][0,i][0])
        else:
            ColUnits[i]=' '
    #ColCell[:] = str(test['ColCell'][0, :][0])
    #ColUnits[np.shape(test['ColUnits'][0,:])[0]>0] = (test['ColUnits'][0, np.shape(test['ColUnits'][0,:])[0]>0][0])
    #ColUnits[np.shape(test['ColUnits'][0, :])[0] <= 0] = ' '
    return cat_onlycone,ColCell, ColUnits

def search_htm_ind(Filename,Long,Lat,Radius,path,VarName=None,CatDir=None):
    """Description: wrapper of htm_search_cone, which select from the vector outputed by htm_search_cone only the
    triangles where there are actually sources.
            Input  : - Filename: the name of the index_file, e.g. FIRST_htm.hdf5
            Output :
            By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018

                    """
    if VarName==None:
        cat_name=Filename.split('_')[0]
        VarName=cat_name+'_HTM'

    filenamex = path+'/'+CatDir+'/'+Filename
    f = h5py.File(filenamex, 'r')#open file for reading #OPTIMIZE
    DataHTM_indexfile = np.array(f[VarName])

    #Son_index=np.arange(2, 6)
    #PolesLong_index=np.arange(6,11,2)
    #PolesLat_index=np.arange(7,12,2)
    ID=htm_search_cone(DataHTM_indexfile,Long,Lat,Radius)#,Son_index=Son_index,PolesLong_index=PolesLong_index,PolesLat_index=PolesLat_index) # returns a list of the ID of the winners mesh, i.e. the meshes that intercept the circle

    ID_array=np.array(ID)
    ID_w_sources=ID_array[DataHTM_indexfile[12,ID]>0] #ou l inverse?
    return ID_w_sources

def htm_search_cone(IndexFile_data,Long,Lat,Radius,Ind=None,Son_index=np.arange(2,6),PolesLong_index=np.arange(6,11,2),PolesLat_index=np.arange(7,12,2)):
    #print('I am running htm_search_cone')
    """Description: Search for all HTM leafs intersecting a small circles
        Input  :-Either a table of HTM data or an open HDF5 object in which the HTM data is stored
                -Longitude (radians) to search
                -Latitutde (radians) to search
                -Radius of the small circle
        Output : a vector of indexes of the winner(s):the "adress" in the indexfile of the smallest leaf(s) intercepting the cone
        By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018

                """
    if Ind is None:
        Sons=np.arange(8)#np.linspace(0,7,8).astype(int) #we start with 8 trixels
    else:
        Sons=Ind.astype(int)
        #print(len(Sons))
        #pdb.set_trace()
    ID=[]#np.empty(1)
    Nsons=len(Sons)
    PolesLong=np.zeros((3,Nsons)) #3 lines, Nsons colomns, on veut mettre a chaque colomne les longitudes des poles du mesh
    PolesLat=np.zeros((3, Nsons)) #3 lignes, Nsons colomnes

    #print(Nsons)
    #print(np.shape(Sons))
    #print(np.shape(PolesLong))
    #print(PolesLong_index)
    #print(np.shape(DataHTM))
    #pdb.set_trace()
    for i in range(Nsons):#OPTIMIZE
        PolesLong[:,i]=IndexFile_data[PolesLong_index[:],Sons[i]] # array where each colomn is the 3 poles longitudes of a son mesh HERE: THIS? OR INVERSE?
        PolesLat[:,i]=IndexFile_data[PolesLat_index[:],Sons[i]] # array where each colomn is the 3 poles latitude of a son mesh HERE: THIS? OR INVERSE?

    #pdb.set_trace()
    Flag=celestial.cone_in_polysphere(PolesLong,PolesLat,Long,Lat,Radius) #check if the cone intercept any of the sons meshes


    for i in range(Nsons): #OPTIMIZE
        if Flag[i]==1: #the cone overlap the son with index i
            if np.isnan(IndexFile_data[Son_index[:],Sons[i]]).all()==True:# there are nans in the index_file at the son's index, which means the data is where you are and you cannot go further in the tree
                ID.append(Sons[i])
            else:
                Ind = IndexFile_data[Son_index[:], Sons[i]] - 1.
                #RECURION IS HERE
                ID.extend(htm_search_cone(IndexFile_data,Long,Lat,Radius,Ind=Ind))

    return ID


#pool = mp.Pool(processes=4)
#results = pool.map(cube, range(1,7))
#print(results)