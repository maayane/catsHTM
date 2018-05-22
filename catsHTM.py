
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
import class_HDF5
import pdb
#import time

d=dict() #this dictionnary containes the names of the index files loaded in search_htm_ind, and allowes us to avoid loading twice the same index file, which can be time consuming e.g. in a loop

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
    elif CatName=='GAIADR2':
        CatDir='GAIA/DR2'
    elif CatName=='GALEX':
        CatDir='GALEX/DR6Plus7'
    elif CatName=='HSCv2':
        CatDir = 'HST/HSCv2'
    elif CatName=='IPHAS':
        CatDir='IPHAS/DR2'
    elif CatName=='NEDz':
        CatDir='NED/20180502'
    elif CatName=='SDSSDR10':
        CatDir='SDSS/DR10'
    elif CatName=='SDSSoffset':
        CatDir='SDSS/DR14offset'
    elif CatName=='SpecSDSS':
        CatDir='SpecSDSS/DR14'
    elif CatName=='SAGE':
        CatDir = 'Spitzer/SAGE'
    elif CatName=='IRACgc':
        CatDir = 'Spitzer/IRACgc'
    elif CatName=='UKIDSS':
        CatDir='UKIDSS/DR10'
    elif CatName=='VISTAviking':
        CatDir='VISTA/Viking/DR2'
    elif CatName=='VSTatlas':
        CatDir='VST/ATLAS/DR3'
    elif CatName=='VSTkids':
        CatDir='VST/KiDS/DR3'
    elif CatName not in ['AKARI','APASS','Cosmos','FIRST','NVSS','PS1','PTFpc','ROSATfsc','SkyMapper','UCAC4','WISE','XMM']:
        print('ERROR: you need to specify a valid name for the catalog (see README file for list of names)')
        sys.exit()
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
    ID=search_htm_ind(IndexFilename,RA,Dec,Radius,catalogs_dir,VarName=IndexVarname,CatDir=CatDir,verbose=verbose) #list of IDs of winners leaf
    ### computes the catalog with the sources located in those trixels
    ID_matlab=ID+1
    FileID=np.floor(ID_matlab/NcatinFile)*NcatinFile
    Nid=len(ID_matlab) #number of leaf intercepting the circle

    if Nid==0:#if none of the catalog's trixel intercept the cone
        if verbose==True:
            print('INFO: the cone does not intercept the catalog')
        cat_onlycone=np.array([])
    else:
        FileName_0 = CatFileTemplate % (CatName, FileID[0])
        DataName_0 = htmTemplate % ID_matlab[0]
        cat = class_HDF5.HDF5(root_to_data + CatDir + '/' + FileName_0).load(DataName_0, numpy_array=True).T
        for Iid in range(Nid)[1:]:
            FileName=CatFileTemplate % (CatName, FileID[Iid])
            DataName=htmTemplate % ID_matlab[Iid]

            cat=np.vstack((cat, class_HDF5.HDF5(root_to_data + CatDir + '/' + FileName).load(DataName, numpy_array=True).T))
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

    return cat_onlycone,ColCell, ColUnits

def search_htm_ind(Filename,Long,Lat,Radius,path,VarName=None,CatDir=None,verbose=False):
    """Description: wrapper of htm_search_cone, which select from the vector outputed by htm_search_cone only the
    triangles where there are actually sources.
            Input  : - Filename: the name of the index_file, e.g. FIRST_htm.hdf5
            Output :
            By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018

                    """
    if VarName==None:
        cat_name=Filename.split('_')[0]
        VarName=cat_name+'_HTM'

    if VarName not in list(d.values()):
        if verbose==True:
            print('I have not seen the index file corresponding to {0} yet'.format(VarName))
        DataHTM_indexfile = class_HDF5.HDF5(path + '/' + CatDir + '/' + Filename).load(VarName,
                                                                                       numpy_array=True)  # load the indexfile content
        d[str(VarName)+'_name']=VarName
        d[str(VarName)+'_array']= DataHTM_indexfile
    else:
        if verbose==True:
            print('I have already loaded the index file corresponding to {0}'.format(VarName))
        DataHTM_indexfile = d[str(VarName) + '_array']

    '''
    #A working alternative to the dictionnay d, with globals()
    if VarName not in list(globals().values()):
        if verbose==True:
            print('I have not see the index file corresponding to {0} yet'.format(VarName))
        print(path + '/' + CatDir + '/' + Filename)
        print(VarName)
        DataHTM_indexfile = class_HDF5.HDF5(path + '/' + CatDir + '/' + Filename).load(VarName,
                                                                                   numpy_array=True)  # load the indexfile content

        globals()[str(VarName)+'_name'] = VarName
        globals()[str(VarName)+'_array']= DataHTM_indexfile

    else:
        if verbose==True:
            print('I have already loaded the index file corresponding to {0}'.format(VarName))
        DataHTM_indexfile = globals()[str(VarName)+'_array']
    '''
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
        Sons=np.arange(8)
    else:
        Sons=Ind.astype(int)
    ID=[]
    Nsons=len(Sons)
    PolesLong=np.zeros((3,Nsons)) #3 lines, Nsons colomns, on veut mettre a chaque colomne les longitudes des poles du mesh
    PolesLat=np.zeros((3, Nsons)) #3 lignes, Nsons colomnes

    for i in range(Nsons):#OPTIMIZE
        PolesLong[:,i]=IndexFile_data[PolesLong_index[:],Sons[i]] # array where each colomn is the 3 poles longitudes of a son mesh HERE: THIS? OR INVERSE?
        PolesLat[:,i]=IndexFile_data[PolesLat_index[:],Sons[i]] # array where each colomn is the 3 poles latitude of a son mesh HERE: THIS? OR INVERSE?

    Flag=celestial.cone_in_polysphere(PolesLong,PolesLat,Long,Lat,Radius) #check if the cone intercept any of the sons meshes

    for i in range(Nsons): #OPTIMIZABLE?
        if Flag[i]==1: #i.e. if the cone overlap the son with index i
            if np.isnan(IndexFile_data[Son_index[:],Sons[i]]).all()==True:# there are nans in the index_file at the son's index, which means the data is where you are and you cannot go further in the tree
                ID.append(Sons[i])
            else:
                Ind = IndexFile_data[Son_index[:], Sons[i]] - 1.
                #RECURION IS HERE
                ID.extend(htm_search_cone(IndexFile_data,Long,Lat,Radius,Ind=Ind))

    return ID
