
"""*******************************************************
A python implementation of catsHTM.m
******************************************************"""
#print __doc__

import math
import tqdm
import numpy as np
from . import celestial
import scipy.io as sio
from . import params
import os.path
import h5py
from . import class_HDF5
import time
import pdb
#import mat73
import hdf5storage
#import time

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

# define FileNotFoundError for Python 2.7
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

__all__=['cone_search','search_htm_ind','search_sortedlat','sources_match','htm_search_cone','xmatch_2cats','load_trix_by_ind','simplify_list','load_colcell','mfind_bin','match_cats','simplify2','simplify3','Example_QueryAllFun','read_ztf_HDF_matched'] #redefinition of '*' for import *

def get_CatDir(CatName):
    if CatName == 'TMASS':
        CatDir = '2MASS'
    elif CatName == 'TMASSxsc':
        CatDir = '2MASSxsc'
    elif CatName == 'DECaLS':
        CatDir = 'DECaLS/DR5'
    elif CatName == 'GAIADR1':
        CatDir = 'GAIA/DR1'
    elif CatName == 'GAIADR2':
        CatDir = 'GAIA/DR2'
    elif CatName == 'GAIAEDR3':
        CatDir = 'GAIA/DRE3'
    elif CatName == 'GALEX':
        CatDir = 'GALEX/DR6Plus7'
    elif CatName == 'HSCv2':
        CatDir = 'HST/HSCv2'
    elif CatName == 'IPHAS':
        CatDir = 'IPHAS/DR2'
    elif CatName == 'NEDz':
        CatDir = 'NED/20180502'
    elif CatName == 'SDSSDR10':
        CatDir = 'SDSS/DR10'
    elif CatName == 'SDSSoffset':
        CatDir = 'SDSS/DR14offset'
    elif CatName == 'SpecSDSS':
        CatDir = 'SpecSDSS/DR14'
    elif CatName == 'SAGE':
        CatDir = 'Spitzer/SAGE'
    elif CatName == 'IRACgc':
        CatDir = 'Spitzer/IRACgc'
    elif CatName == 'UKIDSS':
        CatDir = 'UKIDSS/DR10'
    elif CatName == 'VISTAviking':
        CatDir = 'VISTA/Viking/DR2'
    elif CatName == 'VSTatlas':
        CatDir = 'VST/ATLAS/DR3'
    elif CatName == 'VSTkids':
        CatDir = 'VST/KiDS/DR3'
    elif CatName =='ztfSrcLCDR1':
        CatDir='./ZTF/SrcLCDR1'
    elif CatName=='ztfDR1var':
        CatDir = './ZTF/ztfDR1var'
    elif CatName not in ['AKARI', 'APASS', 'Cosmos', 'FIRST', 'NVSS', 'PS1', 'PTFpc', 'PGC', 'ROSATfsc', 'SkyMapper', 'UCAC4',
                         'WISE', 'XMM','NOAO']:
        raise ValueError('you need to specify a valid name for the catalog (see README file for list of available catalogs and names)')
    else:
        CatDir = CatName
    return CatDir

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
    CatDir=get_CatDir(CatName)

    Rad = 180. / math.pi

    #if RadiusUnits=='arcsec':
    Radius=Radius/(Rad*3600) #converts arcsec radius into radians radius

    ColCelFile=ColCelFile % CatName
    IndexFilename=IndexFileTemplate % CatName
    #print(root_to_data+CatDir+'/'+ColCelFile)
    if os.path.isfile(root_to_data+CatDir+'/'+ColCelFile)==True:
        try:
            test = sio.loadmat(root_to_data+CatDir+'/'+ColCelFile)
            #test = mat73.loadmat(root_to_data + CatDir + '/' + ColCelFile)
        except:
            if verbose==True:
                print('sio.loadmat did not work for opening {0}, I am trying hdf5storage.loadmat'.format(root_to_data+CatDir+'/'+ColCelFile))
            #test = sio.loadmat(root_to_data + CatDir + '/' + ColCelFile)
            test= hdf5storage.loadmat(root_to_data+CatDir+'/'+ColCelFile)
        #print(test1)
        #print(test2)
        #pdb.set_trace()
        if np.shape(test['ColCell'])[1]<np.shape(test['ColCell'])[0]:
            #test=test.transpose()
            Ncol=np.shape(test['ColCell'])[0]
        else:
            Ncol=np.shape(test['ColCell'])[1]
    else:
        raise FileNotFoundError("you need to specify a valid path for the HDF5 catalogs location")
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
        #print('shape of cat is',np.shape(cat))
        for Iid in range(Nid)[1:]:
            FileName=CatFileTemplate % (CatName, FileID[Iid])
            DataName=htmTemplate % ID_matlab[Iid]

            cat=np.vstack((cat, class_HDF5.HDF5(root_to_data + CatDir + '/' + FileName).load(DataName, numpy_array=True).T))
        #if OnlyCone==True:
        D=celestial.sphere_distance_fast(RA,Dec,cat[:,ColRa],cat[:,ColDec])#[0]
        cat_onlycone=cat[D<Radius,:]

    ### a colomne with the cell names:
    if cat_onlycone.ndim>1:
        ColCell=np.empty((np.shape(cat_onlycone)[1]),dtype=object)
        ColUnits=np.empty((np.shape(cat_onlycone)[1]),dtype=object)
    else:
        ColCell=np.empty((Ncol),dtype=object)
        ColUnits=np.empty((Ncol),dtype=object)

    #print(np.shape(test['ColCell']))
    #print(np.shape(ColCell))
    #print(np.shape(cat_onlycone))

    if np.shape(test['ColCell'])[1]>np.shape(test['ColCell'])[0]:
        for i,j in enumerate(test['ColCell'][0,:]):
            #print(test['ColCell'][0,i][0])
            ColCell[i]=str(test['ColCell'][0,i][0])
        for i,j in enumerate(test['ColUnits'][0,:]):
            if len(test['ColUnits'][0,i])>0:
                ColUnits[i]=str(test['ColUnits'][0,i][0])
            else:
                ColUnits[i]=' '

    else: #rare cases: Cosmos and TMASSxsc
        for i,j in enumerate(test['ColCell'][:,0]):
            #print(str(test['ColCell'][i][0][0]))
            ColCell[i]=str(test['ColCell'][i][0][0])
        for i,j in enumerate(test['ColUnits'][0,:]):
            if len(test['ColUnits'][0,i])>0:
                ColUnits[i]=str(test['ColUnits'][0,i][0])
            else:
                ColUnits[i]=' '

    return cat_onlycone,ColCell, ColUnits

def search_sortedlat(Cat,Long,Lat,Radius):
    """Description: .
        Input  : - Catalog name (e.g., 'GAIADR1').
                -
                -
                -
                - Optionnal:-Radius [arcsec]: default is 2 arcsec
        By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Oct 2020
        Output  : """
    Ncat=np.shape(Cat)[0]
    Inear=mfind_bin(Cat[:,1],[Lat-Radius,Lat+Radius])
    Ilow=float(Inear[0])
    Ihigh=min(Ncat,float(Inear[1]+1))#add 1 because of the way mfind_bin works
    #print('Ihigh',Ihigh)#ok
    #print('Ilow',Ilow)#ok
    Dist=celestial.sphere_dist_fast(Long,Lat,Cat[int(Ilow-1):int(Ihigh),0],Cat[int(Ilow-1):int(Ihigh),1])[0]
    #print('Dist is',Dist)
    Ind=Ilow-1+np.argwhere(Dist<=Radius)
    return Ind

def sources_match(CatName,Cat,SearchRadius_arcs=2,catalog_dir='./data'):
    """Description: .
        Input  : - Catalog name (e.g., 'GAIADR1').
                -
                -
                -
                - Optionnal:-Radius [arcsec]: default is 2 arcsec
        By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Oct 2020
        Output  : """

    Rad = 180. / math.pi
    SearchRadius_rad=SearchRadius_arcs/(Rad*3600) #converts arcsec radius into radians radius

    Ra=Cat[:,0] # in rad!
    Dec=Cat[:,1] # in rad!
    MedRa=np.nanmedian(Ra)
    MedDec=np.nanmedian(Dec)
    D=celestial.sphere_dist_fast(MedRa,MedDec,Ra,Dec)[0]
    #print('D is',D)
    Radiusi=np.max(D)*(1+10*np.spacing(1))
    Radius=Radiusi*Rad*3600 #converts to arcsec
    #print('Radius',Radius)
    #print('MedRa',MedRa)
    #print('MedDec',MedDec)
    CatHunsorted,ColCelH,ColUnitsH=cone_search(CatName,MedRa,MedDec,Radius,catalogs_dir=catalog_dir)
    #print('CatHunsorted:',CatHunsorted)
    if CatHunsorted.ndim<2:
        CatH = CatHunsorted
    else:
        CatH=CatHunsorted[np.argsort(CatHunsorted[:, 1])]#sort by declination
    #print('ColCelH',ColCelH)
    #pdb.set_trace()
    Nsrc = np.shape(Cat)[0]
    CatM={}
    CatM['Match'] = np.empty((Nsrc,len(ColCelH)))
    Cat=np.empty((Nsrc,len(ColCelH)+2))
    Cat[:,:]=np.nan
    CatM['Match'][:,:]=np.nan
    CatM['Dist'] = np.empty((Nsrc),dtype=object)
    CatM['Dist'][:]=np.nan
    CatM['Nmatch'] = np.zeros((Nsrc, 1))
    if (len(CatH)!=0):
        for Isrc in range(Nsrc):
            #print('I am running sortedlat')
            Ind = search_sortedlat(CatH, Ra[Isrc], Dec[Isrc], SearchRadius_rad).astype(int)
            #print('I am done running sortedlat')
            if (len(Ind)!=0):
                #print('Isrc',Isrc)
                #print('Ind',Ind)
                Dist = celestial.sphere_dist_fast(Ra[Isrc], Dec[Isrc], CatH[Ind, 0],CatH[Ind, 1])[0]
                Distmin=Dist
                Nmatch = len(Ind)
                if (Nmatch > 1):
                #    print('there is more than one match')
                #    print('Dist',Dist)
                #    print('np.min(Dist)',np.min(Dist))
                    Distmin = np.min(Dist)
                #    print('type(Dist)',type(Dist))
                #    print('np.shape(Dist)',np.shape(Dist))
                #    print('np.argmin(Dist)',np.argmin(Dist))
                    MinInd=np.argmin(Dist)
                    Ind = Ind[MinInd]

                CatM['Match'][Isrc,:] = CatH[Ind,:]
                CatM['Dist'][Isrc] = Distmin
                CatM['Nmatch'][Isrc] = Nmatch
                Cat[Isrc,:-2]=CatH[Ind,:]
                Cat[Isrc,-1]=Distmin
                Cat[Isrc,-2]=Nmatch

    CatM['ColCell'] = ColCelH

    return CatM, Cat

@lru_cache(16)
def get_index_data(path, CatDir, Filename, VarName):
    return class_HDF5.HDF5(path + '/' + CatDir + '/' + Filename).load(VarName, numpy_array=True)  # load the indexfile content

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

    DataHTM_indexfile = get_index_data(path, CatDir, Filename, VarName)

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

def get_index_filename(CatName):
    """Description: gets the name of the index file for
            Input  :- Catalog basename (e.g. 'PS1')
            Output :-name of the index filename : <CatBaseName>_htm.hdf5 (carefull! in the paper we wrote this as  <CatBaseName>_htm_ind.hdf5) (e.g. 'PS1_htm.hdf5')
                    -a string <CatBaseName>_HTM (e.g. 'PS1_HTM'), which is the key of the dataset, in the HDF5 file, that contains the 2 columns of the index file
            example: [IndexFileName,IndexVarName]=catsHTM.get_index_filename('PS1')
            By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018
                    """
    IndexFileName=CatName+'_htm.hdf5'
    IndexVarName=CatName+'_HTM'
    return IndexFileName,IndexVarName

def load_HTM_ind(Filename,VarName,catalogs_dir='./data',CatDir=None):
    """Description: load the content of the catalog index file into a dictionnary
            Input  :- index file: an HDF5 file which exists per catalog, containing a 2D array with as many columns as trixels (the index=the column indixe+1: index1 is in columns 0)and each line being:
            [level,Father index,son1 index,son2 index,son3 index,son4 index, Pole1 long, Pole1 lat,Pole2 long, Pole2 lat,Pole3 long, Pole3 lat, either Nan or the data].
                    - The name of the dataset with the actual 2D array stored in the index file. Default is '<CatName>_HTM'
            Output :- A list of N_trixels dictionnaries containing the 2D matrix info
            example:
            By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018"""
    #print('I am looking for the data in',catalogs_dir + '/' + CatDir + '/' +Filename)
    Data=class_HDF5.HDF5(catalogs_dir + '/' + CatDir + '/' +Filename).load(VarName,numpy_array=True)#as many columns as trixels, 13 lines with:
    #  [index,Father index,son1 index,son2 index,son3 index,son4 index, Pole1 long, Pole1 lat,Pole2 long, Pole2 lat,Pole3 long, Pole3 lat, either Nan or the data]
    N_trixels=np.shape(Data)[1]
    #print('there are {0} trixels'.format(N_trixels))
    #load this data into a dictionnaries
    #each trixel is a dictionnary
    HTM_list=[]#will end up being a list of N_trixels dictionnaries
    for i in range(N_trixels):
        trixel = dict()
        trixel['level']=Data[0,i]#line 1 of column 0
        if np.isnan(np.array(Data[1,i])).all() == True:
            trixel['father']=[]
        else:
            trixel['father']=Data[1,i]
        if np.isnan(np.array(Data[2,i])).all() == True:
            trixel['son']=[]
        else:
            trixel['son']=Data[2:6,i]
        trixel['PolesCoo'] = np.zeros((3, 2))
        trixel['PolesCoo'][0, 0] = Data[6,i]
        trixel['PolesCoo'][0, 1] = Data[7,i]
        trixel['PolesCoo'][1, 0] = Data[8,i]
        trixel['PolesCoo'][1, 1] = Data[9,i]
        trixel['PolesCoo'][2, 0] = Data[10,i]
        trixel['PolesCoo'][2, 1] = Data[11,i]
        trixel['Nsrc']=Data[12,i]
        HTM_list.append(trixel)
    return HTM_list,Data

def load_colcell(CatDir,CatName):
    ColCelFile = CatDir+'/'+CatName + '_htmColCell.mat'
    test = sio.loadmat(ColCelFile)
    if np.shape(test['ColCell'])[1] < np.shape(test['ColCell'])[0]:
        # test=test.transpose()
        Ncol = np.shape(test['ColCell'])[0]
    else:
        Ncol = np.shape(test['ColCell'])[1]
    ColCell = np.empty((Ncol), dtype=object)
    ColUnits = np.empty((Ncol), dtype=object)
    if np.shape(test['ColCell'])[1] < np.shape(test['ColCell'])[0]:
        # test=test.transpose()
        Ncol = np.shape(test['ColCell'])[0]
        for i, j in enumerate(test['ColCell'][:, 0]):
            # print(str(test['ColCell'][i][0][0]))
            ColCell[i] = str(test['ColCell'][i][0][0])
        for i, j in enumerate(test['ColUnits'][0, :]):
            if len(test['ColUnits'][0, i]) > 0:
                ColUnits[i] = str(test['ColUnits'][0, i][0])
            else:
                ColUnits[i] = ' '
    else:
        Ncol = np.shape(test['ColCell'])[1]
        for i, j in enumerate(test['ColCell'][0, :]):
            # print(test['ColCell'][0,i][0])
            ColCell[i] = str(test['ColCell'][0, i][0])
        for i, j in enumerate(test['ColUnits'][0, :]):
            if len(test['ColUnits'][0, i]) > 0:
                ColUnits[i] = str(test['ColUnits'][0, i][0])
            else:
                ColUnits[i] = ' '
    return ColCell, ColUnits

def load_trix_by_ind(CatName,index,SearchParValue=None,num=100,catalogs_dir='./data',Ncol=None,Verbose=True):#load_cat in Eran's library
    """Description: given a catalog basename and the index of a trixel, load the content of the corresponding trixel dataset to a numpy array
                Input  :- CatName
                        - trixel index, or a a dataset name
                        - A two element vector of lower and upper value. Only lines in which the sorted parameter is between the low and high value will be retrieved.
                        If empty, retrieve all lines. Default is empty.
                        -number of columns in the catalog.
                Output :-a numpy array with the content of the trixel, Ind ?
                example:
                By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018"""

    if isinstance(index,str)==False:
        names=get_file_dataset_from_trixel_id(CatName,index,NfilesinHDF=num,Verbose=Verbose)
        Filename=names[0]
        Data_set_name=names[1]
    CatDir=get_CatDir(CatName)

    if SearchParValue is None:
        trixel_data=class_HDF5.HDF5(catalogs_dir + '/'+ CatDir + '/' + Filename).load(Data_set_name, numpy_array=True).T
        Ind=1
    else:
        #load the index file
        VarIndStr=Data_set_name+'_Ind' #the name of the index file
        if Verbose==True:
            print('Filename is',Filename)
        DataInd=class_HDF5.HDF5(catalogs_dir+'/'+CatDir+'/'+Filename).load(VarIndStr,numpy_array=True,Verbose=Verbose).T#the content f the index file
        if len(DataInd)>0:
            Ndi=np.shape(DataInd)[0]
            I1=bin_sear(DataInd[:,1],SearchParValue[0])
            I2=bin_sear(DataInd[:,1],SearchParValue[1])
            #print('before the if, I1 is {0} and I2 is {1}'.format(I1,I2))
            Ind=DataInd[I1,0] #the
            Offset=np.append(DataInd[I1,0]-1,0)
            if I1==I2:
                I2=I2+1
            I2=min(I2,Ndi-1)
            Block=[1+DataInd[I2,0]-DataInd[I1,0],Ncol]
            #print('Block is',Block)
            trixel_data=class_HDF5.HDF5(catalogs_dir+'/'+CatDir+'/'+Filename).load(Data_set_name,Offset=Offset,Block=Block,numpy_array=True,Verbose=Verbose).T
            #seach the indexes of the
        else:
            trixel_data=np.array([])
            Ind=None
    return trixel_data,Ind

def bin_sear(X,Val): #Util.find.of eran
    """Description:
                    Input  :- sorted vector (ascending)
                            - Value to search
                    Output :- Index of closest value
                    example:
                    By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018"""
    N=len(X)
    if N==1:
        IndVal=1
    else:
        Ind1=0
        Ind2=N-1
        IndM=math.floor(0.5*N)
        Y1=X[Ind1]
        Y2=X[Ind2]
        Ym=X[IndM]
        Found=0
        while Found==0:
            if Val>Ym:
                Ind1=IndM
                Y1=X[Ind1]
                if Ind2-Ind1>=2:
                    IndM= math.floor(0.5*(Ind2+Ind1))
                else:
                    Found=1
                    if abs(Val-Y1)<abs(Val-Y2):
                        IndVal=Ind1
                    else:
                        IndVal=Ind2
                Ym=X[IndM]
            elif Val<Ym:
                Ind2=IndM
                Y2=X[Ind2]
                if Ind2-Ind1>=2:
                    IndM=math.floor(0.5*(Ind1+Ind2))
                else:
                    Found=1
                    if abs(Val-Y1)<abs(Val-Y2):
                        IndVal=Ind1
                    else:
                        IndVal=Ind2
                Ym=X[IndM]
            else:
                Found=1
                IndVal=IndM
        return IndVal

def mfind_bin(X,Vals):
    """Description: Binary search on a vector running simolutnously on
              multiple values. A feature of this program is that it
              you need to add 1 to the index in order to make sure
              the found value is larger than the searched value.
                    Input  :- Sorted column vector.
                            - Row vector of values to search.
                    Output :- Indices of nearest values.
                    example:
                    By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018"""
    Nvals=len(Vals)
    N=len(X)
    I1=np.ones(Nvals)
    I2=N*np.ones(Nvals)
    Im=np.floor(0.5*(I1+I2)).astype(int)
    #print('Im is',Im)
    PrevIm=np.zeros(np.shape(Im)[0]).astype(int)
    #print('PrevIm is', PrevIm)
    #pdb.set_trace()
    if np.shape(X)[0]<2:
        if X.size==0:
            Im=[]
        else:
            Im=np.ones(Nvals).astype(int)
    else:
        while np.all(Im==PrevIm)==False:
            #print(np.all(Im==PrevIm))
            #print('X[Im-1] is',X[Im-1])
            FlagU=Vals>X[Im-1]
            #print('FlagU is',FlagU)
            FlagD=np.invert(FlagU)
            #print('FlagD is',FlagD)
            I1[FlagU]=Im[FlagU]
            I2[FlagD]=Im[FlagD]
            PrevIm=Im
            Im=np.floor(0.5*(I1+I2)).astype(int)
            #print('Im is',Im)
            #print('PrevIm is',PrevIm)
    return Im

def get_file_dataset_from_trixel_id(CatName,index,NfilesinHDF,Verbose=True):#get_file_var_from_htmid in Eran's library
    """Description: given a catalog basename and the index of a trixel and the number of trixels in an HDF5 file,
                create the trixel dataset name
                Input  :- CatName
                        - index
                        - NfilesinHDF: number of datasets in an HDF5 files (default is 100)
                Output :- Filename: name of the HDF5 file where the trixel_dataset is stored
                        - Datasetname: name of the trixel_dataset
                example:
                By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018"""


    if Verbose==True:
        print('index is',index)
    num_file=math.floor(index/NfilesinHDF)*NfilesinHDF #equivalent to index//Nfiles*Nfiles
    Filename='%s_htm_%06d.hdf5' % (CatName, num_file)
    DatasetName='htm_%06d' % index
    return Filename,DatasetName

def Number_of_trixels(Catname,catalogs_dir='./data',CatDir=None):
    """Description: finds the number of trixels for a given catalod
                Input  :- catalog basename
                Output :- number of trixels for this catalog
                example:
                By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018"""

    IndexFileName = get_index_filename(Catname)[0] # name of the index file associated with Catname
    IndexVarName=get_index_filename(Catname)[1] # name of the data set containing the index filename content
    List_of_dict=load_HTM_ind(IndexFileName,IndexVarName,catalogs_dir=catalogs_dir,CatDir=CatDir)[0]
    Number_of_trixels_in_cat=len(List_of_dict)
    return Number_of_trixels_in_cat

def simplify_list(val):
    if isinstance(val, list) == False:
        return val
    else:
        if len(val) > 1:
            return val
        else:
            return simplify_list(val[0])

def simplify2(x):
    IDc=[]
    for i in x:
        if isinstance(i, (list, tuple, np.ndarray)) == True:
            for j in i:
                IDc.append(j)
        else:
            IDc.append(i)
    return IDc
    #return simplify2(IDc)

def simplify3(x):
    if isinstance(x[0],(list, tuple, np.ndarray)) == False:
        return x
    else:
        y=simplify2(x)
        #print(y)
        return simplify3(y)

def match_cats(Cat,Refcat,Radius=2,RadiusUnits='arcsec'):
        """Description: translation of VO.search.match_cats of Eran. Given two spherical coordinate catalogs. - for each entry
              in the reference catalog (second input argument), search for all nearby sources in the catalog (first input).
                   Input  :- A catalog sorted by declination. Ra and Dec in Rad
                           - A reference catalog. Ra and Dec in rad
                           - 'Radius' - Search radius. This is either a scalar or a vector which length is identical to that of the reference
                 catalog (second input). If a vector than each source in the reference catalog may have a different search radius.
                 Default is 2 (arcsec).
                           - 'RadiusUnits' - Search radius units. See convert.angular for options. Default is 'arcsec'.
                   Output :-Vec: a dictionnary with the following keys
                    Vec['Nfound']= A vector, the size of RefCat, with the number of sources found in the catalog Cat that are within the search radius from the source with same indice in refcat. in the reference catalog.
                    Vec['MinDist']=A vector, the size of RefCat, with the minimum distance (radians) of matched sources in Cat to the source of same indice in RefCat. NaN if not found.
                         - Res: a list of dictionnaries (one item per *matched* refernce source! this list is not the size of cat1, it is the size of the
                         number of objects in cat1 that DO have at least one cross-matched object in cat2):
                            Res['IndRef']=Index of source in reference catalog.
                            Res['IndCat']=List of indices in the catalog that are matched to
%                       the 'IndRef' source of the reference catalog.
                            Res['Dist']= Vecor of angular distances (radians) for each one
%                       of the sources indicated in 'IndCat'.
                            Res['Num']=Number of sources within search radius
                            - IndCatMinDist: vector, the size of Refcat, with the indice of the cat2 nearest sources to the cat1 source of indice Res[Indref]. NaN if no source was found
                   example:
                   By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018"""
        #print('I am running matchcats')
        if RadiusUnits=='rad':
            Radius=Radius
        if RadiusUnits=='arcsec':
            Radius=math.pi*Radius/(180.*3600.)

        Ncat=np.shape(Cat)[0]
        #print('Ncat is',Ncat)#ok
        #print('Refcat is',Refcat)
        Nref=np.shape(Refcat)[0]
        #print('Nref is', Nref)#ok
        Radius=Radius*np.ones(Nref)
        Res=[]
        Iuppx=mfind_bin(Cat[:,1],Refcat[:,1]+Radius) #only if second column is dec!
        Ilowx=mfind_bin(Cat[:,1],Refcat[:,1]-Radius) #only if second column is dec!
        #print('Iupx is',Iuppx)#ok
        #print('Ilowx is',Ilowx)#ok
        Ilow=np.zeros(np.shape(Ilowx)[0])
        for r,s in enumerate(Ilowx):
            Ilow[r]=max(1,Ilowx[r])
        #Ilow=np.max(1,Ilowx)

        Iupp=np.zeros(np.shape(Iuppx)[0])
        for r,s in enumerate(Iuppx):
            Iupp[r]=min(Ncat,Iuppx[r]+1)
        #print('Iup is',Iupp)#ok
        #print('Ilow is',Ilow)#ok
        Ncand=Iupp-Ilow
        Ic=np.array(np.where(Ncand>=1))[0]
        #print('Ic is',Ic)
        #print(np.shape(Ic))
        #print('Ic is',Ic)#index where condition verified, same as matlab one -1
        Nc=np.shape(Ic)[0]
        #print('Nc is',Nc)
        #pdb.set_trace()

        Vec=dict()
        Vec['Nfound']=np.zeros(Nref)
        #vectornan=np.empty(Nref)
        #vectornan[:]=np.nan
        Vec['MinDist']=np.full(Nref, np.nan)#vectornan
        Vec['MinPa']=np.full(Nref, np.nan)#vectornan
        K=0
        IndCatMinDist=np.full(Nref, np.nan)#vectornan

        for Icr in range(Nc):
            #print("Vec['MinDist']5 is", Vec['MinDist'])
            #print('Nc is',Nc)
            Iref=Ic[Icr]
            #print('Iref is',Iref)#ok
            #pdb.set_trace()
            #print('Iref, type:{0},{1}'.format(Iref,type(Iref)))
            #print('Ilow[Iref], type:{0},{1}'.format(Ilow[Iref], type(Ilow[Iref])))
            #print('Iup[Iref], type:{0},{1}'.format(Iupp[Iref], type(Iupp[Iref])))
            #print('Iupp[Iref]-Ilow[Iref]+1, type:{0},{1}'.format(Iupp[Iref]-Ilow[Iref]+1, type(Iupp[Iref]-Ilow[Iref]+1)))
            #print(Ilow[Iref].is_integer())
            if (Ilow[Iref].is_integer()==False) | (Iupp[Iref].is_integer()==False) :
                raise Exception('Ilow[Iref] and Iupp[Iref] should be integers')
            Icat=np.linspace(int(Ilow[Iref]),int(Iupp[Iref]),int(Iupp[Iref]-Ilow[Iref]+1)).astype(int)
            #print('Icat is',Icat)#ok
            #print('Cat[Icat-1,0] is',Cat[Icat-1,0])#ok
            #print('Cat[Icat-1,1] is',Cat[Icat-1,1])#ok
            #print('Refcat[Iref,0]',Refcat[Iref,0])#ok
            #print( 'Refcat[Iref,1]) is',Refcat[Iref,1])#ok
            Dist=celestial.sphere_dist_fast(Cat[Icat-1,0],Cat[Icat-1,1],Refcat[Iref,0],Refcat[Iref,1])[0]
            #print('In matchcats, Dist is',Dist)
            #print('Radius[Iref] is',Radius[Iref])
            IndRelative=np.where(Dist<=Radius[Iref])[0]
            IndCat=Ilow[Icr]-1+IndRelative
            #print('IndRelative is',IndRelative)#ok
            #print('IndCat is',IndCat)#ok
            Vec['Nfound'][Iref]=np.shape(IndCat)[0]#ok
            #print("Vec['Nfound'][Iref] is",Vec['Nfound'][Iref])#ok
            #pdb.set_trace()
            if Vec['Nfound'][Iref]>0:
                Vec['MinDist'][Iref]=np.min(Dist[IndRelative])
                MinInd=np.argmin(Dist[IndRelative])
                Resi=dict()
                K=K+1
                Resi['IndCat']=IndCat
                Resi['IndRef']=Iref
                Resi['Num']=np.shape(IndCat)[0]
                Resi['Dist']=Dist[IndRelative]
                #print("Resi['Dist']",Resi['Dist'])
                Res.append(Resi)
                #print("Vec['MinDist'] 1.5 is", Vec['MinDist'])
                IndCatMinDist[Iref]=IndCat[MinInd]
                ##print('IndCatMinDist[Iref] is {0} and p.min(Dist[IndRelative]) is {1}'.format(IndCatMinDist[Iref],np.min(Dist[IndRelative])) )
               # #print("Vec['MinDist'] 1.8 is", Vec['MinDist'])# ca met IndCatMinDist[Iref] dans Vec['MinDist'][Iref]
               # print("Vec['MinDist'] 2 is", Vec['MinDist'])
            #print("Vec['MinDist'] 3 is", Vec['MinDist'])
            #pdb.set_trace()
        #print("Vec['MinDist'] 4 is", Vec['MinDist'])
        #pdb.set_trace()
        #print('In matchcat, Ind is',Res)
        return Vec,Res,IndCatMinDist #Match,Ind,IndCatMinDist

def Save_cross_matched_catalogs(Cat1,Cat2Matched,output_dir=None):
    """Description: save the outputs of xmatch_2cats, in a directory with
            Input  :- Catalog 1 basename
                    - Catalog 2 basename
                    -Search_radius: default is 2
                    -Search_radius_units: default is arcsec
                    -QueryFun: function to be applied to the catalog
                    -QUeryFunPar: parameters for QueryFun
            Output :
            example:
            By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018
                            """

'''
def Example_QueryAllFun(Cat1,Ind,Cat2,IndCatMinDist,i):
    print('I am running Example_QueryAllFun')
    print('Cat1 is',Cat1)
    print('Ind is',Ind)
    print('Cat2 is',Cat2)
    print('IndCatMinDist is',IndCatMinDist)
    np.save("./Cat1_"+str(i)+'.txt',Cat1)
    return Cat1
'''

def Example_QueryAllFun(Cat1,Ind,Cat2,IndCatMinDist,i,additionnal_args=None):
    print('****** I am running Example_QueryAllFun *******')
    print("Cat1, the content of the catalog_1's trixel is",Cat1)
    print("Cat2, the content of a catalog_2' trixel overlapping with Cat1 is", Cat2)
    print("Ind is a list of dictionnaries, with one dictionnary per Cat1's object having one or more counterparts in Cat2; ")
    print("""Ind[i]["IndRef"]=Index of the Cat1's source having one or more counterpart in Cat2""")
    print("""Ind[i]["IndCat"]=List of indixes of the Cat2's counterparts.""")
    print("""Ind[i]["Dist"]= Vecor of angular distances (radians) between the Cat1's source and its counterparts in Cat2""")
    print('Ind:',Ind)
    print("IndCatMinDist is a vector, with as many elements as lines in Cat1, with 'nan' at lines where there is no counterpart in Cat2, and at line where there is, the catalog_2's index of the closest counterpart")
    print('IndCatMinDist:',IndCatMinDist)
    if additionnal_args is not None:
        np.savetxt(additionnal_args[0]+"/Cat1_"+str(i)+'.txt',Cat1)
    else:
        np.savetxt("./Cat1_" + str(i) + '.txt', Cat1)
    print('***********************************************')
    print('press "c" to continue, "q" to quit')
    pdb.set_trace()
    return Cat1

def xmatch_2cats(Catname1,Catname2,Search_radius=2,QueryAllFun=None,QueryAllFunPar=None,
                 catalogs_dir='./data',Verbose=False,save_results=False,save_in_one_file=True,
                 save_in_separate_files=True,output='./cross-matching_results',time_it=True,Debug=False):
        """Description: cross match two HDF5/HTM catalogs: for each source in the first catalog, the index of the nearest source in the second catalog
        (nearest within some specified distance) is saved.
                Input  :- Catalog 1 basename
                        - Catalog 2 basename
                        -Search_radius: default is 2 (in arcsec)
                        -QueryFun: function to be applied to the catalog
                        -QUeryFunPar: parameters for QueryFun
                        -Verbose: set to True if yu want the code to tell you what it is doing at each step and output intermediate outputs
                        -save_results: if True the the cross-matching pieces of catalog_1 and catalog_2 will be saved. Beware: only on object of catalog 2 (the closest)
                        is saved per object of catalog 1 having a counterpart.
                        -save_in_one_file: if True the results will be saved in one file, of which the first columns are of catalog1 (only those for which
                        cross matching entries in catalog_2 were found), and then come the columns of catalog2
                        -save_in_two_files: if True the results will be saved in two separate files. One has the entries of catalog_1 having at least one counterpart in catalog2
                        and the second is the entries of catalog 2 for the closest counterparts of catalog_2
                        -catalogs_dir: the directory where the HDF5 catalogs are stored
                Output : if save_results=True, the cross-matching pieces of catalog_1 and catalog_2 are stored in the output directory given as the "output" key.
                example: catsHTM.xmatch_2cats('FIRST','NVSS',Verbose=False,save_in_one_file=True,save_results=True,save_in_separate_files=True)
                By : Maayane Soumagnac (original Matlab function by Eran Ofek)            August 2018
                                """

        #Converts search_radius into radians

        Rad = 180. / math.pi
        Search_radius=Search_radius/(Rad*3600) #converts arcsec radius into radians radius

        ###### find the max level between the level of each catalog #####

        CatDir1=get_CatDir(Catname1) #le catalog 1 sous forme de numpy array
        CatDir2=get_CatDir(Catname2) #le catalog 1 sous forme de numpy array
        ##if Verbose==True:


        IndexFileName1 = get_index_filename(Catname1)[0]  # name of the index file associated with Catname
        IndexVarName1 = get_index_filename(Catname1)[1]  # name of the data set containing the index filename content
        HTM1 = load_HTM_ind(IndexFileName1, IndexVarName1, catalogs_dir=catalogs_dir, CatDir=CatDir1)[0]#content of the catalog index file into a dictionnary

        IndexFileName2 = get_index_filename(Catname2)[0]  # name of the index file associated with Catname
        IndexVarName2 = get_index_filename(Catname2)[1]  # name of the data set containing the index filename content
        HTM2 = load_HTM_ind(IndexFileName2, IndexVarName2, catalogs_dir=catalogs_dir,
                                                        CatDir=CatDir2)[0]

        N_trixels_1=Number_of_trixels(Catname1,catalogs_dir=catalogs_dir,CatDir=CatDir1) # number of trixels in catalog 1
        N_trixels_2=Number_of_trixels(Catname2,catalogs_dir=catalogs_dir,CatDir=CatDir2) # number of trixels in catalog 2
        #if Verbose==True:
        print('Catalog_1 is {0} ({1} trixels)'.format(Catname1,N_trixels_1))
        print('Catalog_2 is {0} ({1} trixels)'.format(Catname2, N_trixels_2))

            #print('Catalog_2 is', CatDir2)
            #print('The number of trixels in {0} is {1}'.format(CatDir1,N_trixels_1))
            #print('The number of trixels in {0} is {1}'.format(CatDir2,N_trixels_2))

        L1=celestial.number_of_trixels_to_level(N_trixels_1)[0] #number of levels in catalog 1
        L2=celestial.number_of_trixels_to_level(N_trixels_2)[0] #number of levels in catalog 2

        if Verbose==True:
            print('The level of {0} is {1}'.format(Catname1,L1))
            print('The level of {0} is {1}'.format(Catname2,L2))

        Lmax=max(L1,L2)
        if Verbose==True:
            print('Lmax is',Lmax)#ok compared with Eran; maximum level between cat1 and cat2

        ####### Create the list of trixel's indexes associated with each level #########
        print('************** I am building all the trixels relevant to our search **************')

        built_array = celestial.htm_build(Lmax,Verbose=Verbose)
        HTM=built_array[0]
        Level=built_array[1] #une liste de Lmax dictionnaires, tels que dic['level']=un nombre designant le level (0 pour le level1) et dic['ptr']=un np array des indices des rixels a ce level
        #print(HTM[0].coo())
        #pdb.set_trace()
        #print('HTM[0] is',HTM[0])#ok compared with eran
        #print('HTM[1] is', HTM[1])#ok compared with Eran
        #print('HTM[0][coo] is',HTM[0]['coo'])#ok w Eran
        #print('HTM[8][coo] is', HTM[8]['coo'])# ok
        #print('HTM[9][coo] is', HTM[9]['coo'])#ok
        #print('HTM[10920][coo] is',HTM[10920]['coo'])#ok

        #('HTM[10920] is',HTM[10920])
        #pdb.set_trace()

        Level1=Level[L1-1] # le dictionnaire de Level correspondant au level L1: Level1['Level']=L1-1 et Level1['ptr']= le unumpy array des index des trixesls a ce niveau
        Level2=Level[L2-1]
        if Verbose==True:
            print('Level1:',Level1)
            print('Level2:',Level2)

        Nh1=len(Level1['ptr'])#the number of trixels in the highest level
        print('The number of trixels in the highest level, for {0} is {1}'.format(Catname1,Nh1))#ok
        #pdb.set_trace()
        Nh2=len(Level2['ptr'])
        print('The number of trixels in the highest level, for {0} is {1}'.format(Catname2, Nh2)) #ok
        #pdb.set_trace()

        ColCell2=load_colcell(catalogs_dir+'/'+CatDir2,Catname2)[0]
        ColUnits2=load_colcell(catalogs_dir+'/'+CatDir2,Catname2)[1]
        Ncol2=np.shape(ColCell2)[0]

        ColCell1=load_colcell(catalogs_dir+'/'+CatDir1,Catname1)[0]
        ColUnits1=load_colcell(catalogs_dir+'/'+CatDir1,Catname1)[1]
        Ncol1=np.shape(ColCell1)[0]
        if Verbose==True:
            print('{0} has the following fields {1}'.format(CatDir1,ColCell1))
            print('in units',ColUnits1)
            print('{0} has the following fields {1}'.format(CatDir2, ColCell2))
            print('in units', ColUnits2)
        #At this stage, we have 2 Level dictionnaries, one per each catalog, such that LevelX['level'] is the number of the highest level (0 for level 1)
        # and LevelX['ptr'] is the list of trixels indexes at the highest level
        #Next, we go through all the highest level trixels of Catalog 1, and for each trixel, if it contains sources, we check if there are some overlapping trixels in catalog 2

        if save_results == True:
            if os.path.exists(output):
                print('the output directory, ' + output + ' exists already')
            else:
                os.mkdir(output)
        header1 = ",".join([Catname1+':'+ColCell1[i] + ' (' + ColUnits1[i] + ')' for i in range(np.shape(ColCell1)[0])])
        header2 = ",".join([Catname2+':'+ColCell2[i] + ' (' + ColUnits2[i] + ')' for i in range(np.shape(ColCell2)[0])])
        cross_matching_result = np.empty((1, np.shape(ColCell1)[0] + np.shape(ColCell2)[0]))
        #print(np.shape(cross_matching_result))
        #print('header1 is',header1)
        #print('header2 is', header2)
        #print(header1+','+header2)
        if save_results==True:
            if save_in_one_file==True:
                if os.path.exists(output + '/cross-matching_result_full.txt'):
                    print('the txt file exists already, I am removing it')
                    os.remove(output + '/cross-matching_result_full.txt')
            if save_in_separate_files==True:
                if os.path.exists(output + '/cross-matching_result_{0}.txt'.format(Catname1)):
                    print('the txt file for {0} exists already, I am removing it'.format(Catname1))
                    os.remove(output + '/cross-matching_result_{0}.txt'.format(Catname1))
                if os.path.exists(output + '/cross-matching_result_{0}.txt'.format(Catname2)):
                    print('the txt file for {0} exists already, I am removing it'.format(Catname2))
                    os.remove(output + '/cross-matching_result_{0}.txt'.format(Catname2))

        #print("Level1['ptr'] is", Level1['ptr'])
        #np.savetxt('indexes.txt',Level1['ptr'])
        print('************** I am looking for overlapping trixels **************')
        start = time.time()
        if Debug == True:
            print('I will stop at the following indexes, if the trixels exists, to debug, ok? press c to continue',
                  [Nh1//1000,Nh1//200,Nh1 // 100, Nh1 //10, Nh1 //4, Nh1 //3, Nh1 // 2, Nh1 // 1.5])
            pdb.set_trace()
        for i in tqdm.tqdm(range(Nh1)): #for each trixels in the highest level of Cat1
            #print("Level1['ptr'][Nh1-1] is",Level1['ptr'][Nh1-1])
            #print("Level1['ptr'][i] is",Level1['ptr'][i])

            index_cat1=Level1['ptr'][i]# takes the index of this trixel and check if this trixel contains sources:
            #print('I am looking for Catalog_2 ({0}) trixels overlapping with the trixel #{2} of Catalog_1 ({1})'.format(Catname2,Catname1,index_cat1))
            if HTM1[index_cat1-1]['Nsrc']>0:#if the trixel contains sources:
                #if index_cat1==27305:
                if Verbose==True:
                    print('I am looking for Catalog_2 ({0}) trixels overlapping with the non-empty trixel #{2} ({3}/{4}) of Catalog_1 ({1})'.format(Catname2, Catname1, index_cat1,i,Nh1))
                    print('there are {0} sources in this trixel'.format(HTM1[index_cat1-1]['Nsrc']))
                #print('not empty')
                #print('I am looking for Catalog_2 ({0}) trixels overlapping with the trixel #{2} of Catalog_1 ({1})'.format(Catname2,Catname1,index_cat1))
                #print('the file with index {0} has {1} sources'.format(index_cat1,HTM1[index_cat1]['Nsrc']))
                #start = time.time()
                Cat1=load_trix_by_ind(Catname1,index_cat1,num=100,catalogs_dir=catalogs_dir,Verbose=Verbose)[0]#load the content of that trixel (in the form of a numpy array)
                #ongoing1=time.time()
                #print(Cat1)#ok
                #Cat 1 is a numpy array with the content of a trixel that contains sources, at the highest level of Catalog1
                #PolesCoo ok
                #print("HTM[index_cat1-1]['coo'] is",HTM[index_cat1-1]['coo'])#ok


                MeanRa=np.mean(HTM[index_cat1-1]['coo'][:,0]) # le meam Ra de ce trixel
                MeanDec=np.mean(HTM[index_cat1-1]['coo'][:,1]) # le mean Dec de ce trixel

                MinDec=np.min(HTM[index_cat1-1]['coo'][:,1])-Search_radius
                MaxDec = np.max(HTM[index_cat1 - 1]['coo'][:, 1]) + Search_radius
                #print('MeanRa is', MeanRa) #ok
                #print('MeanDec is',MeanDec)#ok

                D=celestial.sphere_dist_fast(MeanRa,MeanDec,HTM[index_cat1-1]['coo'][:,0],HTM[index_cat1-1]['coo'][:,1])[0]
                #print('D is',D)
                CircRadius=np.max(D)+Search_radius
                #print('CircRadius is',CircRadius)
                ID2=celestial.htm_search_cone(HTM2,MeanRa,MeanDec,CircRadius,Ind=[])
                #if Verbose==True:
                ID2w=simplify3(ID2)
                ongoing2 = time.time()
                if Verbose==True:
                    print('there are {0} trixel overlapping with it'.format(len(ID2w)))#ok
                    #pdb.set_trace()
                    print('the list of trixels indexes of Catalog_2({0}) overlapping with the trixel #{2} of Catalog_1({1}) is {3}'.format(
                    Catname2, Catname1, index_cat1,ID2w))
                #print('the list without brakets is',ID2w)# a list of indexes of cat2 trixels, which overlap with the cat1 trixel

                #load all the data corresponding to ID2w
                Nid2=len(ID2w) #the number of trixels of cat 2 overlapping with the given trixel of cat1 which we are examining.

                for s in range(Nid2):#for all trixels of catalog 2 overlapping with the given trixel of catalog1
                    if s==0:
                        [Cat2,Ind2]=load_trix_by_ind(Catname2,ID2w[s],[MinDec,MaxDec],catalogs_dir=catalogs_dir,Ncol=Ncol2,Verbose=Verbose)
                        N2=np.shape(Cat2)[0]
                        #Cat2ID=np.array(list(zip(ID2w[i]*np.ones(N2),Ind2+np.array(range(N2)))))#MAYBE Ind2-1?
                        #print('len(Cat2) after i=0 is',len(Cat2))
                        #pdb.set_trace()
                    else:
                        if Verbose==True:
                            print('**********')
                            print("(catalog_2) {0}'s trixel (overlapping with (catalog_1) {1}'s trixel) of index {2}:".format(Catname2,Catname1,index_cat1))
                        [Cat2tmp,Ind2]=load_trix_by_ind(Catname2,ID2w[s],[MinDec,MaxDec],catalogs_dir=catalogs_dir,Ncol=Ncol2,Verbose=Verbose)
                        #print('i={0},shape(Cat2) and shape(Cat2tmp) are {1} and {2}'.format(i,np.shape(Cat2),np.shape(Cat2tmp)))
                        #pdb.set_trace()
                        #ongoing3 = time.time()
                        if len(Cat2)>0:
                            #print('at this (1) stage len(Cat2) is',len(Cat2))
                            #print('Cat2tmp (1) is',Cat2tmp)
                            if len(Cat2tmp)>0:
                                Cat2=np.vstack((Cat2,Cat2tmp))
                                N2 = np.shape(Cat2)[0]
                            #else:

                                #Cat2ID=np.vstack((Cat2ID,np.array(list(zip(ID2w[i]*np.ones(N2),Ind2+np.array(range(N2)))))))#MAYBE Ind2-1?
                            #else: Cat2 reste tel quel
                        else:#si Cat2 etait vide
                            #print('at this (2) stage len(Cat2) is',len(Cat2))
                            #print('Cat2 was empty?')
                            if len(Cat2tmp)>0:#si Cat2tmp n'est pas vide, Cat2 devient lui
                                #print('Cat2tnp.argwhere(np.isnan(x))mp (2) is',Cat2tmp)
                                Cat2=np.copy(Cat2tmp)
                                N2 = np.shape(Cat2)[0]
                    #print('Cat 2 is', Cat2)
                    #pdb.set_trace()
                                #Cat2ID=np.vstack((Cat2ID,np.array(list(zip(ID2w[i]*np.ones(N2),Ind2+np.array(range(N2)))))))#MAYBE Ind2-1?
                            #else: Cat2 reste vide
                #print('Cat2 is',Cat2)
                #print('len(Cat2) is',len(Cat2))
                #print('np.shape(Cat1) is',np.shape(Cat1))
                #print('np.shape(Cat2) is', np.shape(Cat2))
                        #ongoing4 = time.time()
                # C'est quoi Cat2? Cat2 is a catalog with the content of *all the Catalogue 2 trixels overlapping with the given trixel of cat1
                # C'est quoi Cat2ID?*

                #print('Cat2 before sorting is',Cat2)
                #print('Cat2[:, 1] is',Cat2[:,1] )
                #pdb.set_trace()
                #print('len(Cat2) after the loop is',len(Cat2))
                #pdb.set_trace()
                if len(Cat2)>0:
                    cat2=Cat2[Cat2[:, 1].argsort(),] #cat2 est Cat2 -l'ensemble des trixels qui overlappent cat1 -trié par Dec croissant. On a besoin de ca pour applyer match_cats.
                    #np.savetxt('cat2.txt', cat2)
                    #SI=Cat2[:, 1].argsort() #SI est les indexes de Dec croissants de Cat2
                    #print('SI is',SI)# ok, verifie avec matlab
                    #probleme: cat 2 c est toutes les sources des overlapping trixels. Nous on veut que les sources reelelemt overlapping. donc on run match_cat
                    #ongoing5 = time.time()
                    [Match,Ind,IndCatMinDist]=match_cats(cat2,Cat1,Radius=Search_radius,RadiusUnits='rad')
                    #print('Ind',Ind)
                    #pdb.set_trace()

                    if QueryAllFun is not None:
                        if Verbose==True:
                            print('I am applying a function on the cross-match result!')
                        #if i==0:
                        #    Data=np.array([])
                        #else:
                        Data = QueryAllFun(Cat1, Ind, cat2, IndCatMinDist, i, additionnal_args=QueryAllFunPar)
                        #CAREFULLLLLLL it s cat2 not Cat2!!!!
                        #print(Data)
                        #pdb.set_trace()

                    #ongoing6 = time.time()

                    #Match:a dictionnary with the following keys
                    #Match['Nfound']= a vector, the length of cat1, with the number of sources found in the cat2 that are within the search radius from the source in the reference catalog Cat1.
                    #Match['MinDist']=a vector, the size of cat1, wiht the Minimum distance (radians) of sources in cat2 to the source in cat1. NaN if not found
                    #Ind: a list of dictionnaries (as many as sources in Cat1 THAT HAVE CROSS-MTACHED SOURCES in cat2)
                    # Ind[i]['IndRef']=Indice of source in cat1
                    # Ind[i]['IndCat']=List of indices in cat2 that are matched to the 'IndRef' source of Cat1.
                    # Ind[i]['Dist']= Vecor of angular distances (radians) for each one of the sources indicated in 'IndCat'.
                    # Ind[i]['Num']=Number of sources within search radius
                    # IndCatMinDist:  a vector of indices of cat2 objects which are the closest to the source in cat1. NaN if not found ??

                    #print("Match['Nfound'] is",Match['Nfound']) #ok, verifie avec matlab
                    #print("Match['MinDist'] is", Match['MinDist'])  #ok, verifie avec matlab
                    #print("Match['MinPA'] is", Match['MinPa']) #ok, verifie avec matlab
                    #print("Ind is",Ind)
                    #print("the Ind['Num'] are:",[Ind[i]['Num'] for i in range(len(Ind))]) # ok
                    #print("the Ind['IndCat'] are:", [Ind[i]['IndCat'] for i in range(len(Ind))])  # ok, moi=matlab-1, normal
                    #print("the Ind['IndRef'] are:", [Ind[i]['IndRef'] for i in range(len(Ind))])  # ok, moi=matlab-1, normal
                    #print("the Ind['Dist'] are:", [Ind[i]['Dist'] for i in range(len(Ind))]) # ok
                    #pdb.set_trace()
                    #print('IndCatMinDist is',IndCatMinDist)#ok, moi=matlab-1, normal
                    #print('the shape of IndCatMinDist is',np.shape(IndCatMinDist)[0]) #ok
                    """ if (~isempty(InPar.QueryAllFun))
                                % execute InPar.QueryAllFun
                                %  QueryAllFun(Cat1,Ind,Cat2,varargin)
                                if (Ih1==Istart)
                                    Data = [];
                                end

                                Data = InPar.QueryAllFun(Cat1,Ind,Cat2,IndCatMinDist,InPar.QueryAllFunPar{:},'Data',Data,'Ih1',Ih1,'Nh1',Nh1,'SearchRadius',InPar.SearchRadius);
                            end"""
                    IsN=np.isnan(IndCatMinDist)# un tableau de booleans qui est True la ou il y a zero sources cross-matched, et False la ou il y en a
                    #print('IsN is',IsN)
                    #print('IsN is',IsN) ok, mais moi c est des True et False et matlab c est des 0 et 1
                    #print('the shape of IsN is',np.shape(IsN)) ok
                    IndCatMinDist[IsN]=True #
                    #if V
                    #print('IndCatMinDist is now',IndCatMinDist) # un tableau de la taille de cat1 avec : la ou il y a pas de cross-matched dans cat2: 1, et la ou il y en a: l'indice de l'objet de cat2 le plus proche
                    """
                    ceci: pas clair a quoi ca sert dans le code de matlab. Je laisse tomber.
                    print("Cat2ID is",Cat2ID) #ok mais pas sur qu'il dooivent etre identiques
                    print("SI[IndCatMinDist.astype(int)] is",SI[IndCatMinDist.astype(int)]) #pas ok
                    pdb.set_trace()
                    DataInd=Cat2ID[SI[IndCatMinDist.astype(int)],:]
                    DataInd[IsN,:]=np.nan
                    print('DataInd is', DataInd)  # pareil que matlab mais pas sur que c est bien
                    """
                    #print("IndCatMinDist.astype(int) is",IndCatMinDist.astype(int))
                    #print("np.shape(cat2)",np.shape(cat2))
                    #print("np.shape(IndCatMinDist)",np.shape(IndCatMinDist))
                    #print("np.shape(IndCatMinDist.astype(int))",np.shape(IndCatMinDist.astype(int)))

                    #print("cat2[IndCatMinDist.astype(int)-1,:] is",cat2[IndCatMinDist.astype(int)-1,:])
                    #print('IndCatMinDist.astype(int)-1 is',IndCatMinDist.astype(int)-1)
                    #print("cat2[IndCatMinDist.astype(int),:] is", cat2[IndCatMinDist.astype(int), :])
                    #print('IndCatMinDist.astype(int) is', IndCatMinDist.astype(int))
                    indexes_analog_to_matlab=np.zeros(np.shape(IndCatMinDist))
                    indexes_analog_to_matlab[IndCatMinDist!=1]=IndCatMinDist[IndCatMinDist!=1]
                    #THIS CHECK IS CRUCIAL! DON'T ARAISE
                    #if Verbose==True:
                    #    print('i (or matlab Ih1-1)={0},indexes_analog_to_matlab must be matlab Indcatmindist-1 everywhere, check if this is the case: {1}'.format(i,indexes_analog_to_matlab))#ok
                    #
                    Cat2matched = cat2[indexes_analog_to_matlab.astype(int), :]#ok
                    #Cat2matched=cat2[IndCatMinDist.astype(int),:]
                    #Cat2matched=cat2[IndCatMinDist.astype(int),:]
                    #print('cat2 is,',cat2)
                    # Cat2matched est un tableau, de la longueur de cat1 avec:
                    #  -la ligne 0 de cat2 si la ligne correspond a un indice de cat1 qui a pas de cross-match
                    #  -s'il y a un cross-matched dans cat2: la ligne de cat2
                    #print("np.shape(Cat2matched)",np.shape(Cat2matched))
                    #print("np.shape(IsN)",np.shape(IsN))
                    Cat2matched[IsN,:]=np.nan #
                    #print('Cat2matched is', Cat2matched)
                    if Debug==True:
                        if i in [Nh1//1000,Nh1//200,Nh1 // 100, Nh1 //10, Nh1 //4, Nh1 //3, Nh1 // 2, Nh1 // 1.5]:
                            print('******** i={0} ********'.format(i))
                            print('I am saving Cat2matched')
                            np.savetxt(output+'Cat2matched_{0}_4debug.txt'.format(i),Cat2matched) #pas ok
                            pdb.set_trace()
                    #print('Cat2matched at the index of IndCatMinDist is',Cat2matched[IndCatMinDist!=1])
                    #print('IndCatMinDist', IndCatMinDist)
                    #pdb.set_trace()

                    #print('Cat2matched is', Cat2matched)
                    # un tableau, avec le meme nombre de lignes que cat1 et le nombre de colomnes de cat2 avec:
                    #  -NaN si cette ligne de cat1 a pas de cross-match
                    #  -s'il y a un cross-matched dans cat2: la ligne de cat2 correspondant a l objet le plus proche
                    # print("np.shape(Cat2matched)",np.shape(Cat2matched))

                    #print('cat1',cat1)
                    #print('Cat2matched is',Cat2matched)#ok avec matlab
                    #pdb.set_trace()
                    #print('np.shape(Cat2matched is)',np.shape(Cat2matched)) #ok avec matlab

                    #from here it is added by me
                    #create a numpy array with: columns of cat1, columns of Cat2matched
                    #print('let us just make sure that Cat1 and Cat2matched have same number of lines.')ok
                    #print('np.shape(Cat1) is',np.shape(Cat1))
                    #print('np.shape(Cat2matched) is', np.shape(Cat2matched))


                    #if save_results==True:
                    #    if os.path.exists(output):
                    #        print('the output directory, ' + output+ ' exists already')
                    #    else:
                    #        os.mkdir(output)
                    #    if os.path.exists(output+'/trixel_'+str(index_cat1)+'_'+Catname1):
                    #        print('the output directory, ' + output+'/trixel_'+str(index_cat1)+'_'+Catname1 + ' exists already')
                    #    else:
                    #        os.mkdir(output+'/trixel_'+str(index_cat1)+'_'+Catname1)
                    if save_results==True:
                        cross_matching_result_w_nans=np.hstack((Cat1,Cat2matched))
                        #cross_matching_result_intermediate = np.empty((1,np.shape(Cat1)[1]+np.shape(cat2)[1]))
                        cross_matching_result_intermediate = np.zeros((1, np.shape(Cat1)[1] + np.shape(cat2)[1]))
                        for i,j in enumerate(cross_matching_result_w_nans[:,0]): #for all lines,remove the lines where no cross-matched object
                            if np.all(np.isnan(cross_matching_result_w_nans[i, np.shape(Cat1)[1]:])) == False:
                                if Verbose==True:
                                    print('At line {0} of Cat1, there is a cross-matched object in cat2'.format(i))
                                    #print('Cat2matched[i,:] is',Cat2matched[i,:])
                                    #pdb.set_trace()
                                if np.shape(cross_matching_result_intermediate)[0]<2:
                                    #print('np.shape(cross_matching_result_intermediate)[0] is',np.shape(cross_matching_result_intermediate)[0])
                                    cross_matching_result_intermediate=cross_matching_result_w_nans[i,:]
                                    cross_matching_result_intermediate_cat1 = cross_matching_result_w_nans[i, :np.shape(Cat1)[1]]
                                    cross_matching_result_intermediate_cat2 = cross_matching_result_w_nans[i,np.shape(Cat1)[1]:np.shape(Cat1)[1]+np.shape(Cat2matched)[1]]
                                else:
                                    #print('else')
                                    cross_matching_result_intermediate=np.vstack((cross_matching_result_intermediate,cross_matching_result_w_nans[i,:]))
                                    cross_matching_result_intermediate_cat1 = cross_matching_result_intermediate[:, :np.shape(Cat1)[1]]
                                    cross_matching_result_intermediate_cat2 = cross_matching_result_intermediate[:,np.shape(Cat1)[1]:np.shape(Cat1)[1]+np.shape(Cat2matched)[1]]

                            #else:
                                #print('there are no counterparts in cat2')

                        all_zeros = not np.any(cross_matching_result_intermediate)
                        if all_zeros==True:
                            print('There are no counterpart at all in cat 2 for this tri1xel')
                            #pdb.set_trace()
                        else:
                        #print('the shape of cross_matching_result_intermediate_cat1 is',np.shape(cross_matching_result_intermediate_cat1))
                        #print('the shape of cross_matching_result_intermediate_cat2 is',
                        #      np.shape(cross_matching_result_intermediate_cat2))
                        #print('the shape of cross_matching_result_intermediate is',
                        #      np.shape(cross_matching_result_intermediate))
                        #print('ndim of cross_matching_result_intermediate_cat1) is 1?',
                        #      cross_matching_result_intermediate_cat1.ndim)
                        #print('the len of cross_matching_result_intermediate_cat1 is',np.shape(cross_matching_result_intermediate_cat1)[0])
                        #print('the len of cross_matching_result_intermediate_cat2 is',np.shape(cross_matching_result_intermediate_cat2)[0])
                        #print('the len of cross_matching_result_intermediate is',np.shape(cross_matching_result_intermediate)[0])
                        #if np.shape(cross_matching_result_intermediate_cat1)[0]!=np.shape(cross_matching_result_intermediate_cat2)[0]:
                        #    print('ndim of cross_matching_result_intermediate_cat1) is 1?',cross_matching_result_intermediate_cat1.ndim)
                        #    print('the shapes are not the same, probleme!')
                        #    print(cross_matching_result_intermediate_cat1)
                        #    print(cross_matching_result_intermediate_cat2)
                        #    print('np.shape(cross_matching_result_intermediate)[0] is',np.shape(cross_matching_result_intermediate)[0])
                        #    print('cross_matching_result_intermediate is',cross_matching_result_intermediate)
                        #    print('np.shape(cross_matching_result_w_nans[i,:]))',np.shape(cross_matching_result_w_nans[i,:]))
                        #    pdb.set_trace()
                            #if Verbose is True:
                            print('The entries from catalog_1 ({0}) :{1}, cross-matched in catalog_2 ({2}) are {3}'.format(Catname1,cross_matching_result_intermediate_cat1,Catname2,cross_matching_result_intermediate_cat2))
                            pdb.set_trace()
                            #print('cross_matching_result is',cross_matching_result)
                            #print('Is the cross_matching_result the size of Ind?')#yes
                            #print(np.shape(cross_matching_result))
                            #print(len(Ind))
                            #print('Is the number of columns of cross_matching_result the sum of the number of columns of cat1 and cat2?')#yes
                            #print(np.shape(cross_matching_result))
                            #print(np.shape(Cat1))
                            #print(np.shape(cat2))
                            """
                                if (~isempty(InPar.QueryFun))
                                    % execute InPar.QueryFun
                                    % QueryFun can select specific sources (by some
                                    % attributes) from the matched Cat1 and Cat2

                                    FlagSelected       = InPar.QueryFun(Cat1,Cat2matched,InPar.QueryFunPar{:});
                                    % what to do with FlagSelected?
                                    Cat1        = Cat1(FlagSelected,:);
                                    Cat2matched = Cat2matched(FlagSelected,:);

                                end

                                if (~isempty(InPar.SaveFun))
                                    % execute InPar.SaveFun
                                    % Fun(Cat1,Cat2matched)
                                    InPar.SaveFun(Cat1,Cat2matched,InPar.SaveFunPar{:});
                                end
                            """

                            #print('np.shape(cross_matching_result_intermediate) is ',np.shape(cross_matching_result_intermediate))
                            #print(
                            #'np.shape(cross_matching_result_intermediate_cat1) is ', np.shape(cross_matching_result_intermediate_cat1))
                            #print(
                            #'np.shape(cross_matching_result_intermediate_cat2) is ', np.shape(cross_matching_result_intermediate_cat2))
                            #if np.shape(cross_matching_result_intermediate_cat1)[0]!=np.shape(cross_matching_result_intermediate)[0]:
                            #    print('pb!')
                            #    print('cross_matching_result_intermediate is',cross_matching_result_intermediate)
                            #    print('cross_matching_result_intermediate_cat1 is',cross_matching_result_intermediate_cat1)
                            #    print('cross_matching_result_intermediate_cat2 is', cross_matching_result_intermediate_cat2)
                            #    pdb.set_trace()

                            if save_in_one_file==True:
                                if os.path.exists(output +'/cross-matching_result_full.txt')==False:
                                    with open(output +'/cross-matching_result_full.txt', 'ab') as f:
                                        if cross_matching_result_intermediate.ndim>1:
                                            np.savetxt(f, cross_matching_result_intermediate, delimiter=",",header=header1+','+header2)
                                        else:
                                            np.savetxt(f, cross_matching_result_intermediate[None], delimiter=',',header=header1+','+header2)
                                else:
                                    with open(output +'/cross-matching_result_full.txt', 'ab') as f:
                                        if cross_matching_result_intermediate.ndim > 1:
                                            np.savetxt(f, cross_matching_result_intermediate, delimiter=",")
                                        else:
                                            np.savetxt(f, cross_matching_result_intermediate[None], delimiter=",")
                            if save_in_separate_files==True:
                                if os.path.exists(output +'/cross-matching_result_{0}.txt'.format(Catname1))==False:
                                    with open(output +'/cross-matching_result_{0}.txt'.format(Catname1), 'ab') as f:
                                        if cross_matching_result_intermediate_cat1.ndim>1:
                                            np.savetxt(f, cross_matching_result_intermediate_cat1, delimiter=",",header=header1)
                                        else:
                                            np.savetxt(f, cross_matching_result_intermediate_cat1[None], delimiter=",",
                                                       header=header1)
                                else:
                                    with open(output + '/cross-matching_result_{0}.txt'.format(Catname1), 'ab') as f:
                                        if cross_matching_result_intermediate_cat1.ndim>1:
                                            np.savetxt(f, cross_matching_result_intermediate_cat1,
                                                       delimiter=",")
                                        else:
                                            np.savetxt(f, cross_matching_result_intermediate_cat1[None],
                                                       delimiter=",")
                                if os.path.exists(output + '/cross-matching_result_{0}.txt'.format(Catname2)) == False:
                                    with open(output + '/cross-matching_result_{0}.txt'.format(Catname2), 'ab') as f:
                                        if cross_matching_result_intermediate_cat2.ndim>1:
                                            np.savetxt(f, cross_matching_result_intermediate_cat2,
                                                   delimiter=",",header=header2)
                                        else:
                                            np.savetxt(f, cross_matching_result_intermediate_cat2[None],
                                                       delimiter=",", header=header2)
                                else:
                                    with open(output + '/cross-matching_result_{0}.txt'.format(Catname2), 'ab') as f:
                                        if cross_matching_result_intermediate_cat2.ndim>1:
                                            np.savetxt(f, cross_matching_result_intermediate_cat2,
                                                   delimiter=",")
                                        else:
                                            np.savetxt(f, cross_matching_result_intermediate_cat2[None],
                                                       delimiter=",")
                    #time checker:
                    #ongoing7 = time.time()
                    #print(ongoing7 - ongoing6)
                    #print(ongoing6 - ongoing5)#bcp
                    #print(ongoing5 - ongoing4)
                    #print(ongoing4-ongoing3)
                    #print(ongoing3-ongoing2)#bcp
                    #print(ongoing2-ongoing1)
                    #print(ongoing1-start)
                    #print(ongoing7-start)
                    #pdb.set_trace()
                else:
                    print('None of the trixels of catalog_2 ({0}) overlapping with trixel #{1} of catalog_1 ({2}) has sources in it'.format(Catname2,index_cat1,Catname1))
                    #pdb.set_trace()
            else:
                if Verbose == True:
                    print('trixel #{0} of Catalog_1 ({1}) is empty'.format(index_cat1,Catname1))
        if time_it==True:
            ongoing7 = time.time()
            print('it took {0} seconds for the process to run'.format(ongoing7 - start))


def read_ztf_HDF_matched(FieldID,Lines,ColCell=None,path=None):
    """
    Description: Read ZTF matched light curves from local HDF5 light curve files. The HDF5 files are distributed as part of the catsHTM catalogs.
     Input  : - ZTF field number.
              - [start end] lines to read. The lines for a given source are
                available in I1 and I2 in the 'ztfSrcLCDR1' catsHTM catalog.
              - ColCell'  - Column names for catalog.
                           Default is {'HMJD','Mag','MagErr','ColorCoef','Flags'}.
              - path to the data directory. Default is "."
     Output : - Catalog
              - ColCell
     By : Maayane Soumagnac. Trnslated from Eran O. Ofek's matlab routine with the same name
     URL : https://github.com/maayane/catsHTM; http://weizmann.ac.il/home/eofek/matlab/
     Example: Cat,ColCel=catsHTM.read_ztf_HDF_matched(815,[10,25],ColCell=None,path=path)
    """
    if ColCell is None:
        ColCell = np.array(['HMJD','Mag','MagErr','ColorCoef','Flags'])
    if path is None:
        path='.'

    FieldIDstring="{number:06}".format(number=FieldID)#'ztfLCDR1_%06d.hdf5'
    FileName = 'ztfLCDR1_'+FieldIDstring+'.hdf5'

    Cati = class_HDF5.HDF5(path+'/'+FileName).load(dataset_name='/AllLC',numpy_array=True)#,Offset=[Lines[0],Lines[1]-Lines[0]+1])#,Block=[Lines[1]-Lines[0], Ncol-1])
    Cat=Cati.T
    Cat_cut=Cat[Lines[0]-1:Lines[1],:]

    return Cat_cut,ColCell















