
"""*******************************************************
A python implementation of catsHTM.m
******************************************************"""
#print __doc__

import math
import numpy as np
import celestial
import class_HDF5
import scipy.io as sio
import params

#class params(object):
#    def __init__(self,path_catalogs,IndexFileTemplate,CatFileTemplate,htmTemplate,NcatinFile,ColCelFile):
#	self.path_catalogs='/Users/maayanesoumagnac/PostDoc/projects/catsHTM/data/'
#	self.IndexFileTemplate='%s_htm.hdf5'
#	self.CatFileTemplate='%s_htm_%06d.hdf5'
#	self.htmTemplate='htm_%06d'
#	self.NcatinFile=100.
#	self.ColCelFile = '%s_htmColCell.mat'
#root_to_data=params.path_catalogs

def cone_search(CatName,RA,Dec,Radius,catalogs_dir='./data',RadiusUnits='arcsec',IndexFileTemplate=params.IndexFileTemplate,CatFileTemplate=params.CatFileTemplate
                ,htmTemplate=params.htmTemplate,NcatinFile=params.NcatinFile,IndexVarname=None,ColRa = 0,ColDec=1,OnlyCone=True,
                ColCelFile = params.ColCelFile,OutType= 'np_array'):
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
                        IndexVarName - Default is [].
                        ColRA       - Default is 1.
                        ColDec      - Default is2.
                        OnlyCone    - Return only sources within cone. If false will return also some objects outside cone. Default is true.
                        ColCellFile - Default is '%s_htmColCell.mat'.
    By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018
    Output  : a numpy array where each line is the catalog line for the sources inside the cone """

    print '*************'
    print 'Catalog: {0}; cone radius: {1} arcsec; cone center: (RA,DEC)=({2},{3})'.format(CatName,Radius,RA,Dec)
    print '*************'

    root_to_data=catalogs_dir+'/'
    #catdir definition
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

    if RadiusUnits=='arcsec':
        Radius=Radius/(Rad*3600) #converts arcsec radius into radians radius
    #TO DO: the else

    ColCelFile=ColCelFile % CatName
    IndexFilename=IndexFileTemplate % CatName

    test = sio.loadmat(root_to_data+CatDir+'/'+ColCelFile) #CHANGE THIS
    Ncol=np.shape(test['ColCell'])[1]


    ID=search_htm_ind(IndexFilename,RA,Dec,Radius,catalogs_dir,VarName=IndexVarname) #list of IDs of winners leaf
    ID_matlab=ID+1
    FileID=np.floor(ID_matlab/NcatinFile)*NcatinFile
    Nid=np.shape(ID_matlab)[0] #number of leaf intercepting the circle
    cat=[]
    if Nid==0:
	print 'INFO: the cone does not intercept the catalog'
	cat_onlycone=np.array(cat)
    else:
	    for Iid in range(Nid):
		FileName=CatFileTemplate % (CatName, FileID[Iid])

		DataName=htmTemplate % ID_matlab[Iid]

		if Iid==0:
		    cat=class_HDF5.HDF5(root_to_data+CatDir+'/'+FileName).load(DataName,numpy_array=True).T
		else:
		    cat = np.vstack(
			(cat, class_HDF5.HDF5(root_to_data + CatDir + '/' + FileName).load(DataName, numpy_array=True).T))

	    if OnlyCone==True:
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
        if np.shape(test['ColUnits'][0,i])[0]>0:
            ColUnits[i]=(test['ColUnits'][0,i][0])
        else:
            ColUnits[i]=' '

    return cat_onlycone,ColCell, ColUnits


def search_htm_ind(Filename,Long,Lat,Radius,path,VarName=None):
    """Description: wrapper of htm_search_cone, which select from the vector outputed by htm_search_cone only the
    triangles where there are actually sources.
            Input  : - Filename: the name of the index_file, e.g. FIRST_htm.hdf5
            Output :
            By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018

                    """
    #HDF5_file = class_HDF5.HDF5(Filename) #e.g '/Users/maayanesoumagnac/PostDoc/projects/catsHTM/data/FIRST/FIRST_htm.hdf5'
    if VarName==None:
        cat_name=Filename.split('_')[0]
        VarName=cat_name+'_HTM'
    if cat_name=='TMASS':
        CatDir='2MASS'
    elif cat_name=='TMASSxsc':
        CatDir='2MASSxsc'
    elif cat_name=='DECaLS':
        CatDir='DECaLS/DR5'
    elif cat_name=='GAIADR1':
        CatDir='GAIA/DR1'
    elif cat_name=='GALEX':
        CatDir='GALEX/DR6Plus7'
    elif cat_name=='PS1':
        CatDir='/PS1'
    elif cat_name=='SDSSDR10':
        CatDir='SDSS/DR10'
    elif cat_name=='SDSSoffset':
        CatDir='SDSS/DR14offset'
    elif cat_name=='UKIDSS':
        CatDir='UKIDSS/DR10'
    elif cat_name=='VISTAviking':
        CatDir='VISTA/Viking/DR2'
    elif cat_name=='VSTatlas':
        CatDir='VST/ATLAS/DR3'
    elif cat_name=='VSTkids':
        CatDir='VST/KiDS/DR3'
    else:
        CatDir=cat_name

    DataHTM_indexfile= class_HDF5.HDF5(path+'/'+CatDir+'/'+Filename).load(VarName,numpy_array=True) #load the indexfile content
    ID=htm_search_cone(DataHTM_indexfile,Long,Lat,Radius) # ID of the winners mesh, i.e. the meshes that intercept the circle
    ID_array=np.array(ID)
    ID_w_sources=ID_array[DataHTM_indexfile[12,ID]>0] #ou l inverse?
    return ID_w_sources

def htm_search_cone(DataHTM,Long,Lat,Radius,Ind=None):
    """Description: Search for all HTM leafs intersecting a small circles
        Input  : - Either a table of HTM data od an open HDF5 object in which the HTM data is stored
                - Longitude (radians) to search
                -Latitutde (raidans) to search
                -Radius of the small circle
        Output : a vector of indexes of the winner(s):the "adress" in the indexfile of the smallest leaf(s) intercepting the cone
        By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018

                """

    Father_index = 1
    Son_index = np.asarray([2, 3, 4, 5])
    PolesLong_index = np.asarray([6, 8, 10])
    PolesLat_index = np.asarray([7, 9, 11])
    if Ind is None:
        Sons=np.linspace(0,7,8).astype(int) #on debute avec 8 mesh
    else:
        Sons=Ind.astype(int)

    ID=[]#np.empty(1)
    Nsons=np.shape(Sons)[0] 
    PolesLong=np.zeros((3,Nsons)) #3 lignes, Nsons colomnes, on veut mettre a chaque colomne les longitudes des poles du mesh
    PolesLat=np.zeros((3, Nsons)) #3 lignes, Nsons colomnes
    for i in range(Nsons):
        PolesLong[:,i]=DataHTM[PolesLong_index[:].astype(int),Sons[i]] # array where each colomn is the 3 poles longitudes of a son mesh HERE: THIS? OR INVERSE?
        PolesLat[:,i]=DataHTM[PolesLat_index[:].astype(int),Sons[i]] # array where each colomn is the 3 poles latitude of a son mesh HERE: THIS? OR INVERSE?

    Flag=celestial.cone_in_polysphere(PolesLong,PolesLat,Long,Lat,Radius) #check if the cone intercept any of the sons meshes
    for i in range(Nsons):
        if Flag[i]==1: #the cone overlap the son with index i
            CSon=Sons[i] # l'index du son vainqueur
            if np.isnan(DataHTM[Son_index[:],CSon]).all()==True:# there are nans in the index_file at the son's index, which means the data is where you are and you cannot go further in the tree
               ID.append(CSon)
            else:
                Ind=DataHTM[Son_index[:],CSon]-1.
                ID.extend(htm_search_cone(DataHTM,Long,Lat,Radius,Ind=Ind))
    return ID


