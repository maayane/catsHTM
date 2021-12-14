
"""*******************************************************
A python implementation of the celestial functions
******************************************************"""
#print __doc__
import math
import numpy as np
import pdb
import sys
#sys.setrecursionlimit(100000)
def cone_in_polysphere(PolesLong,PolesLat,Long,Lat,Radius):
    """
     Check if a cone (small circle) is within a convex spherical polygon
     Package:
     Description: Check if a cone (small circle) is within a convex spherical
                  polygon which sides are great circles.
     Input  : - Numpy array 3-by-N, in which each column represent the longitude of the
                poles of the half-spaces of a spherical polygin, where the
                pole is directed into the polygon center of mass [rad].
              - Numpy array in which each column represent the latitude of the
                poles of the half-spaces of a spherical polygin, where the
                pole is directed into the polygon center of mass [rad].
              - Vector of longitudes of the cones center [rad].
                The size is either 1 or like the number of columns of the
                first and second input arguments.
              - Vector of latitudes of the cones center [rad].
                The size is either 1 or like the number of columns of the
                first and second input arguments.
              - Vector of radii of the cones [rad].
                The size is either 1 or like the number of columns of the
                first and second input arguments.
     Output : - Flag of logical indicating if cone is in polygon.
         By : Maayane Soumagnac (original Matlab function by Eran Ofek)            Feb 2018
        URL : http://weizmann.ac.il/home/eofek/matlab/
     Example: HTM=celestial.htm.htm_build(4);
              Flag=celestial.htm.cone_in_polysphere(HTM(end).PolesCoo(:,1),HTM(end).PolesCoo(:,2),5.5,-0.6,0.01);
              PLong=rand(3,1000); PLat=rand(3,1000); Long=rand(1,1000); Lat=rand(1000,1); Radius=0.01.*rand(1000,1);
              Flag=celestial.htm.cone_in_polysphere(PLong,PLat,Long,Lat,Radius);
     Reliable:
    """

    #Longitudes_circle=Long # N longitudes de cercles
    #Latitudes_circle=Lat # N latitudes de cercles
    #Radius_circles=Radius # N radius de cercles

    Dist=np.arccos(np.multiply(np.sin(PolesLat),np.sin(Lat))+np.multiply(np.cos(PolesLat),np.cos(Lat))*np.cos(PolesLong-Long))
    Flag=np.zeros(np.shape(Dist)[1])
    for i in range(np.shape(Dist)[1]):#optimize
        Flag[i]=all(Dist[:,i]<=0.5*math.pi+Radius) #1 if all distances are smaller than..
    return Flag

def sphere_distance_fast(RA_1,Dec_1,RA_2,Dec_2):#RADIANS!

    Dist = np.arccos(np.sin(Dec_1)*np.sin(Dec_2) + np.cos(Dec_1)* np.cos(Dec_2)* np.cos(RA_1 - RA_2))

    return Dist

def sphere_dist_fast(RA_1,Dec_1,RA_2,Dec_2):#used by xmatch_2cats and match_cats
    """Description: Names after the function with the same name by Eran.
    Calculate the angular distance between two points on the celestial sphere. Works only with radians and calculate only the distance.
            Input: - np.array of longitudes for the first point [radians]
                   - np.array of latitudes for the first point [radians]
                   - np.array of longitudes for the second point [radians]
                   - np.array of latitudes for the second point [radians]
            Output: - np.array of distances between points [radians]"""

    Dist = np.arccos(np.sin(Dec_1)*np.sin(Dec_2) + np.cos(Dec_1)* np.cos(Dec_2)* np.cos(RA_1 - RA_2))
    dRA = RA_1 - RA_2
    SinPA = np.sin(dRA)* np.cos(Dec_2)/np.sin(Dist)
    CosPA = (np.sin(Dec_2)* np.cos(Dec_1) - np.cos(Dec_2)* np.sin(Dec_1) * np.cos(dRA))/ np.sin(Dist)
    PA = np.arctan2(SinPA, CosPA)
    #print(PA)
    if type(PA) is np.ndarray:
        #I = find(PA < 0);
        PA[(PA<0)] = 2. * math.pi + PA[(PA<0)]
    else:
        if PA<0:
            PA=2*math.pi+PA
    #print('Dist before nan_to_num',Dist)
    #print('RA1:',RA_1)
    #print('DEC1:',Dec_1)
    #print('RA2',RA_2)
    #print('DEC2',Dec_2)
    Distx=np.nan_to_num(Dist)
    #print('Dist after nan_to_num',Distx)
    return Distx,PA

def number_of_trixels_to_level(N_trixels):
    """Description: in Eran's library, this is called
        Input: - N_trixels: number of trixels
        Output: - number of levels
                - number of trixels in the lowest level"""
    number_of_levels=math.floor(math.log(N_trixels/2.)/math.log(4))
    number_of_trixels_in_highest_level=2*4**number_of_levels
    return number_of_levels,number_of_trixels_in_highest_level

def coo2cosined(Long,Lat):#TESTED compared with Eran's, ok
    """Description: Convert coordinates to cosine directions in the same reference frame.
    Input:-np.array of longitudes [radians]
          -np.array of lattitudes [radians]
    Output: - np.array of first cosine directions
            - np.array of second cosine directions
            - np.array of third cosine directions"""
    CosLat=np.cos(Lat)
    CD1=np.cos(Long)*CosLat
    CD2=np.sin(Long)*CosLat
    CD3=np.sin(Lat)
    return CD1,CD2,CD3

def cosined2coo(CD1,CD2,CD3):#TESTED compared with Eran's, ok
    """Description: Convert cosine directions into coordinated in the same reference frame.
    Input: - np.array of first cosine directions
           - np.array of second cosine directions
           - np.array of third cosine directions
    Output:-np.array of longitudes [radians]
           -np.array of lattitudes [radians]
           example: [RA,Dec]=cosined2coo(0.1,0,1)"""
    if type(CD1) is np.ndarray:
        #print(type(CD1))
        #print(type(CD2))
        #print(CD2[0])
        #print(CD1[0])
        #print(np.shape(CD1))
        #print(np.shape(CD2))
        Long=np.arctan2(CD2,CD1)
        SLL=np.sqrt(np.power(CD1,2)+np.power(CD2,2))
        Lat=np.zeros((np.shape(Long)))
        Lat[(SLL!=0)]=np.arctan(CD3[(SLL!=0)]/SLL[(SLL!=0)])
        Lat[(SLL==0)]=np.sign(CD3[SLL==0])*math.pi/2.
        Long[(Long<0)]=2*math.pi+Long[(Long<0)]
    else:
        Long=math.atan2(CD2,CD1)
        #print('Long is',Long)
        SLL=np.sqrt(CD1**2+CD2**2)
        #print('SLL is',SLL)
        #Lat=np.zeros((np.shape(Long)))
        if SLL!=0:
            Lat=np.arctan(CD3/SLL)
            #print('Lat is',Lat)
        else:
            Lat=np.sign(CD3)*math.pi/2.
        if Long<0:
            Long=2*math.pi+Long
    return Long,Lat

def cross_fast(A,B):#TESTED compared with Eran's, ok
    """Description: named after Eran's function with same name in his Util/math library.
     Performs cross product of two 3-columns matrices.
            Input: -First 3-elements vector
                   -Second 3-elements vector
            Output:-vector of cross product between the two input vectors"""
    C=np.array([A[:,1]*B[:,2]-A[:,2]*B[:,1],A[:,2]*B[:,0]-A[:,0]*B[:,2],A[:,0]*B[:,1]-A[:,1]*B[:,0]])
    #print('A is',A)
    #print('B is',B)
    #print('A[:,1]*B[:,2]-A[:,2]*B[:,1] is',A[:,1]*B[:,2]-A[:,2]*B[:,1])
    #print('C is',C)
    #C=C.flatten()
    #print('C is', C)
    #print(np.shape(C.T))
    #print(np.shape(C.T[0]))
    return C.T

def polysphere_poles(VertLong,VertLat):
    """Description: written after the function of same name by Eran. Given a spherical polygone (lattitude and longitude of the vertices),
    finds the poles (center) of each great circles
                Input: - Vector of the longitudes of the vertices
                       - Vector of the lattitudes of the vertices
                Output:- Vector of the poles longitudes
                       - Vector of the poles lattitudes """
    [CD1, CD2, CD3] = coo2cosined(VertLong, VertLat)
    [CenLong, CenLat] = cosined2coo(np.mean(CD1), np.mean(CD2), np.mean(CD3))
    PA=sphere_distance_fast(CenLong,CenLat,VertLong,VertLat)[1]
    SI=np.argsort(PA)
    VertLong_sortPA=VertLong[SI]
    VertLat_sortPA=VertLat[SI]
    CD1_sortPA=CD1[SI]
    CD2_sortPA=CD2[SI]
    CD3_sortPA=CD3[SI]
    #print(CD1_sortPA)
    #print(CD2_sortPA)
    #print(CD3_sortPA)
    Corners=np.zeros((np.shape(CD1_sortPA)[0]+1,3))
    Corners[:-1,0]=CD1_sortPA
    Corners[-1,0]=CD1_sortPA[0]
    Corners[:-1,1]=CD2_sortPA
    Corners[-1,1]=CD2_sortPA[0]
    Corners[:-1,2]=CD3_sortPA
    Corners[-1,2]=CD3_sortPA[0]
    #Corners=np.array([[CD1_sortPA,CD2_sortPA,CD3_sortPA],[CD1_sortPA[0],CD2_sortPA[0],CD3_sortPA[0]]])
    #print('Corners is',Corners)
    #print('the shape of Corners is',np.shape(Corners))
    #print('the shape of Corners[:-1,:] is',np.shape(Corners[:-1,:]))
    #pdb.set_trace()
    Polevec=cross_fast(Corners[:-1,:],Corners[1:,:])
    [PolesLong,PolesLat]=cosined2coo(Polevec[:,0],Polevec[:,1],Polevec[:,2])
    return PolesLong,PolesLat

def htm_build(Level,Verbose=True):
    """Description: build a hieratchical triangular mesh (HTM) class
            Input: - The number of levels in the HTM structure.
            Output: - HTM_list:A list of dictionnaries (one per trixel), with the following keys:
                    trixel['level']=level depth (0 for first level)
                    trixel['father']=index of the father. [] if no father
                    trixel['son']=np array of indexes of the sons. [] if no sons
                    - LevList: A list of dictionnaries (one per level), with the following keys:
                    LevList[i]['level']: level depth index (0 for the first level)
                    LevList[i]['ptr']: a np.array of the indexes of all the trixels in this level"""

    Ind=0
    HTM_list=[]
    #construct 0 level northern and southern hemisphere
    for i in range(4):
        Ind=Ind+1
        trixel = dict()
        #trixel['index']=Data[0,i]#line 1 of column 0
        trixel['index']=Ind
        trixel['level']=0
        trixel['father']=[]
        trixel['son']=[]
        trixel['coo'] = np.array([[0, 0], [math.pi / 2., 0], [0, math.pi / 2.]])#ok
        trixel['coo'][:, 0] = trixel['coo'][:, 0] + i * math.pi / 2.
        [CD1,CD2,CD3]=coo2cosined(trixel['coo'][:,0],trixel['coo'][:,1])
        trixel['cosd']=np.zeros((np.shape(CD1)[0],3))
        trixel['cosd'][:,0]=CD1
        trixel['cosd'][:,1]=CD2
        trixel['cosd'][:,2]=CD3
        #print(trixel['coo'])#ok same as Eran
        #print(trixel['coo'][:,0])#ok same as Eran
        #print(trixel['coo'][:,1])#ok same as Eran
        #print('at i={0}, coo[:,0] is {1}'.format(i,trixel['coo'][:,0])) #ok same as Eran
        #print('at i={0}, coo[:,1] is {1}'.format(i, trixel['coo'][:, 1])) #ok same as Eran
        [Poleslong,Poleslat]=polysphere_poles(trixel['coo'][:,0],trixel['coo'][:,1])
        #print('at i={0}, Poleslong is {1} and polesLat is {2}'.format(i,Poleslong,Poleslat))
        #trixel['PolesCoo']=list(zip(Poleslong,Poleslat))
        trixel['PolesCoo'] = np.zeros((np.shape(Poleslong)[0], 2))
        trixel['PolesCoo'][:, 0] = Poleslong
        trixel['PolesCoo'][:, 1] = Poleslat
        #print("trixel['PolesCoo'] is", trixel['PolesCoo'])
        #print(trixel)
        #print('at i={0}, PolesCoo is {1}'.format(i,trixel['PolesCoo']))
        #pdb.set_trace()
        HTM_list.append(trixel)
        #print(HTM_list)
        #pdb.set_trace()
    for i in range(4):
        Ind=Ind+1
        trixel = dict()
        #trixel['index']=Data[0,i]#line 1 of column 0
        trixel['index']=Ind
        trixel['level']=0
        trixel['father']=[]
        trixel['son']=[]
        trixel['coo'] = np.array([[0, 0], [0,-math.pi / 2.], [math.pi / 2.,0]])
        trixel['coo'][:, 0] = trixel['coo'][:, 0] + i * math.pi / 2.#ok
        [CD1,CD2,CD3]=coo2cosined(trixel['coo'][:,0],trixel['coo'][:,1])
        trixel['cosd']=np.zeros((np.shape(CD1)[0],3))
        trixel['cosd'][:,0]=CD1
        trixel['cosd'][:,1]=CD2
        trixel['cosd'][:,2]=CD3
        [Poleslong,Poleslat]=polysphere_poles(trixel['coo'][:,0],trixel['coo'][:,1])
        trixel['PolesCoo'] = np.zeros((np.shape(Poleslong)[0], 2))
        trixel['PolesCoo'][:, 0] = Poleslong
        trixel['PolesCoo'][:, 1] = Poleslat
        #print("trixel['PolesCoo'] is", trixel['PolesCoo'])
        HTM_list.append(trixel)

    Levlist=[]
    lev_list_dic=dict()
    lev_list_dic['level']=0
    lev_list_dic['ptr']=np.arange(8)+1
    Levlist.append(lev_list_dic)
    #print('Levlist')
    #pdb.set_trace()

    if Level!=0:
        total_HTM_list=htm_build_son(HTM_list,Levlist,Level,Ind)[0]
        total_Levlist = htm_build_son(HTM_list, Levlist, Level, Ind)[1]
    else:
        total_HTM_list=HTM_list
        total_Levlist=Levlist

    if Verbose == True:
        print("lev_list_dic[-1]['level'] is {0}".format(total_Levlist[-1]['level']))
        print("lev_list_dic[-1]['ptr'] is {0}".format(total_Levlist[-1]['ptr']))
        print("HT_list[-1]['level'] is {0}".format(total_HTM_list[-1]['level']))
        print("HT_list[-1]['father'] is {0}".format(total_HTM_list[-1]['father']))
        print("HT_list[-1]['son'] is {0}".format(total_HTM_list[-1]['son']))
        #print("the total Levlist is",total_Levlist)

    return total_HTM_list,total_Levlist

def htm_build_son(HTM_list, Levlist, Level, Ind,Verbose=False):
    """Description: function used by htm_build, to generate an HTM tree of trixels
             Input: - the HTM_list list of dictionnaries: A list of dictionnaries (one per trixel), with the following keys:
                        HTM_list[i]['level']=level depth (0 for first level)
                        HTM_list[i]['father']=index of the father. [] if no father
                        HTM_list[i]['son']=np array of indexes of the sons. [] if no sons
                    - the Levlist list of dictionnaries: A list of dictionnaries (one per level), with the following keys:
                        LevList[i]['level']: level depth index (0 for the first level)
                        LevList[i]['ptr']: a np.array of the indexes of all the trixels in this level
                    - the number of levels required
                    - the last populated index
             Output: - the new HTM_list list of dictionnaries
                     - the new Levlist list of dictionnaries
                    """
    LevelDepth=len(Levlist)#the number of levels so far
    if Verbose==True:
        print('the number of levels so far is {0}'.format(LevelDepth))
        print(type(LevelDepth))
        print('Levlist is {0}'.format(Levlist))
        #pdb.set_trace()
    if Level > LevelDepth: #if there are more Levels to build than the number of levels so far

        #creation du nouveau dictionnaire dans la liste de levels
        lev_list_dic=dict() #je creee un nouveau dictionnaire "level" a mettre dans Levlist
        lev_list_dic['level']=Levlist[LevelDepth-1]['level']+1 # son level index (0 pour 1) est le level index du dernier dictionnaire de LevList +1
        lev_list_dic['ptr']=np.zeros((8*4**(LevelDepth)),dtype=int) #son nombre de trixels est donne par cette formule


        #creation du nouveau dictionnaire dans la liste des trixels
        N_trixels=len(Levlist[LevelDepth-1]['ptr']) #le nombre de trixels dans le plus haut level (le level "pere" de celui qu'on est sur le point de creer)
        K=0
        for i in range(N_trixels): #parcourons ce nombre. Ppour chaque trixel:
            FatherPtr= Levlist[LevelDepth-1]['ptr'][i] # prenons l'index du trixel
            FatherLevel=HTM_list[FatherPtr-1]['level'] # demandons quel est le level de ce trixel
            Vert1=HTM_list[FatherPtr-1]['cosd'][0,:]
            Vert2=HTM_list[FatherPtr-1]['cosd'][1,:]
            Vert3=HTM_list[FatherPtr-1]['cosd'][2,:]
            Cen=gc_mid_section(np.array([Vert1,Vert2,Vert3]),np.array([Vert2,Vert3,Vert1]))
            #print(Cen)
            #pdb.set_trace()
            #maintenant pour chaque trixel, creons les nouveaux trixels: les fils de ceux que l'on est en train de parcourir
            for j in range(4):
                Ind=Ind+1 #on augmente le plus haut index de 1
                trixel=dict() #le dictionnaire qu'on va ajouter
                trixel['level']=FatherLevel+1
                trixel['cosd'] = np.zeros((3, np.shape(Vert1)[0]))
                if j==0:
                    trixel['cosd'][0,:]=Vert1
                    trixel['cosd'][1,:]=Cen[0,:]
                    trixel['cosd'][2,:]=Cen[2,:]
                if j == 1:
                    trixel['cosd'][0, :] = Vert2
                    trixel['cosd'][1, :] = Cen[1, :]
                    trixel['cosd'][2, :] = Cen[0, :]
                if j == 2:
                    trixel['cosd'][0, :] = Vert3
                    trixel['cosd'][1, :] = Cen[2, :]
                    trixel['cosd'][2, :] = Cen[1, :]
                if j == 3:
                    trixel['cosd'][0, :] = Cen[0, :]
                    trixel['cosd'][1, :] = Cen[1, :]
                    trixel['cosd'][2, :] = Cen[2, :]
                trixel['coo']=np.zeros((3,2))
                trixel['coo'][:,0]=cosined2coo(trixel['cosd'][:,0],trixel['cosd'][:,1],trixel['cosd'][:,2])[0]
                trixel['coo'][:,1] = cosined2coo(trixel['cosd'][:, 0], trixel['cosd'][:, 1], trixel['cosd'][:, 2])[1]
                #print(trixel['coo'])
                #print(np.shape(trixel['coo']))
                #pdb.set_trace()
                #trixel['coo']=np.array([[0,0],[math.pi/2.,0],[0,math.pi/2.]])
                #print('FatherPtr is',FatherPtr)
                #print(type(FatherPtr))
                trixel['father']=FatherPtr
                trixel['son']=[]
                HTM_list.append(trixel)
                # populates the sons keys of the fathers with the current IDs
                HTM_list[FatherPtr-1]['son']=np.append(HTM_list[FatherPtr-1]['son'],[Ind])
                K=K+1
                lev_list_dic['ptr'][K-1]=Ind
                #print(trixel['coo'])
                [PoleLong,PoleLat]=polysphere_poles(trixel['coo'][:,0],trixel['coo'][:,1])
                trixel['PolesCoo']=np.zeros((np.shape(PoleLong)[0],2))
                trixel['PolesCoo'][:,0]=PoleLong
                trixel['PolesCoo'][:,1]=PoleLat



                #print("trixel['PolesCoo'] is",trixel['PolesCoo'])

        Levlist.append(lev_list_dic)

    if Level >= LevelDepth+1:#s'il reste des level a creer: continuer en faisant la meme procedure
        HTM_list=htm_build_son(HTM_list,Levlist,Level,Ind)[0]
        Levlist=htm_build_son(HTM_list,Levlist,Level,Ind)[1]

    return HTM_list,Levlist

def htm_search_cone(HTM,Long,Lat,Radius,Ind):#NOT TESTED
    """Description: Names after the function with the same name by Eran.
    Search for all trixels intersecting a small circle (cone search)
            Input: - HTM structure.
                   - Longitude [radians] of the center of the circle.
                   - Latitude [radians] of the center of the circle.
                   - Search radius [radians].
            Output: - list of indexes of the intersecting trixeels"""

    #print('at the begginging of the code, Ind is',Ind)
    if Ind==[]:
        Sons = np.arange(8)+1
    else:
        Sons=Ind


    ID=[]
    Nsons=np.shape(Sons)[0]
    Poleslong=np.zeros((3,Nsons))
    Poleslat=np.zeros((3,Nsons))

    for i in range(Nsons):
        #print(i)
        #print(HTM[Sons[i]])
        #print(HTM[Sons[i]]['PolesCoo'])
        #print(np.shape(HTM[Sons[i]]['PolesCoo']))
        #print('Sons[i] is',Sons[i])
        Poleslong[:,i]=HTM[int(Sons[i])-1]['PolesCoo'][:,0]
        Poleslat[:,i]=HTM[int(Sons[i])-1]['PolesCoo'][:, 1]

    Flag=cone_in_polysphere(Poleslong,Poleslat,Long,Lat,Radius)
    #print('Radius is',Radius) #pas ok
    #print('Lat is',Lat) #ok
    #print('Long is',Long) #ok
    #print('HTM[0] is',HTM[0])#ok compared with Eran
    #print('Poleslat is',Poleslat) #ok compared with Eran
    #print('Poleslong is',Poleslong)#ok
    #print('Flag is',Flag) #pas ok

    #pdb.set_trace()
    #print(Sons)
    for i in range(Nsons):
        #print('at the begginging of the loop, Ind is', Ind)
        #print('at the beggining of the loop, Ind is of type',type(Ind))
        #print('i is {0} and Flag[i] is {1}'.format(i, Flag[i]))
        if Flag[i]==1.:
            CSons=Sons[i]
            #print('HTM[int(CSons)-1] is',HTM[int(CSons)-1])
            #print('Cson is',CSons)
            #print('Ind is',Ind)
            #print('ID is',ID)
            if HTM[int(CSons)-1]['son']==[]:
                #print('ID before appending',ID)
                ID.append(CSons)
                #print('ID after appending',ID)
                #pdb.set_trace()
                #print("HTM[int(CSons)-1]['son']==[]")
                #if ID==[]:
                #    print('ID is []?',ID)
                #    ID=[CSons]
                #else:
                #    #print('I am going to append {0} to ID'.format(CSons))
                #    #print(ID)
                #    #print(CSons)
                #    print('ID right before appending:',ID)
                #    ID.append(CSons)
                #    print('ID right after appending:', ID)
                #   #print('now ID is',ID)
            else:
                #print("HTM[int(CSons)-1]['son']is not[]")
                Ind=HTM[int(CSons)-1]['son']
                #print('Ind is now',Ind)
                #print('ID right before appending #2:', ID)
                ID.append(htm_search_cone(HTM, Long, Lat, Radius, Ind))
                #print('ID right after appending #2:', ID)
                #if isinstance(ID[0],list):
                #    ID=ID[0]
                #if ID==[]:
                #    #print('noooooooo')
                #    print('ID is [] #2?', ID)
                #    ID=[htm_search_cone(HTM,Long,Lat,Radius,Ind)]
                #    #print('ID is now',ID)
                #else:
                #    #print('yesssss')
                #    print('ID right before appending #2:', ID)
                #    ID.append(htm_search_cone(HTM,Long,Lat,Radius,Ind))
                #    print('ID right after appending #2:', ID)
            #print(i)
            #print(Flag[i])
            #print(ID)
            #print(Ind)
        '''
        else:
            print('Flag was 0, i,Flag, ID and Ind are:')
            print(i)
            print(Flag[i])
            print(ID)
            print(Ind)
        '''
        #pdb.set_trace()
    return ID

def gc_mid_section(Pos1,Pos2):
    """Description: Named after Eran's function of same name. Given two points on a sphere, find the central point lying
    on the the shortest great circle section connecting the two points.
            Input: - two numpy arrays, each row is a point, the columns are unit vectors.
%          - A list of the second points (similar the the first point).
%          - Dimension along to operate.
%            If 1, then assume the position vectors in (R) and (C) are
%            in rows, while if 2, then assume they are in columns.
%            Default is 1.
%          - Output type: unit vectors.
% Output : - The coordinates, either unit vectors or [Long, Lat]
%            in radians, of the mid point on the great circle
%            connecting the two points.
            Output: - """

    #[P1,Q1]=np.shape(Pos1)
    #[P2,Q2]=np.shape(Pos2)
    #if Q1==2:
    #    [CD1_1,CD2_1,CD3_1]=coo2cosined(Pos1[0],Pos1[1])
    #    Pos1_tr=np.array([CD1_1,CD2_1,CD3_1])
    #else:
    #    Pos1_tr=Pos1

    #if Q2==2:
    #    [CD1_2, CD2_2, CD3_2] = coo2cosined(Pos2[0], Pos2[1])
    #    Pos2_tr = np.array([CD1_2, CD2_2, CD3_2])
    #else:
    #    Pos2_tr=Pos2

    C=Pos1+Pos2
    #print(C)
    #print(np.sum(C,axis=1))
    #print(np.power(C,2))
    #print(np.sqrt(np.sum(np.power(C,2),axis=1)))
    #print np.sum()

    Cx=np.zeros(np.shape(C))
    if C.ndim >1:
        for i in range(np.shape(C)[0]):
            Cx[i,:]=C[i,:]/np.sqrt(np.sum(np.power(C,2),axis=1)[i])
        return Cx
    else:
        Cx=C/np.sqrt(np.sum(np.power(C,2)))


    #print(np.sqrt(np.sum(np.power(C,2)),axis=1))
    return Cx
