
"""*******************************************************
A python implementation of the celestial functions
******************************************************"""
#print __doc__
import math
import numpy as np

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

def sphere_distance_fast(RA_1,Dec_1,RA_2,Dec_2):
    """Description:
    Input  : - """

    Dist = np.arccos(np.sin(Dec_1)*np.sin(Dec_2) + np.cos(Dec_1)* np.cos(Dec_2)* np.cos(RA_1 - RA_2))

    '''
    dRA = RA_1 - RA_2
    SinPA = np.sin(dRA)* np.cos(Dec_2)/np.sin(Dist)
    CosPA = (np.sin(Dec_2)* np.cos(Dec_1) - np.cos(Dec_2)* np.sin(Dec_1) * np.cos(dRA))/ np.sin(Dist)
    PA = np.atan2(SinPA, CosPA)

        I = find(PA < 0);
        PA(I) = 2. * pi + PA(I);

    end
    '''
    return Dist



