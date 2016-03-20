'''
This file takes the graphs (digitized from the pdf) for a simulated trajectory of Pegasus
    and turns them into a nominal trajectory
'''

# Points to the python scripts needed from Francisco
friscoFiles = '../../../../Prop3Dof/FriscoDebris/pythonFiles/'
# Points to the binaries for propagating trajectories
debrisPropPATH = '../../../../Prop3Dof/FriscoDebris/'
# Points to the files that I've written
tjcFiles = '../../../'

import os
import sys
sys.path.append(friscoFiles)
# sys.path.append(debrisPropPATH)
sys.path.append(tjcFiles)

import numpy as np
from numpy import sin, cos
from numpy.linalg import norm

import matplotlib.pyplot as plt
from scipy import interpolate

import orbitTools as ot

import LaunchSites
import LaunchProviders

# from scipy.interpolate import interp1d

'''
Define the parameters that we want the graphs interpolated to
'''

vehicleName     = LaunchProviders.Pegasus
launchSiteName  = LaunchSites.PegMARS
targetLocation  = LaunchSites.siteDict[launchSiteName]     #this is a custom one that'll be out over the atlantic near MARS

initLon     = 0.
initLat     = 0.

targetLat   = targetLocation['lat']     # Degrees, obviously
targetLon   = targetLocation['lon']

deltaTsec   = 1
timeStart   = 0
timeEnd     = 286   # It looked like this was about the end
timeVec     = np.arange(timeStart, timeEnd + deltaTsec, deltaTsec)

rEarth      = 6371  #km
omegaE      = 7.2921158494529352e-05    #rad/s
planetModel = 0                             # 0 means spherical, 1 means elliptical

doPlot      = True


'''
Process local Velocity [m/s].  Do this one with a linear interp
'''

# These will store the x and y components of the velocity graph (time, velocity magnitude)
xData = []
yData = []

fileData = open('_localVelocity.txt','r')
for line in fileData:
    key = line.split()
    if (len(key) > 0 ) and not (key[0][0] == '#'):
        xData.append(float(key[0]))
        yData.append(float(key[1]))
fileData.close()

# Do the interpolation
f               = interpolate.interp1d(xData, yData)
localVelocity   = f(timeVec)

if doPlot:
    # Plot it
    plt.figure()
    plt.plot(xData, yData, 'rx', timeVec, localVelocity, 'b')
    plt.title('Cubic-spline interpolation')
    plt.show()


'''
Process Altitude [km].  Using cubic.
'''

# These will store the x and y components of the velocity graph (time, velocity magnitude)
xData = []
yData = []

fileData = open('_altitude.txt','r')
for line in fileData:
    key = line.split()
    if (len(key) > 0 ) and not (key[0][0] == '#'):
        xData.append(float(key[0]))
        yData.append(float(key[1]))
fileData.close()

# Do the interpolation
f           = interpolate.interp1d(xData, yData, 'cubic')
altitude    = f(timeVec)

if doPlot:
    # Plot it
    plt.figure()
    plt.plot(xData, yData, 'rx', timeVec, altitude, 'b')
    plt.title('Cubic-spline interpolation')
    plt.show()



'''
Process Flight Path Angle [Deg].  Do this one with a linear interp.
'''

# These will store the x and y components of the velocity graph (time, velocity magnitude)
xData = []
yData = []

fileData = open('_flightPathAngle.txt','r')
for line in fileData:
    key = line.split()
    if (len(key) > 0 ) and not (key[0][0] == '#'):
        xData.append(float(key[0]))
        yData.append(float(key[1]))
fileData.close()

# Do the interpolation
f           = interpolate.interp1d(xData, yData)
FPA         = f(timeVec)

if doPlot:
    # Plot it
    plt.figure()
    plt.plot(xData, yData, 'rx', timeVec, FPA, 'b')
    plt.title('Cubic-spline interpolation')
    plt.show()



'''
Now translate these values into state vectors
'''
v0 = rEarth * omegaE    # the speed of the rotating earth at time of launch (should add in initial altitude)

# rStorage = np.array([[]])
# vStorage = []
altStorage = []

tx = 0
curAlt = altitude[tx]
rPeg = np.array([rEarth + curAlt, 0 , 0])    #[km]

curFPA = FPA[tx] * np.pi/180.
curVel = localVelocity[tx]*1e-3
vPeg = np.array([sin(curFPA), cos(curFPA), 0]) * (v0 + curVel) #[km/s]

thetag = timeVec[tx]*omegaE     # [rad] how far has the earth rotated


# Initialize storage
rStorage = np.array([rPeg])
vStorage = np.array([vPeg])
stateVecStorage = np.array([np.concatenate((rPeg*1.e3, vPeg*1.e3, np.array([thetag])))])


# rStorage = np.concatenate((rStorage, np.array([rPeg])))
# vStorage = np.concatenate((vStorage, vPeg))
altStorage.append(norm(rPeg) - rEarth)

# If i angle everything down, then the trajectories are close enough
# degFudge = -5   # This best approximates the trajectory from the paper
# degFudge = -3   # 200km alt
# degFudge = 0    # 250km alt
degFudge = 1.0    # 400km alt     This seems way more like the youtube video!

# Now run through all the times and get the rPegs
for tx in range(1,len(timeVec)):
    # Get the last rVector
    rLast   = rStorage[-1]

    # Find the current inertial velocity
    curFPA  = (FPA[tx] + degFudge) * np.pi/180.
    curVel  = localVelocity[tx]*1e-3
    vPeg    = np.array([sin(curFPA), cos(curFPA), 0]) * (v0 + curVel) #[km/s]

    # Update the ECI rvector
    rPeg    = rLast + vPeg * deltaTsec

    # Get the thetag for later
    thetag = timeVec[tx]*omegaE     # [rad] how far has the earth rotated
    curStateVec = np.concatenate((rPeg*1.e3, vPeg*1.e3, np.array([thetag])))

    # Save
    rStorage = np.concatenate((rStorage, np.array([rPeg])))
    vStorage = np.concatenate((vStorage, np.array([vPeg])))
    stateVecStorage = np.concatenate((stateVecStorage, np.array([curStateVec])))
    # thetagStorage.append(thetag)

    # Since zion and I use old versions of numpy, i can't norm by axis, so save it here
    altStorage.append(norm(rPeg) - rEarth)

# sys.exit()

if doPlot:
    # Now check by plotting the alititudes
    plt.figure()
    plt.plot(timeVec, altitude, 'rx', timeVec, altStorage, 'b')
    plt.title('Cubic-spline interpolation')
    plt.show()


'''
Rotate the state vectors
'''

# Rename
tfailStorage = timeVec


print '==================== GENERATING A NEW MAIN TRAJECTORY ===================='
fo = open(vehicleName + launchSiteName + "NominalTrajectory.txt", "w")
outString = "FORMAT LATLON METERS"
print str(outString)
fo.write( outString + "\n")

outString = "{0:6}{1:20}{2:20}{3:20}".format("#Time", "LonDeg", "LatDeg","AltM")
print str(outString)
fo.write( outString + "\n")

for ix in range(len(tfailStorage)):

    curVec = stateVecStorage[ix]
    thetag = curVec[-1]
    ecef_C_eci = ot.ecef_C_eci(thetag)
    r_ecef = np.dot(ecef_C_eci,curVec[:3])  # multiply them
    # Now rotate over by lon
    rotated = np.dot(ot.cRnMatrix_Zaxis(initLon * np.pi/180.), r_ecef)
    # Now adjust the latitude
    rotated = np.dot(ot.cRnMatrix_Yaxis( -(targetLat-initLat) * np.pi/180.), rotated)
    # Now that latitude is good, move to target Lon
    rotated = np.dot(ot.cRnMatrix_Zaxis(-targetLon * np.pi/180.), rotated)

    # Check em
    [curLat, curLon, curAlt] = ot.ECEF2latlonalt(rotated,planetModel)
    # print 'current location is lat = {0}, lon = {1}, alt = {2}'.format(curLat, curLon, curAlt)

    outString = "{0:6}{1:20}{2:20}{3:20}".format(str(tfailStorage[ix]), str(curLon), str(curLat), str(curAlt))
    print str(outString)

    fo.write( outString + "\n")
fo.close()







'''
Dump the state vectors
'''









# altitude[1]









# vEarth = np.array([sin(thetag), cos(thetag), 0]) * rEarth * omegaE  #[km/s]
