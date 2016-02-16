#***************************************************************************
# File : AtmosProfile.py
# Reads the atmospheric profile obtained from a text file (using GRAM) and formatted as 
# time height(km) lat lon density Umean Vmean Wmean density SD USD VSD WSD
#
# Created by: Francisco Capristan
#****************************************************************************

import sys
import numpy as np


def readGramAtmos(fileName):

	try:
        
		#inputFile = open(sys.argv[1],'r')
		inputFile = open(fileName,'r')
	except :
		print '\n Error: could not load the Atmospheric profile input file \n'

		exit(1)


	altitudeList = []
	densityMeanList = []
	uMeanList = []
	vMeanList = []
	wMeanList = []
	densitySDList = []
	uSDList = []
	vSDList =[]
	wSDList = []
	lineNumber = 0
		
	for line in inputFile:
	    lineNumber = lineNumber +1
	    if lineNumber > 1:
		key = line.split()

		altitude = [1000.*float(key[1])] #converting from km to meters
		densityMean = [float(key[4])]
		uMean = [float(key[5])]
		vMean = [float(key[6])]
		wMean = [float(key[7])]
		densitySD = [(densityMean[0]*float(key[8]))/100.]
		uSD = [float(key[9])]
		vSD = [float(key[10])]
		wSD = [float(key[11])]

		altitudeList.append(altitude)
		densityMeanList.append(densityMean)
		uMeanList.append(uMean)
		vMeanList.append(vMean)
		wMeanList.append(wMean)
		densitySDList.append(densitySD)
		uSDList.append(uSD)
		vSDList.append(vSD)
		wSDList.append(wSD)

	inputFile.close()
	altitudeList = np.matrix(altitudeList)
	densityMeanList = np.matrix(densityMeanList)
	uMeanList = np.matrix(uMeanList)
	vMeanList = np.matrix(vMeanList)
	wMeanList = np.matrix(wMeanList)
	densitySDList = np.matrix(densitySDList)
	uSDList = np.matrix(uSDList)
	vSDList =np.matrix(vSDList)
	wSDList = np.matrix(wSDList)

	row,col = np.shape(altitudeList)
	
	return(altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,row)


def generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList):
	import random

	density = random.normalvariate(densityMeanList,densitySDList)
	u = random.gauss(uMeanList,uSDList)
	v = random.gauss(vMeanList,vSDList)
	w = random.gauss(wMeanList,wSDList)

	return(density,u,v,w)



def maxAltitudeAtmosFix(altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,row):
# function to ensure that altitude and the other quantities of interest are a one to one function, based on a parabolic type trajectory
# this is not necessary if atmospheric profile is already organized as decreasing in altitude

    newalt = []
    newdenMean =[]
    newuMean = []
    newvMean =[]
    newwMean =[]
    newdenSD = []
    newuSD = []
    newvSD = []
    newwSD = []
    
    
    oldalt = altitudeList[1]
    #newalt.append(oldalt)
    for index in range(0,row):
        if (np.abs(oldalt - altitudeList[index])> 2.0) or (altitudeList[index,0]<0.0): # 2 meter difference in alttitude
            newalt.append([altitudeList[index,0]])
            newdenMean.append([densityMeanList[index,0]])
            newuMean.append([uMeanList[index,0]])
            newvMean.append([vMeanList[index,0]])
            newwMean.append([wMeanList[index,0]])
            newdenSD.append([densitySDList[index,0]])
            newuSD.append([uSDList[index,0]])
            newvSD.append([vSDList[index,0]])
            newwSD.append([wSDList[index,0]])
            oldalt = altitudeList[index,0]
        
    
    
    altitudeList = np.matrix(newalt)
    densityMeanList = np.matrix(newdenMean)
    uMeanList = np.matrix(newuMean)
    vMeanList = np.matrix(newvMean)
    wMeanList = np.matrix(newwMean)
    densitySDList = np.matrix(newdenSD)
    uSDList = np.matrix(newuSD)
    vSDList = np.matrix(newvSD)
    wSDList = np.matrix(newwSD)
    
#print wSDList
        
    for index in range(0,len(altitudeList)-1):
        if altitudeList[index,0] > altitudeList[index+1,0]:
            altitudeList = altitudeList[index:row,0]
            densityMeanList = densityMeanList[index:row,0]
            uMeanList = uMeanList[index:row,0]
            vMeanList = vMeanList[index:row,0]
            wMeanList = wMeanList[index:row,0]
            densitySDList = densitySDList[index:row,0]
            uSDList = uSDList[index:row,0]
            vSDList = vSDList[index:row,0]
            wSDList = wSDList[index:row,0]
            break
    rowNew,col = np.shape(altitudeList)
                #print rowNew,len(altitudeList),col,len(wSDList)

    if rowNew==row:
        print 'Atmospheric Profile already a ONE to ONE function'
    else :
        print 'Adjusting Atmospheric Profile to a ONE to ONE function'
        row=rowNew
		
    return(altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,row)


def getRegion(latarray,lonarray,altarray,Vmagarray,betaaray,gammaarray,mainVelRelOption):
# writes the text files from GRAM for a rectangular area
# it runs all the debris cases in vaccuum to get the max min lat and lon. The desired region is selected from this

    import sys
    gramPath = "../SourceCode/"
    sys.path.append(gramPath)
    import debrisPropagation as DP
    
    maxLonTotal = -9999.9
    minLonTotal = 9999.9
    minLatTotal = 9999.9
    maxLatTotal = -9999.9
        

    altitudeList = 1
    atmosoption = -1
    ncd = 1
    ncl = 1
    nlist = 1
    densitylist = 0.0
    ulist = 0.0
    vlist = 0.0
    wlist = 0.0
    planetModel = 0
    geoptions = 1
    filename = "none"

    debrisRelVel = np.array([0.0,0.0,0.0])

    for index in range(0,len(latarray)):
        initialposition = np.array([latarray[index],lonarray[index],altarray[index]])
        mainVel = np.array([Vmagarray[index],gammaarray[index],betaaray[index]])
        loverd = 0.0
        cd = 0.0
        cl = 0.0
        mass = 1.0
        aref = 1.0
        cloption = 1.0
        minfcd = 1.
        minfcl = 1.


        finalConditions = DP.debrispropagation(initialposition, mainVel,mainVelRelOption,
                                               debrisRelVel, mass,aref, # debris piece state vector wrt to main vehicle
                                               minfcd,cd,cloption,minfcl,cl,loverd, # aerodynamic inputs
                                               atmosoption,altitudeList, # dt for RK4, angrate for tumbling debris. Atmospheric options
                                               densitylist,ulist,vlist,wlist, # atmospheric profile
                                               geoptions,filename,planetModel,[ncd,ncl,nlist]) # geoptions & filename used for Google Earth visualiation. ncd => length of CD array
        


        maxLat = finalConditions[4]
        minLat = finalConditions[5]
        maxLon = finalConditions[6]
        minLon = finalConditions[7]


        if (maxLatTotal<maxLat):
            maxLatTotal = maxLat
        if (minLatTotal>minLat):
            minLatTotal = minLat
        if (maxLonTotal < maxLon):
            maxLonTotal = maxLon
        if (minLonTotal > minLon):
            minLonTotal = minLon
    return (maxLatTotal,minLatTotal,maxLonTotal,minLonTotal)

def writeTXTfiles(maxLat,minLat,maxLon,minLon,ds,alt,deltaalt,gramPath):
# ds is the square spacing in degrees
    import os
    #gramPath = '../SourceCode/debrisGRAM_coupling/'

    latRange = np.arange(minLat,maxLat+ds,ds)
    lonRange = np.arange(minLon,maxLon+ds,ds)

    counter = 0
    for indexlat in range(0,len(latRange)):
        for indexlon in range(0,len(lonRange)):
            lat = latRange[indexlat]
            lon = lonRange[indexlon]
            os.system(gramPath+'simpleGram '+str(lat)+' '+str(lon)+' '+str(alt)+' '+str(deltaalt))
            os.system('mv special.txt '+'lat'+str(counter))
            counter = counter + 1

    

