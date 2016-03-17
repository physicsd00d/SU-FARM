from initVec import * # running initVec...getting r_eci,v_eci. Check the trajectory folder to see what is going on here
import sys

# print "GenerateMainPieceTrajectory hasnt yet been adapted to SUFARM.  Exiting..."
# sys.exit()

# # Points to the python scripts needed from Francisco
# friscoFiles = '../../../../Prop3Dof/FriscoDebris/pythonFiles/'
# # Points to the binaries for propagating trajectories
# debrisPropPATH = '../../../../Prop3Dof/FriscoDebris/'
# sys.path.append(friscoFiles)
# sys.path.append(debrisPropPATH)

from FriscoLegacy import debrisPropagation as dp
import numpy as np

# To read in...
import pickle

def Generate(curMission):
	print "YO"
	# [altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist] = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))

	profiles = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))
	atmStorage = profiles['atmStorage']
	windIX = 0
	altitudeList	= atmStorage[windIX][0]
	densityMeanList	= atmStorage[windIX][1]
	uMeanList	= atmStorage[windIX][2]
	vMeanList	= atmStorage[windIX][3]
	wMeanList	= atmStorage[windIX][4]

	VdebInit = [0,0,0] # impulse velocity for debris in ECI frame

	r0V0 = np.array([r_eci[0],r_eci[1],r_eci[2],v_eci[0],v_eci[1],v_eci[2]]) # vehicle initial state vector.

	# setting parameters for propagation
	massCol = 362.8736
	ArefCol = 0.412861083336890
	cdCol = [0.6]
	clCol = [0.04*cdCol[0]]
	minfcd = [1]
	minfcl = [1]

	dt = 1 # time interval for trajectory calculation

	# ignoring the wind profile
	uListCol = 0.0*np.array(uMeanList)
	vListCol = 0.0*np.array(vMeanList)
	wListCol = 0.0*np.array(wMeanList)


	# propagate main piece of debris all the way to the ground. Check RSAT/source/TrajectoryPropagation/sourceCode/debrisPropAllTime.f90 for more info
	debrisResults, numFinalSteps = dp.debrispropagation(initialstate=r0V0,debrisvel=VdebInit,mass=massCol,sref=ArefCol,minfcd=1,cd=cdCol,cloption=cloption,minfcl=1,cl=clCol,loverd=loverd,atmosoption=2,altitudelist=altitudeList,densitylist=densityMeanList,ulist=uListCol,vlist=vListCol,wlist=wListCol,geoptions=0,filename='colTraj',planetmodel=planetModel,dtinterval=dt,ndtinterval=1000000,thetag0=thetagRef,ncd=len(minfcd),ncl=len(minfcl),nlist=len(altitudeList))

	sizeResults = 7
	sizeKeeping = 6
	debrisResults = debrisResults[0:(numFinalSteps)*sizeResults].reshape(numFinalSteps,sizeResults)[:,:sizeKeeping]

	print '==================== GENERATING A NEW MAIN TRAJECTORY ===================='
	fo = open("ColumbiaMainPieceTrajectory.txt", "w")
	outString = "FORMAT LATLON METERS"
	# print str(outString)
	fo.write( outString + "\n")

	outString = "{0:6}{1:20}{2:20}{3:20}".format("#Time", "LonDeg", "LatDeg","AltM")
	# print str(outString)
	fo.write( outString + "\n")

	for ix in range(len(debrisResults)):
	    outString = "{0:6}{1:20}{2:20}{3:20}".format(str(ix), str(debrisResults[ix][1]), str(debrisResults[ix][0]), str(debrisResults[ix][2]))
	    # print str(outString)
	    fo.write( outString + "\n")
	fo.close()


# sys.exit()


# # Do the propagation but get state vectors in return
# from FriscoLegacy import orbitProp as op # trajectory propagation routine

# propOption = 4  # Propagate until altitude drops below threshold
# propCond = 0    # This is that threshold
# cd=cdCol
# # cloption=1
# cl=clCol
# loverd=0.04
# atmosoption=2
# Tmat = np.zeros((1,3))
# timelist = np.array([1])
# ntime = 1
# isp =1  # zero will blow things up
# geoptions=0
# filename='colTraj'
# # planetmodel=0
# dtinterval=dt
# ndtinterval=1000000
# thetag0=thetagRef
# ThrustOffsetAngDeg = np.zeros((1,2))

# #finalConditions(index1,:) = (/y(1),y(2),y(3),y(5),y(6),y(7),latlonalt(2),latlonalt(1),latlonalt(3),y(4),tout,thetag + tout*omega/)
# finalconditions,finalderivs = op.propagate(r0V0, massCol, propCond,ArefCol,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudeList,densityMeanList,uListCol,vListCol ,wListCol, Tmat,timelist,isp,geoptions,filename,planetModel,dtinterval,ndtinterval,thetag0,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudeList))

# newIndex = finalconditions[:,8]>=0.0

# xvec        = finalconditions[newIndex,0]
# yvec        = finalconditions[newIndex,1]
# zvec        = finalconditions[newIndex,2]
# Vxvec       = finalconditions[newIndex,3]
# Vyvec       = finalconditions[newIndex,4]
# Vzvec       = finalconditions[newIndex,5]
# thetavec    = finalconditions[newIndex,-1]
# tvec        = finalconditions[newIndex,-2]


# fo = open("ColumbiaMainPieceTrajectoryStateVec.txt", "w")

# outString = "FORMAT STATEVEC"
# # print str(outString)
# fo.write( outString + "\n")

# outString = "{0:10}{1:20}{2:20}{3:20}{4:20}{5:20}{6:20}{7:20}".format("#Time", "thetag", "x","y","z","Vx","Vy","Vz")
# # print str(outString)
# fo.write( outString + "\n")

# for ix in range(len(tvec)):

#     outString = "{0:<10}{1:<20}{2:<20}{3:<20}{4:<20}{5:<20}{6:<20}{7:<20}"\
#         .format(tvec[ix], thetavec[ix], xvec[ix], yvec[ix], zvec[ix], Vxvec[ix], Vyvec[ix], Vzvec[ix])

#     # outString = "{0:6}{1:20}{2:20}{3:20}".format(str(ix), str(debrisResults[ix][1]), str(debrisResults[ix][0]), str(debrisResults[ix][2]))
#     # print str(outString)
#     fo.write( outString + "\n")
# fo.close()





