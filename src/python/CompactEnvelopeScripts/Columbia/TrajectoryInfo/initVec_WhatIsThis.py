# Initial State Vector for Columbia Case
import numpy as np
import sys

# Points to the python scripts needed from Francisco
friscoFiles = '../../../../Prop3Dof/FriscoDebris/pythonFiles/'
sys.path.append(friscoFiles)
import orbitTools


#### ==== This has great agreement!!! ==== #### going to alter to change loverd
# These are my updated parameters, needed at a minimum to match the updated map coordinates
planetModel = 0
latitude = 32.95608   + 0.20 + 13. + 0.428669 + 0.003204 + 0.000024 - 0.000002
longitude = -99.04132 - 0.5 - 10. + 0.775222 - 1.144935
altitude = 181193.781    #meters

# Velocity Magnitude.  Increasing this primarily moves the toe forward and only slightly moves the heel
Vmag = 7096.788

# Flight Path Angle tilts the state vec up and down wrt to horizontal.
#   All else being equal, this only really changes the toe-edge of the debris pattern, so it's kind of a length parameter
#   More negative points the state vec more down and decreases the length
FPA = -7.              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed [deg]
beta = -50          # angle from East vector
thetagRef = 0.0
inertial = 0
loverd = 0.3        # Apollo capsule around 0.3
cloption = 0      # if 1, then ignore loverd and use CL profile.  Doesnt change drag profile usage


















Vup = Vmag*np.sin(FPA*np.pi/180.)
Vplane = Vmag*np.cos(FPA*np.pi/180.)
Veast = Vplane*np.cos(beta*np.pi/180.)
Vsouth = -Vplane*np.sin(beta*np.pi/180.)


r_eci = orbitTools.latlonalt2ECI(latitude,longitude,altitude,thetagRef,planetModel)
v_eci = orbitTools.SEZ2ECI(latitude,longitude,altitude,Vsouth,Veast,Vup,inertial,thetagRef,planetModel)

#print r_eci
#print v_eci
