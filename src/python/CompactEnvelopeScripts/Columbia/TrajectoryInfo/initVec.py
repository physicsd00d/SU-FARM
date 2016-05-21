# Initial State Vector for Columbia Case
import numpy as np
from FriscoLegacy import orbitTools

# 8:59:37 -- TJC: the time when the shuttle lost hydraulic pressure and would have spun out of control

# === CAIB says ====
# CAIB V2 D16
# 2.5 DEVELOPMENT OF BREAKUP STATE VECTORS AND ASSOCIATED DEBRIS GROUPS
# The breakup state vectors are based on a progressive
# breakup model that initiated at 13:59:30 GMT and spans
# 120 seconds. 24 debris lists are created for each five-second
# interval. 5 state vectors are used to cover each five-second 150 interval. Each state vector was assigned a failure probabil-
# ity of 0.2. In other words, each debris group is distributed evenly over a five-second span.
# To populate the debris list, a trial run was used to determine the relationship of downrange distance and
# failure time for each ballistic coefficient, beta class (see Figure 2-8). In this study,
# the Loss of Signal (LOS) point at -99.0413E, 32.956N is used as the point of origin.
# Next, the downrange im- pact distance for each fragment is measured from the same point of origin.

# === SpaceFlightNow says ====
# http://spaceflightnow.com/shuttle/sts107/timeline/
#
#  08:59:32 a.m.
# Columbia is approaching Dallas, Texas. Approximate location at initial loss of signal. H=200,700; Mach: 18.1 (+32.9-99.0)
#
#  08:59:37 a.m.
# Begin hypothetical trajectory plot for an object with a ballistic number = 220 psf. This trajectory is projected to compute an approximate impact point.
#
#  09:03:34 a.m.
# Reference trajectory ground impact for hypothetical object with a ballistic number = 220 psf. (+30.781-92.557)





# ACTA estimates from Probabilistic Debris Impact Modeling for Public Risk Analysis
# alt = 61 km
# Earth Rel Vel = 5.4 km/sec
# -.46 deg
# No WINDS and standard 1976 atmosphere

# # This is what Francisco had as his parameters
# planetModel = 0
# latitude = 32.95608
# longitude = -99.04132
# altitude = 61193.781
# Vmag = 5396.788
# FPA = -0.46              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed
# # FPA = -1.5              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed
# beta = -19.3            # angle from East vector
# thetagRef = 0.0
# inertial = 0
# loverd = 0.04
# cloption = 1      # if 1, then ignore loverd


# This is getting close to matching group 1...
# # NO WIND.  Tweaking to bring into better alignment with my better-calibrated debris field image.
# # Actually, should really get the debris model to be more dispersed first.  Come back to this...
# planetModel = 0
# # latitude = 32.95608   - 0.08
# # longitude = -99.04132 + 0.5
# latitude = 32.95608   
# longitude = -99.04132 + 0.25
# altitude = 61193.781
# Vmag = 5396.788
# FPA = -0.46     -1.       # flight path angle -1.5 when assuming winds. -.46 if no winds assumed
# # FPA = -1.5              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed
# beta = -19.3   -0.5         # angle from East vector
# thetagRef = 0.0
# inertial = 0
# loverd = 0.15
# cloption = 0      # if 1, then ignore loverd


# NO WIND.  Tweaking to bring into better alignment with my better-calibrated debris field image.
# Actually, should really get the debris model to be more dispersed first.  Come back to this...
planetModel = 0
# latitude = 32.95608   - 0.08
# longitude = -99.04132 + 0.5
latitude = 32.95608   
longitude = -99.04132 + 0.25
altitude = 61193.781
Vmag = 5396.788
FPA = -0.46     -1.5        # flight path angle -1.5 when assuming winds. -.46 if no winds assumed
# FPA = -1.5              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed
beta = -19.3   -0.5         # angle from East vector
thetagRef = 0.0
inertial = 0
loverd = 0.15
cloption = 0      # if 1, then ignore loverd


# ===== THESE WORK REASONABLY WELL.  A GOOD START BUT ROOM FOR IMPROVEMENT =====
# # These are my updated parameters, needed at a minimum to match the updated map coordinates
# planetModel = 0
# latitude = 32.95608     -0.25
# longitude = -99.04132   + 1.
# altitude = 61193.781
#
# # Velocity Magnitude.  Increasing this primarily moves the toe forward and only slightly moves the heel
# Vmag = 5396.788
#
# # Flight Path Angle tilts the state vec up and down wrt to horizontal.
# #   All else being equal, this only really changes the toe-edge of the debris pattern, so it's kind of a length parameter
# #   More negative points the state vec more down and decreases the length
# FPA = -3.5              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed [deg]
# beta = -19.3    -1.            # angle from East vector
# thetagRef = 0.0
# inertial = 0
# loverd = 0.04
# cloption = 1      # if 1, then ignore loverd


# #### ==== This has great agreement!!! ==== #### going to alter to change loverd
# # These are my updated parameters, needed at a minimum to match the updated map coordinates
# planetmodel = 0
# latitude = 32.95608   - 0.35
# longitude = -99.04132 + 1.25
# altitude = 61193.781
#
# # Velocity Magnitude.  Increasing this primarily moves the toe forward and only slightly moves the heel
# Vmag = 5396.788
#
# # Flight Path Angle tilts the state vec up and down wrt to horizontal.
# #   All else being equal, this only really changes the toe-edge of the debris pattern, so it's kind of a length parameter
# #   More negative points the state vec more down and decreases the length
# FPA = -3.5              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed [deg]
# beta = -19.3    -0.5            # angle from East vector
# thetagRef = 0.0
# inertial = 0
# loverd = 0.04
# cloption = 1      # if 1, then ignore loverd and use CL profile.  Doesnt change drag profile usage


# #### ==== This has great agreement!!! ==== #### going to alter to change loverd
# # These are my updated parameters, needed at a minimum to match the updated map coordinates
# planetModel = 1
# latitude = 32.95608   - 0.35 + 0.06
# longitude = -99.04132 + 1.25
# altitude = 61193.781    #meters
#
# # Velocity Magnitude.  Increasing this primarily moves the toe forward and only slightly moves the heel
# Vmag = 5396.788
#
# # Flight Path Angle tilts the state vec up and down wrt to horizontal.
# #   All else being equal, this only really changes the toe-edge of the debris pattern, so it's kind of a length parameter
# #   More negative points the state vec more down and decreases the length
# FPA = -3.5              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed [deg]
# beta = -19.3    -0.5  -0.8       # angle from East vector
# thetagRef = 0.0
# inertial = 0
# loverd = 0.15
# cloption = 0      # if 1, then ignore loverd and use CL profile.  Doesnt change drag profile usage



# #### ==== This has great agreement!!! ==== #### going to alter to change loverd
# # These are my updated parameters, needed at a minimum to match the updated map coordinates
# planetModel = 1
# latitude = 32.95608   + 0.20
# longitude = -99.04132 - 0.5
# altitude = 61193.781    #meters

# # Velocity Magnitude.  Increasing this primarily moves the toe forward and only slightly moves the heel
# Vmag = 5396.788

# # Flight Path Angle tilts the state vec up and down wrt to horizontal.
# #   All else being equal, this only really changes the toe-edge of the debris pattern, so it's kind of a length parameter
# #   More negative points the state vec more down and decreases the length
# FPA = -1.              # flight path angle -1.5 when assuming winds. -.46 if no winds assumed [deg]
# beta = -19.3          # angle from East vector
# thetagRef = 0.0
# inertial = 0
# loverd = 0.15
# cloption = 0      # if 1, then ignore loverd and use CL profile.  Doesnt change drag profile usage




Vup = Vmag*np.sin(FPA*np.pi/180.)
Vplane = Vmag*np.cos(FPA*np.pi/180.)
Veast = Vplane*np.cos(beta*np.pi/180.)
Vsouth = -Vplane*np.sin(beta*np.pi/180.)


r_eci = orbitTools.latlonalt2ECI(latitude,longitude,altitude,thetagRef,planetModel)
v_eci = orbitTools.SEZ2ECI(latitude,longitude,altitude,Vsouth,Veast,Vup,inertial,thetagRef,planetModel)

#print r_eci
#print v_eci
