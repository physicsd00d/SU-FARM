# This will get passed in to the C++ functions for individual risk calculations

# Concerning the distinction between business jets and commercial aircraft
# (Aircraft Protection for Range Safety Wilde and Draper AIAA-2010-1542-986.pdf)
# The CT aircraft results presented here apply only to commercial transport aircraft class with all the following characteristics:
#   1) Aluminum skin (composite skin aircraft have not been studied),
#   2) Multiple turbofan engines, and
#   3) Governed by the FAA certification requirements of 14 CFR Part 23/25.
# The BJ Class includes all multi-engine, jet propelled aircraft that have the capability to carry
#   no more than 20 passengers for hire. All aircraft within the BJ class primarily exhibit:
#   1) Aluminum skin and structural members,
#   2) Two pilots during operation, and
#   3) Design and maintenance requirements defined by the FAA certification requirements of 14 CFR Part 23/25.
# The BJ class excludes:
#   4) Single pilot versions of otherwise BJ class aircraft,
#   5) Emerging "very light jets" with composite skins or structures, or
#   6) Aircraft that rely on propeller-based propulsion.


AircraftAreaModel = dict()

top2front       = 0.145
m2km            = 1e-3

business        = 0
commercial      = 1

# Starting with just the top few Columbia risky planes
# AircraftFromPaper = [

# 'COA282',
# 'SWA333',
acType      = 'B733'
length      = 33.4      #m
wingArea    = 28.88     #m^2
width       = 3.76      #m
passengers  = 146       #
acClass     = commercial

topArea     = (wingArea + length*width)*(m2km**2)    #km
frontArea   = topArea * top2front               #km

if AircraftAreaModel.has_key(acType):
    print 'ERROR: Youre trying to add the same acType multiple times'
    raise RuntimeError
else:
    AircraftAreaModel[acType] = dict(topArea=topArea, frontArea=frontArea, acClass=acClass)


# 'CAA916',
acType      = 'E120'
length      = 20        #m
wingArea    = 39.4      #m^2
width       = 2.28      #m
passengers  = 30
acClass     = commercial

topArea     = (wingArea + length*width)*(m2km**2)    #km
frontArea   = topArea * top2front               #km

if AircraftAreaModel.has_key(acType):
    print 'ERROR: Youre trying to add the same acType multiple times'
    raise RuntimeError
else:
    AircraftAreaModel[acType] = dict(topArea=topArea, frontArea=frontArea, acClass=acClass)



# 'DAL1055',
acType      = 'MD90'    #
length      = 46.5      #m
wingArea    = 112.3     #m^2    I'm basing this off the MD80-series since they're basically the same size
width       = 3.35      #m      I'm basing this off the MD80-series since they're basically the same size
passengers  = 160
acClass     = commercial

topArea     = (wingArea + length*width)*(m2km**2)    #km
frontArea   = topArea * top2front               #km

if AircraftAreaModel.has_key(acType):
    print 'ERROR: Youre trying to add the same acType multiple times'
    raise RuntimeError
else:
    AircraftAreaModel[acType] = dict(topArea=topArea, frontArea=frontArea, acClass=acClass)


# 'COA1710',
acType      = 'MD82'    #
AircraftAreaModel[acType] = AircraftAreaModel['MD90']


# 'SKW3752',
# 'SKW3825',
acType      = 'CRJ2'
length      = 26.77     #m
wingArea    = 48.35     #m^2
width       = 2.69      #m
passengers  = 50
acClass     = commercial

topArea     = (wingArea + length*width)*(m2km**2)    #km
frontArea   = topArea * top2front               #km

if AircraftAreaModel.has_key(acType):
    print 'ERROR: Youre trying to add the same acType multiple times'
    raise RuntimeError
else:
    AircraftAreaModel[acType] = dict(topArea=topArea, frontArea=frontArea, acClass=acClass)


# 'COA688',
acType      = 'B735'    #
AircraftAreaModel[acType] = AircraftAreaModel['B733']

# 'DAL2137']
acType      = 'B738'    #
AircraftAreaModel[acType] = AircraftAreaModel['B733']




# # Empty Template
# acType      = ''
# length      = #m
# wingArea    = #m^2
# width       = #m
# passengers  =
# acClass     =
#
# topArea     = (wingArea + length*width)*(m2km**2)    #km
# frontArea   = topArea * top2front               #km
#
# if AircraftAreaModel.has_key(acType):
#     print 'ERROR: Youre trying to add the same acType multiple times'
#     raise RuntimeError
# else:
#     AircraftAreaModel[acType] = dict(topArea=topArea, frontArea=frontArea, acClass=acClass)