from collections import namedtuple
from LaunchSites import siteDict
# from LaunchProviders import footprintDict
# Mission = namedtuple('Mission', ['launchDate', 'nominalTime', 'timeZone', 
#                                              'vehicle', 'launchSite', 'azimuth', 'delayAllowed'], verbose=False)

Mission = namedtuple('Mission', ['launchDateTime', 'timeZone', 
                                             'vehicle', 'launchSite', 'azimuth', 'delayAllowed'], verbose=False)

Footprint = namedtuple('Footprint', ['datName'], verbose=False)



def getLatLon(curMission):
    locationName = curMission.launchSite
    try:
        lat = siteDict[locationName]['lat']
        lon = siteDict[locationName]['lon']
    except:
        print "Your launch location " + locationName + " is not recognized"
        lat = 361
        lon = 361
        raise RuntimeError
    return(lat,lon)




# def getMostSimilarFootprint(curMission):
#     datName = footprintDict[curMission.vehicle]
#     return datName

def convertTimeToGMT(curMission):
    curTime = curMission.nominalTime
    print curTime