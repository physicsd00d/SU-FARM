# Read in a scenario file and pump out the FACET folders
import sys
import os

# missionFolder       = "2018Low_SingleDay"
# missionFolder       = "2018Med_SingleDay"
# missionFolder       = "2018High_SingleDay"
# missionFolder       = "2025Low_SingleDay"
# missionFolder       = "2025Med_SingleDay"
missionFolder       = "2025High_SingleDay"

FacetFolder         = "FacetFiles/"
# scenarioFileFolder = "ScenarioFiles/SingleDays/"
missionList = os.listdir(FacetFolder + missionFolder + "/")

aggMissionFileName = "{0}SUA_{1}".format(FacetFolder, missionFolder)
aggMissionFile = open(aggMissionFileName, 'w')


print missionList

missionCounter = 0
for curMission in missionList:
    if curMission[0] == '.':
        continue

    curName =  "_".join(curMission.split("_")[:2]) + "_{0}_".format(missionCounter)
    missionCounter = missionCounter + 1

    curFile = open(FacetFolder + missionFolder + "/" + curMission + "/SUA_GeneratedSUA", 'r')
    for line in curFile:
        if "SUA" in line:
            line = line.replace("SUA", curName)

        aggMissionFile.write(line)
    curFile.close()

aggMissionFile.close()

