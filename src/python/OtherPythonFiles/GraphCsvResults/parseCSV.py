import os
import sys

# Find the path of the current file, then Point to the root of the package so I can run this script from anywhere
curFilePath = os.path.dirname(os.path.abspath(__file__)) + "/"
rootDir =   os.path.abspath(curFilePath + "../../../../") + "/"
outputDir = rootDir + "outputs/" # Where to store results, gitignored
# tempDir =   rootDir + "temp/"   # temp files here, gitignored
csvTradFolder = outputDir + "space2015CsvFiles/"
csvEnvFolder = outputDir + "fullsafeCsvFiles/"
outputPath = outputDir + "fullsafeFigures/"

import glob
import numpy as np
import matplotlib.pyplot as plt

# os.chdir("/Users/marian/Documents/PersonalFinance/aiaaSpace")

# Define some global variables to be used as options
useMean = 1
useMax = 2


csvFolderTrad   = [csvTradFolder + "01_Jan/Traditional/",
                   csvTradFolder + "02_Feb/Traditional/",
                   csvTradFolder + "03_Mar/Traditional/"]
  
csvFolderEnv    = [csvEnvFolder]

# csvFolderEnv    = [csvEnvFolder + "01_Jan/Envelopes/",
#                    csvEnvFolder + "02_Feb/Envelopes/",
#                    csvEnvFolder + "03_Mar/Envelopes/"]

# outputPath="/Users/marian/Desktop/Research Papers/AIAA Space2015/Figures/"
# outputPath="/Users/marian/Desktop/Research Papers/AIAA Space2015/Figures/"

# sys.exit()

#%%
#def loadMissionDict(csvFolder):
#    fileList = glob.glob(csvFolder + "*.csv")   # Get all the csv files
#    fileList = [obj.split('/')[-1] for obj in fileList]          # Strip out the path
#    fileList
#    
#    # Now match them all together by mission and outputName
#    missionDict = dict()
#    for curFile in fileList:
#        # Strip out the info from the filename
#        vehicle, launchLoc, filterLoc, dateStr, scenario = curFile[:-4].split('_')
#        assert launchLoc in filterLoc    
#        outfileName = "{0}_{1}_{2}".format(vehicle, launchLoc, scenario)
#        #print outfileName
#            
#        # Open the file
#        inFile = open(csvFolder + curFile, "r")
#        allLines = inFile.readlines()
#        
#        # Get the names of the quantities of interest, strip endlines and leading hash
#        # All the data we want to aggregate is on the second line of each file.
#        dataKeys = [obj.strip("#").strip() for obj in allLines[0].split(',')]  # First Line of csv
#        dataVals = [float(obj) for obj in allLines[1].strip().split(',')]      # Second line
#        #print dataKeys
#        #print dataVals
#        inFile.close()
#        
#        # Save the results in a dictionary.
#        if not missionDict.has_key(outfileName):
#            missionDict[outfileName] = dict(csvFiles=[], dataVals=[], dataKeys=dataKeys)
#        missionDict[outfileName]["csvFiles"].append(curFile)
#        missionDict[outfileName]["dataVals"].append(dataVals)
#        
#    return missionDict
    
def loadMissionDict(csvFolderList):
    missionDict = dict()

    for csvFolder in csvFolderList:
        fileList = glob.glob(csvFolder + "*.csv")   # Get all the csv files
        fileList = [obj.split('/')[-1] for obj in fileList]          # Strip out the path
        fileList
        
        # Now match them all together by mission and outputName
        for curFile in fileList:
            # Strip out the info from the filename
            vehicle, launchLoc, filterLoc, dateStr, scenario = curFile[:-4].split('_')
            assert launchLoc in filterLoc    
            outfileName = "{0}_{1}_{2}".format(vehicle, launchLoc, scenario)
            #print outfileName
                
            # Open the file
            inFile = open(csvFolder + curFile, "r")
            allLines = inFile.readlines()
            
            # Get the names of the quantities of interest, strip endlines and leading hash
            # All the data we want to aggregate is on the second line of each file.
            dataKeys = [obj.strip("#").strip() for obj in allLines[0].split(',')]  # First Line of csv
            dataVals = [float(obj) for obj in allLines[1].strip().split(',')]      # Second line
            #print dataKeys
            #print dataVals
            inFile.close()
            
            # Save the results in a dictionary.
            if not missionDict.has_key(outfileName):
                missionDict[outfileName] = dict(csvFiles=[], dataVals=[], dataKeys=dataKeys)
            missionDict[outfileName]["csvFiles"].append(csvFolder + curFile)
            missionDict[outfileName]["dataVals"].append(dataVals)
        
    return missionDict

def ProcessAndUpdate(missionDict):
    """This will add fields to the incoming missionDict"""
    for curMission in missionDict:
        curValMat = np.array(missionDict[curMission]["dataVals"])
        curValDict = dict(zip(missionDict[curMission]["dataKeys"], curValMat.T))
                
        sampleSize = len(curValMat)  # Should equal the number of days
        missionDict[curMission]["sampleSize"] = sampleSize

        print "{0}:".format(curMission)
        missionDict[curMission]["average"] = dict()
        missionDict[curMission]["maximum"] = dict()
        missionDict[curMission]["minimum"] = dict()
        missionDict[curMission]["sampleVar"] = dict()
        missionDict[curMission]["confidence"] = dict()
        for curKey, curValVec in curValDict.items():
            print "   {0} ==> {1}".format(curKey, np.mean(curValVec))
            mu = np.mean(curValVec)
            s2 = sum((mu -curValVec)**2)/(sampleSize-1.)
            missionDict[curMission]["average"][curKey] = mu
            missionDict[curMission]["maximum"][curKey] = np.max(curValVec)
            missionDict[curMission]["minimum"][curKey] = np.min(curValVec)
            missionDict[curMission]["sampleVar"][curKey] = s2
            missionDict[curMission]["confidence"][curKey] =  1.96*np.sqrt(s2/sampleSize)
            # WARNING: You can only use 1.96 for 95% if you're sure the means are 
            #  normally distributed!  In general, should use a t-table. 
            # TODO: CHANGE THIS!

def CollapseIntoScenaros(missionDict, missionsInGroup):
    """Collapse the individual missions into the three day-scenarios"""
    scenarioDict = dict()
    for curGroup in missionsInGroup:
        scenarioDict[curGroup] = dict(dataVals=[], csvFiles=[], dataKeys=[])
        missionList = missionsInGroup[curGroup]
        for curMission in missionList:
#            print "{0} : {1}".format(curGroup, curMission)
            subDict = missionDict[curMission]
#            print np.array(subDict["dataVals"])
#            print curMission
            # Add the mission togetherfor this scenario
            if len(scenarioDict[curGroup]["dataVals"]) == 0:
                scenarioDict[curGroup]["dataVals"] = np.array(subDict["dataVals"])
#                print "equals"
            else:
                scenarioDict[curGroup]["dataVals"] += np.array(subDict["dataVals"])
#                print "adds"
                
#            scenarioDict[curGroup]["dataVals"].extend(subDict["dataVals"])
    #        scenarioDict[curGroup]["csvFiles"].append(subDict["csvFiles"])
        scenarioDict[curGroup]["dataKeys"] = subDict["dataKeys"]
    return scenarioDict
    
    
#def PlotNumRerouted(tradDict, envDict, keyToTickMap, topToBottomList, xTrad, xEnv, figsize, outputPath, outputFile):
#    missionsToPlot = list(reversed(topToBottomList))
#    missionsPrettyTicks = [keyToTickMap[name] for name in missionsToPlot]
#    numMissions = len(missionsToPlot)
#    index =  np.arange(numMissions)
#    
#    # Start the plot
#    #plt.figure(1)
#    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=figsize)
#    bar_width = 0.35
#    opacity = 0.5
#    #error_config = {'ecolor': '0.3', 'capsize' : '1'}
#    error_config = {'ecolor': '0.3', 'linewidth' : '3', 'capthick':3}
#    
#    ax1.set_xlabel('Number of Aircraft', fontsize=20)
#    ax1.set_title('Traditional', fontsize=20)
#    ax2.set_xlabel('Number of Aircraft', fontsize=20)
#    ax2.set_title('Compact Envelopes', fontsize=20)
#    plt.yticks(index + bar_width, missionsPrettyTicks)
#    
#    # Do the Traditional boxes first
#    reoutedMeans = [tradDict[name]["average"]["rerouted"] for name in missionsToPlot]
#    reoutedConf = [tradDict[name]["confidence"]["rerouted"] for name in missionsToPlot]
#    numFilterMeans = [tradDict[name]["average"]["numFilter"] for name in missionsToPlot]
#    numFilterConf = [tradDict[name]["confidence"]["numFilter"] for name in missionsToPlot]
#    #sampleSize = tradDict[curMission]["sampleSize"]
#    eps = 0.04
#    
#    print reoutedConf
#    rectsTradRerouted = ax1.barh(index+bar_width+eps, reoutedMeans, bar_width,
#    #                 linewidth=2,  # Looks crappy with opacity not 1
#                     capsize=7,
#                     alpha=opacity,
#                     color='r',
#                     xerr=reoutedConf,
#                     error_kw=error_config,
#                     label='numRerouted')
#                     
#    rectsTradFiltered = ax1.barh(index, numFilterMeans, bar_width,
#    #                 linewidth=2,
#                     alpha=opacity,
#                     color='b',
#                     xerr=numFilterConf,
#                     error_kw=error_config,
#                     label='numFiltered')
#    
#    # Now the Envelopes
#    reoutedMeans = [envDict[name]["average"]["rerouted"] for name in missionsToPlot]
#    reoutedConf = [envDict[name]["confidence"]["rerouted"] for name in missionsToPlot]
#    numFilterMeans = [envDict[name]["average"]["numFilter"] for name in missionsToPlot]
#    numFilterConf = [envDict[name]["confidence"]["numFilter"] for name in missionsToPlot]
#    
#    rectsEnvRerouted = ax2.barh(index+bar_width+eps, reoutedMeans, bar_width,
#                     alpha=opacity,
#                     color='r',
#                     xerr=reoutedConf,
#                     error_kw=error_config,
#                     label='numRerouted')
#                     
#    rectsEnvFiltered = ax2.barh(index, numFilterMeans, bar_width,
#                     alpha=opacity,
#                     color='b',
#                     xerr=numFilterConf,
#                     error_kw=error_config,
#                     label='numFiltered')
#    
#    ax1.set_xbound([0,xTrad])
#    ax2.set_xbound([0,xEnv])
#    #ax1.set_ybound(([0,40]))
#    #ax2.ticklabel_format(size=400)
#    #ax1.tick_params(size=400)
#    
#    plt.setp(ax1.get_yticklabels(), fontweight=800)  # Scenario name font size
#    plt.setp(ax1.get_yticklabels(), fontsize=20)  # Scenario name font size
#                                   # My understanding is that when shared, only first axis is shown, all others set to invisible.
#    plt.setp(ax1.get_xticklabels(), fontsize=20)  # x-axis number size
#    plt.setp(ax2.get_xticklabels(), fontsize=20)  # x-axis number size
#    
#    
#    # NOTE: ticklable objects are just Text objects
#    
#    
#    lg = plt.legend(["Rerouted", "Filtered"])
#    plt.setp(lg.get_texts(), fontsize=20)  # x-axis number size
#    
##    plt.setp(ax1.get_yticklabels(), horizontalalignment='center')
##    ax1.tick_params(axis='y', pad=110) # Moves the ylabels back a litle bit from the edge of the graph.
#    
##    plt.setp(ax1.get_yticklabels(), horizontalalignment='left')
##    ax1.tick_params(axis='y', pad=200) # Moves the ylabels back a litle bit from the edge of the graph.
#    
#    plt.setp(ax1.get_yticklabels(), horizontalalignment='right')
#    ax1.tick_params(axis='y', pad=15) # Moves the ylabels back a litle bit from the edge of the graph.
#        
#    
#    
#    #ax1.yaxis.labelpad = 100
#    
#    plt.tight_layout()
#    #plt.show()
#    
#    #plt.savefig(outputPath+"testOneMonth", dpi=600)
#    plt.savefig(outputPath+outputFile, format='pdf')
#    
#    return ax1, ax2, f


useAll = False

if useAll:
    # Make a dictionary that maps the scenario name into the tick label name I want
    sc2tl = dict()
    sc2tl[ 'SS2_America_2018L'] = "SS2 America \n9:00"
    sc2tl[ 'Pegasus_Wallops_2018L'] = "Pegasus Wallops \n11:20"
    sc2tl[ 'Lynx_Midland_2018L'] = "Lynx Midland \n14:15"

    sc2tl[ 'Atlas5_Vafb_2018M'] = "AtlasV VAFB \n5:30"
    sc2tl[ 'Antares_Wallops_2018M'] = "Antares Wallops \n10:00"
    sc2tl[ 'SS2_America_2018M'] = "SS2 America \n11:45"
    sc2tl[ 'Reentry_PacificReentryLA_2018M'] = "Dragon Pacific \n15:23"

    sc2tl[ 'SS2_Titus_2018H'] = "SS2 Titus \n6:30"
    sc2tl[ 'PM_CornNew_2018H'] = "PM2 Corn \n8:20"
    sc2tl['Lynx_Cecil_2018H'] = "Lynx Cecil \n11:10"
    sc2tl[ 'SS2_America_2018H'] = "SS2 Amerca \n11:45"  
    sc2tl[ 'Sound_America_2018H'] = "Sound America \n12:00"
    sc2tl[ 'Lynx_Midland_2018H'] = "Lynx Midland \n13:50"
    sc2tl[ 'Lynx_FrntRnge_2018H'] = "Lynx FrontRange \n16:45"

    topToBottom2018L = ["SS2_America_2018L",
                         "Pegasus_Wallops_2018L",
                         "Lynx_Midland_2018L"]

    topToBottom2018M = ["Atlas5_Vafb_2018M",
                         "Antares_Wallops_2018M",
                         "SS2_America_2018M",
                         "Reentry_PacificReentryLA_2018M"]

    topToBottom2018H = ["SS2_Titus_2018H",
                         "PM_CornNew_2018H",
                         "Lynx_Cecil_2018H",
    #                     "SS2_America_2018H",
                         "Sound_America_2018H",
                         "Lynx_Midland_2018H",
                         "Lynx_FrntRnge_2018H"]                     
                         
else:

    # For presentations
    sc2tl = dict()
    sc2tl[ 'SS2_America_2018L'] = "SS2 America \n9:00"
    sc2tl[ 'Atlas5_Vafb_2018M'] = "AtlasV VAFB \n5:30"
    sc2tl[ 'Antares_Wallops_2018M'] = "Antares Wallops \n10:00"
    sc2tl[ 'Lynx_FrntRnge_2018H'] = "Lynx FrontRange \n16:45"

    topToBottom2018L = ["SS2_America_2018L",
                         "Lynx_FrntRnge_2018H",
                         "Atlas5_Vafb_2018M",
                         "Antares_Wallops_2018M"]
    topToBottom2018M = []
    topToBottom2018H = []   




# Make a dictionary that maps the mission name into one of the three groups.
missionsInGroup = dict()
missionsInGroup["2018L"] = ["Lynx_Midland_2018L",
                            "Pegasus_Wallops_2018L",
                            "SS2_America_2018L"]

missionsInGroup["2018M"] = ["Reentry_PacificReentryLA_2018M",
                            "Atlas5_Vafb_2018M",
                            "SS2_America_2018M",
                            "Antares_Wallops_2018M"]

missionsInGroup["2018H"] = ["Lynx_Cecil_2018H",
                            "SS2_Titus_2018H",
                            "PM_CornNew_2018H",
                            "Lynx_Midland_2018H",
                            "Lynx_FrntRnge_2018H",
                            "Sound_America_2018H",
                            "SS2_America_2018H"]
scenarioToTickMap=dict()
scenarioToTickMap["2018L"] = "Low Traffic"
scenarioToTickMap["2018M"] = "Medium Traffic"
scenarioToTickMap["2018H"] = "High Traffic"
#%
   
#outputFile      = "testJan.pdf"
#csvFolderTrad   = "csvFiles/01_Jan/Traditional/"
#csvFolderEnv    = "csvFiles/01_Jan/Envelopes/"

#outputFile      = "testFeb.pdf"
#csvFolderTrad   = "csvFiles/02_Feb/Traditional/"
#csvFolderEnv    = "csvFiles/02_Feb/Envelopes/"

#outputFile      = "testMar.pdf"
#csvFolderTrad   = ["csvFiles/03_Mar/Traditional/"]
#csvFolderEnv    = ["csvFiles/03_Mar/Envelopes/"]

# csvFolderTrad   = ["csvFiles/01_Jan/Traditional/",
#                    "csvFiles/02_Feb/Traditional/",
#                    "csvFiles/03_Mar/Traditional/"]
                   
# csvFolderEnv    = ["csvFiles/01_Jan/Envelopes/",
#                    "csvFiles/02_Feb/Envelopes/",
#                    "csvFiles/03_Mar/Envelopes/"]



tradMissions = loadMissionDict(csvFolderTrad)
tradMissions.keys()

envMissions = loadMissionDict(csvFolderEnv)
envMissions.keys()


##%  Broken out by mission
#
# Now process it and take the average
ProcessAndUpdate(tradMissions)
ProcessAndUpdate(envMissions)
#
#outputPath="/Users/marian/Desktop/Research Papers/AIAA Space2015/Figures/"
#outputFile = "testAggregate.pdf"
#figsize=(10, 15)
#xTrad = 100
#xEnv = 2.5
#topToBottomIndividual = []
#topToBottomIndividual.extend(topToBottom2018L)
#topToBottomIndividual.extend(topToBottom2018M)
#topToBottomIndividual.extend(topToBottom2018H)
#ax1, ax2, f = PlotNumRerouted(tradMissions, envMissions, sc2tl, topToBottomIndividual, xTrad, xEnv, figsize, outputPath, outputFile)


#%%  Let's see how closely the other quantities track with numRerouted

def PlotTest(curDict, keyToTickMap, topToBottomList, plotDict, outputPath, outputFile):
    #Unpack
    xlimNumAffected = plotDict["xlimNumAffected"]  
    xlimTime = plotDict["xlimTime"]  
    xlimFuel = plotDict["xlimFuel"]  
    xlimDist = plotDict["xlimDist"]  
    figsize = plotDict["figsize"]  
    legendAnchor = plotDict["legendAnchor"]  
    useLog = plotDict["useLog"]  
    titlePrefix = plotDict["titlePrefix"]  
    plotOption = plotDict["plotOption"]  

    missionsToPlot = list(reversed(topToBottomList))
    missionsPrettyTicks = [keyToTickMap[name] for name in missionsToPlot]
    numMissions = len(missionsToPlot)
    index =  np.arange(numMissions)
    
    # Start the plot
    #plt.figure(1)
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize=figsize)
    bar_width = 0.35
    opacity = 0.5
    #error_config = {'ecolor': '0.3', 'capsize' : '1'}
    error_config = {'ecolor': '0.3', 'linewidth' : '3', 'capthick':3}
    

    ax1.set_title('{0} Aircraft Affected'.format(titlePrefix), fontsize=20, fontweight=800)
    ax2.set_title('{0} Time Disruption'.format(titlePrefix), fontsize=20, fontweight=800)
    ax3.set_title('{0} Fuel Disruption'.format(titlePrefix), fontsize=20, fontweight=800)
    ax4.set_title('{0} Distance Disruption'.format(titlePrefix), fontsize=20, fontweight=800)

    ax1.set_xlabel('Number of Aircraft', fontsize=20, fontweight=800)
    ax2.set_xlabel('Time [min]', fontsize=20, fontweight=800)
    ax3.set_xlabel('Fuel [lbs]', fontsize=20, fontweight=800)
    ax4.set_xlabel('Distance [n.mi]', fontsize=20, fontweight=800)
    plt.yticks(index + bar_width, missionsPrettyTicks)
    
    # Calculate the appropriate bars / errors
    if plotOption == useMean:
        # Number affected
        reroutedBar = [curDict[name]["average"]["rerouted"] for name in missionsToPlot]
        reroutedConf = [curDict[name]["confidence"]["rerouted"] for name in missionsToPlot]
        numFilterBar = [curDict[name]["average"]["numFilter"] for name in missionsToPlot]
        numFilterConf = [curDict[name]["confidence"]["numFilter"] for name in missionsToPlot]
        
        # Time Affected
        deltaTimeBar = [curDict[name]["average"]["DeltaTime"] for name in missionsToPlot]
        deltaTimeConf = [curDict[name]["confidence"]["DeltaTime"] for name in missionsToPlot]
        filterMinutesBar = [curDict[name]["average"]["sumFilterMinutes"] for name in missionsToPlot]
        filterMinutesConf = [curDict[name]["confidence"]["sumFilterMinutes"] for name in missionsToPlot]
    
        # Fuel Affected
        deltaFuelBar = [curDict[name]["average"]["DeltaFuel"] for name in missionsToPlot]
        deltaFuelConf = [curDict[name]["confidence"]["DeltaFuel"] for name in missionsToPlot]
        
        # Dist Affected
        deltaDistBar = [curDict[name]["average"]["DeltaDist"] for name in missionsToPlot]
        deltaDistConf = [curDict[name]["confidence"]["DeltaDist"] for name in missionsToPlot]
    elif plotOption == useMax:
        # Number affected
        reroutedBar = [curDict[name]["maximum"]["rerouted"] for name in missionsToPlot]
        reroutedConf = None
        numFilterBar = [curDict[name]["maximum"]["numFilter"] for name in missionsToPlot]
        numFilterConf = None
        
        # Time Affected
        deltaTimeBar = [curDict[name]["maximum"]["DeltaTime"] for name in missionsToPlot]
        deltaTimeConf = None
        filterMinutesBar = [curDict[name]["maximum"]["sumFilterMinutes"] for name in missionsToPlot]
        filterMinutesConf = None
    
        # Fuel Affected
        deltaFuelBar = [curDict[name]["maximum"]["DeltaFuel"] for name in missionsToPlot]
        deltaFuelConf = None
        
        # Dist Affected
        deltaDistBar = [curDict[name]["maximum"]["DeltaDist"] for name in missionsToPlot]
        deltaDistConf = None        
    else:
        raise NotImplementedError    
    
    #    # Set the number of x-ticks if not making a log plot
    if not useLog:
        ax1.locator_params(axis='x', nbins=4)
        ax2.locator_params(axis='x', nbins=4)
        ax3.locator_params(axis='x', nbins=4)
        ax4.locator_params(axis='x', nbins=4)
    else:
        # If you are using the log plot, don't let any values be negative
        reroutedBar = np.maximum(reroutedBar, 1e-9)
        numFilterBar = np.maximum(numFilterBar, 1e-9)
        deltaTimeBar = np.maximum(deltaTimeBar, 1e-9)
        filterMinutesBar = np.maximum(filterMinutesBar, 1e-9)
        deltaFuelBar = np.maximum(deltaFuelBar, 1e-9)
        deltaDistBar = np.maximum(deltaDistBar, 1e-9)
        
        
    # Do the Number affected first

    #sampleSize = tradDict[curMission]["sampleSize"]
    eps = 0.04
    
    rectsTradRerouted = ax1.barh(index+bar_width+eps, reroutedBar, bar_width,
    #                 linewidth=2,  # Looks crappy with opacity not 1
                     capsize=7,
                     alpha=opacity,
                     color='r',
                     xerr=reroutedConf,
                     error_kw=error_config,
                     log=useLog,
                     label='reroutedBar')
                     
    rectsTradFiltered = ax1.barh(index, numFilterBar, bar_width,
    #                 linewidth=2,
                     alpha=opacity,
                     color='b',
                     xerr=numFilterConf,
                     error_kw=error_config,
                     log=useLog,
                     label='numFilterBar')

    # Now the Times
    rectsTimeRerouted = ax2.barh(index+bar_width+eps, deltaTimeBar, bar_width,
                     alpha=opacity,
                     color='r',
                     xerr=deltaTimeConf,
                     error_kw=error_config,
                     log=useLog,
                     label='deltaTimeBar')
                     
    rectsTimeFiltered = ax2.barh(index, filterMinutesBar, bar_width,
                     alpha=opacity,
                     color='b',
                     xerr=filterMinutesConf,
                     error_kw=error_config,
                     log=useLog,
                     label='filterMinutesBar')
    
    
    # Now the Fuel
    rectsFuelRerouted = ax3.barh(index+bar_width+eps, deltaFuelBar, bar_width,
                     alpha=opacity,
                     color='r',
                     xerr=deltaFuelConf,
                     error_kw=error_config,
                     log=useLog,
                     label='deltaFuelBar') 
    
    # Now the Distance
    rectsDistRerouted = ax4.barh(index+bar_width+eps, deltaDistBar, bar_width,
                     alpha=opacity,
                     color='r',
                     xerr=deltaDistConf,
                     error_kw=error_config,
                     log=useLog,
                     label='deltaDistBar')
              
    ax1.set_xbound(xlimNumAffected)
    ax2.set_xbound(xlimTime)
    ax3.set_xbound(xlimFuel)
    ax4.set_xbound(xlimDist)
    
    # plt.gca().xaxis.grid(True)
    ax1.xaxis.grid(True, linewidth=2)
    ax2.xaxis.grid(True, linewidth=2)
    ax3.xaxis.grid(True, linewidth=2)
    ax4.xaxis.grid(True, linewidth=2)
    
    #ax1.set_ybound(([0,40]))
    #ax2.ticklabel_format(size=400)
    #ax1.tick_params(size=400)
    
    plt.setp(ax1.get_yticklabels(), fontweight=800)  # Scenario name font size
    plt.setp(ax1.get_yticklabels(), fontsize=20)  # Scenario name font size
                                   # My understanding is that when shared, only first axis is shown, all others set to invisible.
    plt.setp(ax1.get_xticklabels(), fontsize=20)  # x-axis number size
    plt.setp(ax1.get_xticklabels(), fontweight=800)  # x-axis number weight
    plt.setp(ax2.get_xticklabels(), fontsize=20)  # x-axis number size
    plt.setp(ax2.get_xticklabels(), fontweight=800)  # x-axis number weight
    plt.setp(ax3.get_xticklabels(), fontsize=20)  # x-axis number size
    plt.setp(ax3.get_xticklabels(), fontweight=800)  # x-axis number weight
    plt.setp(ax4.get_xticklabels(), fontsize=20)  # x-axis number size
    plt.setp(ax4.get_xticklabels(), fontweight=800)  # x-axis number weight


#    ax1.semilogx()
    
    # NOTE: ticklable objects are just Text objects
    
    # bbox(x,y) as percent from origin (lowerleft)
    lg = plt.legend([rectsTradRerouted, rectsTradFiltered], ["Rerouted", "Delayed"], bbox_to_anchor=legendAnchor) 
    plt.setp(lg.get_texts(), fontsize=20)  # x-axis number size
    
#    plt.setp(ax1.get_yticklabels(), horizontalalignment='center')
#    ax1.tick_params(axis='y', pad=110) # Moves the ylabels back a litle bit from the edge of the graph.
    
#    plt.setp(ax1.get_yticklabels(), horizontalalignment='left')
#    ax1.tick_params(axis='y', pad=200) # Moves the ylabels back a litle bit from the edge of the graph.
    
    plt.setp(ax1.get_yticklabels(), horizontalalignment='right')
    ax1.tick_params(axis='y', pad=15) # Moves the ylabels back a litle bit from the edge of the graph.
        
    
    
    #ax1.yaxis.labelpad = 100
    
    plt.tight_layout()
    #plt.show()
    
    #plt.savefig(outputPath+"testOneMonth", dpi=600)
    plt.savefig(outputPath+outputFile, format='pdf')
    
    return ax1, ax2, ax3, ax4, f



outputFile = "allTrad.pdf"
#figsize=(20, 15)  # x,y
#xlimNumAffected = (1e0,1e4)
#xlimTime = (1e0,1e4)
#xlimFuel = (1e1,1e5)
#
#plotDict = dict(figsize=figsize, xlimNumAffected=xlimNumAffected, 
#                xlimTime=xlimTime, xlimFuel=xlimFuel,
#                useLog=True, titlePrefix="Mean", plotOption=useMean)

if useAll:
    # Use this if doing all of the mission
    plotDict = dict(figsize             = (20, 13), 
                    xlimNumAffected     = (1e0,1e4), 
                    xlimTime            = (1e0,1e4), 
                    xlimFuel            = (1e1,1e5),
                    xlimDist            = (1e1,1e5),
                    legendAnchor        = (0.95, 0.95), # bbox(x,y) as percent from origin (lowerleft)
                    useLog              = True, 
                    titlePrefix         = "Mean", 
                    plotOption          = useMean)
else:
    # Use this if only showing 4 missions
    plotDict = dict(figsize             = (20, 4.5), 
                    xlimNumAffected     = (1e-1,1e5), 
                    xlimTime            = (1e-1,1e5), 
                    xlimFuel            = (1e-1,1e5),
                    xlimDist            = (1e-1,1e5),
                    legendAnchor        = (1.05, 1.), # bbox(x,y) as percent from origin (lowerleft)
                    useLog              = True, 
                    titlePrefix         = "Mean", 
                    plotOption          = useMean)



                
topToBottomIndividual = []
topToBottomIndividual.extend(topToBottom2018L)
topToBottomIndividual.extend(topToBottom2018M)
topToBottomIndividual.extend(topToBottom2018H)
ax1, ax2, ax3, ax4, f = PlotTest(tradMissions, sc2tl, topToBottomIndividual, plotDict, outputPath, outputFile)

#%%
outputFile = "allEnv.pdf"

if useAll:
    # Use this if doing all of the mission
    plotDict = dict(figsize             = (20, 13), 
                    xlimNumAffected     = (1e-1,1e4), 
                    xlimTime            = (1e-1,1e4), 
                    xlimFuel            = (1e-1,1e4),
                    xlimDist            = (1e-1,1e4),
                    legendAnchor        = (0.95, 1.), # bbox(x,y) as percent from origin (lowerleft)
                    useLog              = True, 
                    titlePrefix         = "Max", 
                    plotOption          = useMax)
else:
    # Use this if only showing 4 missions
    plotDict = dict(figsize             = (20, 4.5), 
                    xlimNumAffected     = (1e-1,1e5), 
                    xlimTime            = (1e-1,1e5), 
                    xlimFuel            = (1e-1,1e5),
                    xlimDist            = (1e-1,1e5),
                    legendAnchor        = (1.05, 1.), # bbox(x,y) as percent from origin (lowerleft)
                    useLog              = True, 
                    titlePrefix         = "Max", 
                    plotOption          = useMax)


#plotDict = dict(figsize             = (20, 15), 
#                xlimNumAffected     = (0.1,5), 
#                xlimTime            = (0.1,100), 
#                xlimFuel            = (0.1,100),
#                xlimDist            = (0.1,100),
#                legendAnchor        = (0.95, 0.95), # bbox(x,y) as percent from origin (lowerleft)
#                useLog              = False, 
#                titlePrefix         = "Mean", 
#                plotOption          = useMean)
                
                
topToBottomIndividual = []
topToBottomIndividual.extend(topToBottom2018L)
topToBottomIndividual.extend(topToBottom2018M)
topToBottomIndividual.extend(topToBottom2018H)
ax1, ax2, ax3, ax4, f = PlotTest(envMissions, sc2tl, topToBottomIndividual, plotDict, outputPath, outputFile)



# dataKeys
#['rerouted',
# 'DeltaTime',
# 'DeltaDist',
# 'DeltaFuel',
# 'numFilter',
# 'sumFilterMinutes']


#%% Grouped By Scenario / Day

tradScenarios = CollapseIntoScenaros(tradMissions, missionsInGroup)
envScenarios = CollapseIntoScenaros(envMissions, missionsInGroup)

ProcessAndUpdate(tradScenarios)
ProcessAndUpdate(envScenarios)

# Trad
outputFile="tradScenarios.pdf"

plotDict = dict(figsize             = (20, 4.5), 
                xlimNumAffected     = (1e-1,1e4), 
                xlimTime            = (1e-1,1e4), 
                xlimFuel            = (1e-1,1e5),
                xlimDist            = (1e-1,1e5),
                legendAnchor        = (0.95, 0.9), # bbox(x,y) as percent from origin (lowerleft)
                useLog              = True, 
                titlePrefix         = "Mean", 
                plotOption          = useMean)

#plotDict = dict(figsize             = (20, 4), 
#                xlimNumAffected     = (0,250), 
#                xlimTime            = (0,3000), 
#                xlimFuel            = (1,4e4),
#                xlimDist            = (1,4e3),
#                legendAnchor        = (0.95, 0.95), # bbox(x,y) as percent from origin (lowerleft)
#                useLog              = False, 
#                titlePrefix         = "Mean", 
#                plotOption          = useMean)

topToBottom = ["2018L", "2018M", "2018H"]
ax1, ax2, ax3, ax4, f = PlotTest(tradScenarios, scenarioToTickMap, topToBottom, plotDict, outputPath, outputFile)

#%%
# Env
outputFile="envScenarios.pdf"

plotDict = dict(figsize             = (20, 4.5), 
                xlimNumAffected     = (1e-1,1e4), 
                xlimTime            = (1e-1,1e4), 
                xlimFuel            = (1e-1,1e5),
                xlimDist            = (1e-1,1e5),
                legendAnchor        = (0.95, 0.9), # bbox(x,y) as percent from origin (lowerleft)
                useLog              = True, 
                titlePrefix         = "Max", 
                plotOption          = useMax)

#plotDict = dict(figsize             = (20, 4), 
#                xlimNumAffected     = (0,20), 
#                xlimTime            = (0,2000), 
#                xlimFuel            = (0,400),
#                xlimDist            = (0,40),
#                legendAnchor        = (0.95, 0.95), # bbox(x,y) as percent from origin (lowerleft)
#                useLog              = False, 
#                titlePrefix         = "Max", 
#                plotOption          = useMax)

topToBottom = ["2018L", "2018M", "2018H"]
ax1, ax2, ax3, ax4, f = PlotTest(envScenarios, scenarioToTickMap, topToBottom, plotDict, outputPath, outputFile)


#%% Pull outsome values

#['Lynx_Cecil_2018H',
# 'SS2_Titus_2018H',
# 'Lynx_Midland_2018L',
# 'PM_CornNew_2018H',
# 'Reentry_PacificReentryLA_2018M',
# 'Lynx_Midland_2018H',
# 'Atlas5_Vafb_2018M',
# 'Pegasus_Wallops_2018L',
# 'Lynx_FrntRnge_2018H',
# 'Sound_America_2018H',
# 'Antares_Wallops_2018M',
# 'SS2_America_2018H',
# 'SS2_America_2018L',
# 'SS2_America_2018M']


curVehicle = "Lynx_FrntRnge_2018H"
qoi = "maximum"
print "{0} Trad {1}".format(curVehicle, qoi)
for elem in tradMissions[curVehicle][qoi].items():
    print "   {0}".format(elem)
print "{0} Env {1}".format(curVehicle, qoi)
for elem in envMissions[curVehicle][qoi].items():
    print "   {0}".format(elem)


#for elem in zip(np.array(tradMissions[curVehicle]['dataVals'])[:,0], np.array(tradMissions[curVehicle]['dataVals'])[:,4]):
#    print elem

#for elem in zip(np.array(envMissions[curVehicle]['dataVals'])[:,0], np.array(envMissions[curVehicle]['dataVals'])[:,4]):
#    print elem


#%%

#['2018L', '2018M', '2018H']

curVehicle = "2018H"
qoi = "average"
print "{0} Trad {1}".format(curVehicle, qoi)
for elem in tradScenarios[curVehicle][qoi].items():
    print "   {0}".format(elem)
print "{0} Env {1}".format(curVehicle, qoi)
for elem in envScenarios[curVehicle][qoi].items():
    print "   {0}".format(elem)


#tradScenarios[curVehicle]["average"]["rerouted"] + tradScenarios[curVehicle]["average"]["numFilter"]
envScenarios[curVehicle]["average"]["rerouted"] + envScenarios[curVehicle]["average"]["numFilter"]



# topToBottom2018L = ["SS2_America_2018L",
#                      "Lynx_FrntRnge_2018H",
#                      "Atlas5_Vafb_2018M",
#                      "Antares_Wallops_2018M"]
# tradMissions['Lynx_FrntRnge_2018H']
envMissions['SS2_America_2018L']
envMissions['Lynx_FrntRnge_2018H']
envMissions['Atlas5_Vafb_2018M']
envMissions['Antares_Wallops_2018M']

for curMission in envMissions:
    print "{0}  avg rerouted = {1}".format(curMission, envMissions[curMission]['average']['rerouted'])








