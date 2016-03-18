'''
==== General procedure here ====
* I made a spreadsheet of the debris catalog (with my eyeballs) out of the CAIB report
* Read in the results of that spreadsheet (number of pieces per ballistic coefficient range per time)
* Translate the bounds of beta to bounds on mass using correlation from paper
* Aref range is set to be overly inclusive so that in the debris reader, random areas will be sampled and kept
    or thrown out based on if they fit the beta range.  This way we don't over-specify the system and wind up
    creating an effectively smaller-than-it-should-be beta range
* Load this information into debris catalog file that can be run in the main risk code

Author: Thomas J Colvin
'''
import numpy as np
import matplotlib.pyplot as plt
import DebrisCatalogBuilder as dcb

import sys

lbs2kg  = 0.453592
grav    = 9.8          #m/s^2
ft2m    = 0.3048

# Read in the columbia group catalog time-dependent histogram

fileName = "Columbia/Group1.csv"
tag = '#GROUP1'

betaRanges = []
debrisRecord = []
# curCatalog = dict()
file = open(fileName,'r')
for line in file:
    key = line.split(',')

    if (len(key) == 0):
        # This is the first line of the file or is just empty.  Ignore it.
        continue
    elif key[0] == tag:
        print 'key[0]'
    elif key[0] == '#END':
        # Done with reading the file
        break
    elif key[0] == "#Beta":
        # lbs/ft^2
        numRanges = len(key[1:])
        for subkey in key[1:]:
            subkey = subkey.split('-')
            betaRanges.append([float(subkey[0]), float(subkey[1])])

    elif key[0] == "#Group":
        # Figure out how many groups
        # We know the first key is the tag and the last key is newline, so loop through looking for numbers
        numGroups = len(key[1:])
        print 'numGroups = ' + str(numGroups)
        print key
    else:
        # print line
        # THis should be
        lineRecord = [int(key[0])]
        # lineRecord.append(int(key[0]))

        for ix in range(numGroups):
            lineRecord.append(int(key[ix+1]))

        debrisRecord.append(lineRecord)

file.close()

# Time is in reversed order, flip it and make it numpy
debrisRecord = np.flipud(debrisRecord)

# make betaRanges into numpy array and convert units from lbs/ft^2 into kg/m^2
betaRanges = np.array(betaRanges) * lbs2kg / (ft2m**2)

# Having the smallest beta possible be zero will give pieces that take forever and an probably unrealistic
# Find those pieces and change them to something that is a priori more realistic
betaRanges[0][0] = betaRanges[0][1]/10

# Convert the ranges into bins
betaBins = np.append(betaRanges[:,0], betaRanges[-1,1])









# # sample a lognormal distribution
# # mu, sigma = .2, 2 # mean and standard deviation, lbs
# mu, sigma = .2*lbs2kg, 2*lbs2kg # mean and standard deviation, kg
#
# s = np.random.lognormal(mu, sigma, 15470)
# # s = np.random.lognormal(mu, sigma, 75000)
# massBins = np.logspace(-4,3, 100)   #kg
#
# ## NOTE!!!  I am fudging the first bin (0.5->0.66) to make it line up nicely with 300g
# # bins = np.array([0.0, 0.66, 1, 1.78, 3.33, 5.6, 10, 17.78, 30, 56.23, 100, 177, 333]) #beta...fuck
#
# # # Divide the bins in half, append to existing bins, then sort them
# # bins = np.append(bins, bins[:-1] + np.diff(bins)/2)
# # bins.sort()
# print 'massBins = ' + str(massBins)
#
#
#
#
# # bin the data and plot
# plt.figure(1)
# plt.xscale('log')
# count, massBins, ignored = plt.hist(s,bins=massBins)
# plt.axis('tight')
# # plt.show()

# I'm guessing at these parameters to match the graph
corr = 0.5
upperIntercept = 2.9  # These are good for estimating the paper, but perhaps not so good for Columbia
lowerIntercept = 2.2

# # Guessing at values for Columbia
# upperIntercept = 2.9
# lowerIntercept = 1.2

# All pieces in the Columbia catalog effectively have an average Cd=0.6
Cd = 0.6




massHighs   = (betaBins * grav * (10**-lowerIntercept))**(1/corr)
massLows    = (betaBins * grav * (10**-upperIntercept))**(1/corr)

print zip(massLows, massHighs)

debrisCatalog = dict()

for row in debrisRecord:
    # Parse the row
    tstep = row[0]
    curInfo = row[1:]

    debrisCatalogAtTime = []

    # Figure out which bins have debris in them at this tstep
    debrisHereIX = (curInfo > 0)

    curBetaRanges = betaRanges[debrisHereIX]

    # curMassLows = (curBetaRanges * grav * (10**-upperIntercept))**(1/corr)
    # curMassHighs = (curBetaRanges * grav * (10**-lowerIntercept))**(1/corr)


    curMassLowVec = np.mean((curBetaRanges * grav * (10**-upperIntercept))**(1/corr), axis=1)
    curMassHighVec = np.mean((curBetaRanges * grav * (10**-lowerIntercept))**(1/corr), axis=1)

    # curMassRange = np.mean(curMassLows,axis=1)

    # print 'tx = ' + str(tstep)
    # print curBetaRanges
    # print curMassLows
    # print curMassHighs
    # print (curBetaRanges * grav * (10**-lowerIntercept))**(1/corr)
    # print ''

    for ix in range(len(curMassLowVec)):
        numPieces   = curInfo[debrisHereIX][ix]
        curMassLow  = curMassLowVec[ix]
        curMassHigh = curMassHighVec[ix]

        curBetaLow  = curBetaRanges[ix,0]
        curBetaHigh = curBetaRanges[ix,1]

        # Group 1 was not given any impulse velocity?
        curImpVelLow    = 0.
        curImpVelHigh   = 100.

        # Going to allow the debris reader to choose the Arefs such that Betas work out
        #   That said, i want it to sample somewhat close to the right values so multiply these expressions by some margin
        curUpperAref   = 2*curMassHigh/(curBetaLow * Cd)
        curLowerAref   = 0.5*curMassLow/(curBetaHigh * Cd)
        if np.isinf(curUpperAref):
            # Sometimes betaLow = 0, which will make this inf.  Set area to something reasonable
            curUpperAref = curLowerAref*2*100



        curGroup = dict(numPieces   = np.array([numPieces, numPieces]),
                weightIndividual    = np.array([curMassLow, curMassHigh]),
                weightTotal         = numPieces*(curMassLow+curMassHigh)/2,
                aref                = np.array([curLowerAref, curUpperAref]),
                beta                = np.array([curBetaLow,curBetaHigh]),
                impulseVel          = np.array([curImpVelLow,curImpVelHigh]),
                CD                  = 0.6,
                CL                  = 0.04,
                name                = "")

        debrisCatalogAtTime.append(curGroup)

    debrisCatalog[(tstep-1)*5] = debrisCatalogAtTime


from DebrisCatalogBuilder import makeString

fileName = 'debugColumbiaMarch16.txt'

# Print the catalog nicely
outFile = open(fileName, 'w')
outFile.write('units SI\n')

for timeKeys in sorted(debrisCatalog):

    # Get the catalog at this time
    finalCatalog = debrisCatalog[timeKeys]

    outString = 'begin\nvelocity  666.0\ntime   {0}'.format(timeKeys)
    outFile.write(outString + '\n')

    outString = '#{0:14}{1:25}{2:15}{3:25}{4:25}{5:20}{6:10}{7:10}'.\
        format('numPieces', 'weightIndividual', 'weightTotal', 'aref', 'beta', 'impulseVel', 'CD', 'CL')
    print outString

    outFile.write(outString + '\n')

    for topLevelSystem in finalCatalog:

        if type(topLevelSystem) is dict:
            outStr = makeString(topLevelSystem)

            # outFile.write('#' + topLevelSystem['name'] + '\n')
            outFile.write(outStr + '\n')

            print '#' + topLevelSystem['name']
            print outStr

        else:
            # It's an array of subsystems
            print '#' + topLevelSystem[0]['name']
            # outFile.write('#' + topLevelSystem[0]['name'] + '\n')
            for lowLevelSystem in topLevelSystem:
                outStr = makeString(lowLevelSystem)

                outFile.write(outStr + '\n')

                print outStr

    outFile.write('end\n')
    outFile.write('\n')

outFile.close()




sys.exit()




catalogMass = []
catalogBeta = []


for ix in range(len(massBins)-1):

    # Mass of the piece here (middle of bin)
    curMass = ((massBins[ix+1] + massBins[ix])/2.)

    # How many are here (at least one piece)
    numPieces = max(count[ix],1)

    # Only take pieces that are at least one gram and that are smaller than the largest recovered piece
    if (curMass > 0.001) and (curMass < 400.):
        curLowerBeta = ((curMass)**corr) * (10**lowerIntercept) /grav
        curUpperBeta = ((curMass)**corr) * (10**upperIntercept) /grav

        curUpperAref   = curMass/(curLowerBeta * Cd)
        curLowerAref   = curMass/(curUpperBeta * Cd)

        curGroup = dict(numPieces           = np.array([numPieces, numPieces]),
                        weightIndividual    = np.array([curMass, curMass]),
                        weightTotal         = curMass*numPieces,
                        aref                = np.array([curLowerAref, curUpperAref]),
                        beta                = np.array([curLowerBeta,curUpperBeta]),
                        impulseVel          = np.array([0.,100.]),
                        CD                  = 0.6,
                        CL                  = 0.0,
                        name                = "")

        debrisCatalog.append(curGroup)

        # This is non-essential and just for the graph
        catalogMass.append(curMass)
        catalogMass.append(curMass)
        catalogBeta.append(curLowerBeta * grav)
        catalogBeta.append(curUpperBeta * grav)




plt.figure(2)
mass = np.logspace(-1,2,10)
upperBetaLimit = np.power(mass,corr) * (10**upperIntercept)
lowerBetaLimit = np.power(mass,corr) * (10**lowerIntercept)
plt.xscale('log')
plt.yscale('log')
plt.plot(mass,upperBetaLimit,'ro')
plt.plot(mass,lowerBetaLimit,'bo')
plt.plot(catalogMass, catalogBeta,'x')
plt.xlabel('mass (lbs)')
plt.ylabel('ballCoeff (Pa)')
plt.title('Beta correlation imitation of that one plot')


# fileName = '../../Files/Columbia/DebrisCatalog/testFile_Columbia.txt'






plt.show()




# Generate the catalog file



