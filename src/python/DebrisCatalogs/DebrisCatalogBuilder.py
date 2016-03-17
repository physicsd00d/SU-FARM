import numpy as np

# Conversion Factors
lbs2kg  = 0.453592
ft2m    = 0.3048
TBD = 0.


### The ET keys in order of appearance
# Forward Section LOX Tank
# Mid Section Intertank and LH2 Tank
# Aft Section LH2 Tank
#
# SkinEtc
# # LOX Feedline Elbow
# # LOX Feedline Segment
# # LOX Antigeyser Line and Elbow
# # GO2 Pressurization Line
# # GH2 Pressurization Line
# # LOX Vent Valve Actuation Line
# # Nose Cap Purge Duct
# # Helium Injection Line
# # Skin Segment Lox Tank (a)
# # Skin Segment Lox Tank (b)
# # Skin Segment LH2 Tank (a)
# # Skin Segment LH2 Tank (b)
# # Skin Segment LH2 Tank (c)
# # Ring Frame Stabilizer
# # Stability Ring Segment
# # Slosh Baffle Piece
# # Slosh Baffle Piece
# # Slosh Baffle Piece
# # Slosh Baffle Piece
# # Cable Tray Segments
# # LSC Sheath Segments
#
# MiscEtc
# # Misc (Cabling / Bolts / Brackets / etc)


### The ORBITER keys that are left after deletions
# 'Fuselage Mid Section',
# 'Nose Landing Gear',
# 'RCS Nozzles',
# 'RCS Press Tanks',
# 'RCS Fuel Tanks',
# 'Orbiter / ET Attach Struts',
# 'Beams Trusses Pipes',
# 'Payload Bay Doors',
# 'Landing Gear Doors',
# 'Skin Sections',
# 'Misc'



def calcBetaProxy(weightIndividual, aref):
# This is basically assuming Cd = 1
    try:
        minWeight = min(weightIndividual)
        maxWeight = max(weightIndividual)
    except:
        minWeight = weightIndividual
        maxWeight = minWeight

    try:
        minAref = min(aref)
        maxAref = max(aref)
    except:
        minAref = aref
        maxAref = minAref
    return np.array([minWeight/maxAref, maxWeight/minAref])

def readCatalog(fileName):
    import re

    curCatalog = dict()
    file = open(fileName,'r')

    # Figure out what kind of file we're looking at.  There are two kinds:
    #   * CSV, the format of which can be seen in the ShuttleCatalogs
    #   * TXT, this is the format that the propagation code (Debris Reader) wants
    fileSuffix = fileName.split('.')[-1].upper()
    isCSV = fileSuffix == 'CSV'
    isTXT = fileSuffix == 'TXT' # Will need to modify each line to fit it into the default CSV format

    if isCSV:

        for line in file:

            key = line.split(',')

            if (len(key[0]) == 0):
                # This is the first line of the file.  Ignore it.
                continue
            else:
                # If it's possible for the quantity to be a range, then always treat it as a range to be consistent
                # Convert units into kg and m
                name                = key[0]

                subkey              = str.split(key[1],'-')
                numPieces           = np.array([int(subkey[0]), int(subkey[-1])])

                subkey              = str.split(key[2],'-')
                weightIndividual    = np.array([float(subkey[0]),       float(subkey[-1])]) * lbs2kg

                subkey              = str.split(key[3],'-')
                weightTotal         = np.mean([float(subkey[0]),       float(subkey[-1])]) * lbs2kg
                # weightTotal         = float(key[3]) * lbs2kg

                subkey              = str.split(key[4],'-')
                aref                = np.array([float(subkey[0]),  float(subkey[-1])]) * pow(ft2m,2)

                subkey              = str.split(key[5],'-')
                # beta                = np.array([float(subkey[0]),  float(subkey[-1])]) * lbs2kg/pow(ft2m,2)
                beta                = 0.

                subkey              = str.split(key[6],'-')
                impulseVel          = np.array([float(subkey[0]),         float(subkey[-1])]) * ft2m

                CD                  = key[7]
                CL                  = key[8].strip('\n')

                curGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                                weightTotal = weightTotal, aref = aref, beta = beta,
                                impulseVel = impulseVel, CD = CD, CL = CL)
                curCatalog[name] = curGroup


    elif isTXT:
        # Note that the format for this sometimes allows small values to be written in scientific notation,
        #   which is incompatible with the CSV format above because I was splitting ranges based on '-'.  Here,
        #   ranges will have to be divided by a comma.

        categoryCount = 0
        for line in file:
            # We don't know where we are in the file.  Temporarily split by spaces to investigate.
            temp = line.split()

            # Make sure line is not empty
            if len(temp) == 0:
                continue

            # The first line will tell us what the units are
            if (temp[0].upper() == 'UNITS'):
                if (temp[1].upper() == 'SI'):
                    print 'The units are SI.  GOOD!'
                    continue
                else:
                    print '\n\n\nBAD UNITS IN DEBRIS CATALOG.  Exiting'
                    raise RuntimeError

            # I don't care about any of the other non-data lines, so look for them here
            isBegin     = (temp[0].upper() == 'BEGIN')
            isVelocity  = (temp[0].upper() == 'VELOCITY')
            isTime      = (temp[0].upper() == 'TIME')
            isComment   = (temp[0][0] == '#')
            if isBegin or isVelocity or isTime or isComment:
                continue

            isEND       = (temp[0].upper() == 'END')
            if isEND:
                break

            # ALRIGHT!  If you made it here, then it's a data line.  Let's parse it!
            # First, remove any possible spaces after the comma
            line = re.sub(', *',',', line)

            # Now remove any possible spaces between square brackets
            line = re.sub('\[ *','',line)  #left brackets
            line = re.sub(' *\]','',line)  #right brackets

            # Make life a little easier and prepend an empty name
            line = 'UNTITLED_{0}  {1}'.format(categoryCount, line)
            categoryCount += 1

            # Now can split it by spaces
            key = line.split()

            name                = key[0]

            subkey              = str.split(key[1],',')
            numPieces           = np.array([int(subkey[0]), int(subkey[-1])])

            subkey              = str.split(key[2],',')
            weightIndividual    = np.array([float(subkey[0]),       float(subkey[-1])])

            subkey              = str.split(key[3],',')
            weightTotal         = np.mean([float(subkey[0]),       float(subkey[-1])])
            # weightTotal         = float(key[3]) * lbs2kg

            subkey              = str.split(key[4],',')
            aref                = np.array([float(subkey[0]),  float(subkey[-1])])

            subkey              = str.split(key[5],',')
            beta                = np.array([float(subkey[0]),  float(subkey[-1])])
            # beta                = 0.

            subkey              = str.split(key[6],',')
            impulseVel          = np.array([float(subkey[0]),         float(subkey[-1])])

            CD                  = key[7]
            CL                  = key[8].strip('\n')

            curGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                            weightTotal = weightTotal, aref = aref, beta = beta,
                            impulseVel = impulseVel, CD = CD, CL = CL, name = name)
            curCatalog[name] = curGroup




    else:
        print 'THIS IS NOT A RECOGNIZED FILE SUFFIX!!!  DOING NOTHING.'
        # print key[0]
    file.close()

    return curCatalog



def logDistribute(x):
    numPiecesOriginal           = x['numPieces']
    weightIndividualOriginal    = x['weightIndividual']

    # Find the bounds of the range
    try:
        # Usually this will be a tuple
        lo = weightIndividualOriginal[0]
        hi = weightIndividualOriginal[1]
    except:
        # Unless it's not
        lo = weightIndividualOriginal
        hi = lo

    ans = []

    if (hi > lo) and (lo < 0.3):
        # Find how many orders of magnitude there are between them
        loglo = np.log10(lo)
        loghi = np.log10(hi)
        # diff = np.ceil(loghi - loglo)
        # if ((loghi - loglo) == diff):
        #     # Exactly an order(s) of magnitude off
        #     diff += 1

        logloLimit = np.ceil(loglo)
        loghiLimit = np.floor(loghi)
        innerDiff = int(loghiLimit-logloLimit) + 1      # add one for zero

        # Construct the log limits
        innerRange = np.array(range(innerDiff)) + logloLimit
        newLimits = [loglo]
        newLimits.extend(innerRange)
        newLimits.append(loghi)

        # Convert to linear limits
        newLimits = pow(10, np.array(newLimits))

        # Package up the new groups
        for ix in range(len(newLimits)-1):
            numPieces           = np.round(numPiecesOriginal/len(newLimits))
            weightIndividual    = np.array([newLimits[ix], newLimits[ix+1]])
            weightTotal         = numPieces * np.mean(weightIndividual)
            aref                = x['aref']
            # beta                = x['beta']
            beta                = calcBetaProxy(weightIndividual, aref)
            impulseVel          = x['impulseVel']
            CD                  = x['CD']
            CL                  = x['CL']
            name                = x['name'] + '__' + str(ix)
            newGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                            weightTotal = weightTotal, aref = aref, beta = beta,
                            impulseVel = impulseVel, CD = CD, CL = CL, name = name)
            ans.append(newGroup)

    else:
        # There was no range to the individual masses.  Just return the original group.
        ans.append(x)

    return ans


def scaleIndividualWeight(name, GroupsToScale, targetMass, curCatalog, numPiecesIn = -1):
    ansGroup = []

    # First loop through to find the total mass of the pieces in the group
    totalMass = 0.
    for groupStr in GroupsToScale:
        thisGroup = curCatalog[groupStr]
        totalMass += thisGroup['weightTotal']

    # print the scaling factor
    downFactor = targetMass/totalMass

    # Now loop again and scale the masses of the pieces by SkinEtcDownFactor
    for groupStr in GroupsToScale:
        thisGroup           = curCatalog[groupStr]

        # In addition to scaling the weights/areas, allow the user to possibly specify the number of pieces
        if numPiecesIn > 0:
            numPieces = numPiecesIn
        else:
            numPieces           = thisGroup['numPieces']

        weightIndividual    = thisGroup['weightIndividual'] * downFactor
        weightTotal         = numPieces * weightIndividual      # Will eventually output this as the avg, but good to see range for debugging
        aref                = thisGroup['aref'] * downFactor
        # beta                = thisGroup['beta']                 # I scaled mass and area by the same factor, so beta remains unchanged
        beta                = calcBetaProxy(weightIndividual, aref)
        impulseVel          = thisGroup['impulseVel']
        CD                  = thisGroup['CD']
        CL                  = thisGroup['CL']
        curGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                            weightTotal = weightTotal, aref = aref, beta = beta,
                            impulseVel = impulseVel, CD = CD, CL = CL, name = name)

        ansGroup.extend(logDistribute(curGroup))

    return ansGroup

# This scales the number of pieces and sets the total weight equal to the targetMass passed in
def scaleNumPieces(name, groupList, targetMass, curCatalog):
    ansGroup = []

    # First loop through to find the total mass of the pieces in the group
    totalMass = 0.
    for groupStr in groupList:
        thisGroup = curCatalog[groupStr]
        totalMass += thisGroup['weightTotal']

    # print SkinEtcETMass
    downFactor = targetMass/totalMass

    # Now loop again and scale the masses of the pieces by SkinEtcDownFactor
    for groupStr in groupList:
        thisGroup           = curCatalog[groupStr]

        numPieces           = np.round(np.mean(thisGroup['numPieces'] * downFactor))
        weightIndividual    = thisGroup['weightIndividual']
        weightTotal         = targetMass                # It's up to you, later, to distribute the pieces to make this happen
        aref                = thisGroup['aref']
        # beta                = thisGroup['beta']
        beta                = calcBetaProxy(weightIndividual, aref)
        impulseVel          = thisGroup['impulseVel']
        CD                  = thisGroup['CD']
        CL                  = thisGroup['CL']
        curGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                        weightTotal = weightTotal, aref = aref, beta = beta,
                        impulseVel = impulseVel, CD = CD, CL = CL, name = name)

        ansGroup.extend(logDistribute(curGroup))
    return ansGroup


# This scales the number of pieces and also scales the total weight, leaving individual weights intact
def scaleNumPiecesAndTotalWeight(name, groupList, targetMass, curCatalog):
    ansGroup = []

    # First loop through to find the total mass of the pieces in the group
    totalMass = 0.
    for groupStr in groupList:
        thisGroup = curCatalog[groupStr]
        totalMass += thisGroup['weightTotal']

    # print SkinEtcETMass
    downFactor = targetMass/totalMass
    print 'downFactor = {0}'.format(downFactor)

    # Now loop again and scale the masses of the pieces by SkinEtcDownFactor
    for groupStr in groupList:
        thisGroup           = curCatalog[groupStr]

        numPieces           = np.ceil(np.mean(thisGroup['numPieces'] * downFactor))
        weightTotal         = thisGroup['weightTotal'] * downFactor               # It's up to you, later, to distribute the pieces to make this happen
        weightIndividual    = thisGroup['weightIndividual']
        if numPieces == 1:
            weightIndividual = weightTotal
        aref                = thisGroup['aref']
        beta                = thisGroup['beta']
        # beta                = calcBetaProxy(weightIndividual, aref)   #BE CAREFUL!  If you choose an even slightly bad beta, the propagation code could hang.
        impulseVel          = thisGroup['impulseVel']
        CD                  = thisGroup['CD']
        CL                  = thisGroup['CL']
        curGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                        weightTotal = weightTotal, aref = aref, beta = beta,
                        impulseVel = impulseVel, CD = CD, CL = CL, name = name)

        ansGroup.extend([curGroup])
        # ansGroup.extend(logDistribute(curGroup))
    return ansGroup



def noScaleAnalog(name, weightIndividual, numPieces, aref, curCatalog, analogPiece):
    ans = []
    weightTotal         = numPieces * weightIndividual
    # beta                = TBD
    beta                = calcBetaProxy(weightIndividual, aref)
    impulseVel          = curCatalog[analogPiece]['impulseVel']
    CD                  = curCatalog[analogPiece]['CD']
    CL                  = curCatalog[analogPiece]['CL']
    curGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                        weightTotal = weightTotal, aref = aref, beta = beta,
                        impulseVel = impulseVel, CD = CD, CL = CL, name = name)
    ans.extend(logDistribute(curGroup))
    return ans


def totalMassOfCatalog(curCatalog):
    Mass = 0.
    for groupStr in curCatalog:
        thisGroup = curCatalog[groupStr]
        Mass += thisGroup['weightTotal']
    return Mass



def collapseRange(x):
    # If it's not a range, then just output one number
    try:
        # This will work if the value is a range
        if x[0] == x[1]:
            x = '{0:5.4}'.format(x[0])
        else:
            x = '[{0:7.4},{1:7.4}]'.format(x[0], x[1])
    except:
        x = '{0:6.5}'.format(x)

    return x




def makeString(cur):
    numPieces = str(int(np.mean(cur['numPieces'])))

    # If it's not a range, then just output one number
    weightIndividual = collapseRange(cur['weightIndividual'])

    # If it's a range, output the mean.
    weightTotal = '{0:7.5}'.format(np.mean(cur['weightTotal']))

    aref = collapseRange(cur['aref'])
    beta = collapseRange(cur['beta'])


    impulseVel = collapseRange(cur['impulseVel'])

    # impulseVel = cur['impulseVel']
    # if impulseVel[0] == impulseVel[1]:
    #     impulseVel = '{0:4.4}'.format(impulseVel[0])
    # else:
    #     impulseVel = '[{0:4.4},{1:4.4}]'.format(impulseVel[0], impulseVel[1])

    CD = cur['CD']
    CL = cur['CL']


    outString = '{0:15}{1:25}{2:15}{3:25}{4:25}{5:20}{6:10}{7:10}'.\
        format(numPieces, weightIndividual, weightTotal, aref, beta, impulseVel, CD, CL)
    return outString


def getLegend():
    outString = '#{0:14}{1:25}{2:15}{3:25}{4:25}{5:20}{6:10}{7:10}'.\
            format('numPieces', 'weightIndividual', 'weightTotal', 'aref', 'beta', 'impulseVel', 'CD', 'CL')
    return outString


def makeOutputFile(finalCatalog, fileName):
    # Print the catalog nicely
    outFile = open(fileName, 'w')

    outFile.write('units SI\nbegin\nvelocity  666.0\ntime   100.0\n')

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

    outFile.write('end')
    outFile.close()



# =========================================
#
#
#          Testing Group Function
#
#
# =========================================

# Not actually used anymore
def yesOrNo():
    waiting = True
    while waiting:
        isOkay = raw_input('Okay? [y/n]')
        if isOkay == 'y':
            isOkay = True
            waiting = False
        elif isOkay == 'n':
            isOkay = False
            waiting = False
        else:
            print 'Invalid input, try again '

    return isOkay


# This function prompts the user for input, called from groupSimilarPieces
def getGrouping(numHere):

    ans = []

    waiting = True
    while waiting:

        try:
            groupInput = input('all[100], none[-1], or comma separated list:  ')

            waiting = False
            if groupInput == 100:
                ans = range(numHere)
            elif groupInput == -1:
                ans = []
            else:
                try:
                    # This will throw error if only given a single index
                    if len(groupInput) > 1:
                        ans = groupInput
                except:
                    # So catch the error by making the single number into a length-1 array
                    ans = [groupInput]

        except:
            print 'invalid input!  No strings allowed.  Try again.'

    return np.array(ans)




# When you know the subgroups to combine, this function will put them together
def combineGroups(subgroups):
    name = 'Combined_' + str(len(subgroups))
    numPieces = 0
    weightTotal = 0
    minaref = 1e9
    maxaref = 0
    minMass = 1e9
    maxMass = 0

    for curGroup in subgroups:
        numPieces           += curGroup['numPieces']
        weightTotal         += curGroup['weightTotal']
        minaref = min(minaref, curGroup['aref'][0])
        maxaref = max(maxaref, curGroup['aref'][1])
        minMass = min(minMass, curGroup['weightIndividual'][0])
        maxMass = max(maxMass, curGroup['weightIndividual'][1])
        impulseVel           = curGroup['impulseVel']

    aref                = np.array([minaref, maxaref])
    beta                = TBD
    CL                  = 'NONE'
    weightIndividual    = np.array([minMass, maxMass])
    curGroup = dict(numPieces = numPieces, weightIndividual = weightIndividual,
                            weightTotal = weightTotal, aref = aref, beta = beta,
                            impulseVel = impulseVel, CD = CD, CL = CL, name = name)
    return curGroup



# Prompt the user, go through all the CD categories, and determine which groups to combine
def groupSimilarPieces(curCat):
    finalCat = []   # This is where the answers go

    flatCat = []
    CD_Keys = []
    betaVec = []
    minMassVec = []
    maxMassVec = []

    # First, flatten the catalog into a dictionary
    for topLevelSystem in curCat:
        if type(topLevelSystem) is dict:
            outStr = makeString(topLevelSystem)
            flatCat.append(topLevelSystem)
            CD_Keys.append(topLevelSystem['CD'])
            betaVec.append(topLevelSystem['beta'])

            try:
                minMassVec.append(topLevelSystem['weightIndividual'][0])
                maxMassVec.append(topLevelSystem['weightIndividual'][1])
            except:
                minMassVec.append(topLevelSystem['weightIndividual'])
                maxMassVec.append(topLevelSystem['weightIndividual'])

            print '#' + topLevelSystem['name']
            print outStr

        else:
            # It's an array of subsystems
            for lowLevelSystem in topLevelSystem:
                outStr = makeString(lowLevelSystem)
                flatCat.append(lowLevelSystem)
                CD_Keys.append(lowLevelSystem['CD'])
                betaVec.append(lowLevelSystem['beta'])

                minMassVec.append(lowLevelSystem['weightIndividual'][0])
                maxMassVec.append(lowLevelSystem['weightIndividual'][1])

                print outStr

    # Gather up the unique CD values
    CD_Keys = np.array(CD_Keys)
    flatCat = np.array(flatCat)
    betaVec = np.array(betaVec)
    minMassVec = np.array(minMassVec)
    maxMassVec = np.array(maxMassVec)

    for CD in np.unique(CD_Keys):
        # Find the ones that match the unique key
        boolVec = (CD_Keys == CD)

        cdGroup = flatCat[CD_Keys == CD]
        print cdGroup
        cdGroup = sorted(cdGroup, key=lambda k: min(k['weightIndividual']))

        while len(cdGroup) > 0:
            numHere = len(cdGroup)
            if numHere > 1:
                # Print out the members of this group and see how the user wants to group them
                print '\n\n\n' + 'CD = ' + str(CD)
                print 'ID    ' + getLegend()
                for ix in range(len(cdGroup)):
                    print str(ix) + '     ' + makeString(cdGroup[ix])

                groupWhich = getGrouping(numHere)
                print 'groupWhich = ' + str(groupWhich)

                if len(groupWhich) == 0:
                    # Don't group any of them together
                    for ix in range(numHere):
                        finalCat.append(cdGroup[ix])
                    cdGroup = np.delete(cdGroup, range(numHere))    # remove them all from cdGroup
                else:
                    # Group the ones we were told
                    tempGroup = []
                    for ix in groupWhich:
                        tempGroup.append(cdGroup[ix])

                    newGroup = combineGroups(tempGroup)
                    finalCat.append(newGroup)
                    cdGroup = np.delete(cdGroup, groupWhich)    # remove them from cdGroup
            else:
                finalCat.append(cdGroup[0])
                cdGroup = np.delete(cdGroup,0)


    return flatCat, CD_Keys, betaVec, finalCat





















