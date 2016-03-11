

# import os
# import sys
# sys.path.append(friscoFiles)
# sys.path.append(debrisPropPATH)
# sys.path.append(tjcFiles)
from Simulation import TJC

# Note: You have to import all the frisco stuff because TJC demands to know about it.


# Create the atmosphere pickle.  Should only need to be done once ever.
# atmoFile = 'Cape.txt'
# pickleName = 'Cape.pkl'

# atmoFile = 'WSands.txt'
# pickleName = 'WSands.pkl'

# atmoFile = 'FrontRange.txt'
# pickleName = 'FrontRange.pkl'

# atmoFile = 'SpaceportAmerica.txt'
# pickleName = 'SpaceportAmerica.pkl'

atmoFile = 'WestTexas.txt'
pickleName = 'WestTexas.pkl'

atmoFolder = './'
TJC.createAtmoPickle(atmoFolder, atmoFile, pickleName)
