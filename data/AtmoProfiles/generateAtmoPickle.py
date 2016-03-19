import sys

# Error check
if len(sys.argv) != 2:
    print "ERROR: Only one argument allowed.  Must pass in .txt file containing atmo data"
    sys.exit()

atmoFile = sys.argv[1]
baseName = atmoFile.split('.')[0]
pickleName = baseName + ".pkl"

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

# atmoFile = 'WestTexas.txt'
# pickleName = 'WestTexas.pkl'

atmoFolder = './'
TJC.createAtmoPickle(atmoFolder, atmoFile, pickleName)
