# Falcon9 From Environmental Impact
import numpy as np
import sys
from DebrisCatalogs import DebrisCatalogBuilder as dcb

ft2m = 0.3048
sqft2sqm = ft2m**2
lbm2kg = 0.453592


numPieces = [1,9,91,8,175,4,1406,8,977,1,2,1,126,13,12,526,1,2,16,4,1,1081,1]

# lbs mass
wSingle = [7047,
131.8,
[1.2, 213.6],
[10, 175],
[0.5, 32.5],
38,
[0.02, 110.2],
70.3,
[0.3, 40],
840,
30,
653.7,
[0.5, 50.0],
[16, 417.1],
[12.0, 17.2],
[0.02, 110.2],
1124.4,
22.5,
[9.3, 36.3],
70.3,
200,
[0.4, 94.1],
17746]

# psf
beta = [406,
41.7,
[6.7, 89.7],
[11.4, 15.1],
[10.9, 16.2],
7.6,
[7.1, 14.3],
52,
[5.1, 7.4],
75.9,
74,
181.6,
[30.*sqft2sqm/lbm2kg, 60.*sqft2sqm/lbm2kg],     # was [7.1, 7.4] but too narrow for ranges of mass and area
[30.2, 89.7],
[13.0, 13.2] ,
[6.9, 14.2],
52.1,
15.1,
[32.2, 89.4],
52,
19.2,
[6.4, 56.7],
125.2]

# ft/sec
vImpart = [107,
52,
[48.4, 300],
[107, 624],
[58.9, 194.5],
244.6,
[29.6, 1384],
15.3,
[116.9, 770.5],
20.5,
18,
61.2,
[145.2, 365.4],
[14.3, 318.6],
[347.8, 723.2],
[176, 868],
73.7,
106.1,
[54.3, 72.6],
141.5,
173.1,
[13.2, 330.8],
27.5]

# ft^2
aRef = [134.1,
6.87,
[0.4, 5.0],
[2.4, 61.2],
[0.12, 7.7],
13.1,
[0.01, 31.4],
3.1,
[0.16, 21.1],
9,
0.48,
8,
[0.28, 27.6],
[0.5, 4.3],
[3.7, 5.2],
[0.01, 31.6],
48,
3.46,
[0.18, 0.48],
3.1,
12.3,
[0.19, 52.0],
315]

CD = np.ones_like(numPieces)*1.
CL = np.ones_like(numPieces)*0.1

# If the failing ballistic coefficients are falling in a tight range
#   Then need to change CD to adjust.  If in a wide range, then just
#   change the range of allowable beta above.
CD[13] = 10.
CD[14] = 1.29
CD[18] = 2.5
# CD[18] = 2.5

def multiplyList(curList, k):
    # numElem = len(curList)
    ans = []
    for elem in iter(curList):
        if not isinstance(elem, list):
            elem = elem*k
        elif len(elem) == 2:
            elem[0] = elem[0]*k
            elem[1] = elem[1]*k
        else:
            print "ERROR!  Debris catalog is messed up 978234"
            sys.exit()
        ans.append(elem)
    return ans

print wSingle
wConvert = multiplyList(wSingle, 5.)
print wConvert


# Convert the units to SI
wSingle = multiplyList(wSingle, lbm2kg)
beta = multiplyList(beta, lbm2kg/sqft2sqm)
vImpart = multiplyList(vImpart, ft2m)
aRef = multiplyList(aRef, sqft2sqm)

# Find the total weight of each category
wTotal = []
for ix in range(len(numPieces)):
    wTotal.append(np.mean(wSingle[ix]) * numPieces[ix] )

# Now package up so we can print the catalog
# Dict[name] = dict() of properties
catalog = []
for ix in range(len(numPieces)):
    subDict = { 'numPieces':numPieces[ix], \
                    'weightIndividual': wSingle[ix], \
                    'weightTotal': wTotal[ix],\
                    'aref': aRef[ix],\
                    'beta': beta[ix],\
                    'impulseVel': vImpart[ix],\
                    'CD': CD[ix],\
                    'CL': CL[ix],\
                    'name': str(ix)}
    catalog.append(subDict)

# numPieces = str(int(np.mean(cur['numPieces'])))
# # If it's not a range, then just output one number
# weightIndividual = collapseRange(cur['weightIndividual'])
# # If it's a range, output the mean.
# weightTotal = '{0:7.5}'.format(np.mean(cur['weightTotal']))
# aref = collapseRange(cur['aref'])
# beta = collapseRange(cur['beta'])
# impulseVel = collapseRange(cur['impulseVel'])
# CD = cur['CD']
# CL = cur['CL']

fileName="Falcon9_FTS.txt"
dcb.makeOutputFile(catalog, fileName)


# # Print the catalog nicely
# fileName="Falcon9_FTS.txt"
# outFile = open(fileName, 'w')

# outFile.write('units SI\nbegin\nvelocity  666.0\ntime   10.0\n')

# outString = '#{0:14}{1:25}{2:15}{3:25}{4:25}{5:20}{6:10}{7:10}'.\
#     format('numPieces', 'weightIndividual', 'weightTotal', 'aref', 'beta', 'impulseVel', 'CD', 'CL')
# print outString
# outFile.write(outString + '\n')

# for ix in range(len(numPieces)):
#     outString = '{0:<14}{1:<25.3}{2:<15}{3:<25}{4:<25}{5:<20}{6:<10}{7:<10}\n'.\
#     format(numPieces[ix], wSingle[ix], wTotal[ix], aRef[ix], beta[ix], vImpart[ix], '1', '0.1')
#     outFile.write(outString)

# outFile.write('end')
# outFile.close()








