

# pickle.dump([altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist],output,2)

import pickle
profiles = pickle.load(open('WallopsOverOcean_BIG.pkl','rb'))

nlist = profiles[-1]  # original number of points

# Want to scale this down to something small and manageable so it doesn't eat all the machine's memory
targetN = 500   # Not trying to be precise, just get around this many points
delta = nlist/targetN  # int / int = automatically rounded to int.
targetIX = range(0,nlist,delta)
targetIX.append(nlist-1) # Throw in final point too just to make sure we cover the same full range

newProfile = []
for curArray in profiles[:-1]:
    newProfile.append(curArray[targetIX])
newProfile.append(len(targetIX))

output = open("WallopsOverOcean.pkl", 'wb')
pickle.dump(newProfile,output,2)
output.close()







