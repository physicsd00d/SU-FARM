'''
Create a simple failure distribution

Let's say three main modes, one at each staging, say each stage is 60sec

'''
import numpy as np
from scipy.stats import beta
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sqrt
from copy import deepcopy

fig, ax = plt.subplots(1, 1)
# a, b = 2.30984964515, 0.62687954301
# mean, var, skew, kurt = beta.stats(a, b, moments='mvsk')

launchTimeEnd = 400
launchTime = np.array(range(launchTimeEnd + 1)) * 1.

# %% Failure on the pad
# % Use a beta distribution, specify the mean and variance, solve for params
mu = 30./launchTime[-1]
dev = 25./launchTime[-1]

# % Could add in some logic that if alpha is neg or nan, try a new guess
funA = lambda a : ((a**2 *(1/mu - 1))/( (a/mu)**2 * (a/mu + 1)) - dev**2)
alphaTemp = fsolve(funA,10)
betaTemp = alphaTemp*(1/mu - 1)

checkMean = alphaTemp/(alphaTemp + betaTemp)
checkStd = sqrt(  alphaTemp*betaTemp/ ((alphaTemp+betaTemp)**2 * (alphaTemp + betaTemp + 1)))
print "earlier mean {0} and std {1}".format(mu, dev)
print "earlier mean {0} and std {1}".format(checkMean, checkStd)

alphaFirst = alphaTemp
betaFirst = betaTemp

mu = 135./launchTime[-1]  # Assumes a 15 second coast
dev = 20./launchTime[-1]

funA = lambda a : ((a**2 *(1/mu - 1))/( (a/mu)**2 * (a/mu + 1)) - dev**2)
alphaTemp = fsolve(funA,50)
betaTemp = alphaTemp*(1/mu - 1)

checkMean = alphaTemp/(alphaTemp + betaTemp)
checkStd = sqrt(  alphaTemp*betaTemp/ ((alphaTemp+betaTemp)**2 * (alphaTemp + betaTemp + 1)))

alphaSecond = alphaTemp
betaSecond = betaTemp



mu = 200./launchTime[-1]  # Assumes a 15 second coast
dev = 20./launchTime[-1]

funA = lambda a : ((a**2 *(1/mu - 1))/( (a/mu)**2 * (a/mu + 1)) - dev**2)
alphaTemp = fsolve(funA,50)
betaTemp = alphaTemp*(1/mu - 1)

checkMean = alphaTemp/(alphaTemp + betaTemp)
checkStd = sqrt(  alphaTemp*betaTemp/ ((alphaTemp+betaTemp)**2 * (alphaTemp + betaTemp + 1)))

alphaThird  = alphaTemp
betaThird   = betaTemp


x = launchTime / launchTime[-1]
# fullPdf = beta.pdf(x, alphaPad, betaPad) + beta.pdf(x, alphaMaxQ, betaMaxQ)
fullCdf   = (1./3.)*beta.cdf(x, alphaFirst, betaFirst) \
          + (1./3.)*beta.cdf(x, alphaSecond, betaSecond) \
          + (1./3.)*beta.cdf(x, alphaThird, betaThird)

diffCdf = np.hstack((0, fullCdf[1:] - fullCdf[0:-1]))


# x = np.linspace(beta.ppf(0.01, a, b), beta.ppf(0.99, a, b), 100)




totalTimeSec = 865
# reentryTime = np.array(range(totalTimeSec - launchTimeEnd + 1)) * 1. + launchTimeEnd
reentryTime = np.array(range(totalTimeSec - launchTimeEnd + 1)) * 1.

# Explode it for the last two minutes to simulate reentry dispersion
# diffCdf = deepcopy(reentryTime) * 0.

mu = (reentryTime[-1] - 90.)/reentryTime[-1]  # Assumes a 15 second coast
dev = 40./reentryTime[-1]

funA = lambda a : ((a**2 *(1/mu - 1))/( (a/mu)**2 * (a/mu + 1)) - dev**2)
alphaTemp = fsolve(funA,50)
betaTemp = alphaTemp*(1/mu - 1)

checkMean = alphaTemp/(alphaTemp + betaTemp)
checkStd = sqrt(  alphaTemp*betaTemp/ ((alphaTemp+betaTemp)**2 * (alphaTemp + betaTemp + 1)))

alphaReentry  = alphaTemp
betaReentry   = betaTemp

x = reentryTime / reentryTime[-1]
fullCdfReentry   = beta.cdf(x, alphaReentry, betaReentry)

diffCdfReentry = np.hstack((0, fullCdfReentry[1:] - fullCdfReentry[0:-1]))

# Now stitch these into the launch stuff
totalTimeVec = np.hstack((launchTime[:-1], reentryTime + launchTimeEnd))
totalDiffCdf = np.hstack((diffCdf[:-1], diffCdfReentry))



ax.plot(totalTimeVec, totalDiffCdf, 'r-', lw=5, alpha=0.6, label='beta pdf')

plt.axis([0, totalTimeVec[-1], 0, max(totalDiffCdf)])

# Change names for output
time = totalTimeVec
diffCdf = totalDiffCdf








































plt.show()





# # % I have no reason to suspect any time of flight over another, so make it
# # % uniform.  Looking at the flown profile, it seems like there's action
# # % going on for the first 8.5 minutes, after which i assume if you haven't
# # % failed by then then you won't fail at all and will land smoothly.
#
# # For SS2, action happens for about 7 minutes, that's 420 seconds.
# lastTime = 175      # If this is not a multiple deltaTfail, you'll be in trouble
# time = np.linspace(0,lastTime,lastTime+1)
#
# diffCdf = np.ones_like(time) / time[-1]
# # uniformFail = ones(size(time)) / time(end);
# diffCdf[0] = 0.;    # We'll assume it doesn't fail at t=0, this also makes it sum to 1
# print sum(diffCdf)
#


output = open('failProfile.py', 'w')
output.write('import numpy as np\n')
output.write('failProfile = np.array([')
for val in diffCdf[:-1]:
    output.write('{0},\n'.format(val))
output.write('{0}])\n\n'.format(diffCdf[-1]))

output.write('failProfileSeconds = np.array([')
for val in time[:-1]:
    output.write('{0},\n'.format(val))
output.write('{0}])\n\n'.format(time[-1]))

output.close()


