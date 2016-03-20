'''
Create a simple failure distribution

For suborbital vehicles, i treat the launch and reentry as two separate missions.

'''
import numpy as np
from scipy.stats import beta
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sqrt

fig, ax = plt.subplots(1, 1)
# a, b = 2.30984964515, 0.62687954301
# mean, var, skew, kurt = beta.stats(a, b, moments='mvsk')

secondsUntilStaging = 420
time    = np.array(range(secondsUntilStaging + 1)) * 1.
# diffCdf = np.array(range(secondsUntilStaging + 1)) * 0.
#
# # Let's make it uniform over first 90 seconds and last 90 seconds
# secondsUniform = 90
# diffCdf[1:secondsUniform+1] = 1./secondsUniform
# diffCdf[-secondsUniform:]   = 1./secondsUniform
#
# ax.plot(time, diffCdf, 'r-', lw=5, alpha=0.6, label='beta pdf')
#
# plt.axis([0, time[-1], 0, max(diffCdf)])
# plt.show()



# %% Failure on the pad
# % Use a beta distribution, specify the mean and variance, solve for params
mu = 50./time[-1]
dev = 30./time[-1]

# % Could add in some logic that if alpha is neg or nan, try a new guess
funA = lambda a : ((a**2 *(1/mu - 1))/( (a/mu)**2 * (a/mu + 1)) - dev**2)
alphaTemp = fsolve(funA,10)
betaTemp = alphaTemp*(1/mu - 1)

checkMean = alphaTemp/(alphaTemp + betaTemp)
checkStd = sqrt(  alphaTemp*betaTemp/ ((alphaTemp+betaTemp)**2 * (alphaTemp + betaTemp + 1)))
print "earlier mean {0} and std {1}".format(mu, dev)
print "earlier mean {0} and std {1}".format(checkMean, checkStd)

alphaPad = alphaTemp
betaPad = betaTemp

mu = (time[-1]-90.)/time[-1]
dev = 40./time[-1]

funA = lambda a : ((a**2 *(1/mu - 1))/( (a/mu)**2 * (a/mu + 1)) - dev**2)
alphaTemp = fsolve(funA,50)
betaTemp = alphaTemp*(1/mu - 1)

checkMean = alphaTemp/(alphaTemp + betaTemp)
checkStd = sqrt(  alphaTemp*betaTemp/ ((alphaTemp+betaTemp)**2 * (alphaTemp + betaTemp + 1)))

alphaMaxQ = alphaTemp
betaMaxQ = betaTemp
# initialFail = betapdf(time/170,alpha,beta);

x = time / time[-1]
fullPdf = beta.pdf(x, alphaPad, betaPad) + beta.pdf(x, alphaMaxQ, betaMaxQ)
fullCdf = beta.cdf(x, alphaPad, betaPad) + beta.cdf(x, alphaMaxQ, betaMaxQ)

diffCdf = np.hstack((0, fullCdf[1:] - fullCdf[0:-1]))

# x = np.linspace(beta.ppf(0.01, a, b), beta.ppf(0.99, a, b), 100)
ax.plot(time, diffCdf, 'r-', lw=5, alpha=0.6, label='beta pdf')

plt.axis([0, time[-1], 0, max(diffCdf)])
plt.show()


# sys.exit()
#
#
#
#
#
# # # % I have no reason to suspect any time of flight over another, so make it
# # # % uniform.  Looking at the flown profile, it seems like there's action
# # # % going on for the first 8.5 minutes, after which i assume if you haven't
# # # % failed by then then you won't fail at all and will land smoothly.
# #
# # # For SS2, action happens for about 7 minutes, that's 420 seconds.
# # lastTime = 175      # If this is not a multiple deltaTfail, you'll be in trouble
# # time = np.linspace(0,lastTime,lastTime+1)
# #
# # diffCdf = np.ones_like(time) / time[-1]
# # # uniformFail = ones(size(time)) / time(end);
# # diffCdf[0] = 0.;    # We'll assume it doesn't fail at t=0, this also makes it sum to 1
# # print sum(diffCdf)
# #


# output = open('failProfile.py', 'w')
# output.write('import numpy as np\n')
# output.write('failProfile = np.array([')
# for val in diffCdf[:-1]:
#     output.write('{0},\n'.format(val))
# output.write('{0}])\n\n'.format(diffCdf[-1]))
#
# output.write('failProfileSeconds = np.array([')
# for val in time[:-1]:
#     output.write('{0},\n'.format(val))
# output.write('{0}])\n\n'.format(time[-1]))
#
# output.close()


