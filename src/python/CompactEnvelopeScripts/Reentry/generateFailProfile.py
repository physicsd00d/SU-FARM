'''
Create a simple failure distribution
'''
import numpy as np

# % I have no reason to suspect any time of flight over another, so make it
# % uniform.  Looking at the flown profile, it seems like there's action
# % going on for the first 8.5 minutes, after which i assume if you haven't
# % failed by then then you won't fail at all and will land smoothly.

# Actual trajectory here flies until 442, but lets just say that the last 22 seconds have zero probability of explosive failure
lastTime = 420      # If this is not a multiple deltaTfail, you'll be in trouble
time = np.linspace(0,lastTime,lastTime+1)

uniformFail = np.ones_like(time) / time[-1]
# uniformFail = ones(size(time)) / time(end);
uniformFail[0] = 0.;    # We'll assume it doesn't fail at t=0, this also makes it sum to 1
# sum(uniformFail)



output = open('failProfile.py', 'w')
output.write('import numpy as np\n')
output.write('failProfile = np.array([')
for val in uniformFail[:-1]:
    output.write('{0},\n'.format(val))
output.write('{0}])\n\n'.format(uniformFail[-1]))

output.write('failProfileSeconds = np.array([')
for val in time[:-1]:
    output.write('{0},\n'.format(val))
output.write('{0}])\n\n'.format(time[-1]))

output.close()