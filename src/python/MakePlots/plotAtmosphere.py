'''
#
This file creates an envelope for a SpaceX Falcon9 launch
#
'''

import os
import sys

# Want to import some things that are general to all missions
curFilePath = os.path.dirname(os.path.abspath(__file__)) + "/"      # Reference everything from the location of the file
sys.path.insert(0, os.path.abspath(curFilePath+'../'))              # Back up one so we can import CommonThemes
# import CommonThemes as ct

# Find the path of the current file, then Point to the root of the package so I can run this script from anywhere
rootDir =   os.path.abspath(curFilePath + "../../../") + "/"
outputDir = rootDir + "outputs/" # Where to store results, gitignored
# tempDir =   rootDir + ct.tempFolderName   # temp files here, gitignored
debrisPath = rootDir + "src/python/packages/DebrisCatalogs/"

from Simulation import TJC

curMission = dict()
curMission['atmospherePickle']  = rootDir + "data/AtmoProfiles/Cape.pkl"

curMission['numTrajSamples'] = 0
curMission['numWindSamples'] = 100   # Best results if this is a multiple of the number of nodes you're running on.

atmStorage, stateVecStorage, thetagStorage, tfailStorage = \
                        TJC.GenerateWindTrajProfiles(curMission, curMission['numTrajSamples'], curMission['numWindSamples'])






import numpy as np
from bokeh.plotting import figure, output_notebook, output_file, show
from bokeh.palettes import Spectral11
mypalette=Spectral11*10


# output to static HTML file
output_file("lines.html", title="line plot example")

p = figure(plot_width=1000, plot_height=1000, y_range=(0,60), x_range=(0,40))
p.title = "Winds"
p.title_text_font_size = "40pt"

for ix, curProfile in enumerate(atmStorage):
    [altitudeList,densityList,uList,vList,wList] = curProfile

    windMag = np.sqrt(np.multiply(uList, uList) + np.multiply(vList, vList) + np.multiply(wList, wList)) # get the magnitude (as matrix)
    windMag = [x for x in windMag.flat] # Then flatten it out into 1d for bokeh

    p.line(windMag, [x/1000. for x in altitudeList.flat], line_color=mypalette[ix], line_width=3)


p.xaxis.axis_label="Magnitude [m/s]"
p.yaxis.axis_label="Altitude [km]"

p.xaxis.axis_label_text_font_size = "30pt"
p.yaxis.axis_label_text_font_size = "30pt"
p.xaxis.major_label_text_font_size = "30pt"
p.yaxis.major_label_text_font_size = "30pt"

show(p)




# output to static HTML file
output_file("density.html", title="line plot example")

p = figure(plot_width=1000, plot_height=1000, y_range=(0,30), x_range=(0,1.2))
p.title = "Atmospheric Density"
p.title_text_font_size = "40pt"

for ix, curProfile in enumerate(atmStorage):
    [altitudeList,densityList,uList,vList,wList] = curProfile
    densityList = [x for x in densityList.flat] # Then flatten it out into 1d for bokeh

    p.line(densityList, [x/1000. for x in altitudeList.flat], line_color=mypalette[ix], line_width=3)

p.xaxis.axis_label="Magnitude [kg/m^2]"
p.yaxis.axis_label="Altitude [km]"

p.xaxis.axis_label_text_font_size = "30pt"
p.yaxis.axis_label_text_font_size = "30pt"
p.xaxis.major_label_text_font_size = "30pt"
p.yaxis.major_label_text_font_size = "30pt"

show(p)

# Make the titles bigger and then these are probably fine.














