# Set where we want the output files to default to
# cd("/Users/marian/Documents/Research/Prop3Dof/JuliaFiles/Results")

# /Volumes/Storage/Research/SU-FARM/outputs/SpaceportFilters

# Allow me to define modules and "using" them from this folder
# push!(LOAD_PATH, "/Users/marian/Documents/Research/Prop3Dof/JuliaFiles/")

# Juno doesn't ever start with julia running in the same directory as the file, so unfortunately must hardcode some paths.
basePath = "/Volumes/Storage/Research/SU-FARM/"
cd(basePath)

# Allow me to define modules and "using" them from the appropriate folder
push!(LOAD_PATH, abspath("./src/julia/"))

# Get the location of the input file
# suaFileName = abspath("./papers/Space2015/SpaceportFilters/SUA_AmericaFilter")
# suaFolder = "/Volumes/Main Mac/Users/marian/Documents/Research/Prop3Dof/CythonFiles/OtherPythonFiles/ProcessAirTOpHazardAreas/Airtop_2018_Feb2/"
# suaFileName = "SUA_OpeningScheme2025H_CE"

# suaFolder = basePath * "papers/TechCenterReport/HazardAreasSimulated/"
# outputFolder = suaFolder * "KMLs/"

suaFolder = basePath * "temp/SpaceportFilters/Space2015/"
outputFolder = suaFolder * "KMLs/"


# Now move the pwd to where the outputs should go
# cd(abspath("./temp/SpaceportFilters/"))

using Points
include("googleEarth.jl")

function GetSUALatLonAlt(curName, curDict)
  lat = curDict[:lat]
  lon = curDict[:lon]
  numPts = length(lat)
  LLFL      = zeros(Float64, 3, numPts)
#   LLFL[1,:] = lat * pi/180.
#   LLFL[2,:] = lon * pi/180. - π
  LLFL[1,:] = lat
  LLFL[2,:] = lon
  LLFL[3,:] = int(float(curDict[:hiAlt])/100.)

  # Return
  LLFL
end

function GetSUALatLonHiLo(curName, curDict)
  lat = curDict[:lat]
  lon = curDict[:lon]
  numPts = length(lat)
  LLFL      = zeros(Float64, 4, numPts)
#   LLFL[1,:] = lat * pi/180.
#   LLFL[2,:] = lon * pi/180. - π
  LLFL[1,:] = lat
  LLFL[2,:] = lon
  LLFL[3,:] = int(float(curDict[:lowAlt])/100.)
  LLFL[4,:] = int(float(curDict[:hiAlt])/100.)

  # Return
  LLFL
end




# function GetFilterLatLonAlt(curFilter, curDict, filterRadius::Float64)
#   lat0 = float(curDict["lat"]) * pi/180.
#   lon0 = float(curDict["lon"]) * pi/180.
#   arm = filterRadius

#   # Front Range
#   # curFilter = "FrontRange"
#   # lat0 = 39.785057 * pi/180.
#   # lon0 = -104.524242 * pi/180.


#   x0, y0 = projectIntoXY(lat0, lon0)
#   r0 = [x0; y0; 0]

#   numPoints = 36
#   deltaAngle = 2*pi/numPoints
#   XYZ = zeros(Float64, 3, numPoints+1)

#   for ix in 1:(numPoints)
#     # print("$ix \n")
#     XYZ[:,ix] = [arm * cos(deltaAngle * ix); arm * sin(deltaAngle * ix); 0] + r0
#   end
#   XYZ[:,end] = XYZ[:,1]

#   LatLonFL = XYZ2LatLonFL(XYZ)
# end

function ReadIndividualSUA(inFile, suaDict)
  ### ===== SUA Format ===== ###
  # SUA6_2                        # Name
  #
  # 0                             # Dunno, always zero
  # 30001 45002 690 750           # Low Alt (ft), Hi Alt, Start Time (sec since midnight), End Time
  # 0 0 0 1                       # Dunno, 0 0 0 1
  # 45                            # Number of points in the shape
  # 39.3925 255.799
  # 39.3914 255.797
  # 39.3902 255.794
  # 39.3902 255.794
  # ...                           # Chopping out some points for brevity
  # 39.3947 255.804
  # 39.3925 255.799               # SHAPE CLOSES
  # SUA6_3                        # Next name here!  Repeat...

  curName = strip(readline(inFile)) # First line should be the name
  readline(inFile)                  # Empty line
  readline(inFile)                  # A number, probably 0, probably the "probability" of the shape but it doesn't seem to matter
  altTimeLine = split(readline(inFile))  # Contains the low and high altitudes [ft] and start/end times since midnight [sec]
  readline(inFile)                  # Not Important
  numPts = int(readline(inFile))    # Number of points in the shape

  # Now read in all of those points
  latArray = zeros(numPts)
  lonArray = zeros(numPts)
  for ix in 1:numPts
    latArray[ix], lonArray[ix] = [float(obj) for obj in split(readline(inFile))]
  end

  #Finally parse out the altitudes and times
  (lowAlt, hiAlt) = altTimeLine[1:2]
  startTime, endTime = null, null
  if length(altTimeLine) == 4
    (startTime, endTime) = altTimeLine[3:4]
  end

  # Done with current shape
  suaDict[curName] = {:lowAlt=>lowAlt, :hiAlt=>hiAlt, :startTime=>startTime, :endTime=>endTime, :lat=>latArray, :lon=>lonArray}
end

readdir(suaFolder)
suaFileNameVec = Array(String,0)
for curFile in readdir(suaFolder)
  if curFile[1:4] == "SUA_"
    push!(suaFileNameVec, curFile)
  end
end

suaFileNameVec

for suaFileName in suaFileNameVec
  inFile = open(suaFolder*suaFileName, "r")
  suaDict = Dict()
  while !eof(inFile)
    ReadIndividualSUA(inFile, suaDict)
  end
  close(inFile)

  suaDict

  #=
  Now we can use the above functions to create filters for everything
  =#

  # spaceportDict = ReadSpaceportDict(spaceportLocationFile)
  storageDict = Dict()

  filterRadius       = 450.
  lowFL              = 0
  hiFL               = 999

  for (curName, curDict) in suaDict
    println(curName)
    LatLonHiLo = GetSUALatLonHiLo(curName, curDict)
    storageDict[curName] = LatLonHiLo

  #   WriteFacetFilter(LatLonFL, lowFL, hiFL, curFilter)
  #   LatLonFL = GetSUALatLonAlt(curName, curDict)
  #   storageDict[curName] = LatLonFL
  end

  # curKeys = collect(keys(storageDict))
  # storageDict[curKeys[1]][1:3,:]

  # makeLineString(storageDict[curKeys[1]][1:3,:], 10.)

  # Now plot all of the filters
  fileNameGE = outputFolder * suaFileName[5:end] * ".kml"
  staticBoxesV2(fileNameGE, storageDict, 1.)
end


