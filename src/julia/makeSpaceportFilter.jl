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
spaceportLocationFile = abspath("./data/spaceportLocationsSpace2015.csv")

# Now move the pwd to where the outputs should go
cd(abspath("./outputs/SpaceportFilters/"))

using Points
include("googleEarth.jl")


function ReadSpaceportDict(spaceportLocationFile)
  colKeys = []
  spaceportDict = Dict()

  inFile = open(spaceportLocationFile, "r")
  for line in eachline(inFile)
    line = strip(line)        # Strip out any leading or trailing spaces and \n's
    if (length(line) == 0) || (line[1] == '#') || (line[1] == ',')   # skip empty lines and comments
      continue
    elseif line[1] == '@'
      # These are the keys to the columns.  Save them
      colKeys = split(line[2:end], ',')
      continue
    else
      # This should be valid data
#       colVals = split(line,',')
      colVals = [strip(elem) for elem in split(line,',')]  # strip out leading and trailing spaces for each element
      curDict = Dict(zip(colKeys[2:end],colVals[2:end]))
      spaceportDict[colVals[1]] = curDict
    end

  end
  close(inFile)

  spaceportDict
end


function WriteFacetFilter(LatLonFL, lowFL, hiFL, curFilter)
  startTime, endTime  = ("", "")     # Convert times to seconds from midnight.  Empty strings means time is unset

  outFileName = "SUA_" * curFilter * "Filter"
  outFile = open(outFileName, "w")

  # Make sure the shape closes
  isOpen = (LatLonFL[:,end] != LatLonFL[:,1])
  assert(!isOpen)

  numPoints = length(LatLonFL[1,:])
  # preamble = "{0}\n\n0\n{1} {2} {3} {4}\n0 0 0 1\n{5}\n".format(curSUA, lowFL*100, hiFL*100, startTime, endTime, len(lat)+isOpen)
  preamble = "$curFilter\n\n0\n$(lowFL*100) $(hiFL*100) $startTime $endTime\n0 0 0 1\n$(numPoints)\n"

  write(outFile, preamble)

  for ix in 1:(numPoints)
    (curLat, curLon, None) = LatLonFL[:,ix]
    write(outFile, "$curLat $(curLon + 360.) \n")
  end
  close(outFile)
end


function GetFilterLatLonAlt(curFilter, curDict, filterRadius::Float64)
  lat0 = float(curDict["lat"]) * pi/180.
  lon0 = float(curDict["lon"]) * pi/180.
  arm = filterRadius

  # Front Range
  # curFilter = "FrontRange"
  # lat0 = 39.785057 * pi/180.
  # lon0 = -104.524242 * pi/180.


  x0, y0 = projectIntoXY(lat0, lon0)
  r0 = [x0; y0; 0]

  numPoints = 36
  deltaAngle = 2*pi/numPoints
  XYZ = zeros(Float64, 3, numPoints+1)

  for ix in 1:(numPoints)
    # print("$ix \n")
    XYZ[:,ix] = [arm * cos(deltaAngle * ix); arm * sin(deltaAngle * ix); 0] + r0
  end
  XYZ[:,end] = XYZ[:,1]

  LatLonFL = XYZ2LatLonFL(XYZ)
end



#=
Now we can use the above functions to create filters for everything
=#

spaceportDict = ReadSpaceportDict(spaceportLocationFile)
storageDict = Dict()

filterRadius = 450.
lowFL               = 0
hiFL                = 999

for (curFilter, curDict) in spaceportDict
  println(curFilter)
  LatLonFL = GetFilterLatLonAlt(curFilter, curDict, filterRadius)
  WriteFacetFilter(LatLonFL, lowFL, hiFL, curFilter)
  storageDict[curFilter] = LatLonFL

end

# Now plot all of the filters
fileNameGE = "All_Filters.kml"
staticBoxes(fileNameGE, storageDict, 10.)



