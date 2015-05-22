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

  # Return the name of the file
  outFileName
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

LLFL = storageDict["America"]

ECEF = LatLonFL2ECEF(LLFL)

# //ColVec Trajectory::ECEF_To_Geodetic(ColVec Vec) {
# //	//function [gdlat lon alt] = ECEF_To_Geodetic(x, y, z)
# //
# //	//Unpack the ECI vector
# //	double x = Vec(0);
# //	double y = Vec(1);
# //	double z = Vec(2);
# //
# //	//Can pass off to other implementation from here
# //
# //	double rdelta = sqrt(x*x + y*y);
# //	double gdlat = asin(z/sqrt(x*x + y*y + z*z));
# //
# //	// Set the tolerances
# //	double tol = 1e-10;
# //	double err = 50;
# //
# //	double lon;
# //	// I do not remember why this statement is in here.  Something to do with polar orbits it seems
# //	if ((abs(x) <= tol) && (abs(y) <= tol)) {
# //        //		cout << "WOAH!  If you wound up here, then you need to check that everything is okay!" << endl;
# //        //		cout << Vec << endl;
# //		lon = 0; }
# //	else {
# //		lon = atan2(y,x); }   //[deg]
# //
# //
# //	// Iterate to find coordinates
# //	double Nlat;
# //	while( err > tol ){
# //		Nlat = R_equator/sqrt(1-pow(ecc_Earth*sin(gdlat),2));
# //
# //		double gdlat_new = atan((z+Nlat*sin(gdlat)*ecc_Earth*ecc_Earth)/rdelta);
# //
# //		err = abs(gdlat_new - gdlat);
# //		gdlat = gdlat_new;
# //	}
# //
# //	double alt = rdelta/cos(gdlat) - Nlat;
# //
# //	//Pack up the geodetic return vector
# //	ColVec Geodetic(3,1);
# //	Geodetic(0) = gdlat;
# //	Geodetic(1) = lon;
# //	Geodetic(2) = alt;
# //
# //	return Geodetic;
# //}

# Now plot all of the filters
fileNameGE = "All_Filters.kml"
staticBoxes(fileNameGE, storageDict, 10.)






#### Ghetto filter from googleEarth
# GE_string = "-80.48057912537114,26.797443341603,0 -68.86918135730895,35.25457186422343,0 -68.77681870356666,38.25207465413916,0 -71.98312248762959,38.43588281005918,0 -82.34543870737839,28.84395473405241,0 -80.48057912537114,26.797443341603,0"
# GE_string = "-78.00011383794438,37.07946251419159,0 -62.48771305452573,24.23587499232429,0 -58.2448004812892,23.01069864179952,0 -58.26808151047837,26.89621326772835,0 -74.03127049332808,41.04363142366995,0 -78.00011383794438,37.07946251419159,0"
# GE_string = "-72.03216948087703,30.0690786195643,0 -66.18236761588739,24.50791242601199,0 -60.98462632682065,19.57724920013325,0 -52.1364372666887,17.48026404061209,0 -52.87689742292973,23.99766500553687,0 -57.80679469977535,28.88277704380049,0 -65.34175147980434,34.93246663710243,0 -71.94532532654314,39.86817690017847,0 -74.79506031470945,40.3430041108391,0 -77.57217759693916,37.61494410717135,0 -77.27316643677514,35.14037113685587,0 -72.03216948087703,30.0690786195643,0"
# GE_string = "-126.8656063212387,29.91551027948017,0 -128.341950487521,23.25105110503754,0 -129.5814131601475,16.82834379634327,0 -129.8318054463546,10.82932632963048,0 -126.3672915565276,4.722845952146105,0 -120.7243752082968,10.32050228015308,0 -118.9487314646936,16.29722081144956,0 -119.0934186377551,22.21502201234596,0 -118.5558105805195,27.85193540410579,0 -117.9024310410929,33.74146257803556,0 -119.7584896519146,36.01322170701136,0 -122.5519135854599,36.78232055485513,0 -125.4812634705911,35.18492369668617,0 -126.8656063212387,29.91551027948017,0 "

# curName = "PacificReentryLA"
# GE_string = "-125.8111463961022,30.26158040441436,0 -121.5437802182309,26.92137351270492,0 -115.7083759852632,26.19312045186264,0 -115.3428905554021,30.48607119842432,0 -119.197289114417,34.45645285713426,0 -123.3139710561981,37.16253740141109,0 -127.6121879862085,36.7913907103195,0 -128.6559040335125,32.99456941474899,0 -125.8111463961022,30.26158040441436,0"

# curName = "VafbToOrbit225"
# GE_string = "-117.7951706849269,33.73060487147438,0 -119.959565837485,36.16487205303402,0 -122.6328730721677,38.39468768917248,0 -126.0834388651855,37.57386225734376,0 -129.276657821544,31.71867675812874,0 -131.4392472540116,22.38665536172382,0 -131.6091618587042,12.61939397314826,0 -126.2823447267405,4.596714835590557,0 -120.616013833937,9.943151217949483,0 -118.7540503968391,17.3731816880462,0 -117.8275782871151,26.46530578249043,0 -117.7951706849269,33.73060487147438,0"

# curName = "KodiakOrbit120"
# GE_string = "-149.4266507143721,52.2269670958156,0 -144.317828938581,48.66497064025472,0 -139.4435509055136,45.18836983176687,0 -135.0758970599114,42.1061670085459,0 -130.4859239409427,41.78749997179408,0 -129.7837317205255,44.25682396000237,0 -132.4448692886416,48.32156606330273,0 -137.0346553788827,52.22576984186129,0 -142.1894974953659,55.48212860434065,0 -149.215306148733,58.82774467841652,0 -156.1982115684825,58.79805185038602,0 -154.7579337409996,56.01026480034665,0 -149.4266507143721,52.2269670958156,0"

curName = "GeorgiaOrbit90"
GE_string = "-83.57893146201205,31.44109155798375,0 -83.38661310464288,29.42694979137307,0 -79.0241719871466,27.59837472561557,0 -75.94708937128503,25.4097181940247,0 -71.87727255364804,23.35913663811561,0 -67.80012373144125,21.62714447159237,0 -63.6775609401032,22.29468811996032,0 -65.80508576513188,25.86959605202331,0 -67.43450687345323,26.92003619739004,0 -72.68409728689358,29.17036223725201,0 -77.04890312755295,30.91414355048719,0 -81.94140409654086,32.58444667632017,0 -83.57893146201205,31.44109155798375,0"

numPts = length(split(GE_string))  # shape is already closed
LLFL = zeros(Float64, 3, numPts)
for (ix, curLine) in enumerate(split(GE_string))
  lon, lat, alt = [float(obj) for obj in split(strip(curLine), ',')]
  LLFL[:,ix] = [lat, lon + 360, 999]
end

LLFL
outfileName = WriteFacetFilter(LLFL, 0, 999, curName)

outputDir = "/Volumes/Storage/Research/SU-FARM/temp/SpaceportFilters/Space2015/"
mv(outfileName, outputDir*outfileName)

