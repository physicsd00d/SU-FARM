module Points
#=
This is based very loosely on my c++ point class.  I realize now that the projection sucks
and will distort, but it's so fast that I don't care.

Defines a data structure of sorts called LLFL, Lat Lon FlightLevel
Each quantity is a row
What's the convention for Lat / Lon?
=#

export projectIntoXY, projectIntoDeg
export XYZ2LatLonFL, LatLonFL2XYZ, LatLonFL2ECEF
export convertDegMinSecToDecimal

refLat = 38 * (pi/180)
refLon = 0.
refRadius = 6378.1370 # Equator Km
m2ft = 3.28084
ft2m = 0.3048

ecc_Earth = 0.081819221456
R_equator = 6378.1370

function projectIntoDeg(x, y)
  gdlat = (y / refRadius) * (180 / pi)
  lon = refLon + (x / (refRadius * cos(refLat) )) * (180/pi)
  [gdlat, lon]
end

function projectIntoXY(gdLatIN, lonIN)
  # Units radians
  x = refRadius * (lonIN - refLon) * cos(refLat)
  y = refRadius * gdLatIN
  (x,y)
end


function XYZ2LatLonFL(xyz)
  # XYZ in kilometers, is matrix with each column a 3-vector
  LLFL = zeros(xyz)
  for IX in 1:length(xyz[1,:])
    (x, y, z) = xyz[:,IX]
    (lat, lon) = Points.projectIntoDeg(x, y)
    FL = round(z*10*m2ft)
    LLFL[:,IX] = [lat, lon, FL]
  end
  LLFL
end

function LatLonFL2XYZ(LLFL)
  # deg, deg, FL
  XYZ = zeros(LLFL)
  for IX in 1:length(LLFL[1,:])
    (lat, lon, FL) = LLFL[:,IX]
    (x, y) =  projectIntoXY(lat*pi/180, lon*pi/180)
    z = FL/(10*m2ft)
    XYZ[:,IX] = [x, y, z]
  end
  XYZ
end

#### From Python
# def convertDegreesToDegMinSec(val):
#     Deg = int(np.floor(val))
#     Min = int((val - Deg)*60.)
#     Sec = int((val - Deg - Min/60.)*3600.)
#     return [Deg, Min, Sec]

# def convertDegMinSecStringToDegrees(val):
#     if ('.' in val) or ('-' in val):
#         print '# This cannot accomodate decimals in the seconds or negative signs.'
#         raise RuntimeError

#     Seconds      = float(val[-2:])
#     Minutes      = float(val[-4:-2])
#     if len(val) == 6:
#         Degrees  = float(val[:2])
#     elif len(val) == 7:
#         Degrees  = float(val[:3])
#     else:
#         print "BAD LONGITUDE FORMAT: {0}".format(val)
#         raise RuntimeError

#     return Degrees + Minutes/60. + Seconds/3600.

function convertDegMinSecToDecimal(val::String)
  # Converted from Python, but have not yet tested it.  May not work.
  error("Converted from Python, but have not yet tested it.  May not work.")
#   if ('.' in val) or ('-' in val)
#     println("This cannot accomodate decimals in the seconds or negative signs.")
# #     raise RuntimeError
#     error()
#   end

#   Seconds      = float(val[-2:])
#   Minutes      = float(val[-4:-2])
#   if length(val) == 6
#       Degrees  = float(val[:2])
#     elif length(val) == 7
#       Degrees  = float(val[:3])
#   else:
#     print "BAD LONGITUDE FORMAT: {0}".format(val)
#     error()
# #     raise RuntimeError
#   end

  return Degrees + Minutes/60. + Seconds/3600.
end

function convertDegMinSecToDecimal(DMS::Array{Float64,1})
  @assert length(DMS) == 3
  return DMS[1] + DMS[2]/60. + DMS[3]/3600.
end



function LatLonFL2ECEF(LLFL)
  ECEF = zeros(LLFL)
  for ix in 1:length(LLFL[1,:])
    gdlat = LLFL[1,ix] * π/180.
    lon = LLFL[2,ix]   * π/180.
    alt = LLFL[3,ix] * ft2m / 10.  # Convert from FL into km
    Nlat = R_equator/sqrt(1-((ecc_Earth*sin(gdlat))^2))

    ECEF[1,ix] = (Nlat + alt)*cos(gdlat)*cos(lon);
    ECEF[2,ix] = (Nlat + alt)*cos(gdlat)*sin(lon);
    ECEF[3,ix] = (Nlat*(1-ecc_Earth*ecc_Earth) + alt)*sin(gdlat);
  end

  ECEF
end




#### From Python
# def convertDecimalStringToDegMinSecString(LatLonStr):
#     curLat, curLon = LatLonStr

#     if float(curLat) < 0:
#         print "ERROR: CANNOT HAVE NEGATIVE LATITIUDE\n"
#         sys.exit()
#     elif float(curLon) > 0:
#         print "ERROR: CANNOT HAVE POSITIVE LONGITUDE\n"
#         sys.exit()

#     # Do the Latitude first
#     val             = abs(float(curLat))
#     degreesInt      = int(np.floor(val))
#     minutesFloat    = (val - degreesInt)*60
#     minutesInt      = int(np.floor(minutesFloat))
#     secondsInt      = int(np.floor((minutesFloat - minutesInt)*60.))
#     LatDDMMSS       = '{0}{1}{2}N'.format(degreesInt, minutesInt, secondsInt)

#     # Do the Longitude first
#     val             = abs(float(curLon))
#     degreesInt      = int(np.floor(val))
#     minutesFloat    = (val - degreesInt)*60
#     minutesInt      = int(np.floor(minutesFloat))
#     secondsInt      = int(np.floor((minutesFloat - minutesInt)*60.))
#     LonDDMMSS       = '{0}{1}{2}W'.format(degreesInt, minutesInt, secondsInt)

#     return [LatDDMMSS, LonDDMMSS]

end

