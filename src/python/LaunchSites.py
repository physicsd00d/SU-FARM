# Launch Sites
Cape    = 'Cape'
MARS    = 'MARS'
Cecil   = 'Cecil'
OK      = 'OK'
Brown   = 'Brown'
CornNew = 'CornNew'

Vanden  = 'Vanden'       #http://en.wikipedia.org/wiki/Vandenberg_Air_Force_Base#Launch_sites
Odyssey = 'Odyssey'
Atoll   = 'Atoll'
Wallops = 'Wallops'
Kodiak  = 'Kodiak'
Mojave  = 'Mojave'
America = 'America'
McGreg  = 'McGreg'
Poker   = 'Poker'
FrntRnge = 'FrntRnge'
Midland = 'Midland' #The 374 acres is located north of the airport near Farm-to-Market Road 1788 between Highways 191 and 158. The deal is not to exceed $4.01 million.
Georgia = 'Georgia'

WSands  = 'WSands'
WSandsN = 'WSandsN' #north

# Custom
PegMARS = 'PegMARS'
PegMARS2= 'PegMARS2'
Titus   = 'Titus'
Houston = 'Houston'
Hawaii  = 'Hawaii'
PacLA   = 'PacLA'
PegVAFB = 'PegVAFB'

siteDict = {Cape    : {'lat' : 28.445455, 'lon' : -80.564865, 'altkm' : 3e-3, 'timezone' : 'EST'},  #Cape Canveral
            MARS    : {'lat' : 37.8338  , 'lon' : -75.4882  , 'altkm' : 0.0, 'timezone' : 'EST'},   #Mid-Atlantic Regional Spaceport
            Cecil   : {'lat' : 30.218611, 'lon' : -81.876667, 'altkm' : 0.0, 'timezone' : 'EST'},   #Cecil Field Spaceport FL
            OK      : {'lat' : 35.343559, 'lon' : -99.208014, 'altkm' : 0.0, 'timezone' : 'CST'},   #Oklahoma Spaceport
            Brown   : {'lat' : 26.150661, 'lon' : -97.297082, 'altkm' : 0.0, 'timezone' : 'CST'},   #Cameron County TX (SpaceX)
            # Corn    : {'lat' : 31.062125, 'lon' : -104.78185, 'altkm' : 0.0, 'timezone' : 'CST'},   #Corn Ranch TX (in Van Horn) Surely Wrong Coords
            CornNew : {'lat' : 31.423333, 'lon' : -104.75889, 'altkm' : 0.0, 'timezone' : 'CST'},   #Corn Ranch TX (in Van Horn) Surely Wrong Coords
            Vanden  : {'lat' : 34.581111, 'lon' : -120.6275 , 'altkm' : 0.0, 'timezone' : 'PST'},   #SLC6, This is where DeltaIV launches, chekc wikipedia for other pads
            Odyssey : {'lat' : 0.0      , 'lon' : -154.0    , 'altkm' : 0.0, 'timezone' : 'HST'},   #Hawaii Standard Time (no DST)
            Atoll   : {'lat' : 8.716667 , 'lon' : 167.733333, 'altkm' : 0.0, 'timezone' : 'MST'},   #Marshal Islands Time (no DST)
            Wallops : {'lat' : 37.940194, 'lon' : -75.466389, 'altkm' : 0.0, 'timezone' : 'EST'},
            Kodiak  : {'lat' : 57.4353  , 'lon' : -152.3393 -0.15 , 'altkm' : 0.0, 'timezone' : 'AKST'},
            Mojave  : {'lat' : 35.059444, 'lon' : -118.15166, 'altkm' : 0.0, 'timezone' : 'PST'},
            America : {'lat' : 32.990278, 'lon' : -106.96972, 'altkm' : 0.0, 'timezone' : 'MST'},#32.990278, -106.969722
            McGreg  : {'lat' : 31.388104, 'lon' : -97.468994, 'altkm' : 0.0, 'timezone' : 'CST'},#31.38810463880516,-97.46899490911483
            Poker   : {'lat' : 65.117893, 'lon' : -147.43328, 'altkm' : 0.0, 'timezone' : 'AKST'},#Lat: 65.1178935 Lon: -147.4332802
            WSands  : {'lat' : 32.915999, 'lon' : -106.3666 , 'altkm' : 0.0, 'timezone' : 'MST'},#32.563056, -106.57  not actually sure where in WSands this is pointing
            WSandsN  : {'lat' : 33.910837, 'lon' : -106.637063 , 'altkm' : 0.0, 'timezone' : 'MST'},#32.563056, -106.57  not actually sure where in WSands this is pointing
            # WSands  : {'lat' : 32.563056, 'lon' : -106.57   , 'altkm' : 0.0, 'timezone' : 'MST'},#32.563056, -106.57  not actually sure where in WSands this is pointing
            # FrntRnge: {'lat' : 39.793949, 'lon' : -104.56445, 'altkm' : 0.0, 'timezone' : 'MST'},# West end of main runway
            FrntRnge: {'lat' : 39.785057, 'lon' : -104.524242, 'altkm' : 0.0, 'timezone' : 'MST'},# North end of eastern runway?
            # Midland : {'lat' : 32.018138, 'lon' : -102.23713, 'altkm' : 0.0, 'timezone' : 'CST'},#32.018138,-102.237139, SURELEY WRONG, in general spot of proposed site
            Midland : {'lat' : 31.942941, 'lon' : -102.203849, 'altkm' : 0.0, 'timezone' : 'CST'},#32.018138,-102.237139, SURELEY WRONG, in general spot of proposed site
            Georgia : {'lat' : 30.944654 - 0.03, 'lon' : -81.480956 - 0.12 , 'altkm' : 0.0, 'timezone' : 'EST'},   #Georgia Spaceport, gleaned from picture of posterboard haha

            # Custom ones
            PegMARS : {'lat' : 38.201279  , 'lon' : -73.353819  , 'altkm' : 0.0, 'timezone' : 'EST'},   #Pegasus launching over the atlantic from MARS
            PegMARS2: {'lat' : 37.727847  , 'lon' : -73.653940  , 'altkm' : 0.0, 'timezone' : 'EST'},   #Pegasus launching over the atlantic from MARS
            PegVAFB : {'lat' : 35.902172  , 'lon' : -123.055242 , 'altkm' : 0.0, 'timezone' : 'EST'},   #Pegasus launching over the atlantic from MARS
            Titus   : {'lat' : 28.626351  , 'lon' : -80.701945  , 'altkm' : 0.0, 'timezone' : 'EST'},   #I don't think this spaceport is going to happen
                                                                                                        #   so I'm placing this at the shuttle landing
                                                                                                        #   runway at VERY nearby KSC.
            Houston : {'lat' : 29.597021  , 'lon' : -95.163903  , 'altkm' : 0.0, 'timezone' : 'CST'},   # Formerly known as Ellington spaceport, just south of Houston
            Hawaii  : {'lat' : 19.725908  , 'lon' : -156.046778  , 'altkm' : 0.0, 'timezone' : 'HST'},  # South end of Kona Airport runway
            PacLA   : {'lat' : 29.924865  , 'lon' : -119.325785  , 'altkm' : 0.0, 'timezone' : 'PST'},  # Dragon reentering in the pacific by LA

            }


observesDST = ['EST', 'CST', 'MST', 'PST', 'AKST']
# Hawaii does not observe


# print "hey"
# outFile = open("site.csv",'w')
#
#
#
# outFile.write("hey\n")
# outFile.close()


