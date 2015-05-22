ft2m = 0.3048

println("HEY IM HERE!!!")

function makePointString(curLat, curLon, curAlt)
  "        <Point>\n" *
  "          <altitudeMode>relativeToGround</altitudeMode>\n" *
  "          <coordinates>$curLon,$curLat,$curAlt</coordinates>\n" *
  "        </Point>\n"
end

function makeTimeString(curTime, deltaTsec)
  #curTime is sec since midnight UTC
  seconds = int(curTime % 60)
  minutes = int(floor(((curTime) % 3600) / 60))
  hours = int(floor(curTime / 3600.))
  beginString = @sprintf("%02d\:%02d\:%02d",hours, minutes, seconds)

  nextTime = curTime + deltaTsec
  seconds = int(nextTime % 60)
  minutes = int(floor(((nextTime) % 3600) / 60))
  hours = int(floor(nextTime / 3600.))
  endString = @sprintf("%02d\:%02d\:%02d",hours, minutes, seconds)


  "       <TimeSpan>\n" *
  "         <begin>2014-01-01T" * beginString * "Z</begin>\n" *
  "         <end>2014-01-01T" * endString * "Z</end>\n" *
  "       </TimeSpan>\n"
end

function makePlacemarkString(pointStringVec, timeString = "", placemarkName = "")
  collapsedPointString = join(pointStringVec)

  # If a name is provided, use it
  if length(placemarkName) > 0
    placemarkName = "<name>"*placemarkName*"</name>\n"
  end

  "       <Placemark>\n" * placemarkName * "$timeString" *
  "       <styleUrl>#rocketStyleID</styleUrl>\n" *
  "       <MultiGeometry>\n" * "$collapsedPointString" *
  "       </MultiGeometry>\n" *
  "       </Placemark>\n"
end

function googleEarthAircraftAnimation(fileNameGE, acDict, deltaTsec)
  outFileGE = open(fileNameGE, "w")

  if deltaTsec != 60
    error("This function was written with deltaT = 60s and acDict hashed in minutes\n")
  end

  preamble = ["""<?xml version="1.0" encoding="UTF-8"?>\n""",
              """<kml xmlns="http://www.opengis.net/kml/2.2">\n""",
              """  <Document>\n""",
              """    <Style id="styleRed">\n""",
              """      <LineStyle>\n""",
              """        <color>ff0000ff</color>\n""",
              """        <colorMode>normal</colorMode>\n""",
              """        <width>4</width>\n""",
              """      </LineStyle>\n""",
              """    </Style>\n""",
              """    <Style id="styleGreen">\n""",
              """      <LineStyle>\n""",
              """        <color>ff00ff00</color>\n""",
              """        <colorMode>normal</colorMode>\n""",
              """        <width>4</width>\n""",
              """      </LineStyle>\n""",
              """    </Style>\n""",
              """ <name>JuliaAC-Path</name>\n\n"""]

  for line in preamble
      write(outFileGE, line)
  end

  orderedTimes = sort(collect(keys(acDict)))
  for curTime in orderedTimes
    # Make the Google Earth timeString
    timeString = makeTimeString(curTime * deltaTsec, deltaTsec)

    # Run through all AC at this time and make their point strings
    numAcHere = length(acDict[curTime])
    pointStringVec = Array(ASCIIString, numAcHere)
    curIX = 1
    for (curCallsignKey, (curLat,curLon,curFL)) in acDict[curTime]
      curLon = curLon - 360.
      curAlt = round(curFL * 100 * ft2m)

      pointStringVec[curIX] = makePointString(curLat, curLon, curAlt)
      curIX += 1
    end

    # Then dump everything to the output file
    write(outFileGE, makePlacemarkString(pointStringVec, timeString))

  end

  finale = "  </Document>\n" *
           "</kml>"
  write(outFileGE, finale)

  close(outFileGE)
end


#       <Placemark>
#       <TimeSpan>
#         <begin>2014-01-01T02:15:55Z</begin>
#         <end>2014-01-01T02:15:56Z</end>
#       </TimeSpan>
#       <styleUrl>#rocketStyleID</styleUrl>
#       <MultiGeometry>
#         <Point>
#           <altitudeMode>relativeToGround</altitudeMode>
#           <coordinates>
#             -79.7467307331394,29.2377632915929,79357.4676280487
#           </coordinates>
#         </Point>
#         <Point>
#           <altitudeMode>relativeToGround</altitudeMode>
#           <coordinates>
#             -79.743967829952,29.2399467587673,79822.8123209467
#           </coordinates>
#         </Point>


function makeLineString(llfl, k::Float64 = 1.)
  # llfl is lat lon fl matrix
  # k is plot factor, scales altitude
  lineString = ["""    <LineString>\n""",
                """      <altitudeMode>relativeToGround</altitudeMode>\n""",
                """      <coordinates>\n""",
                "", # THIS LINE MUST CONTAIN THE DATA!!!
                """      </coordinates>\n""",
                """    </LineString>\n"""]

#   lineStr = ""
#   curAlt = llfl[3,1]
  for (curLat, curLon, curAlt) in zip(llfl[1,:], ((llfl[2,:]+360.) % 360.), llfl[3,:] * 100 * ft2m * k)
    lineString[4] *= "$curLon,$curLat,$curAlt\n"
  end

  lineString
end



function staticBoxes(fileNameGE, storageDict::Dict, plotFactor)

  # I expect the Lat/Lon/Alt in listOfLists to be in Deg / Deg / M, just what GE expects
#   ft2m = 0.3048

  #fileNameGE = 'GE.kml'
  outFileGE = open(fileNameGE,"w")

  preamble = ["""<?xml version="1.0" encoding="UTF-8"?>\n""",
              """<kml xmlns="http://www.opengis.net/kml/2.2">\n""",
              """ <Document>\n""",
              """  <Style id="rocketStyleID">\n""",
              """   <LineStyle>\n""",
              """    <color>ff0000ff</color>\n""",
              """    <colorMode>normal</colorMode>\n""",
              """    <width>4</width>\n""",
              """   </LineStyle>\n""",
              """  </Style>\n""",
              """ <name>""" * fileNameGE * """</name>\n\n"""]

  for line in preamble
    write(outFileGE, line)
  end

  for (curName, curLatLonFL) in storageDict
#     curLat = curLatLonFL[1,:]
#     curLon = curLatLonFL[2,:]

#     curLatLonFL[3,:] = 10000
    lineStrVec = makeLineString(curLatLonFL, plotFactor)
    placemarkStr = makePlacemarkString(lineStrVec, "", curName)

    write(outFileGE, placemarkStr)

#     curAlt = 100
#     for (curLat, curLon) in zip(curLatLonFL[1,:], curLatLonFL[2,:])
#       print(makePointString(curLat, curLon, curAlt))
#     end
  end

  finale = "  </Document>\n" * "</kml>"
  write(outFileGE, finale)

  close(outFileGE)

end


function staticBoxesV2(fileNameGE, storageDict::Dict, plotFactor)
  # I expect the Lat/Lon/Alt in listOfLists to be in Deg / Deg / FL low, FL hi
  outFileGE = open(fileNameGE,"w")

  preamble = ["""<?xml version="1.0" encoding="UTF-8"?>\n""",
              """<kml xmlns="http://www.opengis.net/kml/2.2">\n""",
              """ <Document>\n""",
              """  <Style id="rocketStyleID">\n""",
              """   <LineStyle>\n""",
              """    <color>ff0000ff</color>\n""",
              """    <colorMode>normal</colorMode>\n""",
              """    <width>4</width>\n""",
              """   </LineStyle>\n""",
              """  </Style>\n""",
              """ <name>""" * fileNameGE * """</name>\n\n"""]

  for line in preamble
    write(outFileGE, line)
  end


  for (curName, curLatLonHiLo) in storageDict
    numFaces = length(curLatLonHiLo[1,:]) - 1  # One less face than number of points for closed shape
    tempLLFL = zeros(Float64, 3, numFaces*5)   # There are five points per face (last one closes)

    curPt = 1
    for curFace in 1:numFaces

      lat1, lon1, minFL, maxFL = curLatLonHiLo[:,curFace]
      lat2, lon2, minFL, maxFL = curLatLonHiLo[:,curFace+1]  #assuming same min/max alts
#       print("$lat1 -> $lat2 \n")

      tempLLFL[:,curPt] = [lat1, lon1, minFL]
      tempLLFL[:,curPt+1] = [lat1, lon1, maxFL]
      tempLLFL[:,curPt+2] = [lat2, lon2, maxFL]
      tempLLFL[:,curPt+3] = [lat2, lon2, minFL]
      tempLLFL[:,curPt+4] = [lat1, lon1, minFL]

      curPt += 5
    end

    lineStrVec = makeLineString(tempLLFL, plotFactor)
    placemarkStr = makePlacemarkString(lineStrVec, "", curName)

    write(outFileGE, placemarkStr)
  end

  finale = "  </Document>\n" * "</kml>"
  write(outFileGE, finale)

  close(outFileGE)
end


# def staticBoxes(fileNameGE, storageDict, plotFactor):
#     ft2m = 0.3048

#     #    # Make this an input!!!!
#     #    numElementsPerLine = 6

#     #fileNameGE = 'GE.kml'
#     outFileGE = open(fileNameGE,'w')

#     line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
#     line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
#     line3GE  = ' <Document>\n'
#     line4GE  = '  <Style id="style1">\n'
#     line5GE  = '   <LineStyle>\n'
#     line6GE  = '    <colorMode>random</colorMode>\n'
#     line7GE  = '    <width>4</width>\n'
#     line8GE  = '   </LineStyle>\n'
#     line9GE  = '  </Style>\n'
#     line10GE = ' <name>' + fileNameGE + '</name>\n\n'
#     outFileGE.write(line1GE)
#     outFileGE.write(line2GE)
#     outFileGE.write(line3GE)
#     outFileGE.write(line4GE)
#     outFileGE.write(line5GE)
#     outFileGE.write(line6GE)
#     outFileGE.write(line7GE)
#     outFileGE.write(line8GE)
#     outFileGE.write(line9GE)
#     outFileGE.write(line10GE)


#     #    lineCounter = 0
#     #    numRuns = len(listOfListsOfLists)
#     #    for curRun in range(0,numRuns):
#     for curItem in sorted(storageDict.keys()):
#         listOfLists = storageDict[curItem]
#         line1GE = '  <Placemark>\n'
#         line2GE = '   <name>'+str(curItem)+'</name>\n'
#         line3GE = '   <styleUrl>#style1</styleUrl>\n'
#         line4GE = '   <MultiGeometry>\n'
#         line5GE = '    <LineString>\n'
#         line6GE = '      <altitudeMode>relativeToGround</altitudeMode>\n'
#         line7GE = '      <coordinates>\n'

#         outFileGE.write(line1GE)
#         outFileGE.write(line2GE)
#         outFileGE.write(line3GE)
#         outFileGE.write(line4GE)
#         outFileGE.write(line5GE)
#         outFileGE.write(line6GE)
#         outFileGE.write(line7GE)

#         numFaces = len(listOfLists) - 1     # One less face than number of points for closed shape

# #        for elem in listOfLists:
#         for curPt in range(numFaces):

#             [lat1, lon1, minAlt, maxAlt] =listOfLists[curPt]
#             [lat2, lon2, minAlt, maxAlt] =listOfLists[curPt+1]  #assuming same min/max alts

#             # Cutoff maxAlt at NAS
# #            maxAlt = min(maxAlt, 60000) * ft2m
#             maxAlt = maxAlt * ft2m
#             minAlt = minAlt * ft2m

#             newLine1 = str(lon1)+','+str(lat1)+','+str(minAlt * plotFactor)+'\n'
#             newLine2 = str(lon1)+','+str(lat1)+','+str(maxAlt * plotFactor)+'\n'
#             newLine3 = str(lon2)+','+str(lat2)+','+str(maxAlt * plotFactor)+'\n'
#             newLine4 = str(lon2)+','+str(lat2)+','+str(minAlt * plotFactor)+'\n'
#             newLine5 = str(lon1)+','+str(lat1)+','+str(minAlt * plotFactor)+'\n'     # close the shape

#             outFileGE.write(newLine1)
#             outFileGE.write(newLine2)
#             outFileGE.write(newLine3)
#             outFileGE.write(newLine4)
#             outFileGE.write(newLine5)


# #            [lat, lon, alt] = elem
#             #        lat = flatArray[lineCounter]
#             #        lon = flatArray[lineCounter+1]
#             #        alt = flatArray[lineCounter+2]
#             #        lineCounter += numElementsPerLine

#             #        if tx > maxTimeSteps:
#             #            continue
#             #        if (alt > 18288) and cutoffNAS:
#             #            continue

# #            newLine = str(lon)+','+str(lat)+','+str(alt)+'\n'
# #            outFileGE.write(newLine)

#         line1GE = '    </coordinates>\n'
#         line2GE = '   </LineString>\n'
#         line3GE = '  </MultiGeometry>\n'
#         line4GE = ' </Placemark>\n\n'

#         outFileGE.write(line1GE)
#         outFileGE.write(line2GE)
#         outFileGE.write(line3GE)
#         outFileGE.write(line4GE)


#     line1GE = '  </Document>\n'
#     line2GE = '</kml>'


#     outFileGE.write(line1GE)
#     outFileGE.write(line2GE)

#     outFileGE.close()






















