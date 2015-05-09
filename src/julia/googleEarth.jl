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
  for (curLat, curLon, curAlt) in zip(llfl[1,:], llfl[2,:]+360., llfl[3,:] * 100 * ft2m * k)
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


























# close(outFileGE)
#   lineCounter = 0
#   for curRun in range(0,numRuns):

#       curStyle = 'styleRed'
#     if not isRed[curRun]:
#         curStyle = 'styleGreen'

#       line1GE = '  <Placemark>\n'
#       line2GE = '   <name>' + callSigns[curRun] + '</name>\n'
#       line3GE = '   <styleUrl>#' + curStyle + '</styleUrl>\n'
#       line4GE = '   <MultiGeometry>\n'
#       line5GE = '    <LineString>\n'
#       line6GE = '      <altitudeMode>relativeToGround</altitudeMode>\n'
#       line7GE = '      <coordinates>\n'

#       outFileGE.write(line1GE)
#       outFileGE.write(line2GE)
#       outFileGE.write(line3GE)
#       outFileGE.write(line4GE)
#       outFileGE.write(line5GE)
#       outFileGE.write(line6GE)
#       outFileGE.write(line7GE)

#       for tx in range(numTimeSteps[curRun]):
#           lat = flatArray[lineCounter]
#         lon = flatArray[lineCounter+1]
#         alt = flatArray[lineCounter+2]
#         lineCounter += numElementsPerLine

#         newLine = str(lon)+','+str(lat)+','+str(alt)+'\n'
#         outFileGE.write(newLine)

#         line1GE = '    </coordinates>\n'
#         line2GE = '   </LineString>\n'
#         line3GE = '  </MultiGeometry>\n'
#         line4GE = ' </Placemark>\n\n'

#         outFileGE.write(line1GE)
#         outFileGE.write(line2GE)
#         outFileGE.write(line3GE)
#         outFileGE.write(line4GE)


#         line1GE = '  </Document>\n'
#         line2GE = '</kml>'


#         outFileGE.write(line1GE)
#         outFileGE.write(line2GE)

#         outFileGE.close()

