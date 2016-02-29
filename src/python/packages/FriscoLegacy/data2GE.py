
# File to get Latitude, Longitude and Altitude from debris.plt
# and then output it to Google Earth

def convert(baseFile,nSamples):

    N = nSamples

    altitude =0
    fileNameGE = 'GE.kml'
    outFileGE = open(fileNameGE,'w')

    line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    line3GE  = ' <Document>\n'
    line4GE  = '  <Style id="style1">\n'
    line5GE  = '   <LineStyle>\n'
    line6GE  = '    <colorMode>random</colorMode>\n'
    line7GE  = '    <width>4</width>\n'
    line8GE  = '   </LineStyle>\n'
    line9GE  = '  </Style>\n'
    line10GE = ' <name>Debris-Path</name>\n\n'  
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    outFileGE.write(line3GE)
    outFileGE.write(line4GE)
    outFileGE.write(line5GE)
    outFileGE.write(line6GE)
    outFileGE.write(line7GE)
    outFileGE.write(line8GE)
    outFileGE.write(line9GE)
    outFileGE.write(line10GE)



    for index in range(1,N+1):

       fileName = baseFile+str(index)
       print fileName
       try:
          inputFile = open(fileName,'r')
    #   outputFile = open(fileNameOut, 'w')
       except:
          print 'Cant open file'

       line1GE = '  <Placemark>\n'
       line2GE = '   <name>Piece'+str(index)+'</name>\n'
       line3GE = '   <styleUrl>#style1</styleUrl>\n'
       line4GE = '   <MultiGeometry>\n'
       line5GE = '    <LineString>\n'
       line6GE = '      <altitudeMode>absolute</altitudeMode>\n'
       line7GE = '      <coordinates>\n'

       outFileGE.write(line1GE)
       outFileGE.write(line2GE)
       outFileGE.write(line3GE)
       outFileGE.write(line4GE)
       outFileGE.write(line5GE)
       outFileGE.write(line6GE)
       outFileGE.write(line7GE)


       lineNumber = 0
       sample = 1
       counter = 0
       altitude = -500

       for line in inputFile:
          counter = counter+1
          lineNumber = lineNumber + 1
          key = line.split() 
        
          if lineNumber > 2 :
             alt = float(key[2])
         #    print alt,index
             if abs(altitude - alt)>2: 
                 # if counter%sample==0:
                newLine = key[0]+','+key[1]+','+key[2]+'\n'
                outFileGE.write(newLine)
                altitude = float(key[2])

       line1GE = '    </coordinates>\n'
       line2GE = '   </LineString>\n'
       line3GE = '  </MultiGeometry>\n'
       line4GE = ' </Placemark>\n\n'
      
       outFileGE.write(line1GE)
       outFileGE.write(line2GE)
       outFileGE.write(line3GE)
       outFileGE.write(line4GE)


    line1GE = '  </Document>\n'
    line2GE = '</kml>'


    outFileGE.write(line1GE)
    outFileGE.write(line2GE)

    inputFile.close()
    outFileGE.close()





def convertTJC(fileNameGE, flatArray, numTimeSteps, numRuns, cutoffNAS = False, maxTimeSteps = 1e10):
    # I expect the Lat/Lon/Alt in flatArray to be in Deg / Deg / M, just what GE expects
    
    # Make this an input!!!!
    numElementsPerLine = 6
    
    #fileNameGE = 'GE.kml'
    outFileGE = open(fileNameGE,'w')
    
    line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    line3GE  = ' <Document>\n'
    line4GE  = '  <Style id="style1">\n'
    line5GE  = '   <LineStyle>\n'
    line6GE  = '    <colorMode>random</colorMode>\n'
    line7GE  = '    <width>4</width>\n'
    line8GE  = '   </LineStyle>\n'
    line9GE  = '  </Style>\n'
    line10GE = ' <name>Debris-Path</name>\n\n'
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    outFileGE.write(line3GE)
    outFileGE.write(line4GE)
    outFileGE.write(line5GE)
    outFileGE.write(line6GE)
    outFileGE.write(line7GE)
    outFileGE.write(line8GE)
    outFileGE.write(line9GE)
    outFileGE.write(line10GE)
    
    
    lineCounter = 0
    for curRun in range(0,numRuns):
        
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>Piece'+str(curRun)+'</name>\n'
        line3GE = '   <styleUrl>#style1</styleUrl>\n'
        line4GE = '   <MultiGeometry>\n'
        line5GE = '    <LineString>\n'
        line6GE = '      <altitudeMode>relativeToGround</altitudeMode>\n'
        line7GE = '      <coordinates>\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
        outFileGE.write(line5GE)
        outFileGE.write(line6GE)
        outFileGE.write(line7GE)
        
        for tx in range(numTimeSteps[curRun]):
            lat = flatArray[lineCounter]
            lon = flatArray[lineCounter+1]
            alt = flatArray[lineCounter+2]
            lineCounter += numElementsPerLine
            
            if tx > maxTimeSteps:
                continue
            if (alt > 18288) and cutoffNAS:
                continue
        
            newLine = str(lon)+','+str(lat)+','+str(alt)+'\n'
            outFileGE.write(newLine)
        
        line1GE = '    </coordinates>\n'
        line2GE = '   </LineString>\n'
        line3GE = '  </MultiGeometry>\n'
        line4GE = ' </Placemark>\n\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
    
    
    line1GE = '  </Document>\n'
    line2GE = '</kml>'
    
    
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    
    outFileGE.close()





def convertTJCAircraft(fileNameGE, flatArray, numTimeSteps, callSigns, numElementsPerLine):
    # I expect the Lat/Lon/Alt in flatArray to be in Deg / Deg / M, just what GE expects
    
    # Make this an input!!!!
#    numElementsPerLine = 6
    
    numRuns = len(numTimeSteps)
    
    #fileNameGE = 'GE.kml'
    outFileGE = open(fileNameGE,'w')
    
    line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    line3GE  = ' <Document>\n'
    line4GE  = '  <Style id="style1">\n'
    line5GE  = '   <LineStyle>\n'
    line6GE  = '    <colorMode>random</colorMode>\n'
    line7GE  = '    <width>4</width>\n'
    line8GE  = '   </LineStyle>\n'
    line9GE  = '  </Style>\n'
    line10GE = ' <name>Debris-Path</name>\n\n'
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    outFileGE.write(line3GE)
    outFileGE.write(line4GE)
    outFileGE.write(line5GE)
    outFileGE.write(line6GE)
    outFileGE.write(line7GE)
    outFileGE.write(line8GE)
    outFileGE.write(line9GE)
    outFileGE.write(line10GE)
    
    
    lineCounter = 0
    for curRun in range(0,numRuns):
        
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>' + callSigns[curRun] + '</name>\n'
        line3GE = '   <styleUrl>#style1</styleUrl>\n'
        line4GE = '   <MultiGeometry>\n'
        line5GE = '    <LineString>\n'
        line6GE = '      <altitudeMode>relativeToGround</altitudeMode>\n'
        line7GE = '      <coordinates>\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
        outFileGE.write(line5GE)
        outFileGE.write(line6GE)
        outFileGE.write(line7GE)
        
        for tx in range(numTimeSteps[curRun]):
            lat = flatArray[lineCounter]
            lon = flatArray[lineCounter+1]
            alt = flatArray[lineCounter+2]
            lineCounter += numElementsPerLine
            
            newLine = str(lon)+','+str(lat)+','+str(alt)+'\n'
            outFileGE.write(newLine)
        
        line1GE = '    </coordinates>\n'
        line2GE = '   </LineString>\n'
        line3GE = '  </MultiGeometry>\n'
        line4GE = ' </Placemark>\n\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
    
    
    line1GE = '  </Document>\n'
    line2GE = '</kml>'
    
    
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    
    outFileGE.close()




def convertTJCAircraftRedGreen(fileNameGE, flatArray, numTimeSteps, callSigns, isRed, numElementsPerLine):
    # I expect the Lat/Lon/Alt in flatArray to be in Deg / Deg / M, just what GE expects
    
    # Make this an input!!!!
    #    numElementsPerLine = 6
    
    numRuns = len(numTimeSteps)
    
    #fileNameGE = 'GE.kml'
    outFileGE = open(fileNameGE,'w')
    
#    line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
#    line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
#    line3GE  = ' <Document>\n'
#    line4GE  = '  <Style id="style1">\n'
#    line5GE  = '   <LineStyle>\n'
#    line6GE  = '    <colorMode>random</colorMode>\n'
#    line7GE  = '    <width>4</width>\n'
#    line8GE  = '   </LineStyle>\n'
#    line9GE  = '  </Style>\n'
#    line10GE = ' <name>Debris-Path</name>\n\n'
#    outFileGE.write(line1GE)
#    outFileGE.write(line2GE)
#    outFileGE.write(line3GE)
#    outFileGE.write(line4GE)
#    outFileGE.write(line5GE)
#    outFileGE.write(line6GE)
#    outFileGE.write(line7GE)
#    outFileGE.write(line8GE)
#    outFileGE.write(line9GE)
#    outFileGE.write(line10GE)
    
    preamble = ['<?xml version="1.0" encoding="UTF-8"?>\n',
                '<kml xmlns="http://www.opengis.net/kml/2.2">\n',
                ' <Document>\n',
                '  <Style id="styleRed">\n',
                '   <LineStyle>\n',
                '    <color>ff0000ff</color>\n',
                '    <colorMode>normal</colorMode>\n',
                '    <width>4</width>\n',
                '   </LineStyle>\n',
                '  </Style>\n',
                '  <Style id="styleGreen">\n',
                '   <LineStyle>\n',
                '    <color>ff00ff00</color>\n',
                '    <colorMode>normal</colorMode>\n',
                '    <width>4</width>\n',
                '   </LineStyle>\n',
                '  </Style>\n',
                ' <name>Debris-Path</name>\n\n']
    
    for line in preamble:
        outFileGE.write(line)
    
    lineCounter = 0
    for curRun in range(0,numRuns):
        
        curStyle = 'styleRed'
        if not isRed[curRun]:
            curStyle = 'styleGreen'
        
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>' + callSigns[curRun] + '</name>\n'
        line3GE = '   <styleUrl>#' + curStyle + '</styleUrl>\n'
        line4GE = '   <MultiGeometry>\n'
        line5GE = '    <LineString>\n'
        line6GE = '      <altitudeMode>relativeToGround</altitudeMode>\n'
        line7GE = '      <coordinates>\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
        outFileGE.write(line5GE)
        outFileGE.write(line6GE)
        outFileGE.write(line7GE)
        
        for tx in range(numTimeSteps[curRun]):
            lat = flatArray[lineCounter]
            lon = flatArray[lineCounter+1]
            alt = flatArray[lineCounter+2]
            lineCounter += numElementsPerLine
            
            newLine = str(lon)+','+str(lat)+','+str(alt)+'\n'
            outFileGE.write(newLine)
        
        line1GE = '    </coordinates>\n'
        line2GE = '   </LineString>\n'
        line3GE = '  </MultiGeometry>\n'
        line4GE = ' </Placemark>\n\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
    
    
    line1GE = '  </Document>\n'
    line2GE = '</kml>'
    
    
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    
    outFileGE.close()








def staticTrajectory(fileNameGE, storageDict):
    # I expect the Lat/Lon/Alt in listOfLists to be in Deg / Deg / M, just what GE expects
    
#    # Make this an input!!!!
#    numElementsPerLine = 6
    
    #fileNameGE = 'GE.kml'
    outFileGE = open(fileNameGE,'w')
    
    line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    line3GE  = ' <Document>\n'
    line4GE  = '  <Style id="style1">\n'
    line5GE  = '   <LineStyle>\n'
    line6GE  = '    <colorMode>random</colorMode>\n'
    line7GE  = '    <width>4</width>\n'
    line8GE  = '   </LineStyle>\n'
    line9GE  = '  </Style>\n'
    line10GE = ' <name>Debris-Path</name>\n\n'
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    outFileGE.write(line3GE)
    outFileGE.write(line4GE)
    outFileGE.write(line5GE)
    outFileGE.write(line6GE)
    outFileGE.write(line7GE)
    outFileGE.write(line8GE)
    outFileGE.write(line9GE)
    outFileGE.write(line10GE)
    
    
#    lineCounter = 0
#    numRuns = len(listOfListsOfLists)
#    for curRun in range(0,numRuns):
    for curItem in sorted(storageDict.keys()):
        listOfLists = storageDict[curItem]
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>'+str(curItem)+'</name>\n'
        line3GE = '   <styleUrl>#style1</styleUrl>\n'
        line4GE = '   <MultiGeometry>\n'
        line5GE = '    <LineString>\n'
        line6GE = '      <altitudeMode>relativeToGround</altitudeMode>\n'
        line7GE = '      <coordinates>\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
        outFileGE.write(line5GE)
        outFileGE.write(line6GE)
        outFileGE.write(line7GE)
        
        for elem in listOfLists:
            [lat, lon, alt] = elem
    #        lat = flatArray[lineCounter]
    #        lon = flatArray[lineCounter+1]
    #        alt = flatArray[lineCounter+2]
    #        lineCounter += numElementsPerLine
            
    #        if tx > maxTimeSteps:
    #            continue
    #        if (alt > 18288) and cutoffNAS:
    #            continue
            
            newLine = str(lon)+','+str(lat)+','+str(alt)+'\n'
            outFileGE.write(newLine)
        
        line1GE = '    </coordinates>\n'
        line2GE = '   </LineString>\n'
        line3GE = '  </MultiGeometry>\n'
        line4GE = ' </Placemark>\n\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
    
    
    line1GE = '  </Document>\n'
    line2GE = '</kml>'
    
    
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    
    outFileGE.close()




def staticBoxes(fileNameGE, storageDict, plotFactor):
    # I expect the Lat/Lon/Alt in listOfLists to be in Deg / Deg / M, just what GE expects
    ft2m = 0.3048

    #    # Make this an input!!!!
    #    numElementsPerLine = 6
    
    #fileNameGE = 'GE.kml'
    outFileGE = open(fileNameGE,'w')
    
    line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    line3GE  = ' <Document>\n'
    line4GE  = '  <Style id="style1">\n'
    line5GE  = '   <LineStyle>\n'
    line6GE  = '    <colorMode>random</colorMode>\n'
    line7GE  = '    <width>4</width>\n'
    line8GE  = '   </LineStyle>\n'
    line9GE  = '  </Style>\n'
    line10GE = ' <name>' + fileNameGE + '</name>\n\n'
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    outFileGE.write(line3GE)
    outFileGE.write(line4GE)
    outFileGE.write(line5GE)
    outFileGE.write(line6GE)
    outFileGE.write(line7GE)
    outFileGE.write(line8GE)
    outFileGE.write(line9GE)
    outFileGE.write(line10GE)
    
    
    #    lineCounter = 0
    #    numRuns = len(listOfListsOfLists)
    #    for curRun in range(0,numRuns):
    for curItem in sorted(storageDict.keys()):
        listOfLists = storageDict[curItem]
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>'+str(curItem)+'</name>\n'
        line3GE = '   <styleUrl>#style1</styleUrl>\n'
        line4GE = '   <MultiGeometry>\n'
        line5GE = '    <LineString>\n'
        line6GE = '      <altitudeMode>relativeToGround</altitudeMode>\n'
        line7GE = '      <coordinates>\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
        outFileGE.write(line5GE)
        outFileGE.write(line6GE)
        outFileGE.write(line7GE)
        
        numFaces = len(listOfLists) - 1     # One less face than number of points for closed shape
        
#        for elem in listOfLists:
        for curPt in range(numFaces):
            
            [lat1, lon1, minAlt, maxAlt] =listOfLists[curPt]
            [lat2, lon2, minAlt, maxAlt] =listOfLists[curPt+1]  #assuming same min/max alts
            
            # Cutoff maxAlt at NAS
#            maxAlt = min(maxAlt, 60000) * ft2m
            maxAlt = maxAlt * ft2m
            minAlt = minAlt * ft2m

            newLine1 = str(lon1)+','+str(lat1)+','+str(minAlt * plotFactor)+'\n'
            newLine2 = str(lon1)+','+str(lat1)+','+str(maxAlt * plotFactor)+'\n'
            newLine3 = str(lon2)+','+str(lat2)+','+str(maxAlt * plotFactor)+'\n'
            newLine4 = str(lon2)+','+str(lat2)+','+str(minAlt * plotFactor)+'\n'
            newLine5 = str(lon1)+','+str(lat1)+','+str(minAlt * plotFactor)+'\n'     # close the shape

            outFileGE.write(newLine1)
            outFileGE.write(newLine2)
            outFileGE.write(newLine3)
            outFileGE.write(newLine4)
            outFileGE.write(newLine5)

        
#            [lat, lon, alt] = elem
            #        lat = flatArray[lineCounter]
            #        lon = flatArray[lineCounter+1]
            #        alt = flatArray[lineCounter+2]
            #        lineCounter += numElementsPerLine
            
            #        if tx > maxTimeSteps:
            #            continue
            #        if (alt > 18288) and cutoffNAS:
            #            continue
            
#            newLine = str(lon)+','+str(lat)+','+str(alt)+'\n'
#            outFileGE.write(newLine)
        
        line1GE = '    </coordinates>\n'
        line2GE = '   </LineString>\n'
        line3GE = '  </MultiGeometry>\n'
        line4GE = ' </Placemark>\n\n'
        
        outFileGE.write(line1GE)
        outFileGE.write(line2GE)
        outFileGE.write(line3GE)
        outFileGE.write(line4GE)
    
    
    line1GE = '  </Document>\n'
    line2GE = '</kml>'
    
    
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    
    outFileGE.close()



