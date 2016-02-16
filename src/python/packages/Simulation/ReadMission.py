import numpy as np
import sys

# This is sloppy, remove it
planetModel =0

# Collection of things from Frisco that I've changed some

def readInput(fileName, basePATH):
    #mission = dictDefinitions.createMain() # defining dictionaries
    
    try:
        inputFile = open(basePATH + fileName,'r')
    except:
        if len(sys.argv) == 1:
            print '\n!!! Error: No input file specified !!!\n' \
            + ' Proper command line usage: $ readInputTrajectory.py [inputFileName] \n \n'
            # print basePATH + fileName
        else:
            print '\n!!! Error: Could not open input file: ' + sys.argv[1] + ' !!!\n'
        # raise
        exit(1)
    
    for line in inputFile:
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            line = line[:iComment]
        if line.find(":") > 0:
            iColon = line.index(":")
            key = (line[:iColon].lower()).split()
            value = (line[iColon+1:]).split()
            if len(value) < 1:
                print '\n!!! Error: Invalid input arguments !!!\n' \
                + ' At line: ' + line.strip() + '\n'
                exit(1)
            
            
            # Vehicle information section
            
            elif key[0] == 'vehicle' and key[1] == 'stages':
                VehicleStages = float(value[0])
                if isint(VehicleStages) == 0 or VehicleStages <= 0:
                    print '\n!!! Error: Invalid number of stages, it must be a positive integer \n'\
                        + ' At line: ' + line.strip() + '\n'
                    exit(1)
                VehicleStages = int(VehicleStages)
            elif key[0] == 'vehicle' and key[1] == 'payload' and key[2]=='mass':
                VehiclePayloadMass = float(value[0])
                if VehiclePayloadMass <= 0:
                    print '\n!!! Error: Invalid Payload Mass \n'\
                        + ' At line: ' + line.strip() + '\n'
                    exit(1)
            
            
            
            elif key[0] == 'vehicle' and key[1] == 'structural' and key[2] == 'mass':
                checkpos = 1 # making sure values are postive.
                checkint = 0# if integer check is needed
                VehicleStructuralMass,err= makeVector(value,VehicleStages,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid Structural Mass \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'vehicle' and key[1] == 'propellant' and key[2] == 'mass':
                checkpos = 1 # making sure values are postive
                checkint = 0# if integer check is needed
                
                VehiclePropellantMass,err= makeVector(value,VehicleStages,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid Propellant Mass \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            
            
            
            
            elif key[0] == 'vehicle' and key[1] == 'isp':
                checkpos = 1 # making sure values are postive
                checkint = 0# if integer check is needed
                
                VehicleIsp,err= makeVector(value,VehicleStages,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid ISP \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            
            
            
            elif key[0] == 'vehicle' and key[1] == 'cd':
                #reading CD text files
                MinfCDlist = []
                CDList = []
                ArefList = []
                if len(value)!= VehicleStages:
                    print '\n!!! ERROR in readInputTrajectory.py. Incorrect number of CD files'
                    exit(1)
                for index in range(len(value)):
                    #print value[index]
                    Minf,CD,Aref = readInputCD(basePATH + value[index])
                    MinfCDlist.append(Minf)
                    CDList.append(CD)
                    ArefList.append(Aref)
            
            
            elif key[0] == 'latitude' and key[1] == 'initial':
                LatitudeInitial = float(value[0])
                if LatitudeInitial < -90 or LatitudeInitial > 90:
                    print '\n!!! Error: Latitude must be within -90 to 90 degrees \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'longitude' and key[1] == 'initial':
                LongitudeInitial = float(value[0])
                if LongitudeInitial < -180 or LongitudeInitial > 180:
                    print '\n!!! Error: Longitude must be within -180 to 180 degrees \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'altitude' and key[1] == 'initial':
                ElevationInitial = float(value[0])
                if ElevationInitial < 0:
                    print '\n!!! Error: Invalid Initial Elevation, it must be positive \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            
            elif key[0] == 'thetag' and key[1] == 'initial':
                thetag0 = float(value[0])
            
            
            elif key[0] == 'thrust' and key[1] == 'file':
                #reading CD text files
                timelist = []
                thrustList = []
                stageList = []
                if len(value)!= VehicleStages:
                    print '\n!!! ERROR in readInputTrajectory.py. Incorrect number of CD files'
                    exit(1)
                for index in range(len(value)):
                    #print value[index]
                    stageVal,timeVal,Tx,Ty,Tz = readInputThrust(basePATH + value[index])
                    timelist.append(timeVal)
                    thrustList.append([Tx,Ty,Tz])
                    stageList.append(stageVal)
            
            
            elif key[0]=='atmospheric' and key[1]=='option':
                atmoOption = float(value[0])
            elif key[0]=='atmospheric' and key[1]=='file':
                atmoFile = value[0]
            elif key[0] == 'launch' and key[1] == 'time':
                checkpos = 1 # making sure values are postive
                checkint = 1# if integer check is needed
                timeLength = 2 # hours and mins
                localTime,err= makeVector(value,timeLength,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid Local Time \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'launch' and key[1] == 'ut':
                UTshift = float(value[0])
            elif key[0] =='dt':
                dtval = float(value[0])
            
            elif key[0] =='launch' and key[1]=='date':
                dateval = value[0].split('/')
                month = float(dateval[0])
                day = float(dateval[1])
                year = float(dateval[2])
                
                if (isint(month)==0 or isint(day)==0 or isint(year)==0):
                    print '\n!!! Error: Invalid date, it must be positive integers \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
        
        elif len(line.strip()) != 0:
            print '\n!!! Error: Unrecognized input parameter !!!!\n' \
            + ' At line: ' + line.strip() + '\n'
            
            exit(1)
    
    inputFile.close()
    
    
    dateList = [month,day,year,UTshift,localTime]
    initialLocation = [LongitudeInitial,LatitudeInitial,ElevationInitial]
    # stage masses
    nend = VehicleStages
    m0 = np.zeros((VehicleStages))
    mf = np.zeros((VehicleStages))
    for index in range(0,VehicleStages):
        m0[index] = np.sum(VehiclePropellantMass[index:nend]) + np.sum(VehicleStructuralMass[index:nend]) + VehiclePayloadMass         # kg
        mf[index] = np.sum(VehiclePropellantMass[index+1:nend]) + np.sum(VehicleStructuralMass[index:nend]) + VehiclePayloadMass      # kg
        print VehiclePropellantMass[index+1:nend]
    massList = [m0,mf]
    
    #outputs = [initialLocation,massList,ArefList,MinfCDlist,CDList,atmoOption,atmoFile,timelist,thrustList,VehicleIsp,thetag0,dateList,dtval]

    outputs = dict(initialLocation = initialLocation, massList = massList, ArefList = ArefList, MinfCDlist = MinfCDlist,
                   CDList = CDList, atmoOption = atmoOption, atmoFile = atmoFile, timelist = timelist,
                   thrustList = thrustList, VehicleIsp = VehicleIsp, thetag0 = thetag0, dateList = dateList, dtval = dtval,
                   VehicleStructuralMass = VehicleStructuralMass)

    '''finalconditions,finalderivs = propagate(initialstate,mass,mass_time_alt_final,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,tlist,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,mass_time_alt_opt,thrustoffangledeg,ncd=len(minfcd),ncl=len(minfcl),ntime=shape(tlist,0),nlist=len(altitudelist))
        '''
    
    return(outputs)


# function to check if number is an integer
# returns 1 if it is an integer, 0 otherwise
def isint(number):
    x = number % 1
    b = 0
    if x == 0:
        b = 1
    return b
# end of integer function

def makeVector(array,stages,checkPos,checkint):
    err = 0
    if len(array)!=int(stages):
        print ' \n  Error: check number of stages agreement'
        err = 1
    ret = []
    for index in range(0,int(stages)):
        val = float(array[index])
        if (checkPos==1 and val<0):
            print ' \n Warning: value must be positive'
            err = 2
        if ((checkint==1) and (isint(val)==0)):
            err = 3
        ret.append(float(array[index]))
    return (np.array(ret),err)
def makeVectorSimple(array,checkPos,checkint):
    err = 0
    
    ret = []
    for index in range(0,len(array)):
        val = float(array[index])
        if (checkPos==1 and val<0):
            print ' \n Warning: value must be positive'
            err = 2
        if ((checkint==1) and (isint(val)==0)):
            err = 3
        ret.append(float(array[index]))
    return (np.array(ret),err)


# From readAeroData.py
def readInputCD(fileName):
    try:
        inputFile = open(fileName,'r')
    except:
        print '\n!!! Error: Could not open input file: ' + str(fileName) + ' !!!\n'
        exit(1)
    arefstring = 'Aref'
    Aref = np.array([])
    CD = np.array([])
    Minf = np.array([])
    for line in inputFile:
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            line = line[:iComment]
        if line.find(":") > 0:
            iColon = line.index(":")
            key = (line[:iColon].lower()).split()
            value = (line[iColon+1:]).split()
            if key[0].lower() == arefstring.lower():
                if len(value) < 1:
                    print '\n!!! Error: Invalid input arguments !!!\n' \
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
                Aref = float(value[0])
        elif len(line.split())==2:
            key = line.split()
            Minf = np.concatenate((Minf,[float(key[0])]),1)
            CD = np.concatenate((CD,[float(key[1])]),1)
    return (Minf,CD,Aref)

# From readThrust.py
def readInputThrust(fileName):
    try:
        inputFile = open(fileName,'r')
    except:
        print '\n!!! Error: Could not open input file: ' + str(fileName) + ' !!!\n'
        exit(1)
    time = np.array([])
    Tx = np.array([])
    Ty = np.array([])
    Tz = np.array([])
    
    for line in inputFile:
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            line = line[:iComment]
        if line.find(":") > 0:
            iColon = line.index(":")
            key = (line[:iColon].lower()).split()
            value = (line[iColon+1:]).split()
            if key[0].lower() == 'stage':
                if len(value) < 1:
                    print '\n!!! Error: Invalid input arguments !!!\n' \
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
                stage = float(value[0])
        elif len(line.split())==4:
            key = line.split()
            time = np.concatenate((time,[float(key[0])]),1)
            Tx = np.concatenate((Tx,[float(key[1])]),1)
            Ty = np.concatenate((Ty,[float(key[2])]),1)
            Tz = np.concatenate((Tz,[float(key[3])]),1)
    time = time - time[0]
    return (stage,time,Tx,Ty,Tz)