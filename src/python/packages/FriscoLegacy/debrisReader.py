# Debris reader, format obtained from Space Shuttle debris catalog
# The number of fields for this catalog is 6
import numpy as np




# routine to help set up sobol indices calculation when running OPENTURNS...
# TODO=> DAKOTA SOBOL

def ProbaModel(debrisCatalog):
    # it sets the distribution to be used in openturns....mainly Aref, Mass, and Vimpulse
    # debris catalog order 
    #=> (num pieces) (weight/piece) (total Weight) (Ref Area) (ball coef) (Vel Imp) (CD) (CL)
    retList = []
    for index in range(len(debrisCatalog)):
        lineVal = debrisCatalog[index]
        for indexVar in range(len(lineVal)):
            # focusing on weight/piece, ref area, and Vel Imp
            varName='None'
            
            if indexVar==1: # weight/piece 
                varName = 'mass'
                val = lineVal[indexVar]
                masscheck,key = isRange(val)
                massVal = key

            elif indexVar==3: #ref area
                varName = 'Sref'
                val = lineVal[indexVar]
                Srefcheck,key = isRange(val)
                SrefVal = key

            elif indexVar==5: #Vel Imp
                varName = 'Vel_imp'
                val = lineVal[indexVar]
                Velcheck,key = isRange(val)     
                velVal = key

        
        retList.append([[masscheck,massVal],[Srefcheck,SrefVal],[Velcheck,velVal]])
    return retList
        




def is_number(s):
    try:
        float(s) # for int, long and float
    except ValueError:
        try:
            complex(s) # for complex
        except ValueError:
            return False
    
    return True




def readDebris(fileName):
    vel = []
    time = []
    delta = []
    inputFile = open(fileName,'r')
    catalogList =[]
    counter = -1
    lineCount = 0
    for line in inputFile:
        lineCount = lineCount+1
        lineComment = ''
        iComment = 0
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            lineComment = line[iComment:-1]
            line = line[:iComment]
        #key = line.split()
    
        import re
        line = re.sub(r'\[ *','[',line)     #strip spaces between left bracket and number
        line = re.sub(r' *\]',']',line)     #strip spaces between number and right bracket
        line = re.sub(r'\, *',',',line)        #strip spaces between comma and right number
        key = line.split()
    
        if len(key)>0:
            if key[0]=='velocity':
                vel.append(float(key[1]))
            elif key[0]=='time':
                time.append(float(key[1]))
            elif key[0]=='begin':
                counter = counter+1             
                catalogList.append([])
            elif key[0].lower()=='units':
                units = key[1]
            elif (len(key)==9 or len(key)==8) and line[0]!='#':
                #print counter,lineCount
                catalogList[counter].append(key)
            elif key[0]!='end' and line[0]!='#':
                print 'ERROR in DebrisCatalog: '+fileName
                print 'Cannot interpret line '+ str(lineCount)
                print line
                exit(1)

    if  units.lower=='si':
        units = 1
    elif units.lower()=='eng':
        units = -1
        print 'Current debris catalog is in ENGLIGH UNITS'
        print 'Convert debris catalog to SI units'
        print 'Error in debrisReader.py'
        exit(3)
        
    if len(catalogList)!=len(vel) or len(catalogList)!=len(time) or counter+1!=len(time):
        print 'Error reading Debris Catalog'
        exit(2)
    return (vel,time,catalogList)



def convertDebrisEng2SI(fileName,outFileName=None):
    
    if outFileName==None:
        outFileName = 'SI_'+fileName

        

    inputFile = open(fileName,'r')
    outFile = open(outFileName,'w')
    catalogList =[]
    counter = -1
    lineCount = 0
    ft2m = 0.3048
    ft2m2 = ft2m**2.0
    lbm2kg = 0.453592
    
    lbperft2_to_kgperm2 = lbm2kg/(ft2m2)
    outFile.write('#Conversion from English to SI Unit\n')
   
    for line in inputFile:
        lineCount = lineCount+1
        lineComment = ''
        iComment = 0
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            lineComment = line[iComment:-1]
            line = line[:iComment]
            lineComment= lineComment.replace('lb','kg')
            lineComment= lineComment.replace('ft','m')
        key = line.split()
       
            
        
        if len(key)>0:
            if key[0]=='velocity':
                key[1]=str(ft2m*float(key[1]))
            elif key[0]=='time':
                key[1]=str(float(key[1]))
            elif key[0]=='begin':
                counter = counter+1             
                catalogList.append([])
            elif key[0].lower()=='units' and key[1].lower=='si':
                print 'Current debris catalog is already in SI UNITS'
                print 'No Conversion required'
            elif  key[0].lower()=='units' and (key[1].lower()=='eng' or key[1].lower()=='english'):
                key = ['units SI']
                    
            elif (len(key)==9 or len(key)==8) and line[0]!='#':
                #mass for each piece
             
                key[1]=convertVar(key[1],lbm2kg)
                #mass for total group
                key[2]=convertVar(key[2],lbm2kg)
                #reference area
                key[3]=convertVar(key[3],ft2m2)
                #ballistic coeff
                key[4]=convertVar(key[4],lbperft2_to_kgperm2)
                #impulse velocity
                key[5]=convertVar(key[5],ft2m)
    
                #print counter,lineCount
                catalogList[counter].append(key)
            elif key[0]!='end' and line[0]!='#':
                print 'ERROR in DebrisCatalog: '+fileName
                print 'Cannot interpret line '+ str(lineCount)
                print line
                exit(1)
        keyOut = ''
        for indexKey in range(len(key)):
            keyOut = keyOut+'    '+key[indexKey]
        outFile.write(keyOut+' '+lineComment+'\n')
    

    outFile.close()
    return 

def convertVar(key,convFactor):
    check,tempkey=isRange(key)
    if check==1:
        key = '['+str(float(convFactor*tempkey[0]))+','+str(convFactor*tempkey[1])+']'
    else:
        key = str(convFactor*float(tempkey))
    return key


            
def desiredList(vel,time,catalogList,desiredVel,desiredTime):
 
        
    #    vel = np.array(vel)
    time = np.array(time)
    #delta = vel*time
    #deltaDesired = desiredVel*desiredTime
    #dmin = abs(delta-deltaDesired)
    dmin = abs(time-desiredTime)
    minval = min(dmin)
    for index in range(0,len(time)):
        if minval==dmin[index]:
            break
    catalog = catalogList[index]
    return catalog

    
## Currently not using this function





def mainReader(fileName):

    inputFile = open(fileName,'r')
    numberOfFields = 6 # this value depends strictly in the debris catalog
    lineNumber  = 0
    lineList=[]

    for line in inputFile:
        key = line.split()
        lineNumber = lineNumber+1
        if len(key)!=numberOfFields and lineNumber>1: # ensuring the number of fields is correct
            print 'check Line Number ',lineNumber
            exit(1)
        elif lineNumber>1:
            lineList.append(key)

    return lineList

            
            
def isRange(variable):
    # Returns 0 if it is not a Range,1 if it is. Also returns the bounds
    key=None
    if variable[0]=='['and variable[-1]==']':
        #we have a range
        key = variable[1:-1].split(',')
 
        if len(key)!=2:
            print 'Error in debrisReader.py deciding if value is a range'
            print 'Check ',variable
            exit(1)
    
    if key!=None:

        check = 1
        
        key[0] = float(key[0])
        key[1] = float(key[1])

                
    else:
        check = 0
        # print 'variable = ' + variable
        key = float(variable)


    return check,key
   


def generateMeanBallisticCoef(catalogList):
# this function returns a catalog with the mean ballistic coefficient to be used for GRAM profile generation
    ballisticCoef = []
    minball = []
    maxball = []
    for index in range(0,len(catalogList)):
        key = catalogList[index]
        ballisticrange,ballisticVal = isRange(key[4])
        if ballisticrange==1:
            minval = np.min(ballisticVal)
            maxval = np.max(ballisticVal)
            value = .5*(ballisticVal[0]+ballisticVal[1])
        else:
            value = ballisticVal
            minval = ballisticVal
            maxval = ballisticVal
        ballisticCoef.append(value)
        minball.append(minval)
        maxball.append(maxval)
    #print ballisticCoef
    return ballisticCoef,np.min(minball),np.max(maxball)

    
def getMassBounds(catalogList):
    minmass = []
    maxmass = []
    totalMass = []
    for index in range(0,len(catalogList)):
        key = catalogList[index]
        massrange,massVal = isRange(key[1])
        if massrange==1:
            minval = np.min(massVal)
            maxval = np.max(massVal)
        #   value = .5*(ballisticVal[0]+ballisticVal[1])
        else:
            value = massVal
            minval = massVal
            maxval = massVal
        #ballisticCoef.append(value)
        minmass.append(minval)
        maxmass.append(maxval)
    #print ballisticCoef
        totalMass.append(float(key[2]))
    return minmass,maxmass,totalMass



def generateRandomCatalog(catalogList):
    #Generate random pieces for each group. 
    # the debris catalog currently expects the element order to be as follows:
    # Number of Pieces-Weight of each piece(lbs)-Total Weight(lbs)-Reference Area(ft^2)-Ballistic Coeff(lb/ft^2)-Velocity Increment(ft/sec)-CD-CL	
    # first calculates number of pieces, then samples from weight. Keeps iterating until the weights in all generated pieces matches the Total Weight
    randomArray = []
    catalogLength = len(catalogList)
    
    for index in range(catalogLength): # step into each subgroup in catalog
        key = catalogList[index]
        #print line
        #key = line.split()
        currentweight = 0
        TotalWeight = float(key[2]) # storing total weight value in subgroup
        piecerange, value = isRange(key[0]) # checking is total number of pieces is a range
        # setting the number of pieces
        if piecerange==1:
            numberOfPieces = np.random.uniform(value[0],value[1])
        else:
            numberOfPieces = value # this value is never really used (weight is actually used to generate number of pieces)
    
        weightValue =0
        
        while (abs(currentweight-TotalWeight)>.01):        
        
            weightrange,weightValue = isRange(key[1]) # looking at weight location
            if weightrange==1:
                weight = (np.random.uniform(weightValue[0],weightValue[1]))
            else:
                weight = weightValue
                    
            arefrange,arefval = isRange(key[3]) #checking if reference area is a range
            if arefrange ==1:# assigning values for reference area
                aref = (np.random.uniform(arefval[0],arefval[1]))
            else:
                aref = arefval
                
            velrange,velval = isRange(key[5])
            if velrange==1:
                vel = (np.random.uniform(velval[0],velval[1]))
            else:
                vel = velval                    
            if is_number(key[6]):
                CD = key[6]
            else:
                CD = key[6]+'CD.txt'
                    
            if is_number(key[7]):
                CL = key[7]
            elif key[7]=='NONE':
                CL = '0'
            else:
                CL = key[7]+'CL.txt'
                

            currentweight = currentweight + weight
            
            if currentweight>=TotalWeight: #ensuring that the generated weights add to the TotalWeight 
                currentweight = currentweight - weight
                weight = TotalWeight - currentweight
                currentweight = currentweight + weight
    
            randomArray.append([weight,aref,vel,CD,CL])
            
    return randomArray

def generateRandomGroupCatalog(catalogList,nSamples):
    #Generate random pieces for each group. 
    # the debris catalog currently expects the element order to be as follows:
    # Number of Pieces-Weight of each piece(lbs)-Total Weight(lbs)-Reference Area(ft^2)-Ballistic Coeff(lb/ft^2)-Velocity Increment(ft/sec)-CD-CL	
    # it generates nSamples from each debris group
    randomArray = []
    catalogLength = len(catalogList)
    
    for index in range(catalogLength): # step into each subgroup in catalog
        key = catalogList[index]
        #print line
        #key = line.split()
        currentweight = 0
        TotalWeight = float(key[2]) # storing total weight value in subgroup
        piecerange, value = isRange(key[0]) # checking is total number of pieces is a range
        # setting the number of pieces
        if piecerange==1:
            numberOfPieces = np.random.uniform(value[0],value[1])
            numberOfPiecesMean = .5*(value[0]+value[1])
        else:
            numberOfPieces = value
            numberOfPiecesMean = 1.*value
        weightValue =0
        
            #while (abs(currentweight-TotalWeight)>.01):        
        for index2 in range(0,nSamples):    
            weightrange,weightValue = isRange(key[1]) # looking at weight location
            if weightrange==1:
                weight = (np.random.uniform(weightValue[0],weightValue[1]))
            else:
                weight = weightValue
            
            arefrange,arefval = isRange(key[3]) #checking if reference area is a range
            if arefrange ==1:# assigning values for reference area
                aref = (np.random.uniform(arefval[0],arefval[1]))
                arefMean = .5*(arefval[0]+arefval[1])
            else:
                aref = arefval
                arefMean = 1.*aref
            
            velrange,velval = isRange(key[5])
            if velrange==1:
                vel = (np.random.uniform(velval[0],velval[1]))
            else:
                vel = velval                    
            if is_number(key[6]):
                CD = key[6]
            else:
                CD = key[6]+'CD.txt'
            
            if is_number(key[7]):
                CL = key[7]
            elif key[7]=='NONE':
                CL = '0'
            else:
                CL = key[7]+'CL.txt'
            
            
            currentweight = currentweight + weight
            
        #   if currentweight>=TotalWeight: #ensuring that the generated weights add to the TotalWeight 
        #       currentweight = currentweight - weight
        #       weight = TotalWeight - currentweight
        #     currentweight = currentweight + weight
            
            randomArray.append([weight,aref,vel,CD,CL,arefMean,numberOfPiecesMean])
    
    return randomArray    

def generateDebris(catalogList,nSamplesList,debrisFilesPATH='None'):
    #Generate random pieces for each group. 
    # the debris catalog currently expects the element order to be as follows:
    # Number of Pieces-Weight of each piece(lbs)-Total Weight(lbs)-Reference Area(ft^2)-Ballistic Coeff(lb/ft^2)-Velocity Increment(ft/sec)-CD-CL	
    # it generates nSamples from each debris group. Returns lists with the desired ballistic coefficient at a given time
    randomArray = []
    catalogLength = len(catalogList)
    weight = [0]*len(catalogList)
    aref = [0]*len(catalogList)
    arefMean = [0]*len(catalogList)
    vel = [0]*len(catalogList)
    CD = [0]*len(catalogList)
    CL = [0]*len(catalogList)
    blast = [0]*len(catalogList)
    minfcdList = [0]*len(catalogList)
    minfclList = [0]*len(catalogList)
    ncdList = [0]*len(catalogList)
    nclList = [0]*len(catalogList)
    numberOfPiecesMean = [0]*len(catalogList)
    #print len(catalogList)
    
    if isinstance(nSamplesList, (int, long, float, complex)):
        nSamplesList = [nSamplesList]*len(catalogList)
    elif len(nSamplesList)==1:
        nSamplesList = [nSamplesList[0]]*len(catalogList)

    #print nSamplesList
    for index in range(catalogLength): # step into each subgroup in catalog
        nSamples=nSamplesList[index]

        key = catalogList[index]
        #print 'key = ' + str(key)
        #print line
        #key = line.split()
        currentweight = 0
        TotalWeight = float(key[2]) # storing total weight value in subgroup
        piecerange, value = isRange(key[0]) # checking is total number of pieces is a range
        # setting the number of pieces
        
        
        if piecerange==1:
            numberOfPieces = np.random.uniform(value[0],value[1])
            numberOfPiecesMeanVal = .5*(value[0]+value[1])
        else:
            numberOfPieces = value
            numberOfPiecesMeanVal = 1.*value
        numberOfPiecesMean[index] = numberOfPiecesMeanVal
        weightrange,weightValue = isRange(key[1]) # looking at weight location
        arefrange,arefval = isRange(key[3]) #checking if reference area is a range
        ballisticrange,ballCoeff = isRange(key[4]) # checking that ballistic coefficient is a range
        weightarray = np.zeros((nSamples))
        arefarray = np.zeros((nSamples))

        #getting CD and CL  
        #initializing values, default to 1 if cd and cl are constants 
        minfcd = 1.
        ncd = 1 
        minfcl =1.
        ncl = 1

        if is_number(key[6]):
            Cd = float(key[6])
        elif key[6]=='NONE':
            Cd = 0.0
        else:
            CDstr = key[6]+'CD.txt'
            minfcd,Cd,ncd =getCLCD(CDstr,debrisFilesPATH)

        if is_number(key[7]):
            Cl = float(key[7])
        elif key[7]=='NONE':
            Cl = 0.0
        else:
            CLstr = key[7]+'CL.txt'       
            minfcl,Cl,ncl =getCLCD(CLstr,debrisFilesPATH)

        CD[index] = Cd#[Cd]*nSamples
        CL[index] = Cl # [Cl]*nSamples
        minfclList[index] = minfcl
        minfcdList[index] = minfcd
        ncdList[index] = ncd
        nclList[index] = ncl       

        if ncd>1:
            Cd = np.array(Cd)
            minfcd = np.array(minfcd)
            #print 'Cd',Cd
            #print 'Minf',minfcd
            #print 'Arange',arefval
            #print 'Mass',weightValue
            #print 'Ball',ballisticrange
            CDaver = 1./(minfcd[-1]-minfcd[0])*np.trapz(Cd,x=minfcd)
            #print 'CD aver is',CDaver
    
        else:
            #print 'Const Cd',Cd
            CDaver = Cd
        for indDebris in range(nSamples):
            validDebris = 0

            while validDebris==0:
                
                if weightrange==1:
                    weightTemp = np.random.uniform(weightValue[0],weightValue[1])
                else:
                    weightTemp = weightValue
                if arefrange ==1:# assigning values for reference area
                    arefTemp = np.random.uniform(arefval[0],arefval[1])
                    arefMean[index] = .5*(arefval[0]+arefval[1])
                else:
                    arefTemp = arefval
                    arefMean[index] = 1.*arefval  

                        
                if ballisticrange==0:
                    validDebris = 1 #if ballistic coeff not a range then accept any value given by mass and aref
                    break
                if arefrange==0 and weightrange==0:
                    validDebris = 1# if area and mass not a range, then accept values
                    break
                ballisticMean = weightTemp/(CDaver*arefTemp)
            
#                print 'index = ' + str(index) + ', indDebris = ' + str(indDebris) + ', mass = ' + str(weightTemp) + ', ballisticMean = ' + str(ballisticMean) \
#                                + ', aref = ' + str(arefTemp) + ', CDaver = ' + str(CDaver)
            
                if ballisticMean<=ballCoeff[1] and ballisticMean>=ballCoeff[0]:
                    validDebris = 1
                    break
            weightarray[indDebris] = weightTemp
            arefarray[indDebris] = arefTemp

        weight[index] = weightarray
        aref[index] = arefarray
        velrange,velval = isRange(key[5])
        if velrange==1:
            velarray = np.random.uniform(velval[0],velval[1],nSamples)
        else:
            velarray = velval*np.ones(nSamples)
        vel[index] = velarray
         

                
        blast[index] = 0
        if len(key)>8:
            if float(key[8]) == 1: # blast on ground
                blast[index] = 1
    

    return (weight,aref,vel,ncdList,minfcdList,CD,nclList,minfclList,CL,arefMean,numberOfPiecesMean,blast)










def pieceRandomGroup(lineList,groupNumber,nSamples):
    #Generate random pieces for each group. Ignore the value created for the number of pieces (currently only the average value is used in Ec)
    # Note that groupNumber is >0
    # the debris catalog currently expects the element order to be as follows:
    # Number of Pieces-Weight of each piece(lbs)-Total Weight(lbs)-Reference Area(ft^2)-Ballistic Coeff(lb/ft^2)-Velocity Increment(ft/sec)-CD-CL	

    randomArray = np.zeros((nSamples,len(lineList[groupNumber-1])))
    line = lineList[groupNumber-1]
    for index1 in range(0,nSamples):
        for index2 in range(0,len(line)):
            value = line[index2]
            val,key=isRange(value)
            if val==1:
                lowerbound = float(key[0])
                upperbound = float(key[1])
                if upperbound<=lowerbound:
                    print 'Check bounds in Debris Catalog'
                    exit(1)
                key= np.random.uniform(float(key[0]),float(key[1]))
            randomArray[index1,index2]=key
    return randomArray
                
                     
def getCLCD(Cxval,filePATH):
    if is_number(Cxval):
        Cx = [float(Cxval)]
        minfcx = 1 # not actually used in this case, but a value needs to be pass to propagator
        ncx = 1
    else:
        minfcx,Cx,ncx = CLCDreader(filePATH+Cxval)

    return minfcx,Cx,ncx


def CLCDreader(filename):
    # Python routine to read CD/CL profile from a text file

    inputFile = open(filename,'r')
    # The format expected is Minf -> first column. Cd or Cl -> second column
    
    M = [] # Minf
    Cx = []# coefficient CL or CD 
    for line in inputFile:
        key = line.split()
        if len(key)>0:
            if (line[0]!='#'):# this way "#" can be used to comment an entire line 
                M.append(float(key[0]))
                Cx.append(float(key[1]))
    n = len(M) 
    inputFile.close()
    
    return (M,Cx,n)

