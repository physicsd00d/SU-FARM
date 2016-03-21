import sys

mainPATH = '../../../source/'
###importing debris catalog scripts
sys.path.append(mainPATH+'DebrisCatalog/')
import debrisReader as dr
fileName = 'Halcon9_1st.dat'
dr.convertDebrisEng2SI(fileName)
fileName = 'Halcon9_2nd.dat'
dr.convertDebrisEng2SI(fileName)


