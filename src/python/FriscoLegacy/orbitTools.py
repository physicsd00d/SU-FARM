# routines that do different coordinate transformation for orbit propagation
import numpy as np

def latlonalt2ECEF(lat,lon,hs,planetModel):
# model can only be used for Earth
# inputs in degrees and in meters.
# planetModel = 0 for circular. 1 for oblate
    Req = 6378137.0 #m
    Rmean = 6371000.8 # m
    ecc = 0.081819221456
    phigc = lat*np.pi/180.
    lam = lon*np.pi/180.
    if planetModel==0:
        x = (Rmean + hs)*np.cos(phigc)*np.cos(lam)
        y = (Rmean + hs)*np.cos(phigc)*np.sin(lam)
        z = (Rmean + hs)*np.sin(phigc)
    elif planetModel==1:
        Nphi = Req/(1.-(ecc*np.sin(phigc))**2)**.5
        x = (Nphi + hs)*np.cos(phigc)*np.cos(lam)
        y = (Nphi + hs)*np.cos(phigc)*np.sin(lam)
        z = (Nphi*(1.-ecc**2)+hs)*np.sin(phigc)
    return [x,y,z] # results in meters ECEF coordinate


def latlonalt2ECI(lat,lon,hs,thetag,planetModel):
    [x,y,z] = latlonalt2ECEF(lat,lon,hs,planetModel)
    rvec = np.array([[x,y,z]]).T
    C = ecef_C_eci(thetag)
    r_eci = np.dot(C.T,rvec)
    x_eci = r_eci[0,0]
    y_eci = r_eci[1,0]
    z_eci = r_eci[2,0]
    return [x_eci,y_eci,z_eci]


def ECI2latlonalt(r_eci,thetag,planetModel):
    dim= np.shape(r_eci)
    C = ecef_C_eci(thetag)
    r_ecef = np.dot(C,r_eci)
    ret = ECEF2latlonalt(r_ecef,planetModel)
    return ret

def ECEF2latlonalt(r,planetModel):
    
    Req = 6378137.0 #m
    Rmean = 6371000.8 # m
    ecc = 0.081819221456
    x = r[0]
    y = r[1]
    z = r[2]
    r_delta= np.sqrt(x**2+y**2)
    lam = np.arctan2(y,x)
    phigd_m1 = np.arcsin(z/(x**2+y**2+z**2)**.5)
    tol = 1.0
    if planetModel==1:
        
        while tol>10.**-6.:
            Nphi = Req/(1.- (ecc*np.sin(phigd_m1))**2.)**.5
            phigd = np.arctan((z+Nphi*np.sin(phigd_m1)*(ecc**2))/(r_delta))
            tol = np.abs(phigd-phigd_m1)
            phigd_m1 = phigd
        he = r_delta/np.cos(phigd) - Nphi
    else:
        phigd = phigd_m1
        he = (x**2+y**2+z**2)**.5 - Rmean
    latitude = phigd*180./np.pi
    longitude = lam*180./np.pi
    return ([latitude,longitude,he])



def SEZ2ECEF(lat,lon,r_ecef,Vsouth,Veast,Vzenith,inertial):
    # Vzenith = Vup, lat,lon in degrees with reference to the vehicle current location in ECEF frame
    # lat lon are calculated from r_ecef (vehicle location), lat lon are inputs to avoid double calculations 
    # inertial = 1. input velocities are inertial velocities in SEZ frame
    # inertial =0. input velocities are local velocities (earth rotation not accounted)
    V = np.array([[Vsouth,Veast,Vzenith]])
    V = V.T
    C = ecef_C_sez(lat*np.pi/180.,lon*np.pi/180)
    Vecef = (np.dot(C,V)).T[0] # getting a nice 3 element array instead of a 3x1 matrix

    if inertial==0:
        omega = 2.0*np.pi/(86164.0906)#7.292115856e-5#rad/sec
        Vecef = Vecef + np.cross([0,0,omega],r_ecef)
    return Vecef
'''
def SEZ2ECI(lat,lon,r_ecef,Vsouth,Veast,Vzenith,inertial,thetag):
    Vecef = SEZ2ECEF(lat,lon,r_ecef,Vsouth,Veast,Vzenith,inertial)
    C = ecef_C_eci(thetag)
    Veci = np.dot(C.T,Vecef)
    return Veci
''' 
def SEZ2ECI(lat,lon,hs,Vsouth,Veast,Vzenith,inertial,thetag,planetModel):
    r_ecef = latlonalt2ECEF(lat,lon,hs,planetModel)
    Vecef = SEZ2ECEF(lat,lon,r_ecef,Vsouth,Veast,Vzenith,inertial)
    C = ecef_C_eci(thetag)
    Veci = np.dot(C.T,Vecef)
    return Veci
def ecef_C_sez(phigc,lam):
    C = np.array([[np.sin(phigc)*np.cos(lam),-np.sin(lam),np.cos(phigc)*np.cos(lam)],
                  [np.sin(phigc)*np.sin(lam),np.cos(lam),np.cos(phigc)*np.sin(lam)],
                  [-np.cos(phigc),0.0,np.sin(phigc)]])
    return C

def ecef_C_eci(thetag):
    C = np.array([[np.cos(thetag),np.sin(thetag),0.0],
                  [-np.sin(thetag),np.cos(thetag),0.0],
                  [0.,0.,1.]])
    return C

def cRnMatrix_Xaxis(theta):
# some rotation matrices ...transform from the n frame to c frame....useful for some frame workarounds
	R = np.array([[1,0,0],
				  [0,np.cos(theta),np.sin(theta)],
				  [0,-np.sin(theta),np.cos(theta)]])
	return R

def cRnMatrix_Yaxis(theta):
    # some rotation matrices ...transform from the n frame to c frame....useful for some frame workarounds
	R = np.array([[np.cos(theta), 0, np.sin(theta)],
                  [0,1,0],
				  [-np.sin(theta), 0, np.cos(theta)]])
	return R

def cRnMatrix_Zaxis(theta):
# some rotation matrices ...transform from the n frame to c frame....useful for some frame workarounds
	R = np.array([[np.cos(theta),np.sin(theta),0],
				  [-np.sin(theta),np.cos(theta),0],
				  [0,0,1]])
	return R
