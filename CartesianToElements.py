import numpy as np
from numpy import rad2deg as deg, cos, sin, sqrt

# Takes in as km, km/s, and km^3/s^2, returns in km and degrees as (Semimajor axis, eccentricity, inclination, argOfPeri, longOfPeri, trueAnomaly).
def cartesianToElements(r, v, SGP):
    # Store norm of position and velocity vectors
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)
    mu = SGP
    
    # Angular Momentum
    h=np.cross(r,v)
    hmag = np.linalg.norm(h)
    
    # Calculate i
    i = np.arccos(h[2]/hmag)
    
    # Calculate e
    eVec = np.cross(v,h)/mu - r/rmag
    e = np.linalg.norm(eVec)
    
    # Calculate a
    a = hmag**2/(mu*(1-e**2))
    
    # Line to AN
    n = np.cross([0,0,1],h)
    
    # Get Omega
    Omega = 0
    if(n[1] >= 0):
        Omega = np.arccos(n[0]/np.linalg.norm(n))
    else:
        Omega = 2*np.pi - np.arccos(n[0]/np.linalg.norm(n))
    
    # Get omega
    omega = 0
    if(eVec[2] >= 0):
        omega = np.arccos(np.dot(n,eVec)/(np.linalg.norm(n)*e))
    else:
        omega = 2*np.pi - np.arccos(np.dot(n,eVec)/(np.linalg.norm(n)*e))
        
    # Get true anomaly
    theta = 0
    if(np.dot(eVec,r) >= 0):
        theta = np.arccos(np.dot(eVec,r)/(e*rmag))
    else:
        theta = 2*np.pi - np.arccos(np.dot(eVec,r)/(e*rmag))
    
    # Convert to degrees and return
    return(a,e,deg(i),deg(omega),deg(Omega),deg(theta))
