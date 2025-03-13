import numpy as np
from numpy import deg2rad as rad, cos, sin, sqrt

# Takes in as km, degrees, and km^3/s^2, and returns km and km/s
def elementsToCartesian(semiMajorAxis, eccentricity, inclination, argOfPeriapsis, longitudeOfPeriapsis, trueAnomaly, SGP):
    # Convert units to km, seconds, and radians, and go unitless
    a = semiMajorAxis
    e = eccentricity
    i = rad(inclination)
    omega = rad(argOfPeriapsis)
    Omega = rad(longitudeOfPeriapsis)
    theta = rad(trueAnomaly)
    mu = SGP
    
    rmag = a*(1-e**2)/(1+e*cos(theta))
   
    r_pqw = np.array([rmag * cos(theta), rmag * sin(theta), 0])
    h = sqrt(a * mu * (1-e**2))
    v_pqw = (mu/h) * np.array([-sin(theta), e+cos(theta),0])
    
    Omega_rot = np.array([[ cos(Omega), -sin(Omega), 0],
                  [ sin(Omega),  cos(Omega), 0],
                  [0, 0, 1]])
    
    i_rot = np.array([[1, 0, 0],
                [0, cos(i), -sin(i)],
                [0, sin(i),  cos(i)]])
    
    longOfPerRot = np.array([[cos(omega), -sin(omega), 0],
                    [sin(omega),  cos(omega), 0],
                    [0, 0, 1]])
    r = (Omega_rot @ i_rot @ longOfPerRot @ r_pqw.T)
    v = (Omega_rot @ i_rot @ longOfPerRot @ v_pqw.T)
    
    return(r,v)