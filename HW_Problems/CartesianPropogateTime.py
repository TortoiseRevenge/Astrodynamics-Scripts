import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, '..')
sys.path.insert(0, parent_dir)

from CartesianToElements import cartesianToElements
from ElementsToCartesian import elementsToCartesian

from scipy.optimize import fsolve
from numpy import deg2rad as rad, rad2deg as deg, cos, sin, sqrt, pi, atan, tan

# Given Orbit, km and km/s
r = [-3381.5, -5885.6, 4374.1]
v = [6.1106, -1.1721, -3.1914]
#r = [-3670, -3870, 4400]
#v = [4.7, -7.4, 1]

# Earth Constants
J2 = 1.08262668e-3  # unitless
SGP = 3.986e5  # km^3/s^2
RE = 6378  # km

# Get Elements from Cartesian
a, e, i, omega, Omega, theta = cartesianToElements(r, v, SGP)
omega = rad(omega)
Omega = rad(Omega)

# Get the period of the orbit
n = sqrt(SGP / a**3)  # rad/s
T = 2 * pi / n  # seconds

# Get Change in Omega
dOmega = -(3 * J2 * n*RE**2*cos(rad(i)))/(2*a**2*(1-e**2)**2)  # rad/s
print("dOmega/day:", dOmega*86400)
domega = ((3*J2*n*RE**2)/(4*a**2*(1-e**2)**2))*(4-5*(sin(rad(i)))**2)  # rad/s
print("domega/day:", domega*86400)

# Find Position and Velocity after 7 days
tDiff = 4*86400  # seconds

# Calculate E 
E = 2*atan(sqrt((1-e)/(1+e))*tan(rad(theta)/2))
if(E < 0):
    E = 2*pi + E
t0 = (E - e*sin(E))/n  # seconds
t = t0 + tDiff  # seconds

t = t % T  # seconds
M = n * t # radians
E=fsolve(lambda E: E - e*sin(E) - M, E)[0]  # Solve for E

# Get Theta
theta = 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
if(theta < 0):
    theta = 2*pi + theta

Omega = Omega + dOmega*tDiff  # radians
omega = omega + domega*tDiff  # radians
# Get new elements
(r,v) = elementsToCartesian(a, e, i, deg(omega), deg(Omega), deg(theta), SGP)

print("New Position:", r)
print("New Velocity:", v)