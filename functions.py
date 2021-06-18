import numpy as np
import random
#Define functions that calculate the energy and the angular momentum of the system
def E(r,vr,vphi):
    '''Takes the radius, radial velocity and angular velocity at a given time and returns the energy.
        r :radius
        vr: radial velocity
        vphi: angular velocity.
        As the problem has axial symmetry, there is no need for the angular position.'''
    T = 0.5 *(vr**2 + (r*vphi)**2)
    U = alpha/r
    return T+U

def L(r,vphi):
    '''Takes the radius and angular velocity at a given time and returns the angular momentum
        r : radius
        vphi: angular velocity'''
    return r**2 * vphi

def r(phi,p,e):
    '''Takes the angular position, eccentricity and semi-latus rectum, returns the radial position
    p : semi-latus rectum
    e: eccenticity
    phi:angular position
    '''
    r = p/(1-e*np.cos(phi))
    return r

#Define how the above variables are calculated.
def dphi(vphi,dt):
    '''Takes the instantaneous rate of change and the step size, returns the total change'''
    return vphi * dt

def r(phi,p,e):
    '''Takes the angular position, eccentricity and semi-latus rectum, returns the radial position
    p : semi-latus rectum
    e: eccenticity
    phi:angular position
    '''
    r = p/(1-e*np.cos(phi))
    return r

def vphi(r,L):
    '''takes  the angular momentum  and radius, uses the angular momentum definition
    to obtain the angular velocity at a given r.
    vphi : angular velocity
    L : angular momentum
    r: radius'''
    return L/r**2

def vr(r,phi,vphi,p,e):
    '''Takes the radial and angular position, semi-latus rectum, eccentricity and the
    angular velocity at a given time. Returns the radial velocity.
    r: radial position
    phi: angular position
    vphi: angular velocity
    p: semi-latus rectum
    e: eccentricity
    '''
    return 1/p * r**2 * e * np.sin(phi) * vphi

#Define the function that updates the state of the system
def updateState(q,dt):
    '''Takes as input an array and a time step. Returns the updated state of the system
    as an array'''
    qUpdated = np.zeros(4)
    #Calculate the new phi
    phiNew = dphi(q[3],dt) + q[1]
    qUpdated[1] = phiNew
    #Calculate the new r
    rNew = r(q[1],p,e)
    qUpdated[0] = rNew
    #New phi dot
    vphiNew = vphi(rNew,L)
    qUpdated[3] = vphiNew
    #New r dot
    vrNew = vr(rNew,phiNew,vphiNew,p,e)
    qUpdated[2] = vrNew
    return qUpdated

#Define how to pass from polar to cartesian
def polarToCartesian(r,theta):
    '''Takes as input polar coordinates and returns cartesian coordinates'''
    x = np.cos(theta) * r
    y = np.cos(theta) * r
    return x,y
