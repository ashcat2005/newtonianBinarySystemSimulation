#simulation conditions:
import numpy as np
import functions.py
G = 4*np.pi**2 # Gravitational constant
c = 63197.8 # Speed of light in units of au/yr

# Masses
m1 = 10. # Solar masses
m2 = 7.  # Solar masses

M = m1 + m2 # Total mass (solar masses)
mu = m1*m2/M # Reduced mass
GM = G*M # gravitational parameter

# Radii
r1 = 2*G*m1/c**2 # Schwarzschild radius for m1
r2 = 2*G*m2/c**2 # Schwarzschild radius for m2
r_merge = r1 + r2 # Separation distance for fusion/collision

# Initial Values
E = -70. # energy
L = 50. # angular momentum
omega = np.pi/3 # argument of the pericenter


# Time grid definition
n = 50000 # number of steps
time = np.zeros(n)
dt = 1E-4 # timestep
stateArray = np.zeros([n,4]) #array with the states
Q = np.zeros([n,3,3]) # Quadrupole tensor

# Initial condition in the grid
Q[0] = qij(position[0,0:3])

# Critical radius to modify dt
r_crit = 1E-1
