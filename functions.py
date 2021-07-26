import numpy as np
import random


#Define the differential equation that describes the newtonian dynamics.
def ODE(q0):
    '''
    ------------------------------------------
    ODE(q0)
    ------------------------------------------
    ODEs system corresponding to the forces.
    ------------------------------------------
    Arguments:
    q0: NumPy array with the coordinates
        defined as.
    q0[0] = x: coordinate x.
    q0[1] = y: coordinate y.
    ------------------------------------------
    Returns:
    a = NumPy array with the components of the
        acceleration.
    '''
    r2 = q0[0]**2 + q0[1]**2
    a = - G*M*q0[0:2]/r2**(3/2)
    return a


#We will use a Velocity-Verlet integrator
def velocityVerletIntegrator(ODE,q,dt):
    '''
    ---------------------------------
    VelocityVerletIntegrator(ODE,q)
    ---------------------------------
    Uses a velocity Verlet integrator
    to obtain the velocity
    and position at the next step.
    ---------------------------------
    Arguments:
    ODE: Force function.
    q: numpy array with the state of
        the system in the order.
    q[0] = x: x coordinate.
    q[1] = y: y coordinate.
    q[2] = vx: velocity in x.
    q[3] = vy: velocity in y.
    ---------------------------------
    Returns:
    qNew: Numpy array with the updated state.
    '''
    qNew = np.zeros((4,)) #here will be the information of the new state
    #Calculate the new acceleration vector at t
    a0 = ODE(q) #Acceleration vector at t
    #Obtain the position at t+dt
    qNew[0:2] = q[0:2] + dt *q[2:4] + 0.5*a0* dt**2
    #Calculate the acceleration at t+dt
    a1 = ODE(qNew)
    #Calculate the velocity at t+dt
    qNew[2:4] = q[2:4] + (a0+a1)/2 * dt
    return qNew

#Define the model for the emision of energy and angular momentum.
def emissionModel(q,dE,dL):
    '''
    -----------------------------------
    EmissionModel(q,dE,dL)
    -----------------------------------
    Models the change in velocity after
    an emission of both energy and
    angular momentum. The details
    of the model and the derivation of the
    emission matrix can be found in the
    document emissionModel.pdf, in the
    theory directory of the repository.
    -----------------------------------
    Arguments:
    q: numpy array with the state of
        the system in the order.
    q[0] = x: x coordinate.
    q[1] = y: y coordinate.
    q[2] = vx: velocity in x.
    q[3] = vy: velocity in y.
    dE: change in energy i.e. energy emited.
    dL: change in angula r momentum
        i.e. angular momentum emited.
    -----------------------------------
    Returns:
    qNew: numpy array with the updated state
    '''
    #Where the information of the updated state will be
    qNew = np.zeros((4,))
    qNew[0:2] = q[0:2]
    #Define the emission matrix
    sigma = x*vx + y*vy
    R = np.array([[[x,-vy],[y,vx]]])
    A = 1/sigma *R #Emission matrix
    b = np.array([dE,dL]) #Emission vector
    #Update the velocities
    qNew[2:4] = np.matmul(A,b)
    return qNew

#Define the quadrupole tensor that describes the radiation
def qij(x):
  '''
  --------------------------------------------------
  qij(x)
  --------------------------------------------------
  Obtains the components of the quadrupole tensor.
  --------------------------------------------------
  Arguments:
  x: numpy array with the position vector.
  --------------------------------------------------
  Returns:
  Q: quadrupole tensor
  '''
  return mu*(np.outer(x,x) - np.identity(2)*np.dot(x,x)/3) #Notice that the dimension of the identity is 2, not 3.

#Taken from: https://github.com/ashcat2005/binary_coalescence/blob/main/main/main.py

#Define the Levi-Civita tensor
ciclicPermutation = [[0,1,2],[1,2,0],[2,0,1]]
leviCivita = np.zeros((3,3,3))
for i in range(3):
    for j in range(3):
        for k in range(3):
            if i!=j and j!=k and i!=k:
                if [i,j,k] in ciclicPermutation:
                    leviCivita[i,j,k]=1
                else: leviCivita[i,j,k]=-1
print(np.shape(leviCivita))
print(leviCivita)

#Functions that calculate second and third derivative.
def second_derivative(f,h):
  '''
  Second derivative of a function
  '''
  return (11*f[0] - 56*f[1] + 114*f[2] - 104*f[3] + 35*f[4])/(12*h**2)


def third_derivative(f,h):
  '''
  Third derivative of a function
  '''
  return (3*f[0] - 14*f[1] + 24*f[2] - 18*f[3] + 5*f[4])/(2*h**3)
#Taken from: https://github.com/ashcat2005/binary_coalescence/blob/main/main/main.py

#Define the loss of energy and angular momentum due to the GW emission.
def dEnergy(Q,dt):
    '''
     ------------------------------------
    dEnergy(Q)
    -------------------------------------
    Models the emission of energy in GW
    according to the theoretical results.
    -------------------------------------
    Arguments:
    Q: quadrupole tensor as a numpy array.
    -------------------------------------
    Returns:
    dE: change in energy
    '''
    a = third_derivative(Q,dt)
    b = 0
    for i in range(np.shape(a)[0]): #This can be more elegantly written as an iteration over the array elements
        for j in range(np.shape(a)[1]):
            b = b + a[i,j]**2
    dE = b*G/(5*c**2) * dt
    return dE

def dAngularMomentum(Q,dt):
    '''
     ------------------------------------
    dAngularMomentum(Q)
    -------------------------------------
    Models the emission of angular momentum
    due to GW according to the theoretical
    results.
    -------------------------------------
    Arguments:
    Q: quadrupole tensor as a numpy array.
    -------------------------------------
    Returns:
    dL: change in angular momentum
    '''
    z = 2
    a = second_derivative(Q,dt)
    b = 0
    for j in range(2):
        for k in range(2):
            for m in range(2): #This can be more elegantly written as an iteration over the array
                b = b + leviCivita[z,j,k]*a[j,m]*a[k,m]
    dL = 2*G/(5*c**5) * dt
    return dL
