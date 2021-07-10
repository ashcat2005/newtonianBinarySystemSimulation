#Define the differential equation that describes the newtonian dynamics.
def ODE(q0):
    '''
    ------------------------------------------
    ODE(q0)
    ------------------------------------------
    ODEs system for the motion of a comet
    around the Sun using cartesian coordinates
    in the orbital plane.
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
def VelocityVerletIntegrator(ODE,q):
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
    qNew = np.zeros((1,4)) #here will be the information of the new state
    #Calculate the new acceleration vector at t
    a0 = ODE(q) #Acceleration vector at t
    #Obtain the position at t+dt
    qNew[0:2] = q[0:2] + dt *q[2:4] + 0.5*a*dt**2
    #Calculate the acceleration at t+dt
    a1 = ODE(qNew)
    #Calculate the velocity at t+dt
    qNew[2:4] = q[2:4] + (a0+a1)/2 * dt
    return qNew

#Define the model for the emision of energy and angular momentum.
def EmissionModel(q,dE,dL):
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
    dL: change in angular momentum
        i.e. angular momentum emited.
    -----------------------------------
    Returns:
    qNew: numpy array with the updated state
    '''
    #Where the information of the updated state will be
    qNew = np.zeros((1,4))
    qNew[0:2] = q[0:2]
    #Define the emission matrix
    sigma = x*vx + y*vy
    R = np.array([[[x,-vy],[y,vx]]])
    A = 1/b *R #Emission matrix
    b = np.array([dE,dL]) #Emission vector
    #Update the velocities
    qNew[2:4] = A*b
    return qNew

def qij(q):
  '''
  ------------------------------------------------------
  qij(q)
  ------------------------------------------------------
  Returns the components of the 3x3 quadrupole tensor
  ------------------------------------------------------
  Arguments:
  q: numpy array with the state of
      the system in the order.
  q[0] = x: x coordinate.
  q[1] = y: y coordinate.
  q[2] = vx: velocity in x.
  q[3] = vy: velocity in y.

  '''
  return mu*(np.outer(x,x) - np.identity(3)*np.dot(x,x)/3)
