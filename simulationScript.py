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
