import numpy as np
import math

def Rx(theta):
    return np.array([[1.00, 0.00, 0.00], [0.00, math.cos(theta), -1*math.sin(theta)], [0.00, math.sin(theta), math.cos(theta)]], dtype=float)

def Ry(theta):
    return np.array([[math.cos(theta), 0.00, math.sin(theta)], [0.00, 1.00, 0.00], [-1*math.sin(theta), 0.00, math.cos(theta)]], dtype=float)

def Rz(theta):
    return np.array([[math.cos(theta), -1*math.sin(theta), 0.00], [math.sin(theta), math.cos(theta), 0.00], [0.00, 0.00, 1.00]], dtype=float)

def Ti(x, y, z):
    return np.transpose(np.array([x, y, z], dtype = float))

#x, y, and z are direction vectors

def Tiplus1(x, y, z, omega, gamma, phi):

    rz1 = Rz(omega/2.0 - phi)
    ry2 = Ry(gamma)
    rz3 = Rz(omega/2.0 + phi)

    result = np.dot(rz1, ry2)
    result = np.dot(result, rz3)
    result = np.dot(result, Ti(x, y, z))

    return result

def Tmst(x, y, z, omega, gamma, phi):

    rz1 = Rz(omega/2.0 - phi)
    ry2 = Ry(gamma/2.0)
    rz3 = Rz(phi)

    result = np.dot(rz1, ry2)
    result = np.dot(result, rz3)
    result = np.dot(result, Ti(x, y, z))

    return result

def Xmst(Tm):
    return Tm[:,0]

def Ymst(Tm):
    return Tm[:,1]

def Zmst(Tm):
    return Tm[:,2]

def newOrigin(oldOrigin, Xm, Ym, Zm, Dx, Dy, Dz):
    return np.array(oldOrigin, dtype=float) + Dx*Xm + Dy*Ym + Dz*Zm

def calcNewBP(oldOrigin, x, y, z, Dx, Dy, Dz, omega, rho, tau):
    
    gamma = math.sqrt(rho*rho + tau*tau)

    if rho != 0:
        phi = math.atan(tau/rho)

    else:
        if tau > 0:
            phi = math.pi/2

        else:
            phi = -math.pi/2

    if math.cos(phi)*rho < 0 or math.sin(phi)*tau < 0:
        gamma = -1*gamma

    tiplus1 = Tiplus1(x, y, z, omega, gamma, phi)

    tmst = Tmst(x, y, z, omega, gamma, phi)

    #ti = Ti(x, y, z)

    Xm = tmst[:,0]
    Ym = tmst[:,1]
    Zm = tmst[:,2]

    new_origin = newOrigin(oldOrigin, Xm, Ym, Zm, Dx, Dy, Dz)

    return (new_origin, tiplus1[:,0], tiplus1[:,1], tiplus1[:,2])

#print(calcNewBP([0,0,0], [1,0,0], [0,1,0], [0,0,1], 3, 0, 0, math.pi/2, 0, 0))
