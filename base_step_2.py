#Andy Jiang, Sherrill Group, Georgia Institute of Technology

#Equations are from J. Mol. Biol (1995) 251, 648-664 by M. A. El Hassan and C. R. Calladine

import numpy as np
import math

#Rx, Ry, and Rz represents rotations about the x, y, or z axis, respectively and are obtained from El Hassan
def Rx(theta):
    return np.array([[1.00, 0.00, 0.00], [0.00, math.cos(theta), -1*math.sin(theta)], [0.00, math.sin(theta), math.cos(theta)]], dtype=float)

def Ry(theta):
    return np.array([[math.cos(theta), 0.00, math.sin(theta)], [0.00, 1.00, 0.00], [-1*math.sin(theta), 0.00, math.cos(theta)]], dtype=float)

def Rz(theta):
    return np.array([[math.cos(theta), -1*math.sin(theta), 0.00], [math.sin(theta), math.cos(theta), 0.00], [0.00, 0.00, 1.00]], dtype=float)

#x, y, and z are the direction vectors associated with the nucleobase
def Ti(x, y, z):
    return np.transpose(np.array([x, y, z], dtype = float))

#See El Hassan Equation 9
def Tiplus1(x, y, z, omega, gamma, phi):

    rz1 = Rz(omega/2.0 - phi)
    ry2 = Ry(gamma)
    rz3 = Rz(omega/2.0 + phi)

    result = np.dot(rz1, ry2)
    result = np.dot(result, rz3)
    result = np.dot(result, Ti(x, y, z))

    return result

#See El Hassan Equation 10
def Tmst(x, y, z, omega, gamma, phi):

    rz1 = Rz(omega/2.0 - phi)
    ry2 = Ry(gamma/2.0)
    rz3 = Rz(phi)

    result = np.dot(rz1, ry2)
    result = np.dot(result, rz3)
    result = np.dot(result, Ti(x, y, z))

    return result

#See El Hassan Equation 11
def newOrigin(oldOrigin, Xm, Ym, Zm, Dx, Dy, Dz):
    return np.array(oldOrigin, dtype=float) + Dx*Xm + Dy*Ym + Dz*Zm

#oldOrigin = the position vector of the first nucleobase; x, y, and z are the direction vectors of the first nucleobase
#Dx = Shift, Dy = Slide, Dz = Rise, omega = Twist (z-rotation), rho = Roll (y-rotation), tau = Tilt (x-rotation)
def calcNewBP(oldOrigin, x, y, z, Dx, Dy, Dz, omega, rho, tau):
    
    # The equations for gamma and phi are derived from rho = gamma*cos(phi) and tau = gamma*sin(phi) from El Hassan
    gamma = math.sqrt(rho*rho + tau*tau)

    if rho != 0:
        phi = math.atan(tau/rho)

    else:
        if tau > 0:
            phi = math.pi/2

        else:
            phi = -math.pi/2

    #if math.cos(phi)*rho < 0 or math.sin(phi)*tau < 0:
    #    gamma = -1*gamma

    tiplus1 = Tiplus1(x, y, z, omega, gamma, phi)

    tmst = Tmst(x, y, z, omega, gamma, phi)

    #ti = Ti(x, y, z)

    Xm = tmst[:,0]
    Ym = tmst[:,1]
    Zm = tmst[:,2]

    new_origin = newOrigin(oldOrigin, Xm, Ym, Zm, Dx, Dy, Dz)
    
    #returns a tuple representing the vectors for the origin and x, y, and z direction vectors for the second nucleobase
    return (new_origin, tiplus1[:,0], tiplus1[:,1], tiplus1[:,2])

#print(calcNewBP([0,0,0], [1,0,0], [0,1,0], [0,0,1], 3, 0, 0, math.pi/2, 0, 0))
