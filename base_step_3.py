#Andy Jiang, Sherrill Group, Georgia Institute of Technology

#Equations are from J. Mol. Biol (1995) 251, 648-664 by M. A. El Hassan and C. R. Calladine

import numpy as np
import math

#The general rotation matrix for a rotation of angle (theta) about a 3-dimensional unit vector (u)
def R(theta, u):
    cos = math.cos(theta)
    sin = math.sin(theta)

    Ux = u[0]
    Uy = u[1]
    Uz = u[2]

    r1 = [cos + Ux*Ux*(1-cos), Ux*Uy*(1-cos) - Uz*sin, Ux*Uz*(1-cos) + Uy*sin]
    r2 = [Uy*Ux*(1-cos) + Uz*sin, cos + Uy*Uy*(1-cos), Uy*Uz*(1-cos) - Ux*sin]
    r3 = [Uz*Ux*(1-cos) - Uy*sin, Uz*Uy*(1-cos) + Ux*sin, cos + Uz*Uz*(1-cos)]

    return np.array([r1, r2, r3], dtype=float)

#x, y, and z are the direction vectors associated with the nucleobase
def Ti(x, y, z):
    return np.transpose(np.array([x, y, z], dtype = float))

#See El Hassan Equation 9
def Tiplus1(x, y, z, omega, gamma, phi):

    rz1 = R(omega/2 - phi, z)
    ry2 = R(gamma, y)
    rz3 = R(omega/2 + phi, z)

    result = np.dot(rz1, ry2)
    result = np.dot(result, rz3)
    result = np.dot(result, Ti(x, y, z))

    return result

#See El Hassan Equation 10
def Tmst(x, y, z, omega, gamma, phi):

    rz1 = R(omega/2.0 - phi, z)
    ry2 = R(gamma/2.0, y)
    rz3 = R(phi, z)

    result = np.dot(rz1, ry2)
    result = np.dot(result, rz3)
    result = np.dot(result, Ti(x, y, z))

    return result

#See El Hassan Equation 11
def newOrigin(oldOrigin, Xm, Ym, Zm, Dx, Dy, Dz):
    return np.array(oldOrigin, dtype=float) + Dx*Xm + Dy*Ym + Dz*Zm

#oldOrigin = the position vector of the first nucleobase's origin; x, y, and z are the direction vectors of the first nucleobase
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

    ti = Ti(x, y, z)

    Xm = tmst[:,0]
    Ym = tmst[:,1]
    Zm = tmst[:,2]

    new_origin = newOrigin(oldOrigin, Xm, Ym, Zm, Dx, Dy, Dz)
    
    #returns a tuple representing the vectors for the origin and x, y, and z direction vectors for the new nucleobase
    return (new_origin, tiplus1[:,0], tiplus1[:,1], tiplus1[:,2])

def calcNewTMST(oldOrigin, x, y, z, Dx, Dy, Dz, omega, rho, tau):
    omega = omega
    tau = tau
    rho = rho

    # The equations for gamma and phi are derived from rho = gamma*cos(phi) and tau = gamma*sin(phi) from El Hassan
    gamma = math.sqrt(rho*rho + tau*tau)

    if rho != 0:
        phi = math.atan(tau/rho)

    else:
        if tau > 0:
            phi = math.pi/2

        else:
            phi = -math.pi/2

    omega = omega
    gamma = gamma
    phi = phi

    #if math.cos(phi)*rho < 0 or math.sin(phi)*tau < 0:
    #    gamma = -1*gamma

    tiplus1 = Tiplus1(x, y, z, omega, gamma, phi)

    tmst = Tmst(x, y, z, omega, gamma, phi)

    return (np.array(oldOrigin), tmst[:,0], tmst[:,1], tmst[:,2])
