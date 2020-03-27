#Andy Jiang, Sherrill Group, Georgia Institute of Technology

import sys

from openbabel import openbabel
from openbabel import pybel

import numpy as np
from numpy import linalg
import math

import base_step_2
from base_step_2 import calcNewBP

#This function allows the calculation of the direction vectors (x, y, and z) of a given molecule (OBMol object)
def getOriginAndVectors(mol, origin):

    atom1 = None
    #atom2 = None

    origin = np.array(origin, dtype=float)

    numRings = GetNumRings(mol)

    if numRings == 1:
        for res in openbabel.OBResidueIter(mol):
            for atom in openbabel.OBResidueAtomIter(res):
                name = res.GetAtomID(atom)
                if name.find("C6") >= 0:
                    atom1 = atom
                #elif name.find("C6") >= 0:
                #    atom2 = atom

                

    elif numRings == 2:
        for res in openbabel.OBResidueIter(mol):
            for atom in openbabel.OBResidueAtomIter(res):
                name = res.GetAtomID(atom)
                if name.find("C8") >= 0:
                    atom1 = atom
                #elif name.find("C4") >= 0:
                #    atom2 = atom

    y_vector = [atom1.GetX(), atom1.GetY(), atom1.GetZ()]

    y_vector = np.array(y_vector, dtype=float) - origin

    y_vector = y_vector/(np.linalg.norm(y_vector))

    atom3 = None

    for res in openbabel.OBResidueIter(mol):
        for atom in openbabel.OBResidueAtomIter(res):
            name = res.GetAtomID(atom)
            if name.find("C1'") >= 0:
                atom3 = atom

    temp_vector = [atom3.GetX() - atom1.GetX(), atom3.GetY() - atom1.GetY(), atom3.GetZ() - atom1.GetZ()]

    temp_vector = np.array(temp_vector, dtype=float)

    z_vector = -1*np.cross(temp_vector, y_vector)

    z_vector = z_vector/np.linalg.norm(z_vector)

    x_vector = np.cross(y_vector, z_vector)

    atom1_vector = [atom1.GetX(), atom1.GetY(), atom1.GetZ()]
    #atom2_vector = [atom2.GetX(), atom2.GetY(), atom2.GetZ()]

    atom1_vector = np.array(atom1_vector, dtype=float)
    #atom2_vector = np.array(atom2_vector, dtype=float)

    #origin = (atom1_vector + atom2_vector)/2

    return (origin, x_vector, y_vector, z_vector)

#This function calculates the positions of the atoms in the molecule to be added (newMol)
#oldCoords is a tuple representing the origin, x-direction vector, y-direction vector, and z-direction vector (in that order)
#oldCoords is from the molecule added before this new molecule is to be added
def getNewBPFromMol(oldCoords, newMol, Dx, Dy, Dz, omega, rho, tau):
    
    #newCoords is a tuple containing the origin and direction vectors in the same order as old tuple
    #newCoords will be the direction vectors of the new base pair after the translations and rotations are applied
    newCoords = calcNewBP(oldCoords[0], oldCoords[1], oldCoords[2], oldCoords[3], Dx, Dy, Dz, omega, rho, tau)

    newMol.SetChainsPerceived()
    
    #oldCoords2 represents the origin and direction vectors of the base pair to be added before the translations and rotations are applied
    oldCoords2 = getOriginAndVectors(newMol, np.array([0, 0, 0], dtype=float))

    for atom in openbabel.OBMolAtomIter(newMol):

        atom.SetVector(atom.GetX() - oldCoords2[0][0], atom.GetY() - oldCoords2[0][1], atom.GetZ() - oldCoords2[0][2])

        atomVector = np.transpose(np.array([atom.GetX(), atom.GetY(), atom.GetZ()], dtype=float))

        A = np.dot(atomVector, oldCoords2[1])
        B = np.dot(atomVector, oldCoords2[2])
        C = np.dot(atomVector, oldCoords2[3])

        superVector = A*newCoords[1] + B*newCoords[2] + C*newCoords[3]

        superVector = superVector + newCoords[0]

        atom.SetVector(superVector[0], superVector[1], superVector[2])\
        
    #Returns the nucleobase, after the translations and rotations have been applied
    return newMol

#Returns the number of rings in a molecule, useful for checking if a nucleobase is a pyrimidine or a purine
def GetNumRings(mol):
    count = 0

    for ring in openbabel.OBMolRingIter(mol):
        count += 1

    return count

#Prints the name of each atom in a molecule, useful for debugging purposes
def printAtomNames(mol):
    for res in openbabel.OBResidueIter(mol):
        for atom in openbabel.OBResidueAtomIter(res):
            name = res.GetAtomID(atom)
            print(name)

#This function returns a nucleobase after the translations and rotations have been applied, also sets the chain of the nucleobase
#i = the order of the nucleobase in the stack, oldCoords = the origin and direction vectors of nucleobase i-1
#newMol = the nucleobase to be added
def build(i, oldCoords, newMol, shift, slide, rise, tilt, roll, twist):

    newMol = getNewBPFromMol(oldCoords, newMol, shift, slide, rise, math.radians(twist), math.radians(roll), math.radians(tilt))

    newMol.GetResidue(0).SetChain(chr(ord('A') + i - 1))

    newMol.SetChainsPerceived()

    return newMol

#Returns the name of a nucleobase file based on its letter (A, U, T, G, C)
def letterToNucleotide(letter):
    if letter == 'A':
        return "adenine.pdb"
    elif letter == 'U':
        return "uracil.pdb"
    elif letter == 'T':
        return "thymine.pdb"
    elif letter == 'G':
        return 'guanine.pdb'
    elif letter == 'C':
        return "cytosine.pdb"

#This function stacks nucleotides based on inputs from a file (stackFile is the file name)
def stack(stackFile):

    conv = openbabel.OBConversion()

    ladder = openbabel.OBMol()

    stackFile = open(stackFile, "r")

    stackLines = stackFile.readlines()

    oldCoords = None

    molArr = None
    
    #ladder is the stack of nucleobases that will be returned by this function
    ladder = openbabel.OBMol()

    firstMol = openbabel.OBMol()

    firstLine = stackLines[1].split()

    conv.ReadFile(firstMol, letterToNucleotide(firstLine[0]))

    oldCoords = getOriginAndVectors(firstMol, [0, 0, 0])

    for i in range(1, len(stackLines)):
        line = stackLines[i].split()

        Dx = float(line[1])
        Dy = float(line[2])
        Dz = float(line[3])
        Rx = float(line[4])
        Ry = float(line[5])
        Rz = float(line[6])

        tempMol = openbabel.OBMol()

        conv.ReadFile(tempMol, letterToNucleotide(line[0]))
        
        #In each iteration, the new nucleobase is added to ladder
        ladder += build(i, oldCoords, tempMol, Dx, Dy, Dz, Rx, Ry, Rz)

        ladder.SetChainsPerceived()
        
        #oldCoords becomes updated with the coordinates of the new nucleobase, so the next nucleobase after it can be adjusted
        oldCoords = calcNewBP(oldCoords[0], oldCoords[1], oldCoords[2], oldCoords[3], Dx, Dy, Dz, math.radians(Rz), math.radians(Ry), math.radians(Rx))


    conv.WriteFile(ladder, "ladder.pdb")
    stackFile.close()
