#Andy Jiang, Sherrill Group, Georgia Institute of Technology

import sys

from openbabel import openbabel
from openbabel import pybel

import numpy as np
from numpy import linalg
import math

import base_step_3
from base_step_3 import calcNewBP
from base_step_3 import calcNewTMST

#This function allows the calculation of the direction vectors (x, y, and z) of a given molecule (OBMol object)
def getOriginAndVectors(mol, origin):

    atom1 = None
    atom2 = None

    origin = np.array(origin, dtype=float)

    numRings = GetNumRings(mol)

    if numRings == 1:
        for res in openbabel.OBResidueIter(mol):
            for atom in openbabel.OBResidueAtomIter(res):
                name = res.GetAtomID(atom)
                if name.find("N3") >= 0:
                    atom1 = atom
                elif name.find("C6") >= 0:
                    atom2 = atom

                

    elif numRings == 2:
        for res in openbabel.OBResidueIter(mol):
            for atom in openbabel.OBResidueAtomIter(res):
                name = res.GetAtomID(atom)
                if name.find("N1") >= 0:
                    atom1 = atom
                elif name.find("C4") >= 0:
                    atom2 = atom

    y_vector = [atom1.GetX(), atom1.GetY(), atom1.GetZ()]

    y_vector = np.array(y_vector, dtype=float)

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
    atom2_vector = [atom2.GetX(), atom2.GetY(), atom2.GetZ()]

    atom1_vector = np.array(atom1_vector, dtype=float)
    atom2_vector = np.array(atom2_vector, dtype=float)

    origin = (atom1_vector + atom2_vector)/2

    origin = np.array([0, 0, 0])
    x_vector = np.array([1, 0, 0])
    y_vector = np.array([0, 1, 0])
    z_vector = np.array([0, 0, 1])

    return (origin, x_vector, y_vector, z_vector)

#This function calculates the positions of the atoms in the molecule to be added (newMol)
#oldCoords is a tuple representing the origin, x-direction vector, y-direction vector, and z-direction vector (in that order)
#oldCoords is from the molecule added before this new molecule is to be added
def applyRotations(coordArr, newMol):
    
    #newCoords is a tuple containing the origin and direction vectors in the same order as old tuple
    #newCoords will be the direction vectors of the new base pair after the translations and rotations are applied
    newCoords = coordArr[-1]

    origin = newCoords[0]

    firstCoords = coordArr[0]

    newMol.SetChainsPerceived()

    for atom in openbabel.OBMolAtomIter(newMol):

        atomVector = np.transpose(np.array([atom.GetX(), atom.GetY(), atom.GetZ()], dtype=float)) #- firstCoords[0]

        length = np.linalg.norm(atomVector)

        #preTransMatrix = np.array([firstCoords[1], firstCoords[2], firstCoords[3]], dtype=float)

        oldTransMatrix = np.array([firstCoords[1], firstCoords[2], firstCoords[3]])
        newTransMatrix = np.transpose(np.array([newCoords[1], newCoords[2], newCoords[3]], dtype=float))

        transMatrix = np.matmul(newTransMatrix, oldTransMatrix)

        superVector = np.matmul(transMatrix, atomVector)

        #superVector = np.matmul(preTransMatrix, superVector)

        #translation = np.array([np.dot(origin, newCoords[1]), np.dot(origin, newCoords[2]), np.dot(origin, newCoords[3])])
        #translation = oldCoords[0] + Dx*tmst[1] + Dy*tmst[2] + Dz*tmst[3]
        #translation = np.array([np.dot(translation, newCoords[1]), np.dot(translation, newCoords[2]), np.dot(translation, newCoords[3])])

        #translation = np.matmul(transMatrix, translation)
        #translation = np.matmul(preTransMatrix, translation)

        #superVector += translation

        atom.SetVector(superVector[0], superVector[1], superVector[2])
        
    #Returns the nucleobase, after the translations and rotations have been applied
    return newMol

def applyTranslations(coordArr, newMol):

    newMol.SetChainsPerceived()

    newCoords = coordArr[-1]
    oldCoords = coordArr[0]

    oldTransMatrix = np.array([oldCoords[1], oldCoords[2], oldCoords[3]])
    newTransMatrix = np.transpose(np.array([newCoords[1], newCoords[2], newCoords[3]]))

    transMatrix1 = np.matmul(newTransMatrix, oldTransMatrix)

    translation = newCoords[0] - np.matmul(transMatrix1, oldCoords[0])

    for atom in openbabel.OBMolAtomIter(newMol):
        atom.SetVector(atom.GetX() + translation[0], atom.GetY() + translation[1], atom.GetZ() + translation[2])

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
def build(i, coordArr, newMol):

    newMol = applyRotations(coordArr, newMol)

    newMol = applyTranslations(coordArr, newMol)

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
    stackFile.close()

    coordArr = []

    firstMol = openbabel.OBMol()
    conv.ReadFile(firstMol, letterToNucleotide(stackLines[1].split()[0]))

    oldCoords = getOriginAndVectors(firstMol, [0, 0, 0])

    #ladder is the stack of nucleobases that will be returned by this function
    ladder = openbabel.OBMol()

    for i in range(1, len(stackLines)):
        line = stackLines[i].split()

        Dx = float(line[1])
        Dy = float(line[2])
        Dz = float(line[3])
        Rx = float(line[4])
        Ry = float(line[5])
        Rz = float(line[6])

        coordArr.append(calcNewBP(oldCoords[0], oldCoords[1], oldCoords[2], oldCoords[3], Dx, Dy, Dz, math.radians(Rz), math.radians(Ry), math.radians(Rx)))
        oldCoords = coordArr[-1]

        tempMol = openbabel.OBMol()

        conv.ReadFile(tempMol, letterToNucleotide(line[0]))

        tempMol = build(i, coordArr, tempMol)

        ladder += tempMol
        
        ladder.SetChainsPerceived()

    conv.WriteFile(ladder, "ladder.pdb")
