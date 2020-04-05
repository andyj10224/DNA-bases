#Andy Jiang, Sherrill Group, Georgia Institute of Technology

import sys

from openbabel import openbabel
from openbabel import pybel

import numpy as np
from numpy import linalg
import math

import base_step_3
from base_step_3 import calcNewBP

#This function returns the direction vectors and origin of each nucleobase before their addition to the nucleotide ladder
def getOriginAndVectors():
    
    origin = np.array([0, 0, 0])
    x_vector = np.array([1, 0, 0])
    y_vector = np.array([0, 1, 0])
    z_vector = np.array([0, 0, 1])

    return (origin, x_vector, y_vector, z_vector)

#This function calculates the positions of the atoms in the molecule to be added (newMol) after the rotations are applied
#coordArr consists of all of the origins and direction vectors of each step of the ladder
def applyRotations(coordArr, newMol):
    
    #newCoords will be the origin and direction vectors of the new base pair after the transformations are applied
    newCoords = coordArr[-1]

    origin = newCoords[0]
    
    #firstCoords is the origin and direction vectors of the nucleobase before the transformations are applied
    firstCoords = coordArr[0]

    newMol.SetChainsPerceived()

    for atom in openbabel.OBMolAtomIter(newMol):

        atomVector = np.transpose(np.array([atom.GetX(), atom.GetY(), atom.GetZ()], dtype=float)) #- firstCoords[0]

        length = np.linalg.norm(atomVector)

        oldTransMatrix = np.array([firstCoords[1], firstCoords[2], firstCoords[3]])
        newTransMatrix = np.transpose(np.array([newCoords[1], newCoords[2], newCoords[3]], dtype=float))

        transMatrix = np.matmul(newTransMatrix, oldTransMatrix)

        superVector = np.matmul(transMatrix, atomVector)

        atom.SetVector(superVector[0], superVector[1], superVector[2])
        
    #Returns the nucleobase, after the rotations have been applied
    return newMol

#This function calculates the positions of the atoms in the molecule to be added (newMol) after the translations are applied
#This function will be applied after the rotations are applied
#coordArr consists of all of the origins and direction vectors of each step of the ladder
def applyTranslations(coordArr, newMol):

    newMol.SetChainsPerceived()
    
    #newCoords will be the origin and direction vectors of the new base pair after the transformations are applied
    newCoords = coordArr[-1]
    
    #oldCoords consists of the origin and direction vectors of the nucleobase before the transformations are applied
    oldCoords = coordArr[0]

    oldTransMatrix = np.array([oldCoords[1], oldCoords[2], oldCoords[3]])
    newTransMatrix = np.transpose(np.array([newCoords[1], newCoords[2], newCoords[3]]))

    transMatrix1 = np.matmul(newTransMatrix, oldTransMatrix)

    translation = newCoords[0] - np.matmul(transMatrix1, oldCoords[0])

    for atom in openbabel.OBMolAtomIter(newMol):
        atom.SetVector(atom.GetX() + translation[0], atom.GetY() + translation[1], atom.GetZ() + translation[2])
        
    #Returns the nucleobase, after the translations are applied
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
#i = the order of the nucleobase in the stack, coordArr = the list of origin-direction vector tuples corresponding to each step of the nucleobase ladder
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

    oldCoords = getOriginAndVectors()

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
