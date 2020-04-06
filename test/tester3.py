import os
import get_coordinate_vectors
import random

testResults = open('testResults.txt', 'w+')

def numToNucleotide(num):
    if num == 0:
        return 'A'
    elif num == 1:
        return 'U'
    elif num == 2:
        return 'G'
    elif num == 3:
        return 'C'

for i in range(1000):
    stackFile = open('stack1.txt', 'w+')
    stackFile.write('Base, shift, slide, rise, tilt, roll, twist\n')
    stackFile.write('G 0 0 0 0 0 0\n')
    paramArr = []
    for j in range(4):
        #shift, slide, rise, tilt, roll, twist
        shift = random.randrange(-100, 100, 1)/10
        slide = random.randrange(-100, 100, 1)/10
        rise = random.randrange(-100, 100, 1)/10
        tilt = random.randrange(-180, 181, 1)
        roll = random.randrange(-180, 181, 1)
        twist = random.randrange(-180, 181, 1)
        nt = numToNucleotide(random.randrange(0, 4, 1))

        stackFile.write(nt + " " + str(shift) + " " + str(slide) + " " + str(rise) + " " + str(tilt) + " " + str(roll) + " " + str(twist) + "\n")

        paramArr.append([shift, slide, rise, tilt, roll, twist])
    
    stackFile.close()

    get_coordinate_vectors.stack("stack1.txt")
    os.system("find_pair -s ladder.pdb | analyze")

    outputFile = open("ladder.outs", "r")

    lines = outputFile.readlines()

    outputFile.close()

    for i in range(len(lines)):
        if 'Local base step parameters' in lines[i]:
            lineNo = i + 2

    for i in range(len(paramArr)):
        elements = lines[lineNo + i].split()

        testResults.write(F'{paramArr[i][0] - float(elements[2])}, {paramArr[i][1] - float(elements[3])}, {paramArr[i][2] - float(elements[4])}, {paramArr[i][3] - float(elements[5])}, {paramArr[i][4] - float(elements[6])}, {paramArr[i][5] - float(elements[7])}')
        testResults.write("\n")

testResults.close()
