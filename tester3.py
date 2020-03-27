import os
import get_coordinate_vectors
import random

#shift, slide, rise, tilt, roll, twist
get_coordinate_vectors.stack("stack1.txt")

os.system("find_pair -s ladder.pdb | analyze")
