import numpy as np
import itertools as it


for z in it.product([2,-1,0,1,2], repeat=10):
    A = np.array([[z[0], z[1], z[2],z[3]],[z[1],z[4],z[5],z[6]],[z[2],z[5],z[7],z[8]],[z[3],z[6],z[8],z[9]]])
    print(A)
                 

