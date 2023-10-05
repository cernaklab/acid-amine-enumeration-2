# structure matrices for each possible starting material hybridization combination
# ac3 = sp3 acid, ac2 = sp2 acid, am3 = sp3 amine, am2 = sp2 amine

ac3_am3 = [[0, 1, 0, 0, 0, 0, 0, 0],
           [1, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 1, 0, 0, 0, 0],
           [0, 0, 1, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 2, 1],
           [0, 1, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 2, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0]]


ac2_am3 = [[0, 1, 0, 0, 0, 0, 0, 0],
           [1, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 2, 0, 0, 0, 0],
           [0, 0, 2, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 2, 1],
           [0, 1, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 2, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0]]

ac3_am2 = [[0, 2, 0, 0, 0, 0, 0, 0],
           [2, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 1, 0, 0, 0, 0],
           [0, 0, 1, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 2, 1],
           [0, 1, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 2, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0]]

ac2_am2 = [[0, 2, 0, 0, 0, 0, 0, 0],
           [2, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 2, 0, 0, 0, 0],
           [0, 0, 2, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 2, 1],
           [0, 1, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 2, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0]]