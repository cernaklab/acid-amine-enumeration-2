import itertools
from itertools import product
import numpy as np
import pickle

"""
enumerate, in this case, a 3-carbon, 2-oxygen, 1-nitrogen system
with no atom-level formal charge, and no octet rule violations.
Assume implicit hydrogens.

Matrices are printed to the command line, which will be later converted to .npy files.
"""

def get_tuple_sums(max_row_sum, max_any_bond_order, n_elements):
    """
    generates list of lists, where:
    1. The sum of each list <= max_row_sum,
    2. Each element <= max_any_bond_order, and
    3. Each list is n_elements long.

    """
    return [list(pdt_set) \
           for pdt_set in product(range(max_any_bond_order+1),repeat=n_elements) \
           if sum(pdt_set) <= max_row_sum]

def generate_adj_rows(LHS_elems, max_any_bond_order, row_valency,remaining_valencies):
    """
    generates all adjacency matrix rows, given:
    1. The elements already populated (LHS_elems)
    2. Maximum permitted bond order (max_any_bond_order)
    3. The total permitted valency (row_valency)

    global variables needed:
       a. n_atoms, the total atoms in the system
    """


    # find out how many valencies are left
    valency_remaining = row_valency - np.sum(LHS_elems)

    # find out how many elements left
    len_remaining = n_atoms - len(LHS_elems)


    # get the raw stuff
    sums = get_tuple_sums(max_row_sum=valency_remaining,\
                          max_any_bond_order=max_any_bond_order,\
                          n_elements=len_remaining)

    # append the left side
    sums = [LHS_elems + row for row in sums]

    # check for valency
    sums = [row for row in sums if np.all(np.less_equal(row,remaining_valencies))]

    return sums


# maximum valencies


# first row
# C C C C C C O O
# here's a dictionary: generate_adj_rows(LHS_elems, max_any_bond_order, row_valency,remaining_valencies):
# carbon: generate_adj_rows(LHS_elems, 3, 4, remaining_valencies)
# nitrogen: generate_adj_rows(LHS_elems, 3, 3, remaining_valencies)
# oxygen: generate_adj_rows(LHS_elems, 2, 2, remaining_valencies)
# boron: generate_adj_rows(LHS_elems, 1, 3, remaining_valencies) check
# bromine: generate_adj_rows(LHS_elems, 1, 1, remaining_valencies)
n_atoms = 8

max_valencies = [4,4,4,4,4,3,2,2]
# C
row1_all = generate_adj_rows([0], 3, 4, max_valencies)

for row1 in row1_all[:5]:
    # C
    row2_valencies = np.subtract(max_valencies,row1)
    row2_LHS = [row1[1], 0]
    row2_all = generate_adj_rows(row2_LHS, 3, 4, row2_valencies)

    for row2 in row2_all:
        # C
        row3_valencies = np.subtract(row2_valencies,row2)
        row3_LHS = [row1[2], row2[2], 0]
        row3_all = generate_adj_rows(row3_LHS, 3, 4, row3_valencies)

        for row3 in row3_all:
            # C
            row4_valencies = np.subtract(row3_valencies,row3)
            row4_LHS = [row1[3], row2[3], row3[3], 0]
            row4_all = generate_adj_rows(row4_LHS, 3, 4,row4_valencies)

            for row4 in row4_all:
                # C
                row5_valencies = np.subtract(row4_valencies,row4)
                row5_LHS = [row1[4], row2[4], row3[4], row4[4], 0]
                row5_all = generate_adj_rows(row5_LHS, 3, 4,row5_valencies)

                for row5 in row5_all:
                    # N
                    row6_valencies = np.subtract(row5_valencies,row5)
                    row6_LHS = [row1[5], row2[5], row3[5], row4[5], row5[5], 0]
                    row6_all = generate_adj_rows(row6_LHS, 3, 3, row6_valencies)

                    for row6 in row6_all:
                        # O
                        row7_valencies = np.subtract(row6_valencies,row6)
                        row7_LHS = [row1[6], row2[6], row3[6], row4[6], row5[6], row6[6], 0]
                        row7_all = generate_adj_rows(row7_LHS, 2, 2, row7_valencies)

                        for row7 in row7_all:
                            # O
                            row8 = [row1[7], row2[7], row3[7], row4[7], row5[7], row6[7], row7[7], 0]
                            adj_mat = [row1,row2,row3,row4,row5,row6,row7,row8]
                            print(adj_mat)
