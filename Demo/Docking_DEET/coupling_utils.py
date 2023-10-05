from rdkit import Chem
from rdkit.Chem.AllChem import ReactionFromSmarts

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import numpy as np


# this is very janky to get past multiprocessing not working properly on ipynb. 

am2 = Chem.MolFromSmarts("[CX3]=[CX3][NH2]")
ac2 = Chem.MolFromSmarts("[CX3]=[CX3][CX3](=O)[OX2H1]")
am3 = Chem.MolFromSmarts("[13C][13C][N]")
ac3 = Chem.MolFromSmarts("CCC(=O)O")





group_dict = {"ac2":ac2, "ac3":ac3,"am2":am2,"am3":am3}


# number of entries in the amat
amat_indices = range(8)
# which entry each index should map to
mapping_order = [0,1,3,4,5,2,6,7]
# mapping_order = [3,4,5,0,1,2,6,7]

# amat index -> amine/acid match index
amat2ind = {k:v for k,v in zip(amat_indices, mapping_order)}
ind2amat = {v:k for k,v in zip(amat_indices, mapping_order)}

bond_dict = {1:Chem.BondType.SINGLE, 2:Chem.BondType.DOUBLE, 3:Chem.BondType.TRIPLE}

starting_numHs = np.array([3, 2, 1, 0, 0, 1, 0, 1])

def apply_amat(acid,acid_match,amine,amine_match,rmat):
    """
    needs the smarts of the acid and amine pre-defined. maybe an external dictionary.
    """
    # combine molecules and make a writable version
    both = Chem.CombineMols(acid,amine)
    bothW = Chem.RWMol(both)
    Chem.Kekulize(both,clearAromaticFlags=True)
    Chem.Kekulize(bothW,clearAromaticFlags=True)
    
    # get the indices for the acid and amine matches. 
    amine_atoms = both.GetSubstructMatch(group_dict[amine_match])
    acid_atoms = both.GetSubstructMatch(group_dict[acid_match])
    
    if not (amine_atoms and acid_atoms):
        print("substruct match not found")
        return
    # join the indices
    atom_list = amine_atoms + acid_atoms
#     print(atom_list)
    
    # go through the amat to find places where the change is not 0
    side_len = len(rmat)
    
    if side_len != len(atom_list):
        print("matrix length does not equal number of matched atoms")
        return
    
    for r in range(side_len):
        for c in range(r+1,side_len):
            bond_order_change = rmat[r][c]
            if bond_order_change != 0: 
                atom1 = atom_list[amat2ind[r]]
                atom2 = atom_list[amat2ind[c]]
                
                # get the current bond order
                current_bond = both.GetBondBetweenAtoms(atom1,atom2)
                if current_bond: 
                    current_bond_order = current_bond.GetBondTypeAsDouble()                
                else: 
                    current_bond_order = 0
                    
                new_bond_order = current_bond_order + bond_order_change
                
                if new_bond_order not in [0,1,2,3]:
                    print("invalid new bond order")
                    return None
                # make bond changes
                
                bothW.RemoveBond(atom1,atom2)
                
                if new_bond_order > 0:
                    bothW.AddBond(atom1,atom2, bond_dict[new_bond_order])

#     try:
#         Chem.SanitizeMol(bothW)

#     except:
#         print("illegal structure")
#         return None

    # may need to become smiles if we want to multitarget, as well as have fragments included.
    return bothW



def run_enumeration(file_tag, acid, amine,side):
    index_start = file_tag * 100000
    file_tag = str(file_tag).zfill(2)

    dmats = np.load(f"../rxn_mats/dmats_ac2_am3_{file_tag}.npy")
    out_file = open(f"./DEET_products/side{side}/deet_side{side}_{file_tag}.csv","w")
    out_file.write("rmat_tag,pdt_smiles\n")
    
    for i, dmat in enumerate(dmats):
        BO_increase = sum(dmat)
        if all(BO_increase <= starting_numHs):

            matrix_index = index_start + i
            m = apply_amat(acid,"ac2",amine,"am3",dmat)
                 
            s = Chem.MolToSmiles(m,isomericSmiles=True)
            out_file.write(f"{matrix_index},{s}\n")

        
    out_file.close()
    
    
    
    
    
    
def couple_amide(amine_smiles):
    
    amideN = Chem.MolFromSmarts('[N;$(NC=[O,S,N])]')
    amidecoup = ReactionFromSmarts('[C:1](=[O:2])[OH1].[N!H0:3]>>[C:1](=[O:2])[N:3]')
    deet_acid = Chem.MolFromSmiles("OC(=O)c1cc(C)ccc1")
    
    amine = Chem.MolFromSmiles(amine_smiles)
    
    
    for match in amine.GetSubstructMatches(amideN):
        amine.GetAtomWithIdx(match[0]).SetProp('_protected','1')
        
    amide_pdt = amidecoup.RunReactants((deet_acid,amine))
    
    if amide_pdt:
        for p in amide_pdt:
            pdt_smiles = Chem.MolToSmiles(p[0],isomericSmiles=True)
            return pdt_smiles
        
    else:
        return None