from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Atom, BondType

from tqdm.notebook import tqdm
import pandas as pd
import numpy as np



def save_npy8(i):
    """
    saves matrix files in 8-bit numpy integer arrays
    """
    file_tag = str(i).zfill(2) 
    
    with open(f'./product_amats/output_split_{file_tag}',"r") as f:
        mats = np.array([eval(line.strip()) for line in f]).astype("int8")
        np.save(f"./product_amats/pdt_amat_{file_tag}_int8.npy",mats)
        
        

# change this for the different atoms
atoms = [6,6,6,6,6,7,8,8]

def molFromAdjMat(atoms, amat):
    """Creates a mol object from an adjacency matrix.
    Inputs:
    atoms: list of atomic numbers of atoms, by row
    amat: adjacency matrix. Has to have same length as atoms (obviously)
    Output: mol object
    """

    m = Chem.RWMol()
    # add in the separate atoms
    for a in atoms: m.AddAtom(Atom(a))
    side_len = len(amat)
    for r in range(side_len):
        for c in range(r+1,side_len):
            bond_order = amat[r][c]
            if bond_order > 0:
                if bond_order == 1: m.AddBond(r,c,BondType.SINGLE)
                if bond_order == 2: m.AddBond(r,c,BondType.DOUBLE)
                if bond_order == 3: m.AddBond(r,c,BondType.TRIPLE)
    try:
        Chem.SanitizeMol(m)
    except: 
        m = Chem.MolFromSmiles("C")
    return m



def compile_smiles_dists(file_index):
    
    """
    given a file index for product matrices, for all products, 
    compute graph edit distance from each starting material hybridization combination
    outputs a csv file.    
    """
    
    # load products
    atoms = [6,6,6,6,6,7,8,8]

    file_tag = str(file_index).zfill(2)
    amat_file = f"./product_amats/pdt_amat_{file_tag}_int8.npy"
    amats = np.load(amat_file)

    # make product smiles
    mols = [molFromAdjMat(atoms,amat) for amat in amats]
    smiles = [Chem.MolToSmiles(m) for m in mols]
    mols = []

    data_df_dict = {}
    data_df_dict["smiles"] = smiles
    
    hybrid_combos = ["ac2_am2","ac2_am3","ac3_am2","ac3_am3"]
    for hc in hybrid_combos:
        # load transformation file and take the sum of bond edits
        dmat_file = f"./rxn_mats/dmats_{hc}_{file_tag}.npy"
        dmats = np.load(dmat_file)
        bond_change_sums = [sum(sum(np.abs(dmat)))/2 for dmat in dmats]
        data_df_dict[hc] = bond_change_sums

    out_df = pd.DataFrame(data=data_df_dict)   
    out_df.to_csv(f"./data_files/smiles_with_all_dists/smiles_with_all_dists_{file_tag}.csv",index=False)
    return None
    
    
# the one just for number of lines - used for the entire drugbank
def search_in_database(smiles,database,kekulize=False):
    
    """ search the smiles in "row" through all the structures in "database"
    """
    
    substruct = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(substruct)
    
    # used to Kekulize, but currently not implemented.
    # SanitizeMol will cause aromatic bonds to be encoded as 1.5 bond order instead of alternating single and double bonds.
    if kekulize:
        Chem.Kekulize(substruct,clearAromaticFlags=True)
    
    substruct_matches = 0
    for database_mol in database:
        if database_mol.HasSubstructMatch(substruct):
            substruct_matches += 1
    return substruct_matches