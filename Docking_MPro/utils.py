from rdkit import Chem
import pandas as pd
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.AllChem import ReplaceSubstructs

ptable = Chem.rdchem.GetPeriodicTable()


def rough_filter(filename):
    
    file_suffix = filename.split("_")[1]
    output = open(f"./pubchem/pubchem_filtered_{file_suffix}.txt","w")
    ptable = Chem.rdchem.GetPeriodicTable()
    
    keep_elems = [1,5,6,7,8,9,15,16,17,35,53]
    #             H B C N O F P  S  Cl Br I
    rm_elems = [i for i in range(1,119) if i not in keep_elems]
    rm_elem_strings = [ptable.GetElementSymbol(i) for i in rm_elems]
    rm_elem_strings.extend(['[2H]','[3H]','.','[13C]'])
    
    data = pd.read_csv(filename,delimiter="\t",names=["index","smiles"])
    
    for s in data.smiles:
        if "N" not in s:
            continue
        
        if any([i in s for i in rm_elem_strings]): 
            continue
            
        else:
            try:
                m = Chem.MolFromSmiles(s)
                Chem.SanitizeMol(m)
                output.write(s+"\n")
            except:
                continue
                
    output.close()
    
    
def getFC(x): return Chem.rdmolops.GetFormalCharge(x)

amine_wH = Chem.MolFromSmarts("[Nv3&!H0]")

ammonium_salt = Chem.MolFromSmarts('[N+&!H0].[F-,Cl-,Br-,I-]')

neutral_N = Chem.MolFromSmarts("N")


def size_filter(file_index):
    
    file_tag = str(file_index).zfill(3)
    filename = f"./pubchem/pubchem_filtered_{file_tag}.txt"
    data = pd.read_csv(filename,names=["smiles"])
    
    amine_wH = Chem.MolFromSmarts("[Nv3&!H0]")
    ammonium_salt = Chem.MolFromSmarts('[N+&!H0].[F-,Cl-,Br-,I-]')
    neutral_N = Chem.MolFromSmarts("N")
    
    output = open(f"./pubchem/pubchem_small_{file_tag}.txt","w")
    for s in data.smiles:

        try:
            m = Chem.MolFromSmiles(s)

            msize = m.GetNumAtoms()
            if msize < 14:
                if m.HasSubstructMatch(ammonium_salt):
                    m = Chem.ReplaceSubstructs(m,ammonium_salt,neutral_N)[0]
                if m.HasSubstructMatch(amine_wH):
                    output.write(f"{s},{msize},{getFC(m)}\n")

        except:
            continue

    output.close()
    
    
    
def couple_amide(file_index):
    
    file_tag = str(file_index).zfill(3)
    filename = f"./pubchem/pubchem_small_{file_tag}.txt"
    data = pd.read_csv(filename,names = ["smiles","natoms","FC"])
    
    dacid = Chem.MolFromSmiles("O=C([C@@H](C1=CC(C#N)=CC=C1)CCC(O)=O)NC2=C(C)C=CN=C2")
    amidecoup = AllChem.ReactionFromSmarts('[C:1](=[O:2])[OH1].[N!H0:3]>>[C:1](=[O:2])[N:3]')
    amideN = Chem.MolFromSmarts('[N;$(NC=[O,S,N])]')
    
    
    
    output = open(f"./pubchem/pubchem_amide_{file_tag}.txt","w")
    output.write("smiles,size\n")
    for r in data.itertuples():
        s = r.smiles
        amine = Chem.MolFromSmiles(s)
        for match in amine.GetSubstructMatches(amideN):
            amine.GetAtomWithIdx(match[0]).SetProp('_protected','1')
        amide_pdt = amidecoup.RunReactants((dacid,amine))
        if amide_pdt:
            for p in amide_pdt:
                pdt_smiles = Chem.MolToSmiles(p[0],isomericSmiles=True)
                output.write(f"{pdt_smiles},{r.natoms}\n")
    output.close()

    
