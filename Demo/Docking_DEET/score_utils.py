from rdkit import Chem
from rdkit.Chem.Descriptors import qed

def check_neutral(s):
    m = Chem.MolFromSmiles(s)
    
    for atm in m.GetAtoms():
        fc = atm.GetFormalCharge()
        if fc != 0:
            return False
        
    else:
        return True
    
    
def count_atoms(s):
    m = Chem.MolFromSmiles(s)
    
    if m:
        n_atoms = m.GetNumAtoms()
        
    else:
        n_atoms = None
    
    return n_atoms


def canonize_smiles(s):
    m = Chem.MolFromSmiles(s)
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


def has_finas_base(s):
    ring_base = Chem.MolFromSmarts("O=C(C=CC12)NC2CCC3C1CC[#6R2][#6R2]3")
    m = Chem.MolFromSmiles(s)
    return m.HasSubstructMatch(ring_base)


def has_deet_base(s):
    deet_base = Chem.MolFromSmarts("[#6!R0]~[#6!R0]~[#6!R0]~[#6!R0](~[#6!R0])~[#6]")
    
    m = Chem.MolFromSmiles(s)
    return m.HasSubstructMatch(deet_base)

def has_aromaticity(s):
    m = Chem.MolFromSmiles(s)
    Chem.SanitizeMol(m)
    
    m_smiles = Chem.MolToSmiles(m)
    
    return "c" in m_smiles


def scrub_smiles(s):
    t = Chem.MolFromSmiles("CC1=CC=NC=C1NC(C(C2=CC(C#N)=CC=C2)CCC(NCC)=O)=O")
    m = Chem.MolFromSmiles(s)
    s_scrub = Chem.MolToSmiles(m,isomericSmiles=False)
    m_scrub = Chem.MolFromSmiles(s_scrub)
    
    return m_scrub.HasSubstructMatch(t)


def get_qed(s):
    m = Chem.MolFromSmiles(s)
    return qed(m)

def get_qed(s):
    
    try:
        m = Chem.MolFromSmiles(s)
        m_qed = qed(m)
        
    except:
        m_qed = None
        
    return m_qed