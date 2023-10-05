from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
opts = StereoEnumerationOptions(unique=True)

def getLogP(x): return Chem.rdMolDescriptors.CalcCrippenDescriptors(x)[0]
def getMW(x): return Chem.Descriptors.MolWt(x)
def getHBD(x): return Chem.rdMolDescriptors.CalcNumHBD(x)
def getHBA(x): return Chem.rdMolDescriptors.CalcNumHBA(x)
def getPSA(x): return Chem.rdMolDescriptors.CalcTPSA(x)
def getROTB(x): return Chem.rdMolDescriptors.CalcNumRotatableBonds(x)
def getAROM(x): return Chem.rdMolDescriptors.CalcNumAromaticRings(x)
def getFSP3(x): return Chem.rdMolDescriptors.CalcFractionCSP3(x)
def getFC(x): return Chem.rdmolops.GetFormalCharge(x)
def getQED(x): return Chem.QED.qed(x)
def getSSSR(x): return Chem.GetSSSR(x)

def getallprops(s):
    x = Chem.MolFromSmiles(s)
    return [getLogP(x),getMW(x),getHBD(x),getHBA(x),getPSA(x),getROTB(x),getFSP3(x),getSSSR(x),getQED(x)]


def getSSSR_smiles(s):
    m = Chem.MolFromSmiles(s)
    return Chem.GetSSSR(m)

def cal_pmi(s):
    try:
        m = Chem.MolFromSmiles(s)
        mol = Chem.AddHs(m)
        AllChem.EmbedMolecule(mol,maxAttempts=1000)
        x = Chem.Descriptors3D.NPR1(mol)
        y = Chem.Descriptors3D.NPR2(mol)
        return x, y
    except:
        return None,None
#         return s



def cal_pmi_2(s):
    try:
        m = Chem.MolFromSmiles(s)
        mol = Chem.AddHs(m)
        AllChem.EmbedMolecule(mol)
        x = Chem.Descriptors3D.NPR1(mol)
#         y = Chem.Descriptors3D.NPR2(mol)
        return s
    except:
#         return None,None
        return None
    
    
def try_embedding(s):
    m = Chem.MolFromSmiles(s)
    
    try:
        mol = Chem.AddHs(m)
        output = AllChem.EmbedMolecule(mol,maxAttempts=100)
        
        return output
    
    except:
        return "Embed error"
        
        
    
    
ring_triple_bond = Chem.MolFromSmarts("[*!R0]#[*!R0]")
cyclopropene = Chem.MolFromSmarts("[*!R0]1=[*!R0]~[*]1")
ring_allene = Chem.MolFromSmarts("[*!R0]=[*!R0]=[*!R0]")
ab_211_2 = Chem.MolFromSmarts("*1(~*2)=*~*~*2~*1")
ab_211_1 = Chem.MolFromSmarts("*1(~*2)=*~*2~*~*1")

ab_221_2 = Chem.MolFromSmarts("*1(~*2)~*~*~*2=*~*1")
ab_221_1 = Chem.MolFromSmarts("*12=*~*(~*~*2)~*~*1")
ab_222_2 = Chem.MolFromSmarts("*12~*~*~*(~*~*2)=*~*1")
ab_210_2 = Chem.MolFromSmarts("*12~*~*=*1~*2")
ab_210_0 = Chem.MolFromSmarts("*1(~*2)=*2~*~*1")
ab_210_1 = Chem.MolFromSmarts("*12~*~*~*1=*2")
ab_114_1 = Chem.MolFromSmarts("*1(~*2)=*~*2~*~*~*~*1")
ab_114_4 = Chem.MolFromSmarts("*1(~*2)=*~*~*~*~*2~*1")
ab_123_1 = Chem.MolFromSmarts("*12=*~*(~*~*~*2)~*~*1")
ab_123_2 = Chem.MolFromSmarts("*12=*~*~*(~*~*~*2)~*1")
ab_123_3 = Chem.MolFromSmarts("*1(~*~*2)=*~*~*~*2~*1")

myrules = [ring_triple_bond,cyclopropene,ring_allene,ab_211_2,ab_211_1,        ab_221_2,ab_221_1,ab_222_2,ab_210_2,ab_210_0,ab_210_1,ab_114_1,ab_114_4,ab_123_1,ab_123_2,ab_123_3]


def check_bredt(s):
    m = Chem.MolFromSmiles(s)

    if not any([m.HasSubstructMatch(patt) for patt in myrules]):
        return True
    else:
        return False
    
    
def get_stereoisomers(s):
    
    try:
        m = Chem.MolFromSmiles(s)
        isomers = tuple(EnumerateStereoisomers(m, options=opts))
        isomers_out = [Chem.MolToSmiles(aa, isomericSmiles=True) for aa in isomers]
    except:
        isomers_out = []
    
    return isomers_out
    
#     
# for r in tqdm(data1.itertuples()):
#     m = Chem.MolFromSmiles(r.largest_cleaned)
#     isomers = tuple(EnumerateStereoisomers(m, options=opts))


proline_6R = Chem.MolFromSmarts("[#6]1~[#6](-C(=O)-[OH1])~N~[#6]~[#6]~[#6]~1")
proline_5R = Chem.MolFromSmarts("[#6]1~[#6](-C(=O)-[OH1])~N~[#6]~[#6]~1")
proline_4R = Chem.MolFromSmarts("[#6]1~[#6](-C(=O)-[OH1])~N~[#6]~1")

prolines = [proline_6R,proline_5R,proline_4R]
def find_proline(s):
    m = Chem.MolFromSmiles(s)
    matches = [m.HasSubstructMatch(pr) for pr in prolines]
    return any(matches)