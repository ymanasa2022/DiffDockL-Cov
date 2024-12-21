from rdkit import Chem
import subprocess 
from rdkit.Chem import rdmolops
import datamol as dm 

# fixing format inconsistencies 
obabel_format_fix = 'obabel /home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/6ary/lig-native-pre.mol2 -O /home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/6ary/lig-native-pre-fixed.mol2'
subprocess.run(obabel_format_fix, shell=True, check=True)

# load without sanitizing to inspect structure
mol = Chem.MolFromMol2File('/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/6ary/lig-native-pre-fixed.mol2', sanitize=False, removeHs=False)

if mol:
    print("Molecule loaded successfully!")
    # mol = s.standardize(mol)
    mol = dm.fix_mol(mol)
    mol = dm.sanitize_mol(mol)
    mol = dm.standardize_mol(mol)
    # mol = Chem.RemoveHs(mol, sanitize=False)
    
    # mol = Chem.AddHs(mol, addCoords=True) # add hydrogens 
    # mol.UpdatePropertyCache()
    # writer = Chem.SDWriter('/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/6ary/6ary_ligand_h.sdf')
    # writer.SetKekulize(False)
    # writer.write(mol)
    # writer.close()
else:
    print("Failed to load the molecule.")

# id problematic atoms and bonds 
for atom in mol.GetAtoms():
    print(f"Atom {atom.GetIdx()}: {atom.GetSymbol()}, Aromatic: {atom.GetIsAromatic()}, Valence: {atom.GetDegree()}, implicit valence: {atom.GetImplicitValence()}")
    atom.UpdatePropertyCache()


for bond in mol.GetBonds():
    print(f"Bond {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}: {bond.GetBondType()}, Aromatic: {bond.GetIsAromatic()}")

# resolve aromaticity
rdmolops.Kekulize(mol, clearAromaticFlags=True)
Chem.MolToMolFile(mol, '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/6ary/6ary_ligand_aromatic.sdf')