from rdkit import Chem
from Bio.Seq import Seq

# Your amino acid sequence
sequence = "FSLESERP"

# Use RDKit to create a molecule from the sequence
mol = Chem.MolFromSequence(sequence)

# Generate the SMILES string
smiles = Chem.MolToSmiles(mol)

print(f"SMILES string for the sequence {sequence}:\n{smiles}")
