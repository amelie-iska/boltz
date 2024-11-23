from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.mmcifio import MMCIFIO

def renumber_residues(structure):
    """Renumber residues sequentially starting at 1."""
    for model in structure:
        for chain in model:
            res_id = 1
            for residue in chain:
                residue.id = (residue.id[0], res_id, residue.id[2])
                res_id += 1

def pdb_to_cif(input_pdb, output_cif):
    """Convert PDB file to CIF file with renumbered residues."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb)
    
    # Renumber residues
    renumber_residues(structure)

    # Write to CIF file
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_cif)

# Example usage
input_pdb = '4eb0.pdb'
output_cif = '4eb0-fixed.cif'
pdb_to_cif(input_pdb, output_cif)
