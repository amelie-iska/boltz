from Bio.PDB import MMCIFParser, Superimposer, MMCIFIO, PDBIO
from Bio.PDB.Polypeptide import is_aa, PPBuilder
from Bio.Align import substitution_matrices, PairwiseAligner
from Bio.PDB import Structure
from Bio.PDB.Chain import Chain

def load_structure(filename):
    """Load a structure from a CIF file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('', filename)
    return structure

def extract_sequence(chain):
    """Extract amino acid sequence from a chain."""
    ppb = PPBuilder()
    peptides = ppb.build_peptides(chain)
    sequence = ''.join(str(peptide.get_sequence()) for peptide in peptides)
    return sequence

def get_aligned_residues(fixed_chain, moving_chain):
    """Align sequences and get corresponding CA atoms."""
    # Extract sequences
    fixed_seq = extract_sequence(fixed_chain)
    moving_seq = extract_sequence(moving_chain)
    
    # Perform sequence alignment
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(fixed_seq, moving_seq)
    best_alignment = alignments[0]

    fixed_aln_seq = best_alignment.aligned[0]
    moving_aln_seq = best_alignment.aligned[1]

    fixed_residues = [res for res in fixed_chain.get_residues() if is_aa(res)]
    moving_residues = [res for res in moving_chain.get_residues() if is_aa(res)]

    fixed_aligned_atoms = []
    moving_aligned_atoms = []

    # Iterate through the aligned regions and collect CA atoms from matching residues
    for (fixed_start, fixed_end), (moving_start, moving_end) in zip(fixed_aln_seq, moving_aln_seq):
        for i in range(fixed_start, fixed_end):
            if i < len(fixed_residues) and i - fixed_start < len(moving_residues):
                fixed_residue = fixed_residues[i]
                moving_residue = moving_residues[i - fixed_start]
                if 'CA' in fixed_residue and 'CA' in moving_residue:
                    fixed_aligned_atoms.append(fixed_residue['CA'])
                    moving_aligned_atoms.append(moving_residue['CA'])

    return fixed_aligned_atoms, moving_aligned_atoms

def is_protein_chain(chain):
    """Check if an object is a valid protein chain with amino acid residues."""
    if not isinstance(chain, Chain):
        return False
    try:
        return any(is_aa(residue) for residue in chain.get_residues())
    except Exception as e:
        print(f"Debug: Chain ID '{chain.id}' is not a valid protein chain")
        return False

def get_unique_chain_id(existing_ids, start_id='B'):
    """Generate a unique chain ID that's not in the existing set."""
    current_id = start_id
    while current_id in existing_ids:
        # Move to next uppercase letter, cycling back to 'A' after 'Z'
        current_id = chr((ord(current_id) + 1 - 65) % 26 + 65)
    return current_id

def rename_chains(structure, start_id='B', exclude_ids=None):
    """Rename chains in the structure to ensure unique chain IDs."""
    if exclude_ids is None:
        exclude_ids = set()

    # Collect existing chain IDs to avoid conflicts
    existing_ids = set(exclude_ids)

    # First pass: identify valid protein chains
    protein_chains = []
    for model in structure:
        for chain in model:
            # Check if it's actually a Chain object and contains protein residues
            if isinstance(chain, Chain) and is_protein_chain(chain):
                if chain.id not in exclude_ids:
                    protein_chains.append(chain)

    # Second pass: rename the identified protein chains
    current_id = start_id
    for chain in protein_chains:
        # Get next available ID
        current_id = get_unique_chain_id(existing_ids, current_id)
        print(f"Renaming protein chain {chain.id} to {current_id}")
        chain.id = current_id
        existing_ids.add(current_id)
        # Move to next potential ID
        current_id = chr((ord(current_id) + 1 - 65) % 26 + 65)

def align_structures(fixed_atoms, moving_atoms):
    """Align two sets of atoms and return the Superimposer object."""
    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    return sup

def main():
    # Input file paths
    fixed_file = '/home/lily/amelie/Workspace/boltz_old/4eb0.cif'
    moving_file = 'boltz_results_4eb0-3xPET/predictions/4eb0-3xPET/4eb0-3xPET_model_5.cif'

    # Output file names
    output_fixed_cif = 'fixed_structure.cif'
    output_moving_cif = 'aligned_moving_structure.cif'
    output_combined_cif = 'combined_aligned_structures.cif'

    output_fixed_pdb = 'fixed_structure.pdb'
    output_moving_pdb = 'aligned_moving_structure.pdb'
    output_combined_pdb = 'combined_aligned_structures.pdb'

    # Chain IDs to use for alignment
    fixed_chain_id = 'A'
    moving_chain_id = 'A'

    # Load structures
    try:
        fixed_structure = load_structure(fixed_file)
        moving_structure = load_structure(moving_file)
    except Exception as e:
        print(f"Error loading structures: {str(e)}")
        return

    try:
        # Get the specified chains
        fixed_chain = fixed_structure[0][fixed_chain_id]
        moving_chain = moving_structure[0][moving_chain_id]
    except KeyError as e:
        print(f"Error: Could not find specified chain: {str(e)}")
        return

    # Get aligned CA atoms based on sequence alignment
    fixed_ca_atoms, moving_ca_atoms = get_aligned_residues(fixed_chain, moving_chain)

    # Check if CA atoms were found
    if not fixed_ca_atoms:
        print("No aligned CA atoms found. Alignment cannot proceed.")
        return

    # Align structures
    sup = align_structures(fixed_ca_atoms, moving_ca_atoms)
    sup.apply(moving_structure.get_atoms())

    # Assign unique model IDs to prevent conflicts
    fixed_model = fixed_structure[0]
    moving_model = moving_structure[0]

    fixed_model.id = 0
    moving_model.id = 1

    # Collect chain IDs from fixed structure to avoid conflicts
    fixed_chain_ids = set(chain.id for chain in fixed_model if isinstance(chain, Chain))

    # Rename chains in the moving structure
    rename_chains(moving_model, start_id='B', exclude_ids=fixed_chain_ids)

    try:
        # Save structures
        io_fixed_cif = MMCIFIO()
        io_fixed_cif.set_structure(fixed_structure)
        io_fixed_cif.save(output_fixed_cif)

        io_moving_cif = MMCIFIO()
        io_moving_cif.set_structure(moving_structure)
        io_moving_cif.save(output_moving_cif)

        io_fixed_pdb = PDBIO()
        io_fixed_pdb.set_structure(fixed_structure)
        io_fixed_pdb.save(output_fixed_pdb)

        io_moving_pdb = PDBIO()
        io_moving_pdb.set_structure(moving_structure)
        io_moving_pdb.save(output_moving_pdb)

        # Combine structures
        combined_structure = Structure.Structure('combined')
        combined_structure.add(fixed_model)
        combined_structure.add(moving_model)

        io_combined_cif = MMCIFIO()
        io_combined_cif.set_structure(combined_structure)
        io_combined_cif.save(output_combined_cif)

        io_combined_pdb = PDBIO()
        io_combined_pdb.set_structure(combined_structure)
        io_combined_pdb.save(output_combined_pdb)

    except Exception as e:
        print(f"Error saving structures: {str(e)}")
        return

    # Print alignment results
    print("Alignment completed.")
    print(f"Aligned moving structure saved to '{output_moving_cif}' and '{output_moving_pdb}'.")
    print(f"Fixed structure saved to '{output_fixed_cif}' and '{output_fixed_pdb}'.")
    print(f"Combined structure saved to '{output_combined_cif}' and '{output_combined_pdb}'.")
    print(f"Number of aligned residues: {len(fixed_ca_atoms)}")
    print(f"RMSD: {sup.rms:.4f} Å")

if __name__ == "__main__":
    main()