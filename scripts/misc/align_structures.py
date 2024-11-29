#!/usr/bin/env python3

import argparse
from Bio.PDB import MMCIFParser, Superimposer, MMCIFIO, PDBIO
from Bio.PDB.Polypeptide import is_aa, PPBuilder
from Bio.Align import substitution_matrices, PairwiseAligner
from Bio.PDB import Structure, Model
from Bio.PDB.Chain import Chain
import sys
import copy

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

def get_aligned_residues(fixed_chain, moving_chain, gap_open=-10, gap_extend=-0.5, matrix="BLOSUM62"):
    """Align sequences and get corresponding CA atoms."""
    # Extract sequences
    fixed_seq = extract_sequence(fixed_chain)
    moving_seq = extract_sequence(moving_chain)
    
    # Perform sequence alignment
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    alignments = aligner.align(fixed_seq, moving_seq)
    best_alignment = alignments[0]

    fixed_aln_seq = best_alignment.aligned[0]
    moving_aln_seq = best_alignment.aligned[1]

    fixed_residues = [res for res in fixed_chain.get_residues() if is_aa(res)]
    moving_residues = [res for res in moving_chain.get_residues() if is_aa(res)]

    fixed_aligned_atoms = []
    moving_aligned_atoms = []

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
        print(f"Debug: Chain ID '{chain.id}' is not a valid protein chain", file=sys.stderr)
        return False

def get_unique_chain_id(existing_ids, start_id='B'):
    """Generate a unique chain ID that's not in the existing set."""
    current_id = start_id
    while current_id in existing_ids:
        current_id = chr((ord(current_id) + 1 - 65) % 26 + 65)
    return current_id

def rename_chains(structure, start_id='B', exclude_ids=None, verbose=False):
    """Rename chains in the structure to ensure unique chain IDs."""
    if exclude_ids is None:
        exclude_ids = set()

    existing_ids = set(exclude_ids)
    protein_chains = []
    
    for model in structure:
        for chain in model:
            if isinstance(chain, Chain) and is_protein_chain(chain):
                if chain.id not in exclude_ids:
                    protein_chains.append(chain)

    current_id = start_id
    for chain in protein_chains:
        current_id = get_unique_chain_id(existing_ids, current_id)
        if verbose:
            print(f"Renaming protein chain {chain.id} to {current_id}", file=sys.stderr)
        chain.id = current_id
        existing_ids.add(current_id)
        current_id = chr((ord(current_id) + 1 - 65) % 26 + 65)

def align_structures(fixed_atoms, moving_atoms):
    """Align two sets of atoms and return the Superimposer object."""
    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    return sup

def create_combined_structure(fixed_structure, moving_structure, verbose=False):
    """Create a new structure containing both the fixed and moving structures."""
    combined = Structure.Structure('combined')
    
    # Create a new model for the combined structure
    combined_model = Model.Model(0)
    combined.add(combined_model)
    
    # Track existing chain IDs
    existing_chain_ids = set()
    
    # Add chains from fixed structure first
    for chain in fixed_structure[0]:
        if is_protein_chain(chain):
            new_chain = copy.deepcopy(chain)
            existing_chain_ids.add(new_chain.id)
            combined_model.add(new_chain)
            if verbose:
                print(f"Added fixed chain {new_chain.id}", file=sys.stderr)
    
    # Add chains from moving structure with new IDs
    for chain in moving_structure[0]:
        if is_protein_chain(chain):
            new_chain = copy.deepcopy(chain)
            # Get a new unique chain ID
            new_id = get_unique_chain_id(existing_chain_ids, 'B')
            if verbose:
                print(f"Renaming moving chain {new_chain.id} to {new_id}", file=sys.stderr)
            new_chain.id = new_id
            existing_chain_ids.add(new_id)
            combined_model.add(new_chain)
    
    return combined

def parse_args():
    parser = argparse.ArgumentParser(description='Align protein structures based on sequence alignment.')
    
    # Required arguments
    parser.add_argument('--fixed', required=True, help='Path to the fixed (reference) structure CIF file')
    parser.add_argument('--moving', required=True, help='Path to the moving structure CIF file')
    
    # Optional arguments
    parser.add_argument('--fixed-chain', default='A', help='Chain ID in the fixed structure (default: A)')
    parser.add_argument('--moving-chain', default='A', help='Chain ID in the moving structure (default: A)')
    parser.add_argument('--prefix', default='', help='Prefix for output files (default: none)')
    parser.add_argument('--gap-open', type=float, default=-10, help='Gap opening penalty (default: -10)')
    parser.add_argument('--gap-extend', type=float, default=-0.5, help='Gap extension penalty (default: -0.5)')
    parser.add_argument('--matrix', default='BLOSUM62', help='Substitution matrix for sequence alignment (default: BLOSUM62)')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress information')
    
    # Output format options
    parser.add_argument('--output-cif', action='store_true', help='Output CIF format files')
    parser.add_argument('--output-pdb', action='store_true', help='Output PDB format files')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # If no output format specified, default to both
    if not args.output_cif and not args.output_pdb:
        args.output_cif = True
        args.output_pdb = True

    # Construct output filenames
    prefix = args.prefix + '_' if args.prefix else ''
    output_files = {
        'fixed_cif': f'{prefix}fixed_structure.cif',
        'moving_cif': f'{prefix}aligned_moving_structure.cif',
        'combined_cif': f'{prefix}combined_aligned_structures.cif',
        'fixed_pdb': f'{prefix}fixed_structure.pdb',
        'moving_pdb': f'{prefix}aligned_moving_structure.pdb',
        'combined_pdb': f'{prefix}combined_aligned_structures.pdb'
    }

    try:
        # Load structures
        if args.verbose:
            print(f"Loading structures...", file=sys.stderr)
        fixed_structure = load_structure(args.fixed)
        moving_structure = load_structure(args.moving)

        # Get the specified chains
        fixed_chain = fixed_structure[0][args.fixed_chain]
        moving_chain = moving_structure[0][args.moving_chain]

        # Get aligned CA atoms
        if args.verbose:
            print("Performing sequence alignment and collecting CA atoms...", file=sys.stderr)
        fixed_ca_atoms, moving_ca_atoms = get_aligned_residues(
            fixed_chain, moving_chain,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
            matrix=args.matrix
        )

        if not fixed_ca_atoms:
            print("Error: No aligned CA atoms found. Alignment cannot proceed.", file=sys.stderr)
            sys.exit(1)

        # Align structures
        if args.verbose:
            print("Performing structural alignment...", file=sys.stderr)
        sup = align_structures(fixed_ca_atoms, moving_ca_atoms)
        sup.apply(moving_structure.get_atoms())

        # Create combined structure with both proteins
        if args.verbose:
            print("Creating combined structure...", file=sys.stderr)
        combined_structure = create_combined_structure(fixed_structure, moving_structure, args.verbose)

        # Save structures
        if args.verbose:
            print("Saving aligned structures...", file=sys.stderr)

        if args.output_cif:
            io_fixed_cif = MMCIFIO()
            io_moving_cif = MMCIFIO()
            io_combined_cif = MMCIFIO()
            
            io_fixed_cif.set_structure(fixed_structure)
            io_moving_cif.set_structure(moving_structure)
            io_combined_cif.set_structure(combined_structure)
            
            io_fixed_cif.save(output_files['fixed_cif'])
            io_moving_cif.save(output_files['moving_cif'])
            io_combined_cif.save(output_files['combined_cif'])

        if args.output_pdb:
            io_fixed_pdb = PDBIO()
            io_moving_pdb = PDBIO()
            io_combined_pdb = PDBIO()
            
            io_fixed_pdb.set_structure(fixed_structure)
            io_moving_pdb.set_structure(moving_structure)
            io_combined_pdb.set_structure(combined_structure)
            
            io_fixed_pdb.save(output_files['fixed_pdb'])
            io_moving_pdb.save(output_files['moving_pdb'])
            io_combined_pdb.save(output_files['combined_pdb'])

        # Print results
        print(f"Alignment completed successfully:")
        print(f"Number of aligned residues: {len(fixed_ca_atoms)}")
        print(f"RMSD: {sup.rms:.4f} Å")
        
        if args.output_cif:
            print("\nCIF files saved:")
            print(f"Fixed structure: {output_files['fixed_cif']}")
            print(f"Aligned moving structure: {output_files['moving_cif']}")
            print(f"Combined structures: {output_files['combined_cif']}")
            
        if args.output_pdb:
            print("\nPDB files saved:")
            print(f"Fixed structure: {output_files['fixed_pdb']}")
            print(f"Aligned moving structure: {output_files['moving_pdb']}")
            print(f"Combined structures: {output_files['combined_pdb']}")

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
    
    
# Basic usage with default parameters
# python script.py --fixed reference.cif --moving target.cif

# Advanced usage with custom parameters
# python script.py --fixed reference.cif --moving target.cif \
#     --fixed-chain B --moving-chain C \
#     --prefix test_alignment \
#     --gap-open -12 --gap-extend -1 \
#     --matrix BLOSUM50 \
#     --verbose \
#     --output-cif