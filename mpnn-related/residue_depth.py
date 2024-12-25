#!/home/amax/.conda/envs/mlfold/bin/python
import os
import argparse
import pandas as pd
from Bio.PDB import PDBParser, ResidueDepth

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

def setup_msms_path(msms_dir):
    if not os.path.exists(msms_dir):
        raise FileNotFoundError(f"MSMS directory not found: {msms_dir}")
    os.environ["PATH"] += os.pathsep + msms_dir

def get_residue_one_letter(residue):
    return THREE_TO_ONE.get(residue.resname, '')

def calculate_residue_depth_with_sequence(pdb_file, msms_tool_path):
    setup_msms_path(msms_tool_path)
    
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    model = structure[0]
    
    rd = ResidueDepth(model, 'msms')
    depth_data = []
    for residue, depth in rd:
        resi_num = residue.id[1]  # Residue number
        resi_depth = depth[0]     # Residue depth
        ca_depth = depth[1]       # Cα depth
        one_letter = get_residue_one_letter(residue)  # Residue one-letter code
        depth_data.append([resi_num, resi_depth, ca_depth, one_letter])
    
    df = pd.DataFrame(depth_data, columns=['resi_num', 'resi_depth', 'Ca_depth', 'amino_acid'])
    return df

def filter_residues_by_depth(df, depth_threshold=4):
    filtered = df[df['resi_depth'] > depth_threshold]
    residue_numbers = filtered['resi_num'].tolist()
    return " ".join(map(str, residue_numbers))

def main():
    parser = argparse.ArgumentParser(description="Calculate residue depth and extract amino acid sequence from a PDB file.")
    parser.add_argument("pdb_file", type=str, help="Path to the PDB file.")
    parser.add_argument("--msms_dir", type=str, default="/data/dj/dj_tutorial/msms", 
                        help="Path to the MSMS tool directory.")
    parser.add_argument("--output_csv", type=str, default="residue_depth_with_sequence.csv", 
                        help="Output CSV file name.")
    parser.add_argument("--depth_threshold", type=float, default=4.0, 
                        help="Depth threshold for filtering residues.")
    args = parser.parse_args()
    
    try:
        depth_df = calculate_residue_depth_with_sequence(args.pdb_file, args.msms_dir)
        depth_df.to_csv(args.output_csv, index=False)
        print(f"Residue depth and sequence data saved to {args.output_csv}")
        filtered_residues = filter_residues_by_depth(depth_df, args.depth_threshold)
        print(f"Residue depth < {args.depth_threshold}: {filtered_residues}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()

