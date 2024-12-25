#!/home/amax/.conda/envs/mlfold/bin/python
import argparse
from Bio import SeqIO
from collections import defaultdict
import pandas as pd

def analyze_mutations(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    reference_seq = str(sequences[0].seq)
    mutation_stats = defaultdict(lambda: defaultdict(int))

    for record in sequences[1:]:  # 跳过参考序列
        seq = str(record.seq)
        for i, (ref_base, alt_base) in enumerate(zip(reference_seq, seq)):
            if ref_base != alt_base:
                mutation_key = f"{ref_base}{i+1}{alt_base}"
                mutation_stats[i + 1][mutation_key] += 1

    mutation_summary = []
    for position, mutations in mutation_stats.items():
        total_mutations = sum(mutations.values())
        mutation_types = ", ".join([f"{key}({count})" for key, count in mutations.items()])
        mutation_summary.append({"Position": position, "Total_Mutations": total_mutations, "Mutation_Types": mutation_types})

    mutation_df = pd.DataFrame(mutation_summary).sort_values(by="Position")
    return mutation_df, mutation_stats, reference_seq

def generate_high_freq_sequence(reference_seq, mutation_stats, top_n=10):
    all_mutations = defaultdict(int)
    for mutations in mutation_stats.values():
        for key, count in mutations.items():
            all_mutations[key] += count

    sorted_mutations = sorted(all_mutations.items(), key=lambda x: x[1], reverse=True)
    top_mutations = [mutation[0] for mutation in sorted_mutations[:top_n]]

    ref_seq_list = list(reference_seq)
    for mutation in top_mutations:
        ref_base, pos, alt_base = mutation[0], int(mutation[1:-1]) - 1, mutation[-1]
        ref_seq_list[pos] = alt_base
    mutated_sequence = "".join(ref_seq_list)

    return mutated_sequence, top_mutations

def main():
    parser = argparse.ArgumentParser(description="Analyze mutations and generate high-frequency mutated sequence.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("--output_csv", type=str, default="./mutation_summary.csv", help="Path to save the mutation summary CSV.")
    parser.add_argument("--output_fasta", type=str, default="./high_freq_mutated_sequence.fasta", help="Path to save the high-frequency mutated sequence FASTA.")
    parser.add_argument("--top_n", type=int, default=10, help="Number of top mutations to apply.")
    args = parser.parse_args()

    mutation_df, mutation_stats, reference_seq = analyze_mutations(args.fasta_file)

    mutation_df.to_csv(args.output_csv, index=False)

    mutated_sequence, top_mutations = generate_high_freq_sequence(reference_seq, mutation_stats, args.top_n)

    with open(args.output_fasta, "w") as f:
        f.write(f">High_Freq_Mutated_Sequence | Mutations: {' '.join(top_mutations)}\n")
        f.write(mutated_sequence + "\n")
    print('ana & generate done!')

if __name__ == "__main__":
    main()
