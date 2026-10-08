import os
import random
import numpy as np
import scipy.stats as stats
from collections import Counter
from multiprocessing import Pool
from functools import partial
import argparse

def read_fasta(file):
    seqs = {}
    with open(file, 'r') as f:
        seq = ''
        seq_id = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    seqs[seq_id] = seq
                seq_id = line[1:]
                seq = ''
            else:
                seq += line
        if seq_id:
            seqs[seq_id] = seq
    return seqs

def get_subseq(sequence, start, length):
    if len(sequence) >= (start + length):
        return sequence[start:start+length]
    return ''
def shuffle_sequence(sequence):
    sequence_list = list(sequence)
    random.shuffle(sequence_list)
    return ''.join(sequence_list)
def reverse_sequence(sequence):
    return sequence[::-1]
def quasi_random_shuffle(sequence):
    counts = Counter(sequence)
    nucleotides = []
    for nucleotide, count in counts.items():
        nucleotides.extend([nucleotide] * count)
    random.shuffle(nucleotides)
    return ''.join(nucleotides)
def compute_mismatches(seq1, seq2):
    mismatches = sum([1 for a, b in zip(seq1, seq2) if a != b])
    return mismatches
def compute_true_positives(virus_seqs, hsa_seqs, subseq_start, subseq_length):
    complement_mapping = str.maketrans({'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'})
    true_positives = []
    for virus_id, virus_seq in virus_seqs.items():
        for hsa_id, hsa_seq in hsa_seqs.items():
            virus_subseq = get_subseq(virus_seq, subseq_start, subseq_length)
            hsa_subseq = get_subseq(hsa_seq, subseq_start, subseq_length)
            complement_seq = hsa_seq.translate(complement_mapping)
            hsa_subseq_complement = complement_seq[-2:-8:-1]
            match = 0
            if virus_subseq == hsa_subseq or virus_subseq == hsa_subseq_complement:
                match = 1
            true_positives.append(match)
    return true_positives

def process_simulation(args_tuple):
    virus_seqs, hsa_seqs, subseq_start, subseq_length = args_tuple
    shuffle_count = 0
    reverse_count = 0
    quasi_random_count = 0

    for virus_id, virus_seq in virus_seqs.items():
        for hsa_id, hsa_seq in hsa_seqs.items():
            virus_subseq_original = get_subseq(virus_seq, subseq_start, subseq_length)
            hsa_subseq_original = get_subseq(hsa_seq, subseq_start, subseq_length)
            if not virus_subseq_original or not hsa_subseq_original:
                continue

            # Shuffle
            virus_subseq_shuffle = virus_seq
            hsa_subseq_shuffle = shuffle_sequence(hsa_seq)
            mismatches_shuffle = compute_mismatches(virus_subseq_shuffle[1:7], hsa_subseq_shuffle[1:7])
            if mismatches_shuffle == 0:
                shuffle_count += 1
            # Reverse
            virus_subseq_reverse = virus_seq
            hsa_subseq_reverse = reverse_sequence(hsa_seq)
            mismatches_reverse = compute_mismatches(virus_subseq_reverse[1:7], hsa_subseq_reverse[1:7])
            if mismatches_reverse == 0:
                reverse_count += 1
            # Quasi-Random Shuffle
            virus_subseq_quasi = quasi_random_shuffle(virus_seq)
            hsa_subseq_quasi = quasi_random_shuffle(hsa_seq)
            mismatches_quasi = compute_mismatches(virus_subseq_quasi[1:7], hsa_subseq_quasi[1:7])
            if mismatches_quasi == 0:
                quasi_random_count += 1
    return shuffle_count, reverse_count, quasi_random_count

def permutation_test(observed_count, random_counts, n):
    count = sum(1 for rc in random_counts if rc >= observed_count)
    p_value = count / n
    return p_value

def main(args):
    virus_seqs = read_fasta(args.virus)
    hsa_seqs = read_fasta(args.hsa)
    true_positives = compute_true_positives(virus_seqs, hsa_seqs, args.subseq_start, args.subseq_length)
    true_positives_count = np.sum(true_positives)
    # 准备模拟任务的参数
    simulation_args = [(virus_seqs, hsa_seqs, args.subseq_start, args.subseq_length) for _ in range(args.num_simulations)]
    # 使用多进程执行模拟
    with Pool(processes=args.num_processes) as pool:
        simulation_results = pool.map(process_simulation, simulation_args)
    shuffle_results_all = [result[0] for result in simulation_results]
    reverse_results_all = [result[1] for result in simulation_results]
    quasi_random_results_all = [result[2] for result in simulation_results]
    true_positives_repeated = [true_positives_count] * args.num_simulations
    #print(true_positives_repeated)
    #print(shuffle_results_all)
    #print(reverse_results_all)
    #print(quasi_random_results_all)

    # 进行pvalue计算
    # 进行 p-value 计算
    shuffle_p_value = permutation_test(true_positives_count, shuffle_results_all, args.num_simulations)
    reverse_p_value = permutation_test(true_positives_count, reverse_results_all, args.num_simulations)
    quasi_random_p_value = permutation_test(true_positives_count, quasi_random_results_all, args.num_simulations)

    #shuffle_p_value = permutation_test(true_positives_repeated, shuffle_results_all, args.num_simulations)
    #reverse_p_value = permutation_test(true_positives_repeated, reverse_results_all, args.num_simulations)
    #quasi_random_p_value = permutation_test(true_positives_repeated, quasi_random_results_all, args.num_simulations)
    print(f"Shuffle Strategy p-value: {shuffle_p_value}")
    print(f"Reverse Strategy p-value: {reverse_p_value}")
    print(f"Quasi-Random Strategy p-value: {quasi_random_p_value}")
    with open(args.output, 'w') as f:
        f.write('Average_True_Positives\tAverage_Shuffle\tAverage_Reverse\tAverage_Quasi_Random\tP-value_Shuffle\tP-value_Reverse\tP-value_Quasi\n')
        f.write(f'{np.mean(true_positives_repeated)}\t{np.mean(shuffle_results_all)}\t{np.mean(reverse_results_all)}\t{np.mean(quasi_random_results_all)}\t{shuffle_p_value}\t{reverse_p_value}\t{quasi_random_p_value}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sequence matching analysis with parallel processing.")
    parser.add_argument('--virus', required=True, help="Path to the virus sequences FASTA file")
    parser.add_argument('--hsa', required=True, help="Path to the host sequences FASTA file")
    parser.add_argument('--subseq_start', type=int, required=True, help="Start position of the subsequence")
    parser.add_argument('--subseq_length', type=int, required=True, help="Length of the subsequence")
    parser.add_argument('--num_simulations', type=int, required=True, help="Number of simulations")
    parser.add_argument('--num_processes', type=int, required=True, help="Number of processes for parallel computing")
    parser.add_argument("--output", required=True, help="Output file for results")

    args = parser.parse_args()
    random.seed(27)
    main(args)

