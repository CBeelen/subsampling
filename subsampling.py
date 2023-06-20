import csv
from argparse import ArgumentParser
from collections import defaultdict
import os
import random
from scipy.stats import fisher_exact

DEFECTS_TO_INVESTIGATE = ['intact', '5defect', 'hypermutated']


class Sequence:
    def __init__(self, clone_id, defect):
        self.clone_id = clone_id
        self.defect = defect


class SequenceList:
    def __init__(self):
        self.sequences = []
        self.unique_counter = 0
        self.clone_counter = 0
        self.distinct_counter = 0
        self.clone_sizes = defaultdict(lambda: 0)

    def add_initial_sequence(self, clone_id, defect, frequency):
        if clone_id == 'unique':
            clone_id = 'unique' + str(self.unique_counter)
            self.unique_counter += 1
            assert frequency == 1
        else:
            # this is a clone we haven't seen yet, it should have a unique id
            assert clone_id not in [sequence.clone_id for sequence in self.sequences]
            self.clone_counter += 1
        for i in range(0, frequency):
            self.sequences.append(Sequence(clone_id, defect))
        self.distinct_counter += 1

    def add_many_sequences(self, sequence_list):
        self.sequences = sequence_list
        distinct_sequences = {sequence.clone_id for sequence in sequence_list}
        self.distinct_counter = len(distinct_sequences)
        counted_clones = set()
        for clone_id in distinct_sequences:
            sequences_this_id = [sequence for sequence in sequence_list if sequence.clone_id == clone_id]
            clone_size = len(sequences_this_id)
            self.clone_sizes[clone_size] += 1
            if clone_size == 1:
                self.unique_counter += 1
            elif clone_id not in counted_clones:
                self.clone_counter += 1
                counted_clones.add(clone_id)

    def print_totals(self, defect):
        print(f"Defect {defect}: total {len(self.sequences)}, distinct {self.distinct_counter}, "
              f"unique {self.unique_counter}, distinct clones {self.clone_counter}")


def get_defect_stats(sequences, defect):
    sequences_this_defect = SequenceList()
    sequences_this_defect.add_many_sequences([sequence for sequence in sequences if sequence.defect == defect])
    other_sequences = SequenceList()
    other_sequences.add_many_sequences([sequence for sequence in sequences if sequence.defect != defect])
    sequences_this_defect.print_totals(defect)
    other_sequences.print_totals(defect=f"NOT {defect}")
    return sequences_this_defect, other_sequences


def read_data(datafile):
    reader = csv.DictReader(datafile)
    all_sequences = SequenceList()
    for row in reader:
        all_sequences.add_initial_sequence(clone_id=row['clonality'],
                                           defect=row['genomicIntegrity'],
                                           frequency=int(row['frequency']))
    return all_sequences


def do_subsampling(defect, defect_seqs, sequences, outfolder, num_replicas=100):
    """ Subsample sequences to the same depth as defect_seqs """
    sampling_depth = len(defect_seqs.sequences)
    defect_unique = defect_seqs.unique_counter
    defect_clonal = defect_seqs.clone_counter
    columns = ["iteration", "unique", "clones", "odds_ratio", "p_value"]
    average_odds = 0
    average_p_value = 0
    average_unique = 0
    average_clonal = 0
    with open(os.path.join(outfolder, f"{defect}_subsampling.csv"), 'w') as outfile:
        writer = csv.DictWriter(outfile, columns)
        writer.writeheader()
        for i in range(0, num_replicas):
            sampled_seqs = random.choices(sequences.sequences, k=sampling_depth)
            sampled_sequences = SequenceList()
            sampled_sequences.add_many_sequences(sampled_seqs)
            sampled_unique = sampled_sequences.unique_counter
            sampled_clonal = sampled_sequences.clone_counter
            table = [[defect_clonal, sampled_clonal], [defect_unique, sampled_unique]]
            stats = fisher_exact(table)
            row = {"iteration": i+1,
                   "unique": sampled_unique,
                   "clones": sampled_clonal,
                   "odds_ratio": stats.statistic,
                   "p_value": stats.pvalue}
            writer.writerow(row)
            average_odds += stats.statistic
            average_p_value += stats.pvalue
            average_clonal += sampled_clonal
            average_unique += sampled_unique
        average_odds /= num_replicas
        average_p_value /= num_replicas
        average_unique /= num_replicas
        average_clonal /= num_replicas
        row = {"iteration": 'average',
               "unique": average_unique,
               "clones": average_clonal,
               "odds_ratio": average_odds,
               "p_value": average_p_value}
        writer.writerow(row)


def main():
    parser = ArgumentParser()
    parser.add_argument('datafile', help='File containing the full set of data')
    parser.add_argument('outfolder', help='Folder to write outputs to')
    parser.add_argument('-N', help='Number of replicas to sample', default=100)
    args = parser.parse_args()

    os.mkdir(args.outfolder)

    with open(args.datafile, 'r') as datafile:
        all_sequences = read_data(datafile)

    all_sequences.print_totals(defect='ALL')

    for defect in DEFECTS_TO_INVESTIGATE:
        seq_defect, seq_other = get_defect_stats(all_sequences.sequences, defect)
        do_subsampling(defect, seq_defect, seq_other, args.outfolder, int(args.N))


if __name__ == '__main__':
    main()
