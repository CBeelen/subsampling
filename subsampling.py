import csv
from argparse import ArgumentParser
from collections import defaultdict
import os
import random
from scipy.stats import fisher_exact, mannwhitneyu
from numpy import std
from datetime import date

DEFECTS_TO_INVESTIGATE = ['intact', '5defect', 'hypermutated']


class Sequence:
    def __init__(self, clone_id, defect, date=None):
        self.clone_id = clone_id
        self.defect = defect
        self.date = date


class SequenceList:
    def __init__(self):
        self.sequences = []
        self.unique_counter = 0
        self.clone_counter = 0
        self.distinct_counter = 0
        self.clone_sizes = defaultdict(lambda: 0)

    def add_initial_sequence(self, clone_id, defect, frequency, date=None):
        if clone_id == 'unique':
            clone_id = 'unique' + str(self.unique_counter)
            self.unique_counter += 1
            assert frequency == 1
        else:
            # this is a clone we haven't seen yet, it should have a unique id
            assert clone_id not in [sequence.clone_id for sequence in self.sequences]
            self.clone_counter += 1
        for i in range(0, frequency):
            self.sequences.append(Sequence(clone_id, defect, date))
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

    def add_many_sequences_to_existing(self, sequence_list):
        self.sequences = self.sequences + sequence_list
        distinct_sequences = {sequence.clone_id for sequence in self.sequences}
        self.distinct_counter = len(distinct_sequences)

    def print_totals(self, identifier):
        print(f"{identifier}: total {len(self.sequences)}, distinct {self.distinct_counter}, "
              f"unique {self.unique_counter}, distinct clones {self.clone_counter}")

    def get_median_date(self):
        dates = [sequence.date for sequence in self.sequences]
        dates.sort()
        length = len(dates)
        if len(dates) % 2 == 0:
            left = dates[int(length/2)]  # python rounds down for 0.5
            right = dates[int(length/2) + 1]
            return left + (right - left)/2
        else:
            return dates[int(length/2)]

    def get_median_date_of_distinct_sequences(self):
        dates = self.get_dates()
        dates.sort()
        length = len(dates)
        if len(dates) % 2 == 0:
            left = dates[int(length/2)]  # python rounds down for 0.5
            right = dates[int(length/2) + 1]
            return date.fromordinal(int(left + (right - left)/2))
        else:
            return date.fromordinal(dates[int(length/2)])

    def get_dates(self):
        sampled_ids = set()
        dates = []
        for sequence in self.sequences:
            if sequence.clone_id in sampled_ids:
                continue
            else:
                sampled_ids.add(sequence.clone_id)
                dates.append(sequence.date.toordinal())
        return dates


def get_defect_stats(sequences, defect):
    sequences_this_defect = SequenceList()
    sequences_this_defect.add_many_sequences([sequence for sequence in sequences if sequence.defect == defect])
    other_sequences = SequenceList()
    other_sequences.add_many_sequences([sequence for sequence in sequences if sequence.defect != defect])
    sequences_this_defect.print_totals(defect)
    other_sequences.print_totals(identifier=f"NOT {defect}")
    return sequences_this_defect, other_sequences


def read_data(datafile):
    reader = csv.DictReader(datafile)
    all_sequences = SequenceList()
    for row in reader:
        all_sequences.add_initial_sequence(clone_id=row['clonality'],
                                           defect=row['genomicIntegrity'],
                                           frequency=int(row['frequency']))
    return all_sequences


def read_dates_data(datafile):
    reader = csv.DictReader(datafile)
    all_sequences = defaultdict(lambda: SequenceList())
    for row in reader:
        person = row['comparison']
        all_sequences[person].add_initial_sequence(clone_id=row['clonality'],
                                                   defect=row['query'],
                                                   frequency=int(row['frequency']),
                                                   date=date.fromisoformat(row['date']))
    return all_sequences


def add_all_data(all_stats, num_unique, num_clonal, odds_ratio, p_value):
    all_stats['unique'].append(num_unique)
    all_stats['clonal'].append(num_clonal)
    all_stats['odds_ratio'].append(odds_ratio)
    all_stats['p_value'].append(p_value)


def calculate_stats(all_stats):
    means = []
    standard_deviations = []
    for entry in all_stats.values():
        means.append(sum(entry)/len(entry))
        standard_deviations.append(std(entry))
    return means, standard_deviations


def do_subsampling(defect_seqs, sequences, outfile, num_replicas=100):
    """ Subsample sequences to the same depth as defect_seqs """
    sampling_depth = len(defect_seqs.sequences)
    defect_unique = defect_seqs.unique_counter
    defect_clonal = defect_seqs.clone_counter
    columns = ["iteration", "unique", "clones", "odds_ratio", "p_value"]
    all_stats = defaultdict(lambda: [])
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
        add_all_data(all_stats, sampled_unique, sampled_clonal, stats.statistic, stats.pvalue)
    means, standard_deviations = calculate_stats(all_stats)
    row = {"iteration": 'averages',
           "unique": means[0],
           "clones": means[1],
           "odds_ratio": means[2],
           "p_value": means[3]}
    writer.writerow(row)
    row = {"iteration": 'standard deviations',
           "unique": standard_deviations[0],
           "clones": standard_deviations[1],
           "odds_ratio": standard_deviations[2],
           "p_value": standard_deviations[3]}
    writer.writerow(row)


def do_subsampling_dates(defect_seqs, sequences, outfile, num_replicas=100):
    """ Subsample sequences to the same depth of distinct sequences as defect_seqs """
    sampling_depth = defect_seqs.distinct_counter
    print(f"Sampling to depth {sampling_depth}")
    sampled_median = defect_seqs.get_median_date_of_distinct_sequences()
    sampled_dates = defect_seqs.get_dates()
    columns = ["iteration", "median date", "p_value"]
    writer = csv.DictWriter(outfile, columns)
    writer.writeheader()
    average_median_date = 0
    average_p = 0
    for i in range(num_replicas):
        sampled_sequences = SequenceList()
        while sampled_sequences.distinct_counter < sampling_depth:
            number_missing = sampling_depth - sampled_sequences.distinct_counter
            new_seqs = random.choices(sequences.sequences, k=number_missing)
            sampled_sequences.add_many_sequences_to_existing(new_seqs)
        subsampled_dates = sampled_sequences.get_dates()
        mann_whitney = mannwhitneyu(sampled_dates, subsampled_dates)
        row = {"iteration": i+1,
               "median date": sampled_sequences.get_median_date_of_distinct_sequences(),
               "p_value": mann_whitney.pvalue}
        writer.writerow(row)
        average_median_date += sampled_sequences.get_median_date_of_distinct_sequences().toordinal()
        average_p += mann_whitney.pvalue
    average_median_date = date.fromordinal(int(average_median_date/num_replicas))
    average_p /= num_replicas
    row = {'iteration': 'Average',
           'median date': average_median_date,
           'p_value': average_p}
    writer.writerow(row)
    row = {'iteration': 'Comparison group',
           'median date': sampled_median,
           'p_value': ''}
    writer.writerow(row)


def defect_based_subsampling(file, outfolder, N):
    with open(file, 'r') as datafile:
        all_sequences = read_data(datafile)

    all_sequences.print_totals(identifier='ALL')

    for defect in DEFECTS_TO_INVESTIGATE:
        seq_defect, seq_other = get_defect_stats(all_sequences.sequences, defect)
        with open(os.path.join(outfolder, f"{defect}_subsampling.csv"), 'w') as outfile:
            do_subsampling(seq_defect, seq_other, outfile, N)


def date_based_subsampling(file, outfolder, N):
    with open(file, 'r') as datafile:
        all_sequences = read_dates_data(datafile)

    for person, sequences in all_sequences.items():
        sequences.print_totals(identifier=f'Person {person}')
        seq_og, seq_subsample = get_defect_stats(sequences.sequences, '0')
        with open(os.path.join(outfolder, f"person_{person}_subsampling.csv"), 'w') as outfile:
            do_subsampling_dates(seq_og, seq_subsample, outfile, N)



def main():
    parser = ArgumentParser()
    parser.add_argument('mode', choices=['defect', 'dates'], help='Choose what to subsample from, dates or defect')
    parser.add_argument('datafile', help='File containing the full set of data')
    parser.add_argument('outfolder', help='Folder to write outputs to')
    parser.add_argument('-N', help='Number of replicas to sample', default=100)
    args = parser.parse_args()

    os.mkdir(args.outfolder)

    if args.mode == 'defect':
        defect_based_subsampling(args.datafile, args.outfolder, int(args.N))
    elif args.mode == 'dates':
        date_based_subsampling(args.datafile, args.outfolder, int(args.N))


if __name__ == '__main__':
    main()
