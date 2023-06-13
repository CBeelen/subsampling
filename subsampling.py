import csv
from argparse import ArgumentParser

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
        for clone_id in distinct_sequences:
            sequences_this_id = [sequence for sequence in sequence_list if sequence.clone_id == clone_id]
            if len(sequences_this_id) == 1:
                self.unique_counter += 1
            else:
                # could do clone size stats here
                self.clone_counter += 1

    def print_totals(self, defect):
        print(f"Defect {defect}: total {len(self.sequences)}, distinct {self.distinct_counter}, "
              f"unique {self.unique_counter}, distinct clones {self.clone_counter}")


def print_defect_stats(sequences):
    for defect in DEFECTS_TO_INVESTIGATE:
        sequences_this_defect = SequenceList()
        sequences_this_defect.add_many_sequences([sequence for sequence in sequences if sequence.defect == defect])
        other_sequences = SequenceList()
        other_sequences.add_many_sequences([sequence for sequence in sequences if sequence.defect != defect])
        sequences_this_defect.print_totals(defect)
        other_sequences.print_totals(defect=f"NOT {defect}")


def read_data(datafile):
    reader = csv.DictReader(datafile)
    all_sequences = SequenceList()
    for row in reader:
        all_sequences.add_initial_sequence(clone_id=row['clonality'],
                                           defect=row['genomicIntegrity'],
                                           frequency=int(row['frequency']))
    return all_sequences


def main():
    parser = ArgumentParser()
    parser.add_argument('datafile', help='File containing the full set of data')
    parser.add_argument('outfile', help='Output file to be written')
    args = parser.parse_args()

    with open(args.datafile, 'r') as datafile:
        all_sequences = read_data(datafile)

    all_sequences.print_totals(defect='ALL')

    print_defect_stats(all_sequences.sequences)


if __name__ == '__main__':
    main()
