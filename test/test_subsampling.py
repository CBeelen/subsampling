import unittest
from subsampling import SequenceList, Sequence


def assert_sequences_equal(sequences1, sequences2):
    assert len(sequences1) == len(sequences2)
    for i, seq in enumerate(sequences1):
        assert seq.clone_id == sequences2[i].clone_id
        assert seq.defect == sequences2[i].defect


class AddInitialSeqTests(unittest.TestCase):
    def setUp(self) -> None:
        self.sequences = SequenceList()

    def test_empty_sequences(self):
        seq_list = self.sequences
        assert seq_list.sequences == []
        assert seq_list.clone_counter == 0
        assert seq_list.unique_counter == 0
        assert seq_list.distinct_counter == 0

    def test_unique_sequence(self):
        seq_list = self.sequences
        seq_list.add_initial_sequence('unique', 'intact', 1)
        expected_seq_list = [Sequence('unique0', 'intact')]
        assert_sequences_equal(expected_seq_list, seq_list.sequences)
        assert seq_list.unique_counter == 1
        assert seq_list.clone_counter == 0
        assert seq_list.distinct_counter == 1

    def test_multiple_unique_sequences(self):
        seq_list = self.sequences
        seq_list.add_initial_sequence('unique', 'intact', 1)
        seq_list.add_initial_sequence('unique', '5defect', 1)
        seq_list.add_initial_sequence('unique', 'hypermutated', 1)
        expected_seq_list = [Sequence('unique0', 'intact'),
                             Sequence('unique1', '5defect'),
                             Sequence('unique2', 'hypermutated')]
        assert_sequences_equal(expected_seq_list, seq_list.sequences)
        assert seq_list.unique_counter == 3
        assert seq_list.clone_counter == 0
        assert seq_list.distinct_counter == 3

    def test_unique_frequency(self):
        seq_list = self.sequences
        with self.assertRaises(AssertionError):
            seq_list.add_initial_sequence('unique', 'intact', 3)

    def test_clonal_sequence(self):
        seq_list = self.sequences
        seq_list.add_initial_sequence('clone1', 'intact', 5)
        expected_seq_list = [Sequence('clone1', 'intact')] * 5
        assert_sequences_equal(expected_seq_list, seq_list.sequences)
        assert seq_list.unique_counter == 0
        assert seq_list.clone_counter == 1
        assert seq_list.distinct_counter == 1

    def test_same_clone(self):
        seq_list = self.sequences
        seq_list.add_initial_sequence('clone1', 'intact', 5)
        with self.assertRaises(AssertionError):
            seq_list.add_initial_sequence('clone1', 'intact', 5)

    def test_clonal_and_unique_sequences(self):
        seq_list = self.sequences
        seq_list.add_initial_sequence('unique', '5defect', 1)
        seq_list.add_initial_sequence('clone1', 'intact', 5)
        seq_list.add_initial_sequence('unique', 'hypermutated', 1)
        seq_list.add_initial_sequence('clone2', 'intact', 7)
        expected_seq_list = [Sequence('clone1', 'intact')] * 5
        expected_seq_list = [Sequence('unique0', '5defect')] + \
                            [Sequence('clone1', 'intact')] * 5 + \
                            [Sequence('unique1', 'hypermutated')] + \
                            [Sequence('clone2', 'intact')] * 7
        assert_sequences_equal(expected_seq_list, seq_list.sequences)
        assert seq_list.unique_counter == 2
        assert seq_list.clone_counter == 2
        assert seq_list.distinct_counter == 4


class AddSeqListTests(unittest.TestCase):
    def setUp(self) -> None:
        self.sequences = SequenceList()

    def test_add_all_unique(self):
        seq_list = self.sequences
        sequences = [Sequence('unique1', 'intact'),
                     Sequence('unique2', 'intact'),
                     Sequence('unique3', 'intact'),
                     Sequence('unique4', 'intact')]
        seq_list.add_many_sequences(sequences)
        assert_sequences_equal(sequences, seq_list.sequences)
        assert seq_list.unique_counter == 4
        assert seq_list.clone_counter == 0
        assert seq_list.distinct_counter == 4
        self.assertDictEqual(seq_list.clone_sizes, {1: 4})

    def test_add_one_clone(self):
        seq_list = self.sequences
        sequences = [Sequence('clone1', 'intact'),
                     Sequence('clone1', 'intact'),
                     Sequence('clone1', 'intact'),
                     Sequence('clone1', 'intact')]
        seq_list.add_many_sequences(sequences)
        assert_sequences_equal(sequences, seq_list.sequences)
        assert seq_list.unique_counter == 0
        assert seq_list.clone_counter == 1
        assert seq_list.distinct_counter == 1
        self.assertDictEqual(seq_list.clone_sizes, {4: 1})

    def test_add_unique_and_clonal(self):
        seq_list = self.sequences
        sequences = [Sequence('clone1', 'intact'),
                     Sequence('unique1', 'intact'),
                     Sequence('clone1', 'intact'),
                     Sequence('clone2', 'hypermutated'),
                     Sequence('clone1', 'intact'),
                     Sequence('clone2', 'hypermutated'),
                     Sequence('unique3', '5defect'),
                     Sequence('clone1', 'intact')]
        seq_list.add_many_sequences(sequences)
        assert_sequences_equal(sequences, seq_list.sequences)
        assert seq_list.unique_counter == 2
        assert seq_list.clone_counter == 2
        assert seq_list.distinct_counter == 4
        self.assertDictEqual(seq_list.clone_sizes, {1: 2, 2: 1, 4: 1})

    def test_unique_becomes_clone(self):
        seq_list = self.sequences
        sequences = [Sequence('unique1', 'intact'),
                     Sequence('unique1', 'intact')]
        seq_list.add_many_sequences(sequences)
        assert_sequences_equal(sequences, seq_list.sequences)
        assert seq_list.unique_counter == 0
        assert seq_list.clone_counter == 1
        assert seq_list.distinct_counter == 1
        self.assertDictEqual(seq_list.clone_sizes, {2: 1})

    def test_clone_becomes_unique(self):
        seq_list = self.sequences
        sequences = [Sequence('clone1', 'intact')]
        seq_list.add_many_sequences(sequences)
        assert_sequences_equal(sequences, seq_list.sequences)
        assert seq_list.unique_counter == 1
        assert seq_list.clone_counter == 0
        assert seq_list.distinct_counter == 1
        self.assertDictEqual(seq_list.clone_sizes, {1: 1})


if __name__ == '__main__':
    unittest.main()
