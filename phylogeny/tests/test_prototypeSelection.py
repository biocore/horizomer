from unittest import TestCase, main

from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance._base import (DistanceMatrixError,
                                        DissimilarityMatrixError,
                                        MissingIDError)
from skbio.util import get_data_path
from skbio.io._exception import UnrecognizedFormatError

from phylogeny.prototypeSelection import (prototype_selection_exhaustive,
                                          distance_sum)


class prototypeSelection(TestCase):
    def setUp(self):
        self.dm_nonNull = get_data_path('distMatrix_nonNull.txt')
        self.dm_repIDs = get_data_path('distMatrix_repIDs.txt')
        self.dm_asym = get_data_path('distMatrix_asym.txt')
        # this file must not exists!
        self.dm_noFile = get_data_path('noFile.txt')
        # any file that is present, but not a DistanceMatrix
        self.dm_wrongFormat = get_data_path('wrongFileformat.txt')

        self.dm100 = DistanceMatrix.read(get_data_path('distMatrix_100.txt'))
        self.dm20 = DistanceMatrix.read(get_data_path('distMatrix_20_f5.txt'))

    def test_distance_sum(self):
        # test that no missing IDs can be used
        self.assertRaisesRegex(
            MissingIDError,
            'The ID \'X\' is not in the dissimilarity matrix.',
            distance_sum,
            ['A', 'B', 'X'],
            self.dm20)

        # test that no ID is duplicated
        self.assertRaisesRegex(
            DissimilarityMatrixError,
            'IDs must be unique. Found the following duplicate IDs',
            distance_sum,
            ['A', 'B', 'C', 'D', 'B'],
            self.dm20)

        # test that list of IDs holds at least 1 element
        self.assertRaisesRegex(
            DissimilarityMatrixError,
            'Data must be at least 1x1 in size',
            distance_sum,
            [],
            self.dm20)

        # test for result correctness
        self.assertAlmostEqual(2454.1437464961, distance_sum(self.dm100.ids,
                                                             self.dm100))
        self.assertAlmostEqual(32.9720926186, distance_sum(
            ['550.L1S173.s.1.sequence', '550.L1S141.s.1.sequence',
             '550.L1S18.s.1.sequence', '550.L1S156.s.1.sequence',
             '550.L1S110.s.1.sequence', '550.L1S143.s.1.sequence',
             '550.L1S134.s.1.sequence', '550.L1S103.s.1.sequence',
             '550.L1S185.s.1.sequence', '550.L1S114.s.1.sequence',
             '550.L1S138.s.1.sequence', '550.L1S137.s.1.sequence'],
            self.dm100))

        self.assertAlmostEqual(81.6313, distance_sum(self.dm20.ids,
                                                     self.dm20))
        self.assertAlmostEqual(13.3887, distance_sum(
            ['A', 'C', 'F', 'G', 'M', 'N', 'P', 'T'],
            self.dm20))

    def test_exhaustive(self):
        # check if execution is rejected if number of combination is too high
        self.assertRaisesRegex(
            RuntimeError,
            'Cowardly refuse to test ',
            prototype_selection_exhaustive,
            self.dm20,
            5,
            max_combinations_to_test=1000)

        self.assertRaisesRegex(
            ValueError,
            "must be >= 2, since a single",
            prototype_selection_exhaustive,
            self.dm20,
            1)

        self.assertRaisesRegex(
            ValueError,
            "otherwise no reduction is necessary",
            prototype_selection_exhaustive,
            self.dm20,
            len(self.dm20.ids)+1)

        res = prototype_selection_exhaustive(self.dm20, 3)
        self.assertCountEqual(('A', 'P', 'Q'), res)
        self.assertAlmostEqual(1.841, distance_sum(res, self.dm20))

        res = prototype_selection_exhaustive(self.dm20, 4)
        self.assertCountEqual(('A', 'J', 'P', 'T'), res)
        self.assertAlmostEqual(3.4347, distance_sum(res, self.dm20))

        res = prototype_selection_exhaustive(self.dm20, 5)
        self.assertCountEqual(('A', 'C', 'O', 'P', 'T'), res)
        self.assertAlmostEqual(5.4494, distance_sum(res, self.dm20))

        res = prototype_selection_exhaustive(self.dm20, 18)
        self.assertCountEqual(
            ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'I', 'J', 'K', 'L', 'M', 'N',
             'O', 'P', 'Q', 'R', 'T'),
            res)
        self.assertAlmostEqual(66.94, distance_sum(res, self.dm20))

        res = prototype_selection_exhaustive(self.dm20, 19)
        self.assertCountEqual(
            ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'I', 'J', 'K', 'L', 'M', 'N',
             'O', 'P', 'Q', 'R', 'S', 'T'),
            res)
        self.assertAlmostEqual(74.1234, distance_sum(res, self.dm20))

    def test_wellformed_distance_matrix(self):
        # tests if matrices are rejected that are not symmetric
        self.assertRaisesRegex(
            DistanceMatrixError,
            "Data must be symmetric",
            DistanceMatrix.read,
            self.dm_asym)

        # tests if matrices are rejected that are hollow, i.e. have non zero
        # entries in their main diagonal
        self.assertRaisesRegex(
            DissimilarityMatrixError,
            "Data must be hollow",
            DistanceMatrix.read,
            self.dm_nonNull)

        # tests if matrices are rejected that have duplicate IDs
        self.assertRaisesRegex(
            DissimilarityMatrixError,
            "IDs must be unique. Found the following duplicate IDs",
            DistanceMatrix.read,
            self.dm_repIDs)

        # tests if absence of files is detected
        self.assertRaisesRegex(
            FileNotFoundError,
            "No such file or directory",
            DistanceMatrix.read,
            self.dm_noFile)

        # tests if a wrong format is rejected
        self.assertRaisesRegex(
            UnrecognizedFormatError,
            "Could not detect the format of",
            DistanceMatrix.read,
            self.dm_wrongFormat)


if __name__ == '__main__':
    main()
