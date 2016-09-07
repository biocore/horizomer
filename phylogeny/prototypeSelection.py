"""
Prototype selection

Given a set of n elements (S) and their pairwise distances in terms of a
skbio distance matrix (DM) and a given number (k << n), find the sub-set (s)
consisting of exactly k elements, such that s best represents the full set S.

Here, we define "best represents" as maximizing the sum of distances between
all points, i.e. sum {over a,b in S} DM[a,b]. This is our objective function.

This problem is known to be NP-hard [1], thus we need to resort to heuristics.
In the future, this module will implement several heuristics, whose quality can
be measured by the objective function for each problem instance, since there is
no global winner.
For completeness, the exact but exponential algorithm is implemented, too.
  "prototype_selection_exhaustive"

[1] Gamez, J. Esteban, François Modave, and Olga Kosheleva.
    "Selecting the most representative sample is NP-hard:
     Need for expert (fuzzy) knowledge."
    Fuzzy Systems, 2008. FUZZ-IEEE 2008.
"""

# needed for signature type annotations, but only works for python >= 3.5
# from typing import Sequence, Tuple
from itertools import combinations

import numpy as np
import scipy as sp
from skbio.stats.distance import DistanceMatrix


def distance_sum(elements, dm):
    '''Compute the sum of pairwise distances for the given elements according to
    the given distance matrix.

    Parameters
    ----------
    elements: sequence of str
        list or elements for which the sum of distances is computed
    dm: skbio.stats.distance.DistanceMatrix
        pairwise distance matrix.

    Returns
    -------
    float:
        the sum of all pairwise distances of dm for IDs in elements

    Notes
    -----
    function signature with type annotation for future use with python >= 3.5
    def distance_sum(elements: Sequence[str], dm: DistanceMatrix) -> float:
    '''

    return np.tril(dm.filter(elements).data).sum()


def prototype_selection_exhaustive(dm, num_prototypes,
                                   max_combinations_to_test=200000):
    '''Select k prototypes for given distance matrix

    Parameters
    ----------
    dm: skbio.stats.distance.DistanceMatrix
        pairwise distances for all elements in the full set S.
        Must be symmetric and non-hollow.
    num_prototypes: int
        Number of prototypes to select for distance matrix.
        Must be >= 2, since a single prototype is useless.
        Must be smaller than the number of elements in the distance matrix,
        otherwise no reduction is necessary.
    max_combinations_to_test: int
        The maximal number of combinations to test. If exceeding, the function
        declines execution.

    Returns
    -------
    sequence of str
        A sequence holding selected prototypes, i.e. a sub-set of the
        elements in the distance matrix.

    Raises
    ------
    RuntimeError
        Combinatorics explode even for small instances. To save the user from
        waiting (almost) forever, this function declines execution if the
        number of combinations to test are too high,
        i.e. > max_combinations_to_test
    ValueError
        The number of prototypes to be found should be at least 2 and at most
        one element smaller than elements in the distance matrix. Otherwise, a
        ValueError is raised.

    Notes
    -----
    This is the reference implementation for an exact algorithm for the
    prototype selection problem. It has an exponential runtime and will only
    operate on small instances (< max_combinations_to_test).
    Idea: test all (n over k) combinations of selecting k elements from n with-
          out replacement. Compute the objective for each such combination and
          report the combination with maximal value.

    function signature with type annotation for future use with python >= 3.5:
    def prototype_selection_exhaustive(dm: DistanceMatrix, num_prototypes: int,
    max_combinations_to_test: int=200000) -> Sequence[str]:
    '''
    if num_prototypes < 2:
        raise ValueError(("'num_prototypes' must be >= 2, since a single "
                          "prototype is useless."))
    if num_prototypes >= len(dm.ids):
        raise ValueError(("'num_prototypes' must be smaller than the number of"
                          " elements in the distance matrix, otherwise no "
                          "reduction is necessary."))

    num_combinations = sp.special.binom(len(dm.ids), num_prototypes)
    if num_combinations >= max_combinations_to_test:
        raise RuntimeError(("Cowardly refuse to test %i combinations. Use a "
                            "heuristic implementation for instances with more "
                            "than %i combinations instead!")
                           % (num_combinations, max_combinations_to_test))

    max_dist, max_set = -1 * np.infty, None
    for s in set(combinations(dm.ids, num_prototypes)):
        d = distance_sum(s, dm)
        if d > max_dist:
            max_dist, max_set = d, s
    return max_set
