"""
Prototype selection

Given a set of n elements (S) and their pairwise distances in terms of a
dissimilarity matrix (DM) and a given number (k << n), find the sub-set (s)
consisting of exactly k elements, such that s best represents the full set S.

Here, we define "best represents" as maximizing the sum of distances between
all points, i.e. sum {over a,b in S} DM[a,b]. This is our objective function.

This problem is known to be NP-hard [1], thus we need to resort to heuristics.
In the future, this module will implement several heuristics, whose quality can
be measured by the objective function for each problem instance, since there is
no global winner.
For completeness, the exact but exponential algorithm is implemented, too.
  "prototypeSelection_exhaustive"

[1] Gamez, J. Esteban, François Modave, and Olga Kosheleva.
    "Selecting the most representative sample is NP-hard:
     Need for expert (fuzzy) knowledge."
    Fuzzy Systems, 2008. FUZZ-IEEE 2008.
"""

from skbio.stats.distance import DistanceMatrix
import itertools
import numpy as np
import scipy as sp


def distanceSum(elements, dm: DistanceMatrix):
    '''Compute the sum of pairwise distances for the given elements according to
    the given dissimilarity matrix.

    Parameters
    ----------
    elements: list or set
        list or set of element IDs for which the sum of distances is computed
    dm: skbio.stats.distance.DistanceMatrix
        pairwise dissimilarity matrix.

    Returns
    -------
    float:
        the sum of all pairwise distances of dm for IDs in elements
    '''

    # some assertions
    if type(dm) != DistanceMatrix:
        raise TypeError('dm is not of type "DistanceMatrix"')
    assert(len(elements) > 1)
    assert(len(elements) <= len(dm.ids))
    if len(set(elements)) != len(elements):
        raise Exception('List of elements contain duplicates!')

    # actual computation
    return sum([dm[pair] for pair in set(itertools.combinations(elements, 2))])


def prototypeSelection_exhaustive(distanceMatrix: DistanceMatrix,
                                  numPrototypes: int,
                                  maxCombinationsToTest=200000):
    '''Select k prototypes for given dissimilarity matrix

    Parameters
    ----------
    distanceMatrix: skbio.stats.distance.DistanceMatrix
        pairwise distances for all elements in the full set S.
        Must be symmetric and non-hollow.
    numPrototypes: int
        Number of prototypes to select for dissimilarity matrix.
        Must be >= 2, since a single prototype is useless.
        Must be <= |distanceMatrix|, otherwise no reduction is neccessary.
    maxCombinationsToTest: int
        The maximal number of combinations to test. If exceeding, the function
        declines execution.

    Returns
    -------
    tuple
        The k-tuple holding the IDs of the selected prototypes.

    Notes
    -----
    This is the reference implementation for an exact algorithm for the
    prototype selection problem. It has an exponential runtime and will only
    operate on small instances (< maxCombinationsToTest).
    Idea: test all (n over k) combinations of selecting k elements from n with-
          out replacement. Compute the objective for each such combination and
          report the combination with maximal value.
    '''
    assert numPrototypes >= 2
    assert numPrototypes < len(distanceMatrix.ids)
    combinations = sp.special.binom(len(distanceMatrix.ids), numPrototypes)
    if combinations >= maxCombinationsToTest:
        raise Exception("Cowardly refuse to test %i combinations. Use a\
                        heuristic implementation for instances with more than\
                        %i combinations instead!"
                        % (combinations, maxCombinationsToTest))

    maxDist, maxSet = -1 * np.infty, None
    for s in set(itertools.combinations(distanceMatrix.ids,
                                        numPrototypes)):
        d = distanceSum(s, distanceMatrix)
        if d > maxDist:
            maxDist, maxSet = d, s
    return maxSet
