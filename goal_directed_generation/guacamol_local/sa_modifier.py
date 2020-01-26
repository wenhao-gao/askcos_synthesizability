from abc import abstractmethod
from typing import List

import sys
sys.path.append('/data/gaowh/synthGen/scscore')
import numpy as np
import guacamol_local.utils.sascorer as sascorer
from scscore.standalone_model_numpy import SCScorer
import rdkit.Chem as Chem

scscorer = SCScorer()
scscorer.restore()

class SAModifier:
    """
    Interface for synthesizability score modifiers.
    """

    @abstractmethod
    def __call__(self, smi, x):
        """
        Apply the modifier on x.

        Args:
            smi: The smiles format of molecule.
            x: The original value to be applied to.

        Returns:
            float or np.array (depending on the type of x) after application of the distance function.
        """


class ChainedModifier(SAModifier):
    """
    Calls several modifiers one after the other, for instance:
        score = modifier3(modifier2(modifier1(raw_score)))
    """

    def __init__(self, modifiers: List[SAModifier]) -> None:
        """
        Args:
            modifiers: modifiers to call in sequence.
                The modifier applied last (and delivering the final score) is the last one in the list.
        """
        self.modifiers = modifiers

    def __call__(self, smi, x):
        score = x
        for modifier in self.modifiers:
            score = modifier(smi, score)
        return score


class LinearModifier(SAModifier):
    """
    Score modifier that multiplies the score by a scalar (default: 1, i.e. do nothing).
    """

    def __init__(self, slope=1.0):
        self.slope = slope

    def __call__(self, smi, x):
        return self.slope * x


class SAScoreModifier(SAModifier):
    """
    Score modifier that multiplies the score by a scalar (default: 1, i.e. do nothing).
    """

    def __init__(self, mu: float = 2.230044, sigma: float = 0.6526308):
        self.mu = mu
        self.sigma = sigma

    def __call__(self, smi, x):
        sa_score = sascorer.calculateScore(Chem.MolFromSmiles(smi))
        mod_score = np.maximum(sa_score, self.mu)
        return np.exp(-0.5 * np.power((mod_score - self.mu) / self.sigma, 2.)) * x


class SCScoreModifier(SAModifier):
    """
    Score modifier that multiplies the score by a scalar (default: 1, i.e. do nothing).
    """

    def __init__(self, mu: float = 2.698693, sigma: float = 0.091844318):
        self.mu = mu
        self.sigma = sigma

    def __call__(self, smi, x):
        sc_score = scscorer.apply(scscorer.smi_to_fp(smi))
        mod_score = np.maximum(sc_score, self.mu)
        return np.exp(-0.5 * np.power((mod_score - self.mu) / self.sigma, 2.)) * x


class SmilesModifier(SAModifier):
    """
    Score modifier that multiplies the score by a scalar (default: 1, i.e. do nothing).
    """

    def __init__(self, a: float = 1.126, b: float = 81130999552.44):
        self.a = a
        self.b = b

    def __call__(self, smi, x):
        smiles = len(smi)
        mod_score = 1 - 1 / (1 + np.exp(- self.a * (smiles - self.b)))
        return mod_score * x


# class SquaredModifier(SAScoreModifier):
#     """
#     Score modifier that has a maximum at a given target value, and decreases
#     quadratically with increasing distance from the target value.
#     """
#
#     def __init__(self, target_value: float, coefficient=1.0) -> None:
#         self.target_value = target_value
#         self.coefficient = coefficient
#
#     def __call__(self, smi, x):
#         return 1.0 - self.coefficient * np.square(self.target_value - x)
#
#
# class AbsoluteScoreModifier(SAScoreModifier):
#     """
#     Score modifier that has a maximum at a given target value, and decreases
#     linearly with increasing distance from the target value.
#     """
#
#     def __init__(self, target_value: float) -> None:
#         self.target_value = target_value
#
#     def __call__(self, x):
#         return 1. - np.abs(self.target_value - x)
#
#
# class GaussianModifier(SAScoreModifier):
#     """
#     Score modifier that reproduces a Gaussian bell shape.
#     """
#
#     def __init__(self, mu: float, sigma: float) -> None:
#         self.mu = mu
#         self.sigma = sigma
#
#     def __call__(self, x):
#         return np.exp(-0.5 * np.power((x - self.mu) / self.sigma, 2.))
#
#
# class MinMaxGaussianModifier(SAScoreModifier):
#     """
#     Score modifier that reproduces a half Gaussian bell shape.
#     For minimize==True, the function is 1.0 for x <= mu and decreases to zero for x > mu.
#     For minimize==False, the function is 1.0 for x >= mu and decreases to zero for x < mu.
#     """
#
#     def __init__(self, mu: float, sigma: float, minimize=False) -> None:
#         self.mu = mu
#         self.sigma = sigma
#         self.minimize = minimize
#         self._full_gaussian = GaussianModifier(mu=mu, sigma=sigma)
#
#     def __call__(self, x):
#         if self.minimize:
#             mod_x = np.maximum(x, self.mu)
#         else:
#             mod_x = np.minimum(x, self.mu)
#         return self._full_gaussian(mod_x)
#
#
# MinGaussianModifier = partial(MinMaxGaussianModifier, minimize=True)
# MaxGaussianModifier = partial(MinMaxGaussianModifier, minimize=False)
#
#
# class ClippedScoreModifier(SAScoreModifier):
#     r"""
#     Clips a score between specified low and high scores, and does a linear interpolation in between.
#
#     The function looks like this:
#
#        upper_x < lower_x                 lower_x < upper_x
#     __________                                   ____________
#               \                                 /
#                \                               /
#                 \__________          _________/
#
#     This class works as follows:
#     First the input is mapped onto a linear interpolation between both specified points.
#     Then the generated values are clipped between low and high scores.
#     """
#
#     def __init__(self, upper_x: float, lower_x=0.0, high_score=1.0, low_score=0.0) -> None:
#         """
#         Args:
#             upper_x: x-value from which (or until which if smaller than lower_x) the score is maximal
#             lower_x: x-value until which (or from which if larger than upper_x) the score is minimal
#             high_score: maximal score to clip to
#             low_score: minimal score to clip to
#         """
#         assert low_score < high_score
#
#         self.upper_x = upper_x
#         self.lower_x = lower_x
#         self.high_score = high_score
#         self.low_score = low_score
#
#         self.slope = (high_score - low_score) / (upper_x - lower_x)
#         self.intercept = high_score - self.slope * upper_x
#
#     def __call__(self, x):
#         y = self.slope * x + self.intercept
#         return np.clip(y, self.low_score, self.high_score)
#
#
# class SmoothClippedScoreModifier(SAScoreModifier):
#     """
#     Smooth variant of ClippedScoreModifier.
#
#     Implemented as a logistic function that has the same steepness as ClippedScoreModifier in the
#     center of the logistic function.
#     """
#
#     def __init__(self, upper_x: float, lower_x=0.0, high_score=1.0, low_score=0.0) -> None:
#         """
#         Args:
#             upper_x: x-value from which (or until which if smaller than lower_x) the score approaches high_score
#             lower_x: x-value until which (or from which if larger than upper_x) the score approaches low_score
#             high_score: maximal score (reached at +/- infinity)
#             low_score: minimal score (reached at -/+ infinity)
#         """
#         assert low_score < high_score
#
#         self.upper_x = upper_x
#         self.lower_x = lower_x
#         self.high_score = high_score
#         self.low_score = low_score
#
#         # Slope of a standard logistic function in the middle is 0.25 -> rescale k accordingly
#         self.k = 4.0 / (upper_x - lower_x)
#         self.middle_x = (upper_x + lower_x) / 2
#         self.L = high_score - low_score
#
#     def __call__(self, x):
#         return self.low_score + self.L / (1 + np.exp(-self.k * (x - self.middle_x)))
#
#
# class ThresholdedLinearModifier(SAScoreModifier):
#     """
#     Returns a value of min(input, threshold)/threshold.
#     """
#
#     def __init__(self, threshold: float) -> None:
#         self.threshold = threshold
#
#     def __call__(self, x):
#         return np.minimum(x, self.threshold) / self.threshold
