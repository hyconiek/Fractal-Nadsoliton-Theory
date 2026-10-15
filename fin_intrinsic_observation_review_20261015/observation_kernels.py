"""Exact parity-snapshot probabilities; no trajectory or physical-data validator."""
from fractions import Fraction
from math import comb


def integer(x):
    if type(x) is not int:
        raise TypeError("exact integer required (not bool/float)")


def hypergeometric(N, M, k, h):
    for x in (N, M, k, h):
        integer(x)
    if N < 1 or not 0 <= M <= N or not 0 <= k <= N:
        raise ValueError("invalid population or block")
    if not max(0, k - (N - M)) <= h <= min(k, M):
        return Fraction(0)
    return Fraction(comb(M, h) * comb(N - M, k - h), comb(N, k))


def two_blocks(N, M, k1, k2, h1, h2):
    for x in (N, M, k1, k2, h1, h2):
        integer(x)
    if N < 1 or not 0 <= M <= N or min(k1, k2) < 0 or k1 + k2 > N:
        raise ValueError("invalid disjoint-block contract")
    first = hypergeometric(N, M, k1, h1)
    if first == 0:
        return first
    if k1 == N:
        return first if h2 == 0 else Fraction(0)
    return first * hypergeometric(N - k1, M - h1, k2, h2)
