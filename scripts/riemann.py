"""Arithmetic counts, zero sums, and spectral estimators for the experiments."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_zeros() -> list[str]:
    """Return the cached decimal ordinates, preserving their source precision."""
    return json.loads((ROOT / "data" / "zeros.json").read_text(encoding="utf-8"))


def prime_sieve(limit: int) -> np.ndarray:
    flags = np.ones(limit + 1, dtype=bool)
    flags[:2] = False
    for p in range(2, math.isqrt(limit) + 1):
        if flags[p]:
            flags[p * p::p] = False
    return np.flatnonzero(flags)


class PrimeTable:
    """Arithmetic counts up to limit, right-continuous unless midpoint is requested."""

    def __init__(self, limit: int):
        self.limit = limit
        self.primes = prime_sieve(limit)
        powers = np.zeros(limit + 1)
        for p in map(int, self.primes):
            power = p
            weight = math.log(p)
            while power <= limit:
                powers[power] = weight
                power *= p
        self.powers = np.flatnonzero(powers)
        self.mangoldt = powers[self.powers]
        self.psi_levels = np.r_[0.0, np.cumsum(self.mangoldt)]

    def pi(self, x, *, midpoint=False):
        count = int(np.searchsorted(self.primes, x, side="right"))
        if midpoint and count and self.primes[count - 1] == x:
            return count - 0.5
        return count

    def psi(self, x, *, midpoint=False):
        count = np.searchsorted(self.powers, x, side="right")
        value = self.psi_levels[count]
        if midpoint:
            left = np.searchsorted(self.powers, x, side="left")
            value = (value + self.psi_levels[left]) / 2
        return value

    def j(self, x, *, midpoint=False) -> Fraction:
        # This also supports evaluating at exact roots via j_at_root.
        return self.j_at_root(x, 1, midpoint=midpoint)

    def j_at_root(self, x, root_index: int, *, midpoint=False) -> Fraction:
        """Evaluate J(x**(1/root_index)) with exact jump comparisons."""
        total = Fraction()
        for k in range(1, int(x).bit_length()):
            exponent = k * root_index
            if 2 ** exponent > x:
                break
            root = integer_root(int(x), exponent)
            count = int(np.searchsorted(self.primes, root, side="right"))
            equality = midpoint and count and int(self.primes[count - 1]) ** exponent == x
            total += Fraction(2 * count - int(equality), 2 * k)
        return total

    def pi_from_j(self, x, *, midpoint=False) -> Fraction:
        mu = mobius_values(int(x).bit_length() - 1)
        return sum((Fraction(mu[m], m) * self.j_at_root(x, m, midpoint=midpoint)
                    for m in range(1, len(mu)) if mu[m]), Fraction())


def integer_root(x: int, exponent: int) -> int:
    """Largest integer r satisfying r**exponent <= x."""
    if x < 2 or exponent == 1:
        return x
    lo, hi = 1, 1 << ((x.bit_length() + exponent - 1) // exponent)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** exponent <= x:
            lo = mid
        else:
            hi = mid
    return lo


def mobius_values(limit: int) -> list[int]:
    mu = [0] + [1] * limit
    for p in map(int, prime_sieve(limit)):
        for multiple in range(p, limit + 1, p):
            mu[multiple] *= -1
        for multiple in range(p * p, limit + 1, p * p):
            mu[multiple] = 0
    return mu


def window(name: str, ordinates: np.ndarray) -> np.ndarray:
    if name == "raw" or not ordinates.size:
        return np.ones_like(ordinates)
    r = ordinates / ordinates[-1]
    if name == "fejer":
        return 1 - r
    if name == "gaussian":
        return np.exp(-r * r)
    if name == "hann":
        return (1 + np.cos(np.pi * r)) / 2
    raise ValueError(f"Unknown window {name}")


def psi_baseline(x):
    return x - np.log(2 * np.pi) - np.log1p(-1 / x**2) / 2


def psi_waves(x, ordinates: np.ndarray) -> np.ndarray:
    """Return one column per conjugate pair, 2 Re(x**rho / rho)."""
    x = np.asarray(x, dtype=float)
    phase = np.arctan2(ordinates, 0.5)
    amplitude = 2 * np.sqrt(x)[..., None] / np.hypot(ordinates, 0.5)
    return amplitude * np.cos(np.log(x)[..., None] * ordinates - phase)


def psi_from_zeros(x, ordinates: np.ndarray, taper="raw"):
    return psi_baseline(np.asarray(x, dtype=float)) - np.sum(
        psi_waves(x, ordinates) * window(taper, ordinates), axis=-1
    )


def gibbs_peak(table: PrimeTable, p: int, ordinates: np.ndarray, taper: str) -> dict:
    """Refine right-side extrema of the reconstructed psi in log coordinates."""
    weights = window(taper, ordinates)
    next_power = table.powers[np.searchsorted(table.powers, p, side="right")]
    lo = math.log(p)
    hi = min((lo + math.log(next_power)) / 2, lo + 2 * math.pi / ordinates[-1])

    def slope(u):
        x = np.exp(u)
        waves = np.sum(np.cos(np.asarray(u)[..., None] * ordinates) * weights, axis=-1)
        return x - 1 / (x * x - 1) - 2 * np.sqrt(x) * waves

    def candidates(intervals):
        grid = np.linspace(lo, hi, intervals + 1)
        gradients = slope(grid)
        extrema = [lo, hi]
        for idx in np.flatnonzero((gradients[:-1] > 0) & (gradients[1:] <= 0)):
            left, right = grid[idx:idx + 2]
            for _ in range(45):
                mid = (left + right) / 2
                if slope(mid) > 0:
                    left = mid
                else:
                    right = mid
            extrema.append((left + right) / 2)
        return extrema

    coarse, fine = candidates(48), candidates(96)
    if len(coarse) != len(fine) or not np.allclose(coarse, fine, atol=1e-13, rtol=0):
        raise RuntimeError(f"Gibbs extrema depend on search density at p={p}")
    points = np.exp(fine)
    points[0] = p
    values = psi_from_zeros(points, ordinates, taper)
    peak = np.argmax(values)
    plateau = table.psi(p)
    excess = values[peak] - plateau
    return dict(p=p, x_peak=points[peak], x_end=points[1], max_value=values[peak],
                plateau=plateau, signed_excess=excess,
                percent=100 * max(0.0, excess) / math.log(p))


def cosine_sum(times: np.ndarray, frequencies: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Evaluate a bounded-memory cosine projection with a fixed reduction order."""
    values = np.empty_like(times, dtype=float)
    for start in range(0, len(times), 32):
        stop = start + 32
        values[start:stop] = np.sum(np.cos(times[start:stop, None] * frequencies) * weights, axis=1)
    return values


def local_dip_depth(target: float, ordinates: np.ndarray, weights: np.ndarray):
    """Return F(target) and its signed depth below the upper median of the ring."""
    offsets = np.array([0.0, -0.06, -0.05, -0.04, 0.04, 0.05, 0.06])
    center, *ring = cosine_sum(target + offsets, ordinates, weights)
    return center, sorted(ring)[3] - center


def prime_spectrum(table: PrimeTable, cutoff: int, frequencies: np.ndarray):
    """Return the Hann-windowed prime cosine sum and its smooth-background residual."""
    stop = np.searchsorted(table.powers, cutoff, side="right")
    powers = table.powers[:stop]
    log_n = np.log(powers)
    start = math.log(2)
    span = math.log(cutoff) - start
    taper = (1 - np.cos(2 * np.pi * (log_n - start) / span)) / 2
    raw = cosine_sum(frequencies, log_n, taper * table.mangoldt[:stop] / np.sqrt(powers))

    # Integrate the smooth density exp(u/2) against the same Hann and cosine.
    def exp_integral(angular_frequency):
        z = 0.5 + 1j * angular_frequency
        return np.expm1(z * span) / z

    omega = 2 * np.pi / span
    integral = (0.5 * exp_integral(frequencies)
                - 0.25 * (exp_integral(frequencies + omega) + exp_integral(frequencies - omega)))
    background = np.real(np.exp((0.5 + 1j * frequencies) * start) * integral)
    return raw, raw - background


def local_maxima(values: np.ndarray) -> np.ndarray:
    return 1 + np.flatnonzero((values[1:-1] > values[:-2]) & (values[1:-1] >= values[2:]))


def quadratic_extremum(grid: np.ndarray, values: np.ndarray, index: int):
    """Interpolate the extremum of three equally spaced samples."""
    left, center, right = values[index - 1:index + 2]
    curvature = left - 2 * center + right
    offset = 0.5 * (left - right) / curvature
    return grid[index] + offset * (grid[1] - grid[0]), center - 0.25 * (left - right) * offset


def half_height_width(grid: np.ndarray, values: np.ndarray, index: int):
    half = values[index] / 2
    left = np.flatnonzero(values[:index] <= half)
    right = np.flatnonzero(values[index + 1:] <= half)
    if not len(left) or not len(right):
        return None

    def crossing(i, j):
        return grid[i] + (grid[j] - grid[i]) * (half - values[i]) / (values[j] - values[i])

    a = left[-1]
    b = index + 1 + right[0]
    return crossing(b - 1, b) - crossing(a, a + 1)


def bohr_components(n: int, times: np.ndarray, terms: int) -> np.ndarray:
    """Integrate the absolutely convergent zeta Dirichlet series termwise."""
    m = np.arange(1, terms + 1, dtype=float)
    theta = times[:, None] * np.log(n / m) / 2
    kernel = np.exp(1j * theta) * np.sinc(theta / np.pi)
    return np.sum(kernel / m**2, axis=1)


def li(x):
    return mp.li(x)


def riemann_r(x):
    return mp.riemannr(x)


class ComplexGram:
    """Evaluate R(e^z) with z unwrapped using a stable Möbius split."""

    def __init__(self, split=100, terms=100):
        self.split = split
        self.terms = terms
        mu = mobius_values(split)
        self.nonzero = [(n, mu[n]) for n in range(1, split + 1) if mu[n]]
        with mp.workdps(240):
            self.residual = [
                1 / mp.zeta(k + 1) - mp.fsum(mp.mpf(sign) / mp.mpf(n)**(k + 1)
                                             for n, sign in self.nonzero)
                for k in range(1, terms + 1)
            ]

    def __call__(self, z):
        total = mp.mpf(1)
        if z == 0:
            return total
        for n, sign in self.nonzero:
            w = z / n
            total += (mp.mpf(sign) / n) * (mp.ei(w) - mp.euler - mp.log(w))
        term = mp.mpf(1)
        for k, coefficient in enumerate(self.residual, 1):
            term *= z / k
            total += coefficient * term / k
        return total

    def tail_bound(self, z):
        k = self.terms + 1
        ratio = abs(z) / self.split
        return ratio**k / (k**2 * mp.factorial(k) * (1 - ratio / (k + 1)))


def j_from_zeros(y, zeros, counts):
    """Evaluate the same height-ordered J truncation at several pair counts."""
    log_y = mp.log(y)
    tail = mp.quad(lambda t: 1 / (t * (t * t - 1) * mp.log(t)), [y, mp.inf])
    baseline = li(y) - mp.log(2) + tail
    terms = [2 * mp.re(mp.ei((mp.mpf("0.5") + 1j * gamma) * log_y))
             for gamma in zeros[:max(counts)]]
    return {n: baseline - mp.fsum(terms[:n]) for n in counts}


def pi_from_corrected_j(x, zeros, counts):
    """Apply the finite Möbius inversion after reconstructing J."""
    log_x = mp.log(x)
    cutoff = 0
    while 2**(cutoff + 1) <= x:
        cutoff += 1
    mu = mobius_values(cutoff)
    estimates = {n: mp.mpf(0) for n in counts}
    for m in range(1, cutoff + 1):
        if mu[m]:
            j = j_from_zeros(mp.exp(log_x / m), zeros, counts)
            for n in counts:
                estimates[n] += (mp.mpf(mu[m]) / m) * j[n]
    return estimates
