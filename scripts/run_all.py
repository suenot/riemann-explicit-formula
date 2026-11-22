"""Reproduce all experiments offline and write results/results.json."""

from __future__ import annotations

import argparse
from itertools import accumulate
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import mpmath as mp
import numpy as np

if __package__:
    from . import riemann as rf
else:
    import riemann as rf


COUNTS = (0, 1, 2, 3, 4, 5, 10, 20, 30, 100, 1000)
CUTOFFS = (1000, 10_000, 100_000, 1_000_000)
WINDOWS = ("raw", "fejer", "gaussian", "hann")
FLOAT_PLACES = 10


def error_metrics(observed, exact):
    error = np.abs(observed - exact)
    return dict(l1_error=error.sum(), mean_absolute_error=error.mean(), max_absolute_error=error.max())


def waves(table, zeros):
    xs = np.array([10.5, 19.5, 50.5, 100.5, 1000.5])
    grid = np.arange(30, 1001) / 20
    truth = table.psi(xs, midpoint=True)
    grid_truth = table.psi(grid, midpoint=True)
    log_lcm = [math.log(lcm) for lcm in accumulate(range(1, int(xs[-1]) + 1), math.lcm, initial=1)]
    for x, target in zip(np.r_[xs, grid], np.r_[truth, grid_truth]):
        n = int(x)
        independent = (log_lcm[n - 1] + log_lcm[n]) / 2 if x == n else log_lcm[n]
        if not math.isclose(target, independent, rel_tol=0, abs_tol=1e-11):
            raise ArithmeticError(f"Prime-power and log-lcm counts disagree at x={x}")
    next_jump = np.searchsorted(table.powers, xs)
    distance = np.minimum(xs - table.powers[next_jump - 1], table.powers[next_jump] - xs)

    points = {str(x): dict(x=x, psi=target, nearest_prime_power_distance=d, estimates={})
              for x, target, d in zip(xs, truth, distance)}
    convergence = {}
    for n in COUNTS:
        estimate = rf.psi_from_zeros(xs, zeros[:n])
        cutoff = zeros[n - 1] + 1e-12 if n else None
        bound = [None] * len(xs) if cutoff is None else (
            xs * np.log(xs * cutoff)**2 / cutoff + np.log(xs) * np.minimum(1, xs / (cutoff * distance))
        )
        for point, y, error, scale in zip(points.values(), estimate, np.abs(estimate - truth), bound):
            point["estimates"][str(n)] = dict(psi=y, absolute_error=error, bound_scale=scale)
        convergence[str(n)] = dict(
            cutoff=cutoff, fixed_points=error_metrics(estimate, truth),
            grid=error_metrics(rf.psi_from_zeros(grid, zeros[:n]), grid_truth)
        )

    periods = []
    with mp.workdps(38):
        for g in zeros[:5]:
            rho = mp.mpc(0.5, float(g))
            period = 2 * mp.pi / rho.imag
            phase = mp.arg(rho) / rho.imag

            def normalized_wave(u):
                x = mp.exp(u)
                return 2 * mp.re(x**rho / rho) / mp.sqrt(x)

            a, b, c = (mp.findroot(lambda u: mp.diff(normalized_wave, u), phase + i * period / 2)
                       for i in range(3))
            periods.append(dict(gamma=g, period_log_x=float(period), measured_period=float(c - a),
                                envelope_at_20=float(2 * mp.sqrt(20) / abs(rho)),
                                measured_envelope_at_20=float(mp.sqrt(20) * (normalized_wave(a) - normalized_wave(b)) / 2)))
    example_gamma = 14.1347
    rho = complex(0.5, example_gamma)
    examples = []
    for x in (10, 20, 30):
        phase = (rho.imag * math.log(x) - math.atan2(rho.imag, rho.real)) % (2 * math.pi)
        envelope = 2 * math.sqrt(x) / abs(rho)
        pair = 2 * (x**rho / rho).real
        if not math.isclose(pair, envelope * math.cos(phase), abs_tol=1e-14):
            raise ArithmeticError("Polar and complex wave terms disagree")
        examples.append(dict(x=x, log_x=math.log(x), phase_mod_2pi=phase, envelope=envelope,
                             cosine=math.cos(phase), pair=pair, change_in_psi=-pair))
    return dict(counts=COUNTS, points=points, convergence=convergence,
                grid=dict(start=grid[0], end=grid[-1], step=0.05, count=len(grid)),
                first_wave_pairs=periods, example_gamma=example_gamma, examples=examples)


def gibbs(table, zeros):
    grid = np.arange(30, 1001) / 20
    truth = table.psi(grid, midpoint=True)
    jump_counts = (30, 100, 300, 1000)
    primes = (7, 11, 13)
    metrics = {}
    eligible = []
    for n in sorted((*COUNTS, 300)):
        metrics[str(n)] = {}
        for taper in WINDOWS:
            mae = np.mean(np.abs(rf.psi_from_zeros(grid, zeros[:n], taper) - truth))
            jumps = {str(p): rf.gibbs_peak(table, p, zeros[:n], taper) for p in primes} if n in jump_counts else {}
            worst = max((row["percent"] for row in jumps.values()), default=None)
            metrics[str(n)][taper] = dict(mean_absolute_error=mae, worst_percent=worst, jumps=jumps)
            if worst is not None and mae <= 0.25 and worst <= 2:
                eligible.append((n, mae, taper))
    best_n, _, best_window = min(eligible)
    return dict(jump_counts=jump_counts, primes=primes, metrics=metrics,
                fourier_gibbs_percent=float(100 * (mp.si(mp.pi) / mp.pi - mp.mpf("0.5"))),
                thresholds=dict(mean_absolute_error=0.25, overshoot_percent=2),
                first_qualifying_count=best_n, selected_window=best_window)


def spectrum(table, zeros):
    frequencies = np.arange(3001) / 50
    reference = zeros[:10]
    by_cutoff = {}
    for cutoff in CUTOFFS:
        raw, centered = rf.prime_spectrum(table, cutoff, frequencies)
        amplitude = np.abs(centered)
        peaks = [i for i in rf.local_maxima(amplitude) if 10 <= frequencies[i] <= 50]
        strongest = sorted(peaks, key=lambda i: amplitude[i], reverse=True)[:len(reference)]
        resolution = 2 * np.pi / math.log(cutoff / 2)
        recovered = {str(k): dict(gamma_true=g, peak=None, absolute_error=None, amplitude=None)
                     for k, g in enumerate(reference, 1)}
        unmatched = []
        for i in sorted(strongest):
            peak, height = rf.quadratic_extremum(frequencies, amplitude, i)
            nearest = np.argmin(np.abs(reference - peak))
            error = abs(reference[nearest] - peak)
            row = recovered[str(nearest + 1)]
            if error <= resolution / 2 and row["peak"] is None:
                row.update(peak=peak, absolute_error=error, amplitude=height)
            else:
                unmatched.append(dict(peak=peak, amplitude=height))
        first = max((i for i in peaks if frequencies[i] <= 18), key=lambda i: amplitude[i])
        width = rf.half_height_width(frequencies, amplitude, first)
        by_cutoff[str(cutoff)] = dict(
            recovered=recovered, unmatched=unmatched, raw_cosine_sum=raw, amplitude=amplitude,
            nominal_resolution=2 * np.pi / math.log(cutoff), effective_resolution=resolution,
            first_peak_half_height_width=width, width_to_effective_resolution=width / resolution,
            resolved_count=sum(row["peak"] is not None for row in recovered.values()),
            max_absolute_error=max(row["absolute_error"] for row in recovered.values() if row["peak"] is not None)
        )
    return dict(frequencies=frequencies, frequency_step=0.02, cutoffs=by_cutoff)


def duality(table, zeros):
    times = np.linspace(0.3, 3.5, 1601)
    zeros = zeros[:1000]
    gamma_max = zeros[-1]
    weights = rf.window("hann", zeros)
    values = rf.cosine_sum(times, zeros, weights)
    radius = np.pi / gamma_max
    minima = [rf.quadratic_extremum(times, values, i) for i in rf.local_maxima(-values)]
    normalization = gamma_max / (4 * np.pi)
    pronounced_threshold = 0.3
    powers = {}
    controls = {}
    for n in range(2, 32):
        target = math.log(n)
        if not times[0] <= target <= times[-1]:
            continue
        nearby = [(u, f) for u, f in minima if abs(u - target) <= radius]
        dip_t, dip_value = min(nearby, key=lambda pair: pair[1]) if nearby else (None, None)
        at_target, depth = rf.local_dip_depth(target, zeros, weights)
        normalized_depth = depth * math.sqrt(n) / normalization
        row = dict(log_n=target, at_log_n=at_target, dip_t=dip_t, dip_value=dip_value, depth=depth,
                   location_error=abs(dip_t - target) if dip_t is not None else None,
                   depth_times_sqrt_n=depth * math.sqrt(n), normalized_depth=normalized_depth,
                   pronounced=bool(normalized_depth > pronounced_threshold))
        index = np.searchsorted(table.powers, n)
        if table.powers[index] == n:
            log_p = table.mangoldt[index]
            row.update(log_p=log_p, landau_predicted_depth=normalization * log_p / math.sqrt(n),
                       relative_depth_error=abs(normalized_depth / log_p - 1))
            powers[str(n)] = row
        else:
            controls[str(n)] = row
    return dict(n_pairs=len(zeros), gamma_max=gamma_max, time=times, values=values, search_radius=radius,
                landau_normalization=normalization, pronounced_threshold=pronounced_threshold,
                prime_powers=powers, composites=controls,
                pronounced_prime_powers=sum(row["pronounced"] for row in powers.values()),
                pronounced_composites=sum(row["pronounced"] for row in controls.values()),
                max_prime_power_location_error=max(row["location_error"] for row in powers.values()),
                max_prime_power_relative_depth_error=max(row["relative_depth_error"] for row in powers.values()),
                max_composite_normalized_depth=max(row["normalized_depth"] for row in controls.values()))


def bohr():
    indices = (2, 3, 4, 5, 6, 7, 12)
    times = np.array([100, 400, 1600])
    terms = 100_000
    coefficients = {n: rf.bohr_components(n, times, terms) for n in indices}
    rows = {}
    for j, time in enumerate(times):
        by_n = {}
        for n, measured in coefficients.items():
            value = measured[j]
            expected = n**-2
            by_n[str(n)] = dict(real=value.real, imaginary=value.imag, limit=expected,
                                absolute_error=abs(value - expected))
        product = coefficients[2][j] * coefficients[3][j]
        rows[str(time)] = dict(
            coefficients=by_n, max_error=max(row["absolute_error"] for row in by_n.values()),
            mean_absolute_error=np.mean([row["absolute_error"] for row in by_n.values()]),
            uniform_truncation_bound=2 / (time * terms * math.log(terms / max(indices))),
            product_2_3=dict(real=product.real, imaginary=product.imag,
                             difference_from_6=abs(product - coefficients[6][j]))
        )
    return dict(terms=terms, times=rows, limit_2_times_3=2**-2 * 3**-2, limit_6=6**-2)


def counting(table, zero_strings):
    with mp.workdps(40):
        zeros = [mp.mpf(v) for v in zero_strings]
        gram = rf.ComplexGram()
        counts = (0, 10, 100)
        arithmetic = {}
        for x in (100, 1000, 10_000, 100_000, 1_000_000):
            r, li = rf.riemann_r(x), rf.li(x)
            pi = table.pi(x)
            j = table.j(x)
            inverted = table.pi_from_j(x)
            if inverted != pi:
                raise ArithmeticError(f"Möbius inversion and sieve count disagree at x={x}")
            arithmetic[str(x)] = dict(x=x, pi=pi, j=float(j), j_exact=str(j),
                                      mobius_pi=int(inverted), r=float(r), li=float(li),
                                      r_error=float(r - pi), li_error=float(li - pi))

        corrections = {}
        for x in (100, 1000):
            r = rf.riemann_r(x)
            log_x = mp.log(x)
            terms = [2 * mp.re(gram((mp.mpf("0.5") + 1j * gamma) * log_x)) for gamma in zeros[:max(counts)]]
            corrected = rf.pi_from_corrected_j(x, zeros, counts)
            pi = table.pi(x)
            corrections[str(x)] = {
                str(n): dict(r_minus_zeros=float(r - mp.fsum(terms[:n])),
                             r_minus_zeros_error=float(r - mp.fsum(terms[:n]) - pi),
                             mobius_corrected_j=float(corrected[n]),
                             mobius_corrected_j_error=float(corrected[n] - pi))
                for n in counts
            }

        jump = {}
        for tenths in range(65, 76):
            x = mp.mpf(tenths) / 10
            budget = (*counts, 1000) if x == 7 else counts
            exact = table.pi(float(x), midpoint=True)
            values = rf.pi_from_corrected_j(x, zeros, budget)
            jump[str(float(x))] = dict(x=float(x), pi_midpoint=exact,
                                       estimates={str(n): float(v) for n, v in values.items()})
        roots = [table.j_at_root(10, n) for n in (1, 2, 3)]
        worst_z = (mp.mpf("0.5") + 1j * zeros[99]) * mp.log(1000)
        if abs(gram(worst_z) - rf.ComplexGram(split=120)(worst_z)) >= mp.mpf("1e-30"):
            raise ArithmeticError("Gram series depends on the Möbius split")
        with mp.workdps(90):
            z = (mp.mpf("0.5") + 1j * mp.mpf(zero_strings[0])) * mp.log(100)
            direct = 1 + mp.fsum(z**k / (k * mp.factorial(k) * mp.zeta(k + 1)) for k in range(1, 251))
            if abs(gram(z) - direct) >= mp.mpf("1e-30"):
                raise ArithmeticError("Split and direct Gram series disagree")
        return dict(arithmetic=arithmetic, corrections=corrections, jump_at_7=jump,
                    example_at_10=dict(j=[str(v) for v in roots], pi=int(table.pi_from_j(10))),
                    gram=dict(split=gram.split, residual_terms=gram.terms, precision_dps=mp.mp.dps,
                              residual_tail_bound_at_x1000_zero100=mp.nstr(gram.tail_bound(worst_z), 6)))


def rounded(value):
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [rounded(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = round(float(value), FLOAT_PLACES)
        return 0.0 if number == 0 else number
    return value


def run(output: Path) -> None:
    start = perf_counter()
    zeros = rf.load_zeros()
    ordinates = np.array(zeros, dtype=float)
    table = rf.PrimeTable(max(CUTOFFS))
    experiments = (
        ("E1_waves", waves, (table, ordinates)),
        ("E2_gibbs", gibbs, (table, ordinates)),
        ("E3_spectrum", spectrum, (table, ordinates)),
        ("E4_duality", duality, (table, ordinates)),
        ("E5_bohr", bohr, ()),
        ("E6_counting", counting, (table, zeros)),
    )
    results = {}
    for name, experiment, arguments in experiments:
        before = perf_counter()
        results[name] = experiment(*arguments)
        print(f"{name}: {perf_counter() - before:.2f} s", flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temporary:
        pending = Path(temporary) / output.name
        pending.write_text(json.dumps(rounded(results), indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
        pending.replace(output)
    print(f"Wrote {output} in {perf_counter() - start:.2f} s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/results.json"),
                        help="write to this path relative to the current working directory")
    args = parser.parse_args()
    run(args.output)
