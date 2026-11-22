from fractions import Fraction
from math import log
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from scripts import riemann
from scripts.check_paper_numbers import check
from scripts.run_all import duality, spectrum


ROOT = Path(__file__).resolve().parents[1]


def test_psi_at_19_5_by_primes():
    expected = log(2**4 * 3**2 * 5 * 7 * 11 * 13 * 17 * 19)
    assert riemann.PrimeTable(20).psi(19.5) == pytest.approx(expected, abs=1e-12)


def test_prime_and_prime_power_counts():
    table = riemann.PrimeTable(100)
    assert table.pi(10) == 4
    assert table.j(10) == Fraction(16, 3)
    assert table.pi_from_j(100) == 25
    assert table.pi_from_j(7, midpoint=True) == Fraction(7, 2)
    assert table.pi_from_j(64, midpoint=True) == table.pi(64)


def test_spectrum_recovers_first_two_zeros():
    table = riemann.PrimeTable(1_000_000)
    zeros = np.array(riemann.load_zeros(), dtype=float)
    recovered = spectrum(table, zeros)["cutoffs"]["1000000"]["recovered"]
    assert recovered["1"]["peak"] == pytest.approx(14.1347, abs=0.05)
    assert recovered["2"]["peak"] == pytest.approx(21.0220, abs=0.05)


def test_duality_depth_matches_landau_at_prime_powers():
    table = riemann.PrimeTable(31)
    zeros = np.array(riemann.load_zeros(), dtype=float)
    results = duality(table, zeros)
    assert results["max_prime_power_relative_depth_error"] <= 0.005
    assert results["max_composite_normalized_depth"] < 0.01


def test_full_run_is_deterministic(tmp_path):
    outputs = [tmp_path / f"results-{i}.json" for i in range(2)]
    for output in outputs:
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "run_all.py"), "--output", str(output)],
            cwd=ROOT, check=True, timeout=600,
        )
    assert outputs[0].read_bytes() == outputs[1].read_bytes()


def test_paper_checker_rejects_unmapped_or_changed_claims(tmp_path):
    results = tmp_path / "results.json"
    results.write_text('{"waves": {"psi": 19.26566}}')
    numbers = tmp_path / "numbers.json"
    numbers.write_text('{"/waves/psi": 19.26566}')
    tex = tmp_path / "main.tex"
    tex.write_text(r"\newcommand{\ResultPsi}{\num{/waves/psi}} \ResultPsi % \num{/ignored}")
    assert check(results, numbers, tex)[1:] == (1, 1, [])
    tex.write_text(r"\ResultPsi \num{/wrong}")
    assert len(check(results, numbers, tex)[-1]) == 2
    numbers.write_text('{"/waves/psi": 19.26}')
    assert "paper has 19.26" in check(results, numbers, tex)[-1][0]
