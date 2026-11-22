# Hearing the zeros

[Read the paper online](https://riemann-explicit-formula.marketmaker.cc) · [PDF](https://riemann-explicit-formula.marketmaker.cc/paper.pdf)

Six reproducible numerical experiments explore Riemann's explicit formula in both directions. Zeta zeros reconstruct the Chebyshev staircase, windowed prime-power sums recover low zero ordinates, and a cosine sum over zeros reveals prime powers. The experiments also measure finite-time Bohr coefficients of the zeta function and compare two orders of operations for reconstructing the prime count.

## Headline results

Every number in the headline tables comes from [results/results.json](results/results.json). Values are rounded for display.

### Waves and Gibbs overshoot

| At $x=19.5$ | $\psi$ | Absolute error |
| --- | ---: | ---: |
| Arithmetic target | 19.26566 | — |
| 0 pairs | 17.66344 | 1.60222 |
| 5 pairs | 18.63061 | 0.63505 |
| 30 pairs | 19.48757 | 0.22191 |
| 1,000 pairs | 19.24637 | 0.01929 |

| Window at 1,000 pairs | Grid mean absolute error | Worst overshoot at the tested primes |
| --- | ---: | ---: |
| Raw | 0.02805 | 9.15% |
| Fejér | 0.04495 | 0.00% |
| Gaussian | 0.01427 | 2.21% |
| Hann | 0.01705 | 0.64% |

### Recovering frequencies from primes

| Prime-power cutoff $X$ | First peak | Second peak | Effective spacing | Resolved targets |
| ---: | ---: | ---: | ---: | ---: |
| 1,000 | 14.13457 | 21.02139 | 1.01103 | 9 |
| 10,000 | 14.13480 | 21.02224 | 0.73771 | 10 |
| 100,000 | 14.13472 | 21.02201 | 0.58071 | 10 |
| 1,000,000 | 14.13474 | 21.02204 | 0.47882 | 10 |

### Finding prime powers and measuring Bohr coefficients

| Test of the zero sum | Result |
| --- | ---: |
| Pronounced dips at prime powers | 17 |
| Pronounced dips at other composites | 0 |
| Largest dip-location error at a prime power | 0.00008565 |
| Largest relative error in normalized depth at a prime power | 0.197% |
| Largest normalized depth at other composites | 0.00333 |

E4 subtracts $F(\log n)$ from the upper median of six samples around $\log n$,
then applies the Landau normalization described in
[paper/FORMULAS.md](paper/FORMULAS.md#e4-prime-powers-from-zeros).

| Averaging time $T$ | Largest Bohr coefficient error | Error in the finite-time product for $n=6$ |
| ---: | ---: | ---: |
| 100 | 0.024377 | 0.008630 |
| 400 | 0.004668 | 0.000796 |
| 1,600 | 0.001441 | 0.000748 |

### Counting primes

| $x$ | Sieve count | $J(x)$ | Möbius inversion | $R(x)$ | $\operatorname{Li}(x)$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 25 | 28.53333 | 25 | 25.66163 | 30.12614 |
| 1,000 | 168 | 176.69563 | 168 | 168.35945 | 177.60966 |
| 10,000 | 1,229 | 1247.09799 | 1,229 | 1226.93122 | 1246.13722 |
| 100,000 | 9,592 | 9633.76922 | 9,592 | 9587.43174 | 9629.80900 |
| 1,000,000 | 78,498 | 78597.11166 | 78,498 | 78527.39943 | 78627.54916 |

| At 100 zero pairs | $R$ minus the zero sum | Möbius inversion of corrected $J$ |
| --- | ---: | ---: |
| $x=100$ | 25.02575 | 24.91474 |
| $x=1,000$ | 167.75950 | 167.59506 |

At the jump at 7, Möbius inversion of the zero-corrected $J$ gives 3.50008 with 1,000 pairs. The midpoint is 3.5.

## Reproduce

Use Python 3.9–3.12. After dependency setup, all experiments run offline.

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt pytest -q
.venv/bin/python scripts/run_all.py && shasum results/results.json
.venv/bin/python scripts/run_all.py && shasum results/results.json
.venv/bin/python -m pytest -q
.venv/bin/python scripts/check_paper_numbers.py
```

The runner sorts JSON keys and rounds floats to ten decimal places. It prints wall times to the terminal.
The two runs should produce identical hashes in the same Python and NumPy environment.

The unchanged [zero cache](data/zeros.json) comes from the blog research file
`marketmaker-cc-landing/docs/blog-research/riemann/data/zeros.json`.
It contains the first 1,000 positive ordinates, computed with `mpmath.zetazero`
at 38 decimal places and saved to 35 significant digits.
SHA-256 `ca34ed5458d0956bdb567cba2deafbf575b868303dca7fd61c1a9aeb8159d1f8`.

The paper number map uses paths into the results, such as
`\num{/E1_waves/points/19.5/psi}`. An optional named alias takes the form
`\newcommand{\ResultPsi}{\num{/E1_waves/points/19.5/psi}}`.
After adding such references, generate the TeX definitions with
`.venv/bin/python scripts/check_paper_numbers.py --write-tex`. Build the stub
from `paper/` with `latexmk -pdf main.tex`.

## Repository layout

```text
data/zeros.json                Cached positive zero ordinates
scripts/riemann.py             Arithmetic, explicit formulas, windows and spectra
scripts/run_all.py             Six deterministic experiments
scripts/check_paper_numbers.py Result paths and LaTeX claim checks
results/results.json           Full generated results
tests/test_experiments.py      Arithmetic, spectra, determinism and checker tests
paper/main.tex                 Article stub
paper/numbers.json             Numerical claims keyed by result path
paper/FORMULAS.md              Formulas and numerical conventions
paper/refs.bib                  Bibliography
paper/REFS_CHECK.md             Bibliographic verification
CODEX_REPORT.md                 Results and acceptance log
```

## License

Code is under the [MIT License](LICENSE). Paper text and figures use [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
