"""Check keyed LaTeX claims against a snapshot of the experiment results."""

from __future__ import annotations

import argparse
import json
import re
import sys
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
KEY = re.compile(r"/[A-Za-z0-9_./-]+")
NUM = re.compile(r"\\num(?![A-Za-z])\s*\{([^{}]*)\}")
ALIAS = re.compile(
    r"\\(?:newcommand|providecommand)\s*\{\\(Result[A-Za-z]+)\}\s*\{\s*"
    r"\\num\s*\{([^{}]*)\}\s*\}"
)
ALIAS_USE = re.compile(r"\\(Result[A-Za-z]*)\b")


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)


def result_at(results, key):
    for component in key[1:].split("/"):
        results = results[int(component)] if isinstance(results, list) else results[component]
    return results


def strip_comments(tex: str) -> str:
    # A percent sign starts a comment when preceded by an even number of backslashes.
    return re.sub(r"(?<!\\)((?:\\\\)*)%[^\n]*", r"\1", tex)


def check(results_path: Path, numbers_path: Path, tex_path: Path):
    results = read_json(results_path)
    numbers = read_json(numbers_path)
    if not isinstance(numbers, dict) or not numbers:
        raise ValueError("numbers.json must contain a nonempty mapping")

    errors = []
    for key, value in numbers.items():
        if not KEY.fullmatch(key) or type(value) not in (int, Decimal):
            errors.append(f"Invalid numeric mapping {key!r}: {value!r}")
            continue
        try:
            expected = result_at(results, key)
        except (KeyError, IndexError, ValueError, TypeError):
            errors.append(f"Result path does not exist: {key}")
            continue
        if type(expected) not in (int, Decimal) or value != expected:
            errors.append(f"{key}: paper has {value}, results have {expected}")

    direct_uses = alias_uses = 0
    if tex_path.exists():
        tex = strip_comments(tex_path.read_text(encoding="utf-8"))
        for match in NUM.finditer(tex):
            key = match[1].strip()
            if key not in numbers:
                errors.append(f"Unmapped \\num argument in {tex_path.name}: {key!r}")
            direct_uses += 1
        malformed = re.sub(r"\\(?:newcommand|providecommand)\s*\{\\num\}", "", tex)
        if re.search(r"\\num(?![A-Za-z])", NUM.sub("", malformed)):
            errors.append(f"Malformed \\num reference in {tex_path.name}")
        aliases = {}
        for alias, key in ALIAS.findall(tex):
            if alias in aliases:
                errors.append(f"Duplicate result alias: \\{alias}")
            aliases[alias] = key
        for alias in ALIAS_USE.findall(ALIAS.sub("", tex)):
            if alias not in aliases:
                errors.append(f"Undefined result alias: \\{alias}")
            alias_uses += 1
    return numbers, direct_uses, alias_uses, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-tex", action="store_true",
                        help="write the validated values to paper/numbers.tex for LaTeX")
    args = parser.parse_args()
    try:
        numbers, direct, aliased, errors = check(
            ROOT / "results" / "results.json", PAPER / "numbers.json", PAPER / "main.tex"
        )
    except (OSError, ValueError) as error:
        print(error, file=sys.stderr)
        return 1
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Verified {len(numbers)} mapped values, {direct} \\num references, and {aliased} alias uses.")
    if args.write_tex:
        definitions = [
            rf"\expandafter\def\csname result:{key}\endcsname{{{format(value, 'f') if isinstance(value, Decimal) else value}}}"
            for key, value in sorted(numbers.items())
        ]
        (PAPER / "numbers.tex").write_text("\n".join(definitions) + "\n", encoding="utf-8")
        print("Wrote paper/numbers.tex")
    return 0


if __name__ == "__main__":
    sys.exit(main())
