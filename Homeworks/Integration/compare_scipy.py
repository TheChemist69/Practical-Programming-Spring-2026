#!/usr/bin/env python3
"""Compare benchmark integrals with scipy.integrate.quad if SciPy is available."""

import math
import os


def maybe_import_scipy():
    """Return scipy.integrate if available, otherwise None with a friendly message."""
    try:
        from scipy import integrate  # type: ignore
    except Exception as exc:  # pragma: no cover - runtime environment dependent
        print(f"SciPy comparison skipped: {exc}")
        return None
    return integrate


def benchmark_cases():
    """Return the benchmark case list shared with the C++ executable."""
    pi = math.acos(-1.0)
    return [
        ("int_0_1_sqrt_x", lambda x: math.sqrt(x), 0.0, 1.0, 2.0 / 3.0),
        ("int_0_1_inv_sqrt_x", lambda x: 1.0 / math.sqrt(x), 0.0, 1.0, 2.0),
        ("int_0_1_sqrt_1_minus_x2", lambda x: math.sqrt(1.0 - x * x), 0.0, 1.0, pi / 4.0),
        ("int_0_1_log_x_over_sqrt_x", lambda x: math.log(x) / math.sqrt(x), 0.0, 1.0, -4.0),
    ]


def read_cpp_rows(path):
    """Read the C++ benchmark table for side-by-side call-count comparison."""
    rows = {}
    if not os.path.exists(path):
        return rows

    with open(path, "r", encoding="utf-8") as fobj:
        for line in fobj:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            name = parts[0]
            rows[name] = {
                "value": float(parts[3]),
                "abs_error": float(parts[5]),
                "evaluations": int(parts[7]),
                "status": parts[9],
            }
    return rows


def main():
    """Run SciPy quad on all benchmarks and write comparison output."""
    integrate = maybe_import_scipy()
    if integrate is None:
        return 0

    cpp_rows = read_cpp_rows("benchmark_results.data")

    with open("scipy_benchmarks.data", "w", encoding="utf-8") as fout:
        fout.write("# name value exact abs_error est_error neval\\n")

        # Print and store one comparison row per benchmark.
        print("SciPy benchmark comparison (quad):")
        print("name                           abs_error      est_error      neval   cpp_calls")
        for name, func, a, b, exact in benchmark_cases():
            value, est_error, info = integrate.quad(func, a, b, full_output=1, limit=1000)
            abs_error = abs(value - exact)
            neval = int(info.get("neval", -1))

            cpp_calls = "n/a"
            if name in cpp_rows:
                cpp_calls = str(cpp_rows[name]["evaluations"])

            print(f"{name:28s} {abs_error: .3e}   {est_error: .3e}   {neval:6d}   {cpp_calls}")
            fout.write(f"{name} {value:.17g} {exact:.17g} {abs_error:.17g} {est_error:.17g} {neval}\\n")

    print("Wrote scipy_benchmarks.data")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
