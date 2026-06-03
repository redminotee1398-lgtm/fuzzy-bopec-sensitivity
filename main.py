"""
main.py
=======
Main entry point for BOPEC + Fuzzy sensitivity analysis.

Reproduces all paper tables and generates Figure 1.

Usage:
    python main.py              # full run (tables + figure)
    python main.py --no-plot    # tables only
    python main.py --outdir /path/to/dir
"""

import sys
import io
import os
import time
import argparse
import numpy as np

# Force UTF-8 on Windows
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from bopec_sensitivity import (
    linear_sensitivity, Z_linear,
    nonlinear_sensitivity, Z_nonlinear,
    cournot_sensitivity, theorem2
)
from fuzzy_bopec import compute_alpha_cuts, print_summary
from plots import make_figure1

SEP  = "=" * 65
SEP2 = "-" * 65


def parse_args():
    p = argparse.ArgumentParser(description="BOPEC Sensitivity Analysis")
    p.add_argument("--no-plot", action="store_true",
                   help="Skip figure generation")
    p.add_argument("--outdir", default="results",
                   help="Output directory (default: results/)")
    return p.parse_args()


def run_tables():
    """Reproduce Tables 1-4 from the paper."""

    # ── TABLE 1 ──────────────────────────────────────────────
    print(f"\n{SEP}")
    print("TABLE 1  Linear System  B=1.0, n=3")
    print(SEP)
    t0 = time.perf_counter()
    Sd, Si, St, _, _ = linear_sensitivity(1.0, 3)
    t_ana = time.perf_counter() - t0
    h = 1e-7
    t0 = time.perf_counter()
    St_fd = (Z_linear(1 + h, 3) - Z_linear(1 - h, 3)) / (2 * h)
    t_fd = time.perf_counter() - t0

    tgts = {"Direct": 9.200000, "Indirect": -0.502755, "Total": 8.697245}
    vals = {"Direct": Sd, "Indirect": Si, "Total": St}
    print(f"  {'Component':<14}{'Paper':>12}{'Computed':>14}{'Match':>8}")
    print(f"  {SEP2}")
    for k in tgts:
        ok = "✓" if abs(vals[k] - tgts[k]) < 5e-4 else "✗"
        print(f"  {k:<14}{tgts[k]:>+12.6f}{vals[k]:>+14.6f}{ok:>8}")
    print(f"  FD check{'':6}{St_fd:>+14.6f}  err={abs(St - St_fd):.2e}")
    print(f"  CPU: analytical={t_ana * 1e3:.2f} ms  |  FD={t_fd * 1e3:.2f} ms")

    # ── TABLE 3 ──────────────────────────────────────────────
    print(f"\n{SEP}")
    print("TABLE 3  Linear Scalability  B=1.0")
    print(SEP)
    paper3 = {3: (-0.503, 8.697), 5: (-0.570, 8.630),
               10: (-0.586, 8.614), 20: (-0.500, 8.700), 50: (-0.305, 8.895)}
    print(f"  {'n':>4}{'Direct':>10}{'Indirect':>12}{'Total':>10}{'Match':>8}")
    print(f"  {SEP2}")
    for n in [3, 5, 10, 20, 50]:
        Sd, Si, St, _, _ = linear_sensitivity(1.0, n)
        pi, pt = paper3[n]
        ok = "✓" if abs(Si - pi) < 0.002 and abs(St - pt) < 0.002 else "✗"
        print(f"  {n:>4}{Sd:>10.3f}{Si:>12.3f}{St:>10.3f}{ok:>8}")

    # ── TABLE 4 ──────────────────────────────────────────────
    print(f"\n{SEP}")
    print("TABLE 4  Nonlinear System  B=1.0")
    print(SEP)
    print(f"  {'n':>4}{'Direct':>10}{'Indirect':>12}{'Total':>10}"
          f"{'SignRev':>9}{'FD err':>10}")
    print(f"  {SEP2}")
    for n in [10, 20, 50]:
        Sd, Si, St, _, _ = nonlinear_sensitivity(1.0, n)
        St_fd = (Z_nonlinear(1 + 1e-6, n) - Z_nonlinear(1 - 1e-6, n)) / 2e-6
        rev = "YES" if Sd * St < 0 else "no"
        print(f"  {n:>4}{Sd:>10.3f}{Si:>12.3f}{St:>10.3f}"
              f"{rev:>9}{abs(St - St_fd):>10.2e}")

    # ── COURNOT ──────────────────────────────────────────────
    print(f"\n{SEP}")
    print("COURNOT DUOPOLY  B=1.0")
    print(f"  Target: Direct=-0.500 / Indirect=+1.916 / Total=+1.416")
    print(SEP)
    Sd, Si, St, y = cournot_sensitivity(1.0)
    for label, val, tgt in [
            ("Direct", Sd, -0.500),
            ("Indirect", Si, +1.916),
            ("Total", St, +1.416)]:
        ok = "✓" if abs(val - tgt) < 0.01 else "✗"
        print(f"  {label:<12}{val:>+8.3f}   target {tgt:>+6.3f}  {ok}")

    # ── THEOREM 2 BOUNDS ─────────────────────────────────────
    print(f"\n{SEP}")
    print("THEOREM 2  Cournot Bounds  B in [0.5, 2.0]")
    print(SEP)
    LB, UB, bLB, bUB = theorem2(cournot_sensitivity, 0.5, 2.0, n=2)
    print(f"  LB = {LB:.4f}  at B = {bLB:.3f}")
    print(f"  UB = {UB:.4f}  at B = {bUB:.3f}")
    print(f"  Welfare-improving (LB>0): {LB > 0}")


def run_fuzzy():
    """Run fuzzy alpha-cut analysis and return results."""
    print(f"\n{SEP}")
    print("FUZZY EXTENSION — Alpha-Cut Analysis")
    print(SEP)

    res_lin = compute_alpha_cuts(0.8, 1.0, 1.2, linear_sensitivity, N=51)
    res_cou = compute_alpha_cuts(0.5, 1.0, 2.0, cournot_sensitivity, N=51)

    print_summary(res_lin, label="LINEAR  n=3")
    print_summary(res_cou, label="COURNOT")

    return res_lin, res_cou


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(SEP)
    print("BOPEC Sensitivity Analysis — Full Run")
    print("Paper: Aghaei Ata, Jahromi, Keshavarzfard (2025)")
    print(SEP)

    run_tables()
    res_lin, res_cou = run_fuzzy()

    if not args.no_plot:
        print(f"\n{SEP}")
        print(f"Generating Figure 1 -> {args.outdir}/")
        print(SEP)
        make_figure1(res_lin, res_cou, save_dir=args.outdir)

    print(f"\n{SEP}")
    print("DONE")
    print(SEP)


if __name__ == "__main__":
    main()
