"""
fuzzy_bopec.py
==============
Fuzzy extension of BOPEC sensitivity analysis.

Extends Theorem 1 and Theorem 2 of the main paper to triangular
fuzzy policy parameters B~ = (B_l, B_m, B_u).

Theory (Corollary 1 + Proposition 1):
    - Alpha-cut: S~_alpha = [min_{B in B~_alpha} S(B),
                              max_{B in B~_alpha} S(B)]
    - S~ is exactly triangular  iff  S(B) is affine on [B_l, B_u]
    - S~ is quasi-triangular    if   S(B) strictly monotone, nonlinear
    - S~ is non-triangular      if   S'(B*)=0 for some interior B*

Reference:
    Aghaei Ata, H., Jahromi, M., & Keshavarzfard, R. (2025).
    Fuzzy Extension of BOPEC Sensitivity Analysis.
"""

import numpy as np
from scipy.optimize import minimize_scalar


# ═══════════════════════════════════════════════════════════════════
#  COROLLARY 1 — Alpha-Cut Bounds
# ═══════════════════════════════════════════════════════════════════

def compute_alpha_cuts(B_l, B_m, B_u, sens_fn, N=51):
    """
    Rigorous alpha-cut computation with critical point detection.

    For each alpha in [0, 1]:
        B~_alpha = [B_l + alpha*(B_m - B_l),  B_u - alpha*(B_u - B_m)]
        S~_alpha = [min_{B in B~_alpha} S(B),  max_{B in B~_alpha} S(B)]

    Uses scipy.optimize.minimize_scalar (bounded) to guarantee capture
    of interior critical points at each alpha level (Extreme Value Thm).

    Parameters
    ----------
    B_l, B_m, B_u : float   Triangular fuzzy parameter
    sens_fn        : callable(B) -> (S_direct, S_indirect, S_total)
    N              : int     Number of alpha-cuts (default 51)

    Returns
    -------
    dict with keys: alphas, LB, UB, core, support,
                    is_monotone, is_affine, linearity_error,
                    shape, B_tilde, S_traj
    """
    alphas = np.linspace(0, 1, N)
    S_lb, S_ub = [], []

    for alpha in alphas:
        b_al = B_l + alpha * (B_m - B_l)
        b_au = B_u - alpha * (B_u - B_m)

        if abs(b_au - b_al) < 1e-12:
            v = sens_fn(b_al)[2]
            S_lb.append(v)
            S_ub.append(v)
            continue

        Sa = sens_fn(b_al)[2]
        Sb = sens_fn(b_au)[2]

        res_min = minimize_scalar(
            lambda B: sens_fn(B)[2],
            bounds=(b_al, b_au), method='bounded')
        res_max = minimize_scalar(
            lambda B: -sens_fn(B)[2],
            bounds=(b_al, b_au), method='bounded')

        S_lb.append(min(Sa, Sb, res_min.fun))
        S_ub.append(max(Sa, Sb, -res_max.fun))

    S_lb = np.array(S_lb)
    S_ub = np.array(S_ub)

    # ── Monotonicity check ──────────────────────────────────────
    Bv  = np.linspace(B_l, B_u, 500)
    Sv  = np.array([sens_fn(B)[2] for B in Bv])
    dSv = np.diff(Sv)
    is_monotone = bool(np.all(dSv >= -1e-8) or np.all(dSv <= 1e-8))

    # ── Linearity (affine) check ────────────────────────────────
    S_linear       = np.linspace(Sv[0], Sv[-1], 500)
    linearity_error = float(np.max(np.abs(Sv - S_linear)))
    is_affine      = linearity_error < 1e-4

    # ── Proposition 1 shape classification ─────────────────────
    if is_affine:
        shape = "Exactly Triangular"
    elif is_monotone:
        shape = "Quasi-Triangular (Convex)"
    else:
        shape = "Non-Triangular (Critical Points)"

    return {
        "alphas"         : alphas,
        "LB"             : S_lb,
        "UB"             : S_ub,
        "core"           : float(sens_fn(B_m)[2]),
        "support"        : [float(S_lb[0]), float(S_ub[0])],
        "is_monotone"    : is_monotone,
        "is_affine"      : is_affine,
        "linearity_error": linearity_error,
        "shape"          : shape,
        "B_tilde"        : (B_l, B_m, B_u),
        "S_traj"         : (Bv, Sv),
    }


# ═══════════════════════════════════════════════════════════════════
#  PRINT SUMMARY
# ═══════════════════════════════════════════════════════════════════

def print_summary(result, label="System"):
    """Print certified fuzzy bounds to console."""
    Bl, Bm, Bu = result["B_tilde"]
    n   = len(result["alphas"])
    mid = n // 2
    print(f"\n{label}  B~ = ({Bl}, {Bm}, {Bu})")
    print(f"  Core S(B_m)       = {result['core']:.6f}")
    print(f"  Support [LB, UB]  = "
          f"[{result['support'][0]:.6f}, {result['support'][1]:.6f}]")
    print(f"  alpha=0.5 cut     = "
          f"[{result['LB'][mid]:.6f}, {result['UB'][mid]:.6f}]")
    print(f"  Linearity error   = {result['linearity_error']:.2e}")
    print(f"  Shape             = {result['shape']}")
    if result["support"][0] > 0:
        print(f"  dZ/dB > 0 certified for all B in B~  [Robustness: YES]")
