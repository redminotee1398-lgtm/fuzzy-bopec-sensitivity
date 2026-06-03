"""
plots.py
========
Publication-quality figures for BOPEC fuzzy sensitivity analysis.

Figure 1 (4-panel):
    (A) Linear system sensitivity trajectory S_direct, S_indirect, S_total
    (B) Linear fuzzy membership function (exactly triangular)
    (C) Cournot duopoly sensitivity trajectory
    (D) Cournot fuzzy membership function (quasi-triangular)
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from bopec_sensitivity import linear_sensitivity, cournot_sensitivity

# ── Style constants ──────────────────────────────────────────────
C_DIR   = "#1f77b4"    # blue   — direct effect
C_IND   = "#d62728"    # red    — indirect effect
C_TOT   = "#2ca02c"    # green  — total (linear)
C_TOT_C = "#9467bd"    # purple — total (Cournot)
C_GOLD  = "#FFD700"    # gold   — support region
NAVY    = "#1F3864"    # navy   — core line / title
FONT    = "DejaVu Sans"


def _style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.28, linestyle="--")


def _trajectory_panel(ax, B_range, sens_fn, B_tilde, color_tot,
                       panel_label, title):
    """
    Plot S_direct, S_indirect, S_total over B_range with support shading.
    """
    Bl, Bm, Bu = B_tilde
    Bv  = B_range
    Sd_v = np.array([sens_fn(B)[0] for B in Bv])
    Si_v = np.array([sens_fn(B)[1] for B in Bv])
    St_v = np.array([sens_fn(B)[2] for B in Bv])

    ax.axvspan(Bl, Bu, alpha=0.11, color=C_GOLD, zorder=0,
               label=r"Support $[B^l, B^u]$")
    ax.plot(Bv, Sd_v, color=C_DIR, lw=1.8, ls="--",
            label=r"$S_{\rm direct}$")
    ax.plot(Bv, Si_v, color=C_IND, lw=1.8, ls=":",
            label=r"$S_{\rm indirect}$")
    ax.plot(Bv, St_v, color=color_tot, lw=2.5,
            label=r"$S_{\rm total}$")
    ax.axvline(Bm, color=NAVY, lw=1.3, ls="--", alpha=0.70,
               label=f"$B^m={Bm}$")
    ax.axhline(0, color="black", lw=0.7, ls=":", alpha=0.45)

    # Annotate key S_total values
    for B_, ha in [(Bl, "right"), (Bm, "left"), (Bu, "left")]:
        sv = sens_fn(B_)[2]
        ax.plot(B_, sv, "o", color="#333333", ms=5.5, zorder=6)
        ax.annotate(f"{sv:.3f}", (B_, sv),
                    xytext=(5 if ha == "left" else -5, 6),
                    textcoords="offset points",
                    fontsize=8.5, ha=ha, color="#333333")

    ax.set_xlabel(r"Policy parameter $B$", fontsize=11)
    ax.set_ylabel(r"$dZ/dB$", fontsize=11)
    ax.set_title(f"({panel_label}) {title}", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, loc="best", framealpha=0.85)
    _style_ax(ax)


def _membership_panel(ax, result, panel_label, title,
                      color_lb, color_ub):
    """
    Plot fuzzy membership function from alpha-cut results.
    """
    a    = result["alphas"]
    lb   = result["LB"]
    ub   = result["UB"]
    core = result["core"]
    sup  = result["support"]
    lin_err = result["linearity_error"]
    shape   = result["shape"]

    # Fill
    ax.fill_betweenx(a, lb, ub, alpha=0.18, color=color_lb, zorder=0)

    # Boundaries
    ax.plot(lb, a, color=color_lb, lw=2.5, label="Lower bound")
    ax.plot(ub, a, color=color_ub, lw=2.5, label="Upper bound")

    # Ideal TFN reference (dashed black)
    ax.plot([sup[0], core, sup[1]], [0, 1, 0],
            color="black", lw=1.1, ls="--", alpha=0.45,
            label="Ideal TFN (ref.)")

    # Core line
    ax.axvline(core, color=NAVY, lw=1.7, ls="--",
               label=f"Core = {core:.3f}")

    # LB / UB annotations
    ax.annotate(f"LB={sup[0]:.3f}", (sup[0], 0.04),
                fontsize=8.5, color=color_lb,
                ha="left", fontweight="bold")
    ax.annotate(f"UB={sup[1]:.3f}", (sup[1], 0.04),
                fontsize=8.5, color=color_ub,
                ha="right", fontweight="bold")

    ax.set_xlabel(r"Sensitivity value $dZ/dB$", fontsize=11)
    ax.set_ylabel(r"Membership degree $\alpha$", fontsize=11)
    ax.set_title(
        f"({panel_label}) {title}\n"
        f"Shape: {shape}  |  Lin. err = {lin_err:.1e}",
        fontsize=10, fontweight="bold")
    ax.set_ylim(-0.06, 1.12)
    ax.legend(fontsize=9, framealpha=0.85)
    _style_ax(ax)


def make_figure1(res_lin, res_cou, save_dir="results"):
    """
    Generate Figure 1 (4-panel) and save to save_dir.

    Parameters
    ----------
    res_lin  : dict   Output of compute_alpha_cuts for linear system
    res_cou  : dict   Output of compute_alpha_cuts for Cournot
    save_dir : str    Directory to save PNG and PDF
    """
    os.makedirs(save_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.patch.set_facecolor("white")
    fig.suptitle(
        r"Fuzzy BOPEC Sensitivity Analysis — $\alpha$-Cut Certified Bounds",
        fontsize=13, fontweight="bold", color=NAVY, y=0.995)

    Bl_l, Bm_l, Bu_l = res_lin["B_tilde"]
    Bl_c, Bm_c, Bu_c = res_cou["B_tilde"]

    # (A) Linear trajectory
    _trajectory_panel(
        axes[0, 0],
        np.linspace(0.55, 1.45, 400),
        linear_sensitivity,
        (Bl_l, Bm_l, Bu_l),
        C_TOT,
        "A", r"Sensitivity Trajectory — Linear ($n=3$)")

    # (B) Linear membership
    _membership_panel(
        axes[0, 1], res_lin,
        "B", r"Fuzzy Membership — Linear ($n=3$)",
        C_TOT, C_DIR)

    # (C) Cournot trajectory
    _trajectory_panel(
        axes[1, 0],
        np.linspace(0.30, 2.20, 400),
        cournot_sensitivity,
        (Bl_c, Bm_c, Bu_c),
        C_TOT_C,
        "C", "Sensitivity Trajectory — Cournot Duopoly")

    # (D) Cournot membership
    _membership_panel(
        axes[1, 1], res_cou,
        "D", "Fuzzy Membership — Cournot Duopoly",
        C_TOT_C, "#e377c2")

    plt.tight_layout(rect=[0, 0, 1, 0.985], pad=2.2)

    path_png = os.path.join(save_dir, "fuzzy_bopec_figure1.png")
    path_pdf = os.path.join(save_dir, "fuzzy_bopec_figure1.pdf")
    fig.savefig(path_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(path_pdf, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path_png}")
    print(f"  Saved: {path_pdf}")
    return path_png, path_pdf
