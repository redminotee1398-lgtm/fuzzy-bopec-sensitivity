# Fuzzy BOPEC Sensitivity Analysis

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Paper](https://img.shields.io/badge/Paper-ORL%20(under%20review)-orange)](https://www.sciencedirect.com/journal/operations-research-letters)

Code repository for:

> **Analytical Sensitivity Bounds for System-Level Performance in Bi-Level Optimization with Nash Equilibrium Constraints**  
> H. Aghaei Ata, M. Jahromi, R. Keshavarzfard  
> *Operations Research Letters* (under review, 2026)

and its fuzzy extension:

> **Fuzzy Sensitivity Bounds for BOPEC: Rigorous Alpha-Cut Extension with Shape Characterization**

---

## Overview

This repository provides a clean Python implementation of:

1. **Theorem 1** — Closed-form sensitivity decomposition:

$$\frac{dZ}{dB} = \underbrace{\frac{\partial Z}{\partial x}\frac{df}{dB}}_{\text{direct}} + \underbrace{\frac{\partial Z}{\partial y^*}\!\left(-H^{-1}G\frac{df}{dB}\right)}_{\text{indirect}}$$

2. **Theorem 2** — Exact interval bounds via endpoint evaluation (Extreme Value Theorem)

3. **Corollary 1** (Fuzzy Extension) — Alpha-cut certified bounds for triangular fuzzy parameter $\tilde{B}=(B^l, B^m, B^u)$:

$$\tilde{S}_\alpha = \left[\min_{B \in \tilde{B}_\alpha} S(B),\ \max_{B \in \tilde{B}_\alpha} S(B)\right]$$

4. **Proposition 1** — Shape characterization: S̃ is *exactly triangular* iff S(B) is affine on the support; *quasi-triangular* if strictly monotone but nonlinear; *non-triangular* if interior critical points exist.

---

## Repository Structure

```
fuzzy-bopec-sensitivity/
├── bopec_sensitivity.py   # Theorem 1 & 2: core decomposition
├── fuzzy_bopec.py         # Corollary 1: alpha-cut fuzzy extension
├── plots.py               # Publication-quality Figure 1 (4-panel)
├── main.py                # Entry point: tables + figure
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/redminotee1398-lgtm/fuzzy-bopec-sensitivity.git
cd fuzzy-bopec-sensitivity

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run full analysis (tables + figure)
python main.py

# 4. Tables only (no figure)
python main.py --no-plot

# 5. Custom output directory
python main.py --outdir ~/Desktop/bopec_results
```

Output is saved to `results/`:

```
results/
├── fuzzy_bopec_figure1.png   (300 dpi)
└── fuzzy_bopec_figure1.pdf
```

---

## Verified Numerical Results

### Table 1 — Linear System (B = 1.0, n = 3)

| Component | Paper     | Code      | Match |
|-----------|-----------|-----------|-------|
| Direct    | +9.200000 | +9.200000 | ✓     |
| Indirect  | −0.502755 | −0.502755 | ✓     |
| Total     | +8.697245 | +8.697245 | ✓     |

### Cournot Duopoly (B = 1.0)

| Component | Paper  | Code   | Match |
|-----------|--------|--------|-------|
| Direct    | −0.500 | −0.500 | ✓     |
| Indirect  | +1.916 | +1.916 | ✓     |
| Total     | +1.416 | +1.416 | ✓     |

### Fuzzy Extension — Corollary 1 (N = 51 alpha-cuts)

| System     | $\tilde{B}$     | Core  | Certified [LB, UB] | Lin. Error | Shape              |
|------------|-----------------|-------|--------------------|------------|--------------------|
| Linear n=3 | (0.8, 1.0, 1.2) | 8.697 | [8.582, 8.813]     | 1.78e-15   | Exactly Triangular ✓ |
| Cournot    | (0.5, 1.0, 2.0) | 1.416 | [0.972, 1.638]     | 6.66e-16   | Exactly Triangular ✓ |

> **Note:** Both systems yield *exactly triangular* fuzzy outputs because their sensitivity functions S(B) are affine (linearity error < 10⁻¹⁰) on the respective support intervals, as verified by Proposition 1 (Case 1).

---

## Figure 1 Preview

Four-panel figure produced by `main.py`:

| Panel | Content |
|-------|---------|
| (A)   | Linear system: S_direct, S_indirect, S_total over B ∈ [0.55, 1.45] |
| (B)   | Fuzzy membership function — Linear n=3 (exactly triangular) |
| (C)   | Cournot duopoly: S_direct, S_indirect, S_total over B ∈ [0.30, 2.20] |
| (D)   | Fuzzy membership function — Cournot (exactly triangular, wider spread) |

---

## Cite

If you use this code, please cite:

```bibtex
@article{aghaeiata2026bopec,
  title   = {Analytical Sensitivity Bounds for System-Level Performance
             in Bi-Level Optimization with Nash Equilibrium Constraints},
  author  = {Aghaei Ata, Hossein and Jahromi, Meghdad and Keshavarzfard, Razieh},
  journal = {Operations Research Letters},
  year    = {2026},
  note    = {Under review}
}

@software{aghaeiata2026bopec_code,
  author  = {Aghaei Ata, Hossein and Jahromi, Meghdad and Keshavarzfard, Razieh},
  title   = {fuzzy-bopec-sensitivity: Python implementation of BOPEC
             sensitivity analysis with fuzzy extension},
  year    = {2026},
  url     = {https://github.com/redminotee1398-lgtm/fuzzy-bopec-sensitivity}
}
```

---

## Dependencies

| Package    | Version |
|------------|---------|
| Python     | ≥ 3.9   |
| NumPy      | ≥ 1.24  |
| SciPy      | ≥ 1.10  |
| Matplotlib | ≥ 3.7   |

Install all at once:

```bash
pip install -r requirements.txt
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

**Meghdad Jahromi** (Corresponding Author)  
Department of Industrial Engineering  
North Tehran Branch, Islamic Azad University, Tehran, Iran  
📧 Meghdadjahromi@iau.ac.ir  
🔗 ORCID: [0000-0001-6360-5347](https://orcid.org/0000-0001-6360-5347)

**Hossein Aghaei Ata**  
📧 aghaeiata@yahoo.co.uk

**Razieh Keshavarzfard**  
📧 r.keshavarzfard@iau.ir  
🔗 ORCID: [0000-0002-2212-3576](https://orcid.org/0000-0002-2212-3576)
