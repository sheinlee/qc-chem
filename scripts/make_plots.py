"""Figures from results/*.json  ->  results/figures/*.png (+ results/summary_table.csv).

fig1_ladder:   per-system qubits / Pauli terms / lambda / UCCSD params & CNOTs / QPE-proxy queries, with the
               simulation-regime bands (statevector VQE <=12-14, sparse ED <=16, explicit JW <=24, beyond).
fig2_curves:   N2 and H2O: E - E_CASCI vs R for HF and UCCSD-VQE (+ ADAPT if present), chemical-accuracy line.
fig3_scaling:  Pauli terms and lambda vs active orbitals, power-law fit, extrapolation to 40-60 qubits.
"""
from __future__ import annotations

import csv
import glob
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"
FIG.mkdir(exist_ok=True, parents=True)

# validated categorical palette (fixed order), text tokens
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100", "magenta": "#e87ba4",
     "green": "#008300", "violet": "#6250d6", "red": "#e34948"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
METHOD_COLOR = {"HF": C["blue"], "UCCSD-VQE": C["orange"], "ADAPT-VQE": C["aqua"], "CASCI": INK}
CHEM_ACC = 1.6e-3

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "figure.dpi": 150, "savefig.dpi": 200, "lines.linewidth": 1.6})


def load_all() -> dict:
    out = {}
    for f in sorted(glob.glob(str(RES / "*.json"))):
        d = json.load(open(f, encoding="utf-8"))
        if "points" in d:
            out[d["name"]] = d
    return out


def good_points(d):
    return [p for p in d["points"] if "error" not in p]


# ------------------------------------------------------------------------------------------ table
LADDER_ORDER = ["h2_sto3g", "lih_cas44", "h2o_cas86", "n2_cas66", "ce3_4f_sa", "pr3_4f_sa", "eu3_4f_sa", "gd3_4f",
                "ce3_4f5d", "eu3_4f5d", "cecl3_L1_4f", "cecl3_L2_4f5d", "cecl3_L3_4f5d_sigma",
                "cecl3_L4_4f_cl3p", "cecl3_L5_full"]
PRETTY = {"h2_sto3g": "H$_2$", "lih_cas44": "LiH (4,4)", "h2o_cas86": "H$_2$O (8,6)", "n2_cas66": "N$_2$ (6,6)",
          "ce3_4f_sa": "Ce$^{3+}$ 4f", "pr3_4f_sa": "Pr$^{3+}$ 4f", "eu3_4f_sa": "Eu$^{3+}$ 4f", "gd3_4f": "Gd$^{3+}$ 4f",
          "ce3_4f": "Ce$^{3+}$ 4f (ROHF ref)", "pr3_4f": "Pr$^{3+}$ 4f (ROHF ref)", "eu3_4f": "Eu$^{3+}$ 4f (ROHF ref)",
          "pr3_4f5d": "Pr$^{3+}$ 4f5d", "gd3_4f5d": "Gd$^{3+}$ 4f5d", "cecl3_L1_4f_sa": "CeCl$_3$ 4f (SA-CASSCF)",
          "ce3_4f5d": "Ce$^{3+}$ 4f5d", "eu3_4f5d": "Eu$^{3+}$ 4f5d", "cecl3_L1_4f": "CeCl$_3$ 4f",
          "cecl3_L2_4f5d": "CeCl$_3$ 4f5d", "cecl3_L3_4f5d_sigma": "CeCl$_3$ 4f5d+σ",
          "cecl3_L4_4f_cl3p": "CeCl$_3$ 4f+Cl3p", "cecl3_L5_full": "CeCl$_3$ full"}


def representative(d):
    """One point per system: equilibrium-ish for scans (min CASCI energy), the single point otherwise."""
    pts = good_points(d)
    if not pts:
        return None
    if len(pts) == 1:
        return pts[0]
    return min(pts, key=lambda p: p.get("e_casci") if p.get("e_casci") == p.get("e_casci") else 1e9)


def row_from_point(name, p) -> dict:
    lam_int = p.get("lambda_from_integrals", {})
    jw = p.get("jw", {})
    v = p.get("vqe_uccsd", {})
    tp = p.get("taper", {})
    row = {"system": name, "tag": p.get("tag"), "n_electrons": p.get("n_active_electrons"),
           "n_orbitals": p.get("n_active_orbitals"), "n_qubits": p.get("n_qubits"),
           "n_qubits_tapered": tp.get("n_qubits_tapered"),
           "pauli_terms_jw": jw.get("n_pauli_terms"), "lambda_jw": jw.get("lambda_pauli"),
           "lambda_from_integrals": lam_int.get("lambda_total"), "lambda_kind": lam_int.get("kind"),
           "e_hf": p.get("e_hf"), "e_casci": p.get("e_casci"),
           "ed_err": p.get("ed", {}).get("error_vs_casci"),
           "uccsd_params": v.get("n_params"), "uccsd_cnot": v.get("n_cnot"), "uccsd_depth": v.get("circuit_depth"),
           "uccsd_err": v.get("error_vs_casci"), "uccsd_evals": v.get("n_energy_evals"), "uccsd_t_s": v.get("t_seconds"),
           "qpe_queries": (p.get("qpe_proxy_jw") or p.get("qpe_proxy_from_integrals") or {}).get("queries_qpe"),
           "toffoli_proxy": (p.get("qpe_proxy_jw") or p.get("qpe_proxy_from_integrals") or {}).get("toffoli_proxy_sparse")}
    return row


def write_table(data):
    rows = []
    for name in LADDER_ORDER + [n for n in data if n not in LADDER_ORDER]:
        if name not in data:
            continue
        p = representative(data[name])
        if p:
            r = row_from_point(name, p)
            r["in_ladder_figure"] = name in LADDER_ORDER
            rows.append(r)
    with open(RES / "summary_table.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


# ------------------------------------------------------------------------------------------ fig 2
def fig_curves(data):
    systems = [("n2_cas66", "N$_2$ CAS(6,6), 12 qubits"), ("h2o_cas86", "H$_2$O CAS(8,6), 12 qubits")]
    systems = [(k, t) for k, t in systems if k in data]
    if not systems:
        return
    fig, axes = plt.subplots(2, len(systems), figsize=(4.2 * len(systems), 6.2), squeeze=False, sharex="col")
    for j, (key, title) in enumerate(systems):
        pts = sorted(good_points(data[key]), key=lambda p: float(p["tag"].split("=")[1]))
        R = np.array([float(p["tag"].split("=")[1]) for p in pts])
        e_ref = np.array([p["e_casci"] for p in pts])
        e_hf = np.array([p["e_hf"] for p in pts])
        ax = axes[0, j]
        ax.plot(R, e_ref, color=INK, lw=1.2, label="CASCI (FCI in CAS)")
        ax.plot(R, e_hf, color=METHOD_COLOR["HF"], label="HF")
        if "vqe_uccsd" in pts[0]:
            ax.plot(R, [p["vqe_uccsd"]["energy"] for p in pts], color=METHOD_COLOR["UCCSD-VQE"], ls="--", marker="o", ms=3.5, label="UCCSD-VQE")
        if "vqe_adapt" in pts[0]:
            ax.plot(R, [p["vqe_adapt"]["energy"] for p in pts], color=METHOD_COLOR["ADAPT-VQE"], ls=":", marker="s", ms=3.5, label="ADAPT-VQE")
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_ylabel("E / Ha")
        if j == 0:
            ax.legend(fontsize=8)
        ax = axes[1, j]
        ax.axhspan(0, CHEM_ACC, color=GRID, alpha=0.9, lw=0)
        ax.text(R[-1], CHEM_ACC * 1.15, "chemical accuracy (1.6 mHa)", ha="right", va="bottom", fontsize=7, color=INK2)
        ax.plot(R, e_hf - e_ref, color=METHOD_COLOR["HF"], marker="o", ms=3, label="HF")
        if "vqe_uccsd" in pts[0]:
            err = np.array([p["vqe_uccsd"]["error_vs_casci"] for p in pts])
            ax.plot(R, np.maximum(err, 1e-12), color=METHOD_COLOR["UCCSD-VQE"], marker="o", ms=3, label="UCCSD-VQE")
        if "vqe_adapt" in pts[0]:
            err = np.array([p["vqe_adapt"]["error_vs_casci"] for p in pts])
            ax.plot(R, np.maximum(err, 1e-12), color=METHOD_COLOR["ADAPT-VQE"], marker="s", ms=3, label="ADAPT-VQE")
        ax.set_yscale("log")
        ax.set_ylim(1e-7, 1)
        ax.set_ylabel("E − E$_{CASCI}$ / Ha")
        ax.set_xlabel("R / Å")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_curves.png")
    plt.close(fig)


# ------------------------------------------------------------------------------------------ fig 3
def fig_scaling(data):
    xs, terms, lams, labels, lam_only = [], [], [], [], []
    for name, d in data.items():
        for p in good_points(d):
            n = p.get("n_active_orbitals")
            li = p.get("lambda_from_integrals", {})
            if "jw" in p:
                xs.append(n); terms.append(p["jw"]["n_pauli_terms"]); lams.append(p["jw"]["lambda_pauli"]); labels.append(name)
            elif li:
                lam_only.append((n, li["lambda_total"], name))
    if len(xs) < 3:
        return
    xs, terms, lams = np.array(xs, float), np.array(terms, float), np.array(lams, float)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    # Pauli terms vs N: fit on distinct N (use max per N to avoid scan duplicates dominating)
    ax = axes[0]
    ax.scatter(xs, terms, s=14, color=C["blue"], alpha=0.7, label="explicit JW (≤ 24 qubits)")
    uniq = sorted(set(xs))
    tmax = np.array([terms[xs == u].max() for u in uniq]); uniq = np.array(uniq)
    if len(uniq) >= 3:
        k, b = np.polyfit(np.log(uniq), np.log(tmax), 1)
        xx = np.linspace(min(uniq), 30, 100)
        ax.plot(xx, np.exp(b) * xx ** k, color=C["blue"], lw=1, ls="--", label=f"fit ∝ N$^{{{k:.2f}}}$ (extrapolated to 60 qubits)")
        for n in (16, 21, 30):
            ax.scatter([n], [np.exp(b) * n ** k], marker="x", color=INK2, s=25)
    ax.set_yscale("log"); ax.set_xticks([2, 4, 6, 7, 12, 16, 21, 30])
    ax.set_xlabel("active spatial orbitals N (qubits = 2N)"); ax.set_ylabel("JW Pauli terms")
    ax.legend(fontsize=7.5, loc="upper left")
    ax = axes[1]
    ax.scatter(xs, lams, s=14, color=C["orange"], alpha=0.7, label="λ, explicit JW")
    if lam_only:
        ax.scatter([a for a, _, _ in lam_only], [b for _, b, _ in lam_only], s=22, marker="D", color=C["violet"],
                   label="λ from integrals (exact identity, > 24 qubits)")
        for a, b, nm in lam_only:
            ax.annotate(PRETTY.get(nm, nm), (a, b), fontsize=6.5, color=INK2, xytext=(4, 2), textcoords="offset points")
    ax.set_yscale("log"); ax.set_xticks([2, 4, 6, 7, 12, 16, 17, 21])
    ax.set_xlabel("active spatial orbitals N"); ax.set_ylabel("Pauli 1-norm λ / Ha")
    ax.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_scaling.png")
    plt.close(fig)


# ------------------------------------------------------------------------------------------ fig 1
def fig_ladder(rows):
    rows = [r for r in rows if r["n_qubits"] and r.get("in_ladder_figure")]
    if not rows:
        return
    names = [PRETTY.get(r["system"], r["system"]) for r in rows]
    x = np.arange(len(rows))
    nq = np.array([r["n_qubits"] for r in rows], float)
    fig, axes = plt.subplots(4, 1, figsize=(max(7, 0.62 * len(rows) + 2), 9.5), sharex=True)
    # regime bands drawn on each panel as background by qubit count of the system
    def bands(ax):
        for i, q in enumerate(nq):
            if q <= 14:
                col = "#e8f1fb"   # statevector VQE feasible
            elif q <= 24:
                col = "#fdeee6"   # explicit Pauli statistics / classical CASCI only
            else:
                col = "#f0eef9"   # beyond explicit expansion: lambda from integrals, fault-tolerant target
            ax.axvspan(i - 0.5, i + 0.5, color=col, lw=0, zorder=0)
    ax = axes[0]; bands(ax)
    ax.bar(x, nq, width=0.55, color=C["blue"], label="qubits (JW)")
    nt = [r["n_qubits_tapered"] or np.nan for r in rows]
    ax.scatter(x, nt, color=INK, s=16, zorder=3, label="after Z$_2$ tapering")
    ax.set_ylabel("qubits"); ax.legend(fontsize=7.5, loc="upper left")
    ax = axes[1]; bands(ax)
    terms = np.array([r["pauli_terms_jw"] or np.nan for r in rows], float)
    ax.bar(x, terms, width=0.55, color=C["orange"])
    ax.set_yscale("log"); ax.set_ylabel("JW Pauli terms")
    ax = axes[2]; bands(ax)
    lam = np.array([r["lambda_jw"] if r["lambda_jw"] else r["lambda_from_integrals"] for r in rows], float)
    ax.bar(x, lam, width=0.55, color=C["violet"])
    ax.set_yscale("log"); ax.set_ylabel("λ / Ha")
    ax = axes[3]; bands(ax)
    q = np.array([r["qpe_queries"] or np.nan for r in rows], float)
    ax.bar(x, q, width=0.55, color=C["aqua"])
    ax.set_yscale("log"); ax.set_ylabel("QPE queries\n(πλ/2ε, ε = 1 mHa)")
    params = [r["uccsd_params"] for r in rows]
    for i, (pp, cn) in enumerate(zip(params, [r["uccsd_cnot"] for r in rows])):
        if pp:
            axes[0].text(i, nq[i] + 0.8, f"UCCSD\n{pp} θ\n{cn} CNOT", ha="center", va="bottom", fontsize=6, color=INK2)
    axes[0].set_ylim(0, max(nq) * 1.9)
    axes[-1].set_xticks(x); axes[-1].set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    from matplotlib.patches import Patch
    axes[-1].legend(handles=[Patch(color="#e8f1fb", label="statevector VQE / sparse ED (≤ 14–16 qubits)"),
                             Patch(color="#fdeee6", label="explicit Pauli statistics, classical CASCI (≤ 24)"),
                             Patch(color="#f0eef9", label="λ from integrals only; fault-tolerant target (> 24)")],
                    fontsize=7, loc="upper left")
    fig.suptitle("From H$_2$ to CeCl$_3$: qubit and Hamiltonian resources along the active-space ladder",
                 fontsize=10, color=INK, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG / "fig1_ladder.png")
    plt.close(fig)


if __name__ == "__main__":
    data = load_all()
    rows = write_table(data)
    fig_curves(data)
    fig_scaling(data)
    fig_ladder(rows)
    print("systems:", list(data))
    print("wrote", sorted(p.name for p in FIG.glob("*.png")), "and summary_table.csv")
