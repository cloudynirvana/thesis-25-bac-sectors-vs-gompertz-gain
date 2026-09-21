#!/usr/bin/env python3
"""Shared subsystem toy: BAC sector labels and Gompertz residual tests.

The weight paths, the grounded Laplacian, and the sector classifier are the
objects of Thesis #18. The load×gain hazard and the weighted residual test
are the objects of Thesis #13. Both are evaluated on the same five-index
weight paths. No number is copied from either deposit's results file.

The hazard grid is deterministic. Seed 20260921 is used only for the
monotonicity draws on the Laplacian.

Research only. Not a medical device. Not a dose. Not a rejuvenation claim.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
SEED = 20260921

LABELS = ["molecular", "cellular", "tissue", "organism", "evolutionary"]
N = 5
ORG = 3
GROUND = 4

# Classifier and weight schedules, fixed as in Thesis #18 before any label is read.
DELTA_AGE = 0.04
DELTA_CUT = 0.08
DELTA_AMB_CUT = 0.06
DELTA_AMB_BLOCK = 0.02
SIGMA_HOLD = 0.15
T_END = 80.0
N_FINE = 801
DIAG0 = 0.90

THR_LAM_FALL = 0.5
THR_BOTH_FALL = 0.5
THR_SECTOR_BAND = 1.5
THR_CUT_COLLAPSE = 0.25
THR_BLOCK_HELD = 0.9
THR_SECTOR_RISE = 3.0

# Hazard side. Cuts and the gain list follow the Thesis #13 procedure.
# They are applied here to this five-index toy, not to that deposit's six-node output.
GAMMA = 8.0
MU_B = 2.0e-4
LOG_MU_B = float(np.log(MU_B))
V = np.array([0.010, 0.012, 0.009, 0.011, 0.008], dtype=float)
D0 = np.array([0.20, 0.18, 0.22, 0.19, 0.21], dtype=float)
VBAR = float(np.mean(V))
L0 = float(np.mean(D0))

E0 = 2.0e5
EXPOSURE_DECAY = 0.035
EXPOSURE_FLOOR = 500.0

# Named-path gain scan, the Thesis #13 list, declared before the calls are read.
G_SCAN = np.array([0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08], dtype=float)
# Plane panels. Subset of G_SCAN, not chosen after looking at a Jaccard index.
G_PLANE = np.array([0.0, 0.005, 0.02], dtype=float)
# Hazard-map gain factors. Damage does not depend on γ.
GAMMA_FACTORS = np.array([0.25, 0.5, 1.0, 2.0, 4.0], dtype=float)
GAMMA_AT_G = (0.0, 0.005)

CHI2_95_1 = 3.841
NEAR_LINEAR_R2 = 0.995
NEAR_LINEAR_MAX_REL = 0.02
LOG_MU_CAP = 80.0  # above this, the exponential hazard is left unscored

N_GRID = 13
KINDS = ("held", "aging", "cancer", "ambiguous", "diagonal")


def initial_weights() -> np.ndarray:
    pairs = {
        (0, 1): 0.80,
        (0, 2): 0.35,
        (0, 3): 0.55,
        (0, 4): 0.20,
        (1, 2): 0.75,
        (1, 3): 0.60,
        (1, 4): 0.25,
        (2, 3): 0.70,
        (2, 4): 0.30,
        (3, 4): 0.40,
    }
    w = np.zeros((N, N), dtype=float)
    for (i, j), value in pairs.items():
        w[i, j] = w[j, i] = value
    return w


def symmetrize(weights: np.ndarray) -> np.ndarray:
    a = 0.5 * (np.asarray(weights, dtype=float) + np.asarray(weights, dtype=float).T)
    np.fill_diagonal(a, 0.0)
    return a


def laplacian(weights: np.ndarray) -> np.ndarray:
    a = symmetrize(weights)
    return np.diag(a.sum(axis=1)) - a


def algebraic_connectivity(weights: np.ndarray) -> float:
    ev = np.linalg.eigvalsh(laplacian(weights))
    return float(np.sort(ev)[1])


def grounded_lambda_min(weights: np.ndarray, ground: int = GROUND) -> float:
    ell = laplacian(weights)
    keep = [i for i in range(N) if i != ground]
    return float(np.linalg.eigvalsh(ell[np.ix_(keep, keep)])[0])


def induced_block_lambda2(weights: np.ndarray, drop: int = ORG) -> float:
    keep = [i for i in range(N) if i != drop]
    block = symmetrize(weights)[np.ix_(keep, keep)]
    return float(np.sort(np.linalg.eigvalsh(laplacian(block)))[1])


def cut_and_block_means(weights: np.ndarray) -> tuple[float, float]:
    a = symmetrize(weights)
    cut = []
    block = []
    for i in range(N):
        for j in range(i + 1, N):
            if i == ORG or j == ORG:
                cut.append(a[i, j])
            else:
                block.append(a[i, j])
    return float(np.mean(cut)), float(np.mean(block))


def edge_rates(kind: str) -> np.ndarray:
    rates = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(i + 1, N):
            incident = i == ORG or j == ORG
            if kind == "held":
                rate = 0.0
            elif kind == "aging":
                rate = DELTA_AGE
            elif kind == "cancer":
                rate = DELTA_CUT if incident else 0.0
            elif kind == "ambiguous":
                rate = DELTA_AMB_CUT if incident else DELTA_AMB_BLOCK
            elif kind == "diagonal":
                rate = 0.0
            else:
                raise ValueError(kind)
            rates[i, j] = rates[j, i] = rate
    return rates


def plane_rates(delta_cut: float, delta_block: float) -> np.ndarray:
    rates = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(i + 1, N):
            rate = delta_cut if (i == ORG or j == ORG) else delta_block
            rates[i, j] = rates[j, i] = rate
    return rates


def weights_at(t: float, rates: np.ndarray, w0: np.ndarray) -> np.ndarray:
    w = np.zeros_like(w0)
    for i in range(N):
        for j in range(i + 1, N):
            value = w0[i, j] * np.exp(-rates[i, j] * t)
            w[i, j] = w[j, i] = value
    return w


def classify(lam_ratio: float, block_ratio: float, cut_ratio: float, sector_change: float) -> str:
    aging = (
        lam_ratio < THR_LAM_FALL
        and block_ratio < THR_BOTH_FALL
        and cut_ratio < THR_BOTH_FALL
        and (1.0 / THR_SECTOR_BAND) < sector_change < THR_SECTOR_BAND
    )
    cancer = (
        lam_ratio < THR_LAM_FALL
        and cut_ratio < THR_CUT_COLLAPSE
        and block_ratio > THR_BLOCK_HELD
        and sector_change > THR_SECTOR_RISE
    )
    if aging and cancer:
        return "both"
    if aging:
        return "aging-like"
    if cancer:
        return "cancer-like"
    return "neither"


def first_crossing(t: np.ndarray, lam: np.ndarray, sigma: float) -> float | None:
    v = lam - sigma
    if np.all(v > 0):
        return None
    idx = int(np.argmax(v <= 0))
    if idx == 0:
        return float(t[0])
    t0, t1 = float(t[idx - 1]), float(t[idx])
    v0, v1 = float(v[idx - 1]), float(v[idx])
    if v1 == v0:
        return t1
    return t0 + (0.0 - v0) * (t1 - t0) / (v1 - v0)


def bac_series(rates: np.ndarray, w0: np.ndarray, t: np.ndarray) -> dict:
    lam = np.empty_like(t)
    lam2 = np.empty_like(t)
    mean_cut = np.empty_like(t)
    mean_block = np.empty_like(t)
    for k, time in enumerate(t):
        w = weights_at(float(time), rates, w0)
        lam[k] = grounded_lambda_min(w)
        lam2[k] = algebraic_connectivity(w)
        mean_cut[k], mean_block[k] = cut_and_block_means(w)
    sector0 = mean_block[0] / mean_cut[0]
    sector_t = mean_block[-1] / mean_cut[-1]
    label = classify(
        float(lam[-1] / lam[0]),
        float(mean_block[-1] / mean_block[0]),
        float(mean_cut[-1] / mean_cut[0]),
        float(sector_t / sector0),
    )
    return {
        "t": t,
        "lambda_min": lam,
        "lambda2": lam2,
        "mean_cut": mean_cut,
        "mean_block": mean_block,
        "class": label,
        "lambda_min_ratio": float(lam[-1] / lam[0]),
        "mean_cut_ratio": float(mean_cut[-1] / mean_cut[0]),
        "mean_block_ratio": float(mean_block[-1] / mean_block[0]),
        "sector_change": float(sector_t / sector0),
        "crossing_hold": first_crossing(t, lam, SIGMA_HOLD),
        "lambda_min_0": float(lam[0]),
        "lambda_min_T": float(lam[-1]),
        "lambda2_0": float(lam2[0]),
        "lambda2_T": float(lam2[-1]),
    }


def coupling_operator(w0: np.ndarray) -> tuple[np.ndarray, float]:
    """A_ij = W_ij(0) / mean initial off-diagonal weight, zero diagonal."""
    a = symmetrize(w0)
    off = [a[i, j] for i in range(N) for j in range(i + 1, N)]
    mean_w = float(np.mean(off))
    op = a / mean_w
    np.fill_diagonal(op, 0.0)
    return op, mean_w


def integrate_damage(rates: np.ndarray, g: float, t_out: np.ndarray, op: np.ndarray) -> np.ndarray:
    """dD/dt = v + (g/(n-1)) (W(t)/mean W0) D, with W_ij(t) = W_ij(0) exp(-δ_ij t)."""
    if abs(g) < 1e-15:
        return D0[None, :] + np.outer(t_out, V)
    scale = g / (N - 1)

    def f(_t: float, y: np.ndarray) -> np.ndarray:
        atten = np.exp(-rates * _t)
        return V + scale * ((atten * op) @ y)

    sol = solve_ivp(
        f,
        (float(t_out[0]), float(t_out[-1])),
        D0.copy(),
        t_eval=t_out,
        method="DOP853",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success or sol.y.shape[1] != t_out.size:
        raise RuntimeError(f"damage integration failed: {sol.message}")
    return sol.y.T


def constant_coefficient_damage(g: float, t_out: np.ndarray, op: np.ndarray) -> np.ndarray:
    """Closed form when every rate is zero: dD/dt = v + (g/(n-1)) A D."""
    if abs(g) < 1e-15:
        return D0[None, :] + np.outer(t_out, V)
    mate = (g / (N - 1)) * op
    out = np.empty((t_out.size, N), dtype=float)
    mate_inv = np.linalg.inv(mate)
    eye = np.eye(N)
    for i, ti in enumerate(t_out):
        big = expm(mate * float(ti))
        out[i] = big @ D0 + mate_inv @ ((big - eye) @ V)
    return out


def r2(y: np.ndarray, yhat: np.ndarray) -> float:
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 0.0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def linearity(x: np.ndarray, y: np.ndarray) -> dict:
    design = np.column_stack([np.ones_like(x), x])
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    yhat = design @ coef
    increment = abs(float(coef[1]) * (float(x[-1]) - float(x[0])))
    max_abs = float(np.max(np.abs(y - yhat)))
    return {
        "r2": r2(y, yhat),
        "max_rel_resid": max_abs / max(increment, 1e-12),
    }


def node_linearity(t: np.ndarray, damage: np.ndarray) -> dict:
    stats = [linearity(t, damage[:, i]) for i in range(N)]
    min_r2 = float(min(s["r2"] for s in stats))
    max_rel = float(max(s["max_rel_resid"] for s in stats))
    return {
        "min_r2": min_r2,
        "max_rel_resid": max_rel,
        "near_linear": bool(min_r2 >= NEAR_LINEAR_R2 and max_rel <= NEAR_LINEAR_MAX_REL),
    }


def exposure(t: np.ndarray) -> np.ndarray:
    return np.maximum(E0 * np.exp(-EXPOSURE_DECAY * t), EXPOSURE_FLOOR)


def reference_sigma(t: np.ndarray) -> np.ndarray:
    """Noise budget of the linear load at the declared γ. It does not move with g."""
    load = L0 + VBAR * t
    mu = np.exp(LOG_MU_B + GAMMA * load)
    deaths = np.maximum(mu * exposure(t), 1e-8)
    return 1.0 / np.sqrt(deaths)


def wls_poly(t: np.ndarray, y: np.ndarray, sigma: np.ndarray, deg: int) -> tuple[np.ndarray, float]:
    scale = 1.0 / sigma
    design = np.column_stack([(t**k) * scale for k in range(deg + 1)])
    coef, _, _, _ = np.linalg.lstsq(design, y * scale, rcond=None)
    yhat = sum(float(coef[k]) * t**k for k in range(deg + 1))
    rss = float(np.sum(((y - yhat) / sigma) ** 2))
    return coef, rss


def log_chord(log0: float, log1: float, t: np.ndarray) -> np.ndarray:
    """Log of the straight line in μ that matches the exponential map at the endpoints."""
    span = float(t[-1] - t[0])
    s = (t - float(t[0])) / span
    out = np.empty_like(s, dtype=float)
    out[0] = log0
    out[-1] = log1
    mid = (s > 0.0) & (s < 1.0)
    sm = s[mid]
    out[mid] = np.logaddexp(log0 + np.log1p(-sm), log1 + np.log(sm))
    return out


def gompertz_call(t: np.ndarray, log_mu: np.ndarray, sigma: np.ndarray, crit: float) -> dict:
    finite = bool(np.all(np.isfinite(log_mu)) and np.max(log_mu) <= LOG_MU_CAP and np.min(log_mu) > -LOG_MU_CAP)
    if not finite:
        return {
            "scored": False,
            "alpha": None,
            "beta": None,
            "quadratic_c": None,
            "rss_linear": None,
            "rss_quadratic": None,
            "lr_vs_quadratic": None,
            "gompertz_like": False,
            "curvature_detectable": True,
            "reason": "log hazard left the declared scoring window",
        }
    coef1, rss1 = wls_poly(t, log_mu, sigma, 1)
    coef2, rss2 = wls_poly(t, log_mu, sigma, 2)
    lr = rss1 - rss2
    return {
        "scored": True,
        "alpha": float(coef1[0]),
        "beta": float(coef1[1]),
        "quadratic_c": float(coef2[2]),
        "rss_linear": rss1,
        "rss_quadratic": rss2,
        "lr_vs_quadratic": lr,
        "gompertz_like": bool(rss1 <= crit),
        "curvature_detectable": bool(lr > CHI2_95_1),
        "reason": None,
    }


def hazard_from_damage(t: np.ndarray, damage: np.ndarray, gamma: float, kind: str) -> np.ndarray:
    load = damage.mean(axis=1)
    log_exp = LOG_MU_B + gamma * load
    if kind == "exponential":
        return log_exp
    if kind == "additive":
        return log_chord(float(log_exp[0]), float(log_exp[-1]), t)
    raise ValueError(kind)


def coincidence(labels: list[str], flags: list[bool]) -> dict:
    aging = np.array([lab == "aging-like" for lab in labels], dtype=bool)
    gomp = np.array(flags, dtype=bool)
    both = int(np.sum(aging & gomp))
    aging_only = int(np.sum(aging & ~gomp))
    gomp_only = int(np.sum(~aging & gomp))
    neither = int(np.sum(~aging & ~gomp))
    union = both + aging_only + gomp_only
    return {
        "n": int(aging.size),
        "n_aging": int(np.sum(aging)),
        "n_gompertz": int(np.sum(gomp)),
        "aging_and_gompertz": both,
        "aging_not_gompertz": aging_only,
        "gompertz_not_aging": gomp_only,
        "neither_label": neither,
        "sets_equal": bool(aging_only == 0 and gomp_only == 0),
        "jaccard": None if union == 0 else both / union,
    }


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
        }
    )


COLORS = {
    "held": "#4C5760",
    "aging": "#1F4E79",
    "cancer": "#8C3A3A",
    "ambiguous": "#8A6A12",
    "diagonal": "#2F6B4F",
}
NAMES = {
    "held": "held weights",
    "aging": "global decay",
    "cancer": "organism-scale cut",
    "ambiguous": "uneven decay",
    "diagonal": "diagonal only",
}


def draw(
    bac: dict[str, dict],
    t_year: np.ndarray,
    log_mu: dict[tuple[str, float], np.ndarray],
    scan_rows: list[dict],
    plane: dict[float, dict],
    crit: float,
) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    style()

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for key in KINDS:
        s = bac[key]
        ax.plot(
            s["t"],
            s["lambda_min"],
            color=COLORS[key],
            lw=1.6,
            ls="--" if key == "diagonal" else "-",
            label=f"{NAMES[key]} ({s['class']})",
        )
    ax.axhline(SIGMA_HOLD, color="#666666", ls="--", lw=0.9, label="proxy σ = 0.15")
    ax.set_xlabel("toy time")
    ax.set_ylabel("λ_min of the grounded operator")
    ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIG / "lambda_min_with_classes.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), sharey=False)
    for ax, g, title in (
        (axes[0], 0.0, "exponential map, g = 0"),
        (axes[1], 0.02, "exponential map, g = 0.02"),
    ):
        for key in KINDS:
            y = log_mu[(key, g)] / np.log(10.0)
            ax.plot(
                t_year,
                y,
                color=COLORS[key],
                lw=1.5,
                ls="--" if key == "diagonal" else "-",
                label=NAMES[key],
            )
        ax.set_xlabel("toy time")
        ax.set_title(title)
        ax.set_ylabel("log10 hazard")
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "hazard_same_paths.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for key in KINDS:
        rows = [r for r in scan_rows if r["kind"] == key and r["map"] == "exponential" and r["gamma"] == GAMMA]
        gs = [r["g"] for r in rows]
        rss = [np.nan if r["rss_linear"] is None else r["rss_linear"] for r in rows]
        ax.plot(
            gs,
            rss,
            color=COLORS[key],
            marker="o",
            ms=3.5,
            lw=1.4,
            ls="--" if key == "diagonal" else "-",
            label=f"{NAMES[key]} ({bac[key]['class']})",
        )
    ax.axhline(crit, color="#666666", ls="--", lw=0.9, label=f"χ² budget {crit:.2f}")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("damage-rate gain g")
    ax.set_ylabel("weighted residual of the best Gompertz line")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "residual_versus_gain.png")
    plt.close(fig)

    deltas = plane["deltas"]
    fig, axes = plt.subplots(1, len(G_PLANE), figsize=(9.2, 3.5), sharex=True, sharey=True)
    marker = {"aging-like": "o", "cancer-like": "s", "neither": "^", "both": "X"}
    for ax, g in zip(axes, G_PLANE):
        block = plane["by_g"][float(g)]
        for lab, mk in marker.items():
            pts = [p for p in block if p["class"] == lab]
            if not pts:
                continue
            xs = [p["delta_block"] for p in pts]
            ys = [p["delta_cut"] for p in pts]
            gomp = [p["gompertz_like"] for p in pts]
            ax.scatter(
                [x for x, flag in zip(xs, gomp) if flag],
                [y for y, flag in zip(ys, gomp) if flag],
                marker=mk,
                s=28,
                c="#1F4E79",
                linewidths=0.0,
                label=None,
                zorder=3,
            )
            ax.scatter(
                [x for x, flag in zip(xs, gomp) if not flag],
                [y for y, flag in zip(ys, gomp) if not flag],
                marker=mk,
                s=32,
                facecolors="none",
                edgecolors="#8C3A3A",
                linewidths=1.0,
                label=None,
                zorder=3,
            )
        ax.plot([0.0, float(deltas[-1])], [0.0, float(deltas[-1])], color="#CCCCCC", lw=0.7, zorder=0)
        ax.set_title(f"g = {g:g}")
        ax.set_xlabel("block rate")
        ax.set_xlim(-0.004, float(deltas[-1]) + 0.004)
        ax.set_ylim(-0.004, float(deltas[-1]) + 0.004)
    axes[0].set_ylabel("cut rate")
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#1F4E79", markeredgecolor="none", markersize=6, label="Gompertz-like"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#8C3A3A", markersize=6, label="not Gompertz-like"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#333333", markeredgecolor="none", markersize=6, label="circle: aging-like"),
            Line2D([0], [0], marker="s", color="none", markerfacecolor="#333333", markeredgecolor="none", markersize=6, label="square: cancer-like"),
            Line2D([0], [0], marker="^", color="none", markerfacecolor="#333333", markeredgecolor="none", markersize=6, label="triangle: neither"),
        ],
        frameon=False,
        fontsize=7.5,
        loc="lower center",
        ncol=5,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.24)
    fig.savefig(FIG / "regime_plane.png")
    plt.close(fig)


def round_floats(obj, nd=6):
    if isinstance(obj, float):
        if obj != obj or obj in (float("inf"), float("-inf")):
            return None
        return round(obj, nd)
    if isinstance(obj, dict):
        return {k: round_floats(v, nd) for k, v in obj.items()}
    if isinstance(obj, list):
        return [round_floats(v, nd) for v in obj]
    if isinstance(obj, (np.floating,)):
        return round_floats(float(obj), nd)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def snapshot_bac(series: dict, times: list[float]) -> list[dict]:
    rows = []
    t = series["t"]
    for time in times:
        k = int(np.argmin(np.abs(t - time)))
        rows.append(
            {
                "t": float(t[k]),
                "lambda_min": float(series["lambda_min"][k]),
                "lambda2": float(series["lambda2"][k]),
                "mean_cut": float(series["mean_cut"][k]),
                "mean_block": float(series["mean_block"][k]),
            }
        )
    return rows


def main() -> None:
    w0 = initial_weights()
    op, mean_w = coupling_operator(w0)
    t_fine = np.linspace(0.0, T_END, N_FINE)
    t_year = np.arange(0.0, T_END + 1e-9, 1.0)
    if t_year.size != 81:
        raise SystemExit(f"expected 81 annual nodes, got {t_year.size}")
    df = int(t_year.size - 2)
    crit = float(chi2.ppf(0.95, df))
    sigma = reference_sigma(t_year)

    bac = {kind: bac_series(edge_rates(kind), w0, t_fine) for kind in KINDS}
    # Classifier on the annual nodes must agree: the label uses only the endpoints.
    for kind in KINDS:
        annual = bac_series(edge_rates(kind), w0, t_year)
        if annual["class"] != bac[kind]["class"]:
            raise SystemExit(f"class changed with the grid for {kind}")

    damage_cache: dict[tuple[str, float], np.ndarray] = {}
    for kind in KINDS:
        rates = edge_rates(kind)
        for g in G_SCAN:
            damage_cache[(kind, float(g))] = integrate_damage(rates, float(g), t_year, op)

    scan_rows: list[dict] = []
    log_mu_store: dict[tuple[str, float], np.ndarray] = {}
    for kind in KINDS:
        for g in G_SCAN:
            damage = damage_cache[(kind, float(g))]
            lin = node_linearity(t_year, damage)
            for kind_map in ("exponential", "additive"):
                log_mu = hazard_from_damage(t_year, damage, GAMMA, kind_map)
                if kind_map == "exponential" and float(g) in (0.0, 0.02):
                    log_mu_store[(kind, float(g))] = log_mu
                call = gompertz_call(t_year, log_mu, sigma, crit)
                scan_rows.append(
                    {
                        "kind": kind,
                        "bac_class": bac[kind]["class"],
                        "g": float(g),
                        "gamma": GAMMA,
                        "map": kind_map,
                        "near_linear": lin["near_linear"],
                        "min_r2": lin["min_r2"],
                        "max_rel_resid": lin["max_rel_resid"],
                        "load_0": float(damage[0].mean()),
                        "load_T": float(damage[-1].mean()),
                        **call,
                    }
                )

    gamma_rows: list[dict] = []
    for g in GAMMA_AT_G:
        for kind in KINDS:
            damage = damage_cache[(kind, float(g))]
            for factor in GAMMA_FACTORS:
                gamma = float(GAMMA * factor)
                log_mu = hazard_from_damage(t_year, damage, gamma, "exponential")
                call = gompertz_call(t_year, log_mu, sigma, crit)
                gamma_rows.append(
                    {
                        "kind": kind,
                        "bac_class": bac[kind]["class"],
                        "g": float(g),
                        "gamma": gamma,
                        "factor": float(factor),
                        "map": "exponential",
                        **{k: call[k] for k in ("scored", "rss_linear", "lr_vs_quadratic", "gompertz_like", "curvature_detectable", "beta")},
                    }
                )

    deltas = np.linspace(0.0, DELTA_CUT, N_GRID)
    plane_points: dict[float, list[dict]] = {float(g): [] for g in G_PLANE}
    for g in G_PLANE:
        for delta_block in deltas:
            for delta_cut in deltas:
                rates = plane_rates(float(delta_cut), float(delta_block))
                series = bac_series(rates, w0, np.array([0.0, T_END]))
                damage = integrate_damage(rates, float(g), t_year, op)
                log_mu = hazard_from_damage(t_year, damage, GAMMA, "exponential")
                call = gompertz_call(t_year, log_mu, sigma, crit)
                plane_points[float(g)].append(
                    {
                        "delta_cut": float(delta_cut),
                        "delta_block": float(delta_block),
                        "class": series["class"],
                        "gompertz_like": call["gompertz_like"],
                        "scored": call["scored"],
                        "rss_linear": call["rss_linear"],
                        "sector_change": series["sector_change"],
                        "lambda_min_ratio": series["lambda_min_ratio"],
                    }
                )

    # --- checks that do not depend on a preferred answer ---
    rng = np.random.default_rng(SEED)
    mono_fail = 0
    for _ in range(200):
        i, j = rng.choice(N, size=2, replace=False)
        lifted = w0.copy()
        lifted[i, j] += 0.05
        lifted[j, i] += 0.05
        if grounded_lambda_min(lifted) + 1e-10 < grounded_lambda_min(w0):
            mono_fail += 1

    ray_err = 0.0
    for _ in range(50):
        x = rng.normal(size=N)
        x = x - x.mean()
        x = x / np.linalg.norm(x)
        quad = float(x @ laplacian(w0) @ x)
        acc = 0.0
        a = symmetrize(w0)
        for i in range(N):
            for j in range(i + 1, N):
                acc += a[i, j] * (x[i] - x[j]) ** 2
        ray_err = max(ray_err, abs(quad - acc))

    aging = bac["aging"]
    scale = aging["lambda_min"][0] * np.exp(-DELTA_AGE * aging["t"])
    scale_err = float(np.max(np.abs(aging["lambda_min"] - scale)))
    diag_range = float(np.max(bac["diagonal"]["lambda_min"]) - np.min(bac["diagonal"]["lambda_min"]))
    held_range = float(np.max(bac["held"]["lambda_min"]) - np.min(bac["held"]["lambda_min"]))

    closed = constant_coefficient_damage(0.02, t_year, op)
    ivp_held = damage_cache[("held", 0.02)]
    held_gap = float(np.max(np.abs(closed - ivp_held)))

    g0_ref = damage_cache[("held", 0.0)]
    g0_gap = max(float(np.max(np.abs(damage_cache[(k, 0.0)] - g0_ref))) for k in KINDS)
    linear_gap = float(np.max(np.abs(g0_ref - (D0[None, :] + np.outer(t_year, V)))))

    exp_g0 = [r for r in scan_rows if r["map"] == "exponential" and r["g"] == 0.0 and r["gamma"] == GAMMA]
    rss_g0 = [r["rss_linear"] for r in exp_g0]
    classes_g0 = {r["kind"]: r["bac_class"] for r in exp_g0}

    expected = {
        "held": "neither",
        "aging": "aging-like",
        "cancer": "cancer-like",
        "ambiguous": "neither",
        "diagonal": "neither",
    }
    class_ok = {name: bac[name]["class"] == label for name, label in expected.items()}

    analytic_cross = float(np.log(aging["lambda_min"][0] / SIGMA_HOLD) / DELTA_AGE)
    cross_err = abs(float(aging["crossing_hold"]) - analytic_cross)

    # Diagonal decay does not move off-diagonal W, so it must share the held hazard at every g.
    diag_held_hazard_gap = 0.0
    for g in G_SCAN:
        for kind_map in ("exponential", "additive"):
            d1 = damage_cache[("held", float(g))]
            d2 = damage_cache[("diagonal", float(g))]
            y1 = hazard_from_damage(t_year, d1, GAMMA, kind_map)
            y2 = hazard_from_damage(t_year, d2, GAMMA, kind_map)
            diag_held_hazard_gap = max(diag_held_hazard_gap, float(np.max(np.abs(y1 - y2))))

    add_g0_not = all(
        (not r["gompertz_like"])
        for r in scan_rows
        if r["map"] == "additive" and r["g"] == 0.0
    )
    exp_g0_yes = all(r["gompertz_like"] and r["scored"] for r in exp_g0)

    named_exp = [r for r in scan_rows if r["map"] == "exponential" and r["gamma"] == GAMMA]
    named_coin = coincidence([r["bac_class"] for r in named_exp], [r["gompertz_like"] for r in named_exp])
    plane_coin = {}
    for g, pts in plane_points.items():
        plane_coin[str(g)] = coincidence([p["class"] for p in pts], [p["gompertz_like"] for p in pts])

    # γ moves the residual and does not move the class.
    gamma_classes = {row["kind"]: row["bac_class"] for row in gamma_rows}
    gamma_class_stable = all(gamma_classes[k] == bac[k]["class"] for k in KINDS)

    checks = {
        "monotonicity_grounded_failures": mono_fail,
        "rayleigh_max_abs_err": float(ray_err),
        "aging_scale_max_abs_err": scale_err,
        "diagonal_lambda_min_range": diag_range,
        "held_lambda_min_range": held_range,
        "held_expm_ivp_gap_g002": held_gap,
        "g0_damage_gap_across_paths": g0_gap,
        "g0_linear_gap": linear_gap,
        "classes_match_predeclaration": class_ok,
        "aging_analytic_crossing": analytic_cross,
        "aging_numeric_crossing_abs_err": float(cross_err),
        "diagonal_held_hazard_gap": diag_held_hazard_gap,
        "exponential_g0_all_gompertz": exp_g0_yes,
        "additive_g0_all_not_gompertz": add_g0_not,
        "g0_rss_max": float(max(rss_g0)),
        "named_exponential_sets_equal": named_coin["sets_equal"],
        "plane_sets_equal": {k: v["sets_equal"] for k, v in plane_coin.items()},
        "gamma_class_stable": gamma_class_stable,
        "initial_lambda_min": float(bac["held"]["lambda_min"][0]),
        "initial_lambda2": float(bac["held"]["lambda2"][0]),
        "block_lambda2_0": induced_block_lambda2(w0),
    }

    hard = [
        checks["monotonicity_grounded_failures"] == 0,
        checks["rayleigh_max_abs_err"] < 1e-9,
        checks["aging_scale_max_abs_err"] < 1e-8,
        checks["diagonal_lambda_min_range"] < 1e-10,
        checks["held_lambda_min_range"] < 1e-10,
        checks["held_expm_ivp_gap_g002"] < 1e-6,
        checks["g0_damage_gap_across_paths"] < 1e-10,
        checks["g0_linear_gap"] < 1e-12,
        all(class_ok.values()),
        checks["aging_numeric_crossing_abs_err"] < 1e-3,
        checks["diagonal_held_hazard_gap"] < 1e-10,
        checks["exponential_g0_all_gompertz"] is True,
        checks["additive_g0_all_not_gompertz"] is True,
        checks["g0_rss_max"] < 1e-6,
        checks["named_exponential_sets_equal"] is False,
        checks["plane_sets_equal"].get("0.0") is False,
        checks["gamma_class_stable"] is True,
        abs(checks["initial_lambda_min"] - 0.285102) < 5e-6,
        abs(checks["initial_lambda2"] - 1.413084) < 5e-6,
    ]
    if not all(hard):
        raise SystemExit("toy checks failed:\n" + json.dumps(round_floats(checks), indent=2))

    times = [0.0, 20.0, 40.0, 60.0, 80.0]
    summaries = {}
    for kind in KINDS:
        s = bac[kind]
        summaries[kind] = {
            "class": s["class"],
            "crossing_hold": s["crossing_hold"],
            "lambda_min_ratio": s["lambda_min_ratio"],
            "mean_cut_ratio": s["mean_cut_ratio"],
            "mean_block_ratio": s["mean_block_ratio"],
            "sector_change": s["sector_change"],
            "lambda_min_T": s["lambda_min_T"],
            "lambda2_T": s["lambda2_T"],
        }

    beta_ref = GAMMA * VBAR
    payload = {
        "seed": SEED,
        "labels": LABELS,
        "organism_index": ORG,
        "ground_index": GROUND,
        "depends_on": ["T13", "T18"],
        "shared_horizon": {"t0": 0.0, "t_end": T_END, "n_annual": int(t_year.size), "n_fine_bac": N_FINE},
        "parameters": {
            "delta_age": DELTA_AGE,
            "delta_cut": DELTA_CUT,
            "delta_ambiguous_cut": DELTA_AMB_CUT,
            "delta_ambiguous_block": DELTA_AMB_BLOCK,
            "sigma_hold": SIGMA_HOLD,
            "gamma": GAMMA,
            "mu_b": MU_B,
            "v": V.tolist(),
            "D0": D0.tolist(),
            "vbar": VBAR,
            "L0": L0,
            "beta_linear": beta_ref,
            "doubling_time_linear": float(np.log(2.0) / beta_ref),
            "mean_initial_weight": mean_w,
            "g_scan": G_SCAN.tolist(),
            "g_plane": G_PLANE.tolist(),
            "gamma_factors": GAMMA_FACTORS.tolist(),
            "log_mu_cap": LOG_MU_CAP,
        },
        "thresholds": {
            "lambda_ratio_fall": THR_LAM_FALL,
            "both_sectors_fall": THR_BOTH_FALL,
            "sector_change_band": THR_SECTOR_BAND,
            "cut_collapse": THR_CUT_COLLAPSE,
            "block_held": THR_BLOCK_HELD,
            "sector_rise": THR_SECTOR_RISE,
            "chi2_95_df": crit,
            "df": df,
            "quadratic_lr": CHI2_95_1,
            "near_linear_r2": NEAR_LINEAR_R2,
            "near_linear_max_rel": NEAR_LINEAR_MAX_REL,
        },
        "initial": {
            "lambda_min": checks["initial_lambda_min"],
            "lambda2": checks["initial_lambda2"],
            "lambda2_block": checks["block_lambda2_0"],
            "mean_cut": float(bac["held"]["mean_cut"][0]),
            "mean_block": float(bac["held"]["mean_block"][0]),
        },
        "summaries": summaries,
        "snapshots": {kind: snapshot_bac(bac[kind], times) for kind in KINDS},
        "scan": scan_rows,
        "gamma_scan": gamma_rows,
        "coincidence_named_exponential": named_coin,
        "coincidence_named_by_g": {},
        "coincidence_plane": plane_coin,
        "checks": checks,
        "reference_sigma_ends": {"t0": float(sigma[0]), "tT": float(sigma[-1])},
    }
    for g in G_SCAN:
        rows = [r for r in named_exp if r["g"] == float(g)]
        payload["coincidence_named_by_g"][str(float(g))] = coincidence(
            [r["bac_class"] for r in rows],
            [r["gompertz_like"] for r in rows],
        )

    additive_rows = [r for r in scan_rows if r["map"] == "additive"]
    payload["coincidence_named_additive"] = coincidence(
        [r["bac_class"] for r in additive_rows],
        [r["gompertz_like"] for r in additive_rows],
    )

    draw(
        bac,
        t_year,
        log_mu_store,
        scan_rows,
        {"deltas": deltas, "by_g": plane_points},
        crit,
    )

    # Plane points are numerous. Keep the class and the call; drop nothing the tables need.
    payload["plane"] = {
        "deltas": deltas.tolist(),
        "points": {str(g): pts for g, pts in plane_points.items()},
    }

    out = ROOT / "results.json"
    out.write_text(json.dumps(round_floats(payload), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print(f"chi2_crit df={df} value={crit:.6f} beta={beta_ref:.6f} doubling={np.log(2)/beta_ref:.6f}")
    for kind in KINDS:
        s = summaries[kind]
        print(
            f"{kind:10} {s['class']:12} cross={s['crossing_hold']} "
            f"lmin_ratio={s['lambda_min_ratio']:.6f} sector={s['sector_change']:.4f}"
        )
    print("--- exponential scan (g, class, near, rss, lr, gomp) ---")
    for r in named_exp:
        rss = r["rss_linear"]
        lr = r["lr_vs_quadratic"]
        rss_s = "None" if rss is None else f"{rss:.6g}"
        lr_s = "None" if lr is None else f"{lr:.6g}"
        print(
            f"{r['kind']:10} g={r['g']:<6} {r['bac_class']:12} "
            f"near={str(r['near_linear']):5} rss={rss_s:>12} lr={lr_s:>12} gomp={r['gompertz_like']}"
        )
    print("--- additive g=0 ---")
    for r in scan_rows:
        if r["map"] == "additive" and r["g"] == 0.0:
            print(f"{r['kind']:10} rss={r['rss_linear']:.6g} gomp={r['gompertz_like']} curv={r['curvature_detectable']}")
    print("--- coincidence by g (named, exponential) ---")
    for key, val in payload["coincidence_named_by_g"].items():
        print(key, {k: val[k] for k in ("n_aging", "n_gompertz", "aging_and_gompertz", "aging_not_gompertz", "gompertz_not_aging", "sets_equal", "jaccard")})
    print("--- plane ---")
    for key, val in plane_coin.items():
        print(key, {k: val[k] for k in ("n", "n_aging", "n_gompertz", "aging_and_gompertz", "aging_not_gompertz", "gompertz_not_aging", "sets_equal", "jaccard")})
    print("--- gamma scan where call differs from gamma=8, or all at g=0.005 ---")
    for r in gamma_rows:
        if r["g"] == 0.005 or r["g"] == 0.0:
            print(
                f"{r['kind']:10} g={r['g']:<6} γ={r['gamma']:<6} "
                f"class={r['bac_class']:12} rss={r['rss_linear']} gomp={r['gompertz_like']} curv={r['curvature_detectable']}"
            )


if __name__ == "__main__":
    main()
