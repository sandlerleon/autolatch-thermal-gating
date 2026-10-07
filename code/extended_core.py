"""Core routines for the extended analyses: vectorised RK4 of the coupled T-C-alpha system with a capsule population.
Imported by extended_analyses.py and by test_extended.py (no analysis is run on import)."""
import json
import os
import sys
import time

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.stats import norm, qmc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures") + os.sep
R_GAS = 8.314
A1, EA1 = 2.0e3, 5.0e4
A2, EA2 = 2.5e7, 6.5e4
M_EXP, N_EXP = 1.0, 1.5
A_REL, EA_REL = 6.4e6, 6.0e4
T_MELT, W_DEF = 118.0, 3.0
TAU_TH, DT_AD = 3.0, 220.0
T_REF = 140.0
T_STORE = T_MELT - 6.0         # storage temperature used for the leakage objective (as in v2 Section 4.6)
K2_REF = A2 * np.exp(-EA2 / (R_GAS * (T_REF + 273.15)))
TAU_REL_REF = 1.0 / (A_REL * np.exp(-EA_REL / (R_GAS * (T_REF + 273.15))))
plt.rcParams.update({"font.size": 11, "axes.titlesize": 11, "axes.labelsize": 11, "savefig.dpi": 300})


def col(x, S):
    x = np.asarray(x, dtype=float)
    return np.broadcast_to(x, (S,)).copy() if x.ndim == 0 else x


def gate(T_C, Tm, w):
    z = np.clip(-(T_C - Tm) / w, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def integrate(Tinf, tau_th=TAU_TH, rel_scale=1.0, Tmelt=T_MELT, width=W_DEF, off=(0.0,), mode="auto", gated=True,
              p=None, dt=0.05, tmax=150.0, T0_off=30.0, t_query=25.0):
    """Vectorised fixed-step RK4 of the coupled T-C-alpha system for S scenarios at once.
    C has one component per capsule (K offsets of the melt temperature); the reaction is driven by the capsule-mean C.
    Returns dict of per-scenario arrays: t10, t50, t90 (linear interpolation), a_q (alpha at t_query), a_end, Tmax."""
    p = p or {}
    Tinf = np.atleast_1d(np.asarray(Tinf, dtype=float))
    S = len(Tinf)
    off = np.asarray(off, dtype=float)
    K = len(off)
    tau_th, Tmelt, width = col(tau_th, S), col(Tmelt, S), col(width, S)
    rel_scale = np.asarray(rel_scale, dtype=float)
    rel_scale = np.broadcast_to(rel_scale if rel_scale.ndim == 2 else (rel_scale.reshape(-1, 1) if rel_scale.ndim == 1 else rel_scale.reshape(1, 1)), (S, K)).copy()
    ea1, ea2, a2, m, n = (col(p.get(k, v), S) for k, v in (("ea1", EA1), ("ea2", EA2), ("a2", A2), ("m", M_EXP), ("n", N_EXP)))
    dtad = col(p.get("dtad", DT_AD), S)
    kfirst, nord, eac = col(p.get("kfirst", K2_REF), S), col(p.get("nord", 1.0), S), col(p.get("eac", EA2), S)
    T = Tinf - T0_off
    C = np.zeros((S, K)) if gated else np.ones((S, K))
    a = np.zeros(S)
    Tmax = T.copy()
    t10 = np.full(S, np.nan); t50 = t10.copy(); t90 = t10.copy()
    a_q = np.full(S, np.nan)
    nsteps = int(round(tmax / dt))
    iq = int(round(t_query / dt))
    Tref_K = T_REF + 273.15

    def f(T, C, a):
        TK = np.clip(T + 273.15, 50.0, None)
        ac = np.clip(a, 0.0, 0.999)
        if gated:
            krel = (A_REL * np.exp(-EA_REL / (R_GAS * TK)))[:, None] / rel_scale
            dC = krel * (gate(T[:, None], Tmelt[:, None] + off[None, :], width[:, None]) - C)
        else:
            dC = np.zeros_like(C)
        Ce = C.mean(axis=1)
        if mode == "auto":
            da = Ce * (A1 * np.exp(-ea1 / (R_GAS * TK)) + a2 * np.exp(-ea2 / (R_GAS * TK)) * np.maximum(ac, 1e-6) ** m) * (1 - ac) ** n
        else:
            da = Ce * kfirst * np.exp(-eac / R_GAS * (1.0 / TK - 1.0 / Tref_K)) * (1 - ac) ** nord
        da = np.where(a >= 0.999, 0.0, da)       # reaction complete: no further heat (conserves energy; the v2 clip left a small residual source)
        return (Tinf - T) / tau_th + dtad * da, dC, da

    for i in range(1, nsteps + 1):
        k1 = f(T, C, a)
        k2 = f(T + dt / 2 * k1[0], C + dt / 2 * k1[1], a + dt / 2 * k1[2])
        k3 = f(T + dt / 2 * k2[0], C + dt / 2 * k2[1], a + dt / 2 * k2[2])
        k4 = f(T + dt * k3[0], C + dt * k3[1], a + dt * k3[2])
        a_prev = a
        T = T + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        C = np.clip(C + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]), 0.0, 1.0)
        a = np.clip(a + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]), 0.0, 0.999)
        Tmax = np.maximum(Tmax, T)
        for arr, thr in ((t10, 0.10), (t50, 0.50), (t90, 0.90)):
            new = np.isnan(arr) & (a >= thr)
            if new.any():
                arr[new] = (i - 1) * dt + dt * (thr - a_prev[new]) / np.maximum(a[new] - a_prev[new], 1e-30)
        if i == iq:
            a_q = a.copy()
    return dict(t10=t10, t50=t50, t90=t90, a_q=a_q, a_end=a.copy(), Tmax=Tmax, overshoot=Tmax - Tinf)


def sharp(r):
    return (r["t90"] - r["t10"]) / r["t10"]


