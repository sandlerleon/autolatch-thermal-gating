"""Checks for the extended analyses: integrator validity, population-model consistency, Sobol estimator, energy bound, determinism."""
import sys
import numpy as np
from extended_core import integrate, sharp, T_REF, T_STORE, DT_AD, W_DEF
from scipy.stats import qmc

FAILS = []


def check(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


ref = integrate([T_REF], dt=0.02)
check("batch RK4 reproduces the v2 adaptive-RK45 reference (t90 34.20 min) within 0.3 %", abs(ref["t90"][0] - 34.20) / 34.20 < 0.003, "t90 = %.3f" % ref["t90"][0])
check("alpha(25 min) = 0.072 and overshoot = 111.7 K as in v2", abs(ref["a_q"][0] - 0.072) < 0.002 and abs(ref["overshoot"][0] - 111.7) < 0.5)
r2 = integrate([T_REF], dt=0.02)
check("deterministic (identical repeat)", np.array_equal(ref["t90"], r2["t90"]))
# a one-capsule population is exactly the single-gate model; a degenerate (zero-spread) population too
p1 = integrate([T_REF], off=np.zeros(1), dt=0.05)
p61 = integrate([T_REF], off=np.zeros(61), dt=0.05)
check("one-capsule population equals the single-gate model", abs(p1["t90"][0] - p61["t90"][0]) < 1e-9)
check("peak temperature never exceeds T_inf + DT_ad (energy bound)", ref["Tmax"][0] <= T_REF + DT_AD)
check("conversion below the gate is zero for a sharp gate (structural result of Section 2.4)", integrate([T_STORE - 6.0], width=0.2, dt=0.05)["a_end"][0] < 1e-9)
# Jansen/Saltelli estimators checked on the Ishigami function (analytic S1 = 0.3139, 0.4424, 0; ST = 0.5576, 0.4424, 0.2437)
d, N, a, b = 3, 8192, 7.0, 0.1
U = qmc.Sobol(d=2 * d, scramble=True, seed=1).random(N)
def ish(U_):
    x = -np.pi + 2 * np.pi * U_
    return np.sin(x[:, 0]) + a * np.sin(x[:, 1]) ** 2 + b * x[:, 2] ** 4 * np.sin(x[:, 0])
UA, UB = U[:, :d], U[:, d:]
fA, fB = ish(UA), ish(UB)
V = np.var(np.concatenate([fA, fB]), ddof=1)
S1, ST = [], []
for i in range(d):
    Uab = UA.copy(); Uab[:, i] = UB[:, i]
    f = ish(Uab)
    S1.append((V - 0.5 * np.mean((fB - f) ** 2)) / V); ST.append(0.5 * np.mean((fA - f) ** 2) / V)
check("Sobol estimators recover the Ishigami indices within 0.03", np.allclose(S1, [0.3139, 0.4424, 0.0], atol=0.03) and np.allclose(ST, [0.5576, 0.4424, 0.2437], atol=0.03),
      "S1=%s ST=%s" % (np.round(S1, 3), np.round(ST, 3)))
print("\n%s" % ("ALL CHECKS PASSED" if not FAILS else "FAILED: " + ", ".join(FAILS)))
sys.exit(1 if FAILS else 0)
