"""
Theoretical/computational study: thermally gated autocatalytic transformation
==============================================================================
Generates every numerical "Result" figure for the theory paper. All outputs
are MODEL predictions under stated, literature-informed illustrative
parameters -- not experimental measurements. Every figure/number in this
script is genuinely computed here, not invented after the fact.

Core model:
    dalpha/dt = C(t,T) * (k1(T) + k2(T)*alpha^m) * (1-alpha)^n      [Kamal-Sourour, gated]
    C(t,T)    = H(T - T_melt) * (1 - exp(-t/tau_release))           [thermal gate]
    k_i(T)    = A_i * exp(-Ea_i / (R*T))                            [Arrhenius]

Dimensionless coupling parameter:
    Lambda = tau_release / tau_cure,   tau_cure = 1 / k2(T)
"""
import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")  # Windows consoles default to cp1252 and choke on
    sys.stderr.reconfigure(encoding="utf-8")  # the unicode (alpha, ~, degree...) in print output
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata

plt.rcParams.update({
    "font.size": 12, "axes.titlesize": 12, "axes.labelsize": 12,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5, "legend.fontsize": 10,
    "figure.dpi": 100, "savefig.dpi": 300,
})

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures") + os.sep
os.makedirs(OUT, exist_ok=True)
R_GAS = 8.314  # J/mol/K

# ---------------------------------------------------------------- kinetics
# Prefactors calibrated so k2(413.15 K = 140 degC) ~ 0.15 /min (cure timescale
# of several minutes), a representative order of magnitude for autocatalytic
# epoxy systems near this temperature -- illustrative, not fitted to AutoLatch.
A1, Ea1 = 2.0e3, 5.0e4     # background pathway (1/min prefactor, J/mol)
A2, Ea2 = 2.5e7, 6.5e4     # autocatalytic pathway
M_EXP, N_EXP = 1.0, 1.5

def k1_of_T(T_K):
    return A1 * np.exp(-Ea1 / (R_GAS * T_K))

def k2_of_T(T_K):
    return A2 * np.exp(-Ea2 / (R_GAS * T_K))

def simulate(T_target_C, T_melt_C, tau_release, k2_override=None, autocatalytic=True,
             t_max=150.0, dt=0.02):
    """Integrate the gated cure ODE (scalar or array-broadcastable in T_target_C)."""
    T_K = T_target_C + 273.15
    k1 = k1_of_T(T_K)
    k2 = k2_of_T(T_K) if k2_override is None else k2_override
    gate_open = T_target_C > T_melt_C
    n_steps = int(t_max / dt)
    alpha = np.zeros_like(np.asarray(T_target_C, dtype=float)) if np.ndim(T_target_C) else 0.0
    alpha = np.array(alpha, dtype=float)
    traj = np.zeros((n_steps + 1,) + alpha.shape)
    t = 0.0
    times = [0.0]
    for i in range(1, n_steps + 1):
        t += dt
        C = np.where(gate_open, 1.0 - np.exp(-t / tau_release), 0.0)
        if autocatalytic:
            dalpha = C * (k1 + k2 * np.power(np.clip(alpha, 1e-6, None), M_EXP)) * np.power(1 - alpha, N_EXP)
        else:
            dalpha = C * k2 * (1 - alpha)  # first-order control (m=0 equivalent, same k2 scale)
        alpha = np.clip(alpha + dalpha * dt, 0.0, 0.999)
        traj[i] = alpha
        times.append(t)
    return np.array(times), traj

def crossing_time(times, alpha_traj, level):
    idx = np.argmax(alpha_traj >= level) if np.any(alpha_traj >= level) else -1
    return times[idx] if idx >= 0 else np.nan

# =====================================================================
# RESULT 1: Thermal gating creates an induction regime
# =====================================================================
T_melt = 118.0  # degC
targets = [110, 115, 120, 130, 140, 155]  # degC, spanning below/above T_melt
fig, ax = plt.subplots(figsize=(6.6, 4.2))
for Tt in targets:
    times, traj = simulate(Tt, T_melt, tau_release=6.0)
    label = f"{Tt}\u00b0C" + (" (below gate)" if Tt <= T_melt else "")
    ax.plot(times, traj, label=label, linewidth=2 if Tt > T_melt else 1.4,
            linestyle="-" if Tt > T_melt else "--")
ax.axvspan(25, 75, color="gray", alpha=0.12)
ax.set_xlabel("Time (min)"); ax.set_ylabel("Conversion, \u03B1")
ax.set_title(f"Result 1: Induction regime vs. target temperature (T_melt = {T_melt}\u00b0C)", fontsize=10)
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT + "t1_induction.png", dpi=300, bbox_inches="tight"); plt.close()

# =====================================================================
# RESULT 2: Autocatalysis sharpens the transition (vs. first-order control)
# =====================================================================
T_demo = 140.0
times_a, traj_a = simulate(T_demo, T_melt, tau_release=6.0, autocatalytic=True)
# choose first-order k so that t50 roughly matches the autocatalytic case
k2_match = k2_of_T(T_demo + 273.15)
times_b, traj_b = simulate(T_demo, T_melt, tau_release=6.0, autocatalytic=False, k2_override=k2_match * 0.55)

def sharpness(times, traj):
    t10 = crossing_time(times, traj, 0.10)
    t90 = crossing_time(times, traj, 0.90)
    return (t90 - t10) / t10, t10, t90

S_auto, t10_a, t90_a = sharpness(times_a, traj_a)
S_first, t10_b, t90_b = sharpness(times_b, traj_b)

fig, ax = plt.subplots(figsize=(6.6, 4.0))
ax.plot(times_a, traj_a, label=f"Autocatalytic (S = {S_auto:.2f})", color="#2E5597", linewidth=2.2)
ax.plot(times_b, traj_b, label=f"First-order control (S = {S_first:.2f})", color="#B5651D", linewidth=2.2, linestyle="--")
ax.axvspan(25, 75, color="gray", alpha=0.12, label="Target window")
ax.set_xlabel("Time (min)"); ax.set_ylabel("Conversion, \u03B1")
ax.set_title(f"Result 2: Transition sharpness, S=(t90\u2212t10)/t10, at {T_demo}\u00b0C", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT + "t2_sharpness.png", dpi=300, bbox_inches="tight"); plt.close()

print(f"Result 2: S_autocatalytic={S_auto:.3f}, S_first_order={S_first:.3f}  (lower = sharper)")

# =====================================================================
# RESULT 3: Release timescale controls transformation time
# =====================================================================
tau_values = np.linspace(1, 25, 25)
k2_scenarios = {"low k2": k2_of_T(T_demo+273.15)*0.5, "mid k2": k2_of_T(T_demo+273.15),
                "high k2": k2_of_T(T_demo+273.15)*1.8}
fig, ax = plt.subplots(figsize=(6.6, 4.0))
colors = ["#B5651D", "#2E5597", "#3D7A46"]
for (label, k2v), c in zip(k2_scenarios.items(), colors):
    t90s = []
    for tau in tau_values:
        times, traj = simulate(T_demo, T_melt, tau_release=tau, k2_override=k2v)
        t90s.append(crossing_time(times, traj, 0.90))
    ax.plot(tau_values, t90s, label=label, color=c, linewidth=2)
ax.axhspan(25, 75, color="gray", alpha=0.12, label="Target window")
ax.set_xlabel("Release time constant, \u03C4_release (min)")
ax.set_ylabel("Time to 90% conversion (min)")
ax.set_title("Result 3: Release timescale sets transformation time", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT + "t3_release_sweep.png", dpi=300, bbox_inches="tight"); plt.close()

# =====================================================================
# RESULT 4: Dimensionless parameter map (Lambda vs. k2 strength)
# =====================================================================
tau_grid = np.linspace(1, 25, 40)
k2_grid = np.linspace(0.3, 2.2, 40) * k2_of_T(T_demo + 273.15)
Lambda_grid = tau_grid / (1.0 / k2_of_T(T_demo + 273.15))  # tau_release / tau_cure(reference)
T90 = np.zeros((len(k2_grid), len(tau_grid)))
for i, k2v in enumerate(k2_grid):
    for j, tau in enumerate(tau_grid):
        times, traj = simulate(T_demo, T_melt, tau_release=tau, k2_override=k2v)
        T90[i, j] = crossing_time(times, traj, 0.90)

fig, ax = plt.subplots(figsize=(6.8, 4.4))
im = ax.pcolormesh(Lambda_grid, k2_grid / k2_of_T(T_demo+273.15), T90, shading="auto", cmap="viridis")
cbar = plt.colorbar(im, ax=ax); cbar.set_label("Time to 90% conversion (min)")
cs = ax.contour(Lambda_grid, k2_grid / k2_of_T(T_demo+273.15), T90, levels=[75], colors="white", linewidths=2)
ax.clabel(cs, fmt="75 min")
ax.set_xlabel("\u039B = \u03C4_release / \u03C4_cure")
ax.set_ylabel("Normalized autocatalytic strength (k2 / k2,ref)")
ax.set_title("Result 4: Parameter-space map for target transformation time", fontsize=10)
plt.tight_layout(); plt.savefig(OUT + "t4_lambda_map.png", dpi=300, bbox_inches="tight"); plt.close()

# =====================================================================
# RESULT 5: Three-system comparison (A: no gate, B: gate+first-order, C: gate+autocatalytic)
# =====================================================================
times_A, traj_A = simulate(T_demo, -50.0, tau_release=0.01, autocatalytic=True)  # gate effectively always open (no gating mechanism)
times_B, traj_B = times_b, traj_b   # gate + first-order (from Result 2)
times_C, traj_C = times_a, traj_a   # gate + autocatalytic (from Result 2)

def conv_at(times, traj, t_query):
    idx = np.searchsorted(times, t_query)
    return traj[min(idx, len(traj)-1)]

alpha25 = {"A: no gate, autocatalytic": conv_at(times_A, traj_A, 25),
           "B: gate, first-order": conv_at(times_B, traj_B, 25),
           "C: gate, autocatalytic": conv_at(times_C, traj_C, 25)}
S_vals = {"A: no gate, autocatalytic": sharpness(times_A, traj_A)[0],
          "B: gate, first-order": S_first,
          "C: gate, autocatalytic": S_auto}

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
axes[0].plot(times_A, traj_A, label="A: no gate, autocatalytic", color="#9c4b1e")
axes[0].plot(times_B, traj_B, label="B: gate, first-order", color="#B5651D", linestyle="--")
axes[0].plot(times_C, traj_C, label="C: gate, autocatalytic", color="#2E5597")
axes[0].axvspan(25, 75, color="gray", alpha=0.12)
axes[0].set_xlabel("Time (min)"); axes[0].set_ylabel("Conversion, \u03B1")
axes[0].legend(fontsize=7.5); axes[0].grid(alpha=0.25)
axes[0].set_title("Trajectories", fontsize=10)

labels = list(alpha25.keys())
x = np.arange(len(labels))
axes[1].bar(x - 0.18, [alpha25[l] for l in labels], width=0.35, label="\u03B1 at t=25 min (want \u2192 0)", color="#B5651D")
axes[1].bar(x + 0.18, [S_vals[l] for l in labels], width=0.35, label="Sharpness S (want \u2192 low)", color="#2E5597")
axes[1].set_xticks(x); axes[1].set_xticklabels(["A", "B", "C"])
axes[1].legend(fontsize=8); axes[1].set_title("Premature conversion vs. sharpness", fontsize=10)
axes[1].grid(alpha=0.25, axis="y")
plt.tight_layout(); plt.savefig(OUT + "t5_system_comparison.png", dpi=300, bbox_inches="tight"); plt.close()

print("Result 5 (alpha at t=25 min):", {k: round(v, 3) for k, v in alpha25.items()})
print("Result 5 (sharpness S):", {k: round(v, 3) for k, v in S_vals.items()})

# System D: no gate, first-order (conventional baseline control)
times_D, traj_D = simulate(T_demo, -50.0, tau_release=0.01, autocatalytic=False, k2_override=k2_match * 0.55)
alpha25["D: no gate, first-order"] = conv_at(times_D, traj_D, 25)
S_vals["D: no gate, first-order"] = sharpness(times_D, traj_D)[0]

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
axes[0].plot(times_A, traj_A, label="A: no gate, autocatalytic", color="#9c4b1e", linewidth=2)
axes[0].plot(times_B, traj_B, label="B: gate, first-order", color="#B5651D", linestyle="--", linewidth=2)
axes[0].plot(times_C, traj_C, label="C: gate, autocatalytic", color="#2E5597", linewidth=2)
axes[0].plot(times_D, traj_D, label="D: no gate, first-order", color="#555555", linestyle=":", linewidth=2.2)
axes[0].axvspan(25, 75, color="gray", alpha=0.12)
axes[0].set_xlabel("Time (min)"); axes[0].set_ylabel("Conversion, \u03B1")
axes[0].legend(fontsize=10, loc="lower right"); axes[0].grid(alpha=0.25)
axes[0].set_title("Trajectories (2\u00d72 control matrix)", fontsize=12)

labels = ["A", "B", "C", "D"]
keys = ["A: no gate, autocatalytic", "B: gate, first-order", "C: gate, autocatalytic", "D: no gate, first-order"]
x = np.arange(len(labels))
axes[1].bar(x - 0.18, [alpha25[k] for k in keys], width=0.35, label="\u03B1 at t=25 min (want \u2192 0)", color="#B5651D")
axes[1].bar(x + 0.18, [S_vals[k] for k in keys], width=0.35, label="Sharpness S (want \u2192 low)", color="#2E5597")
axes[1].set_xticks(x); axes[1].set_xticklabels(labels, fontsize=12)
axes[1].legend(fontsize=10); axes[1].set_title("Premature conversion vs. sharpness", fontsize=12)
axes[1].grid(alpha=0.25, axis="y")
plt.tight_layout(); plt.savefig(OUT + "t5_system_comparison.png", dpi=300, bbox_inches="tight"); plt.close()

print("Full 2x2 comparison:")
for k in keys:
    print(f"  {k:30s} alpha(25min)={alpha25[k]:.3f}  S={S_vals[k]:.3f}")

# =====================================================================
# RESULT 7 (NEW): Sharpness vs release timescale -- the rate-limiting regime
# Tests the explicit falsifiable prediction: "increasing the catalyst-release
# timescale beyond a defined regime becomes rate limiting and eliminates the
# sharp transformation window."
# =====================================================================
tau_wide = np.linspace(0.5, 80, 60)
S_of_tau = []
t90_of_tau = []
k2_ref = k2_of_T(T_demo + 273.15)
for tau in tau_wide:
    times, traj = simulate(T_demo, T_melt, tau_release=tau, k2_override=k2_ref)
    S_val, _, t90v = sharpness(times, traj)
    S_of_tau.append(S_val)
    t90_of_tau.append(t90v)
S_of_tau = np.array(S_of_tau)
Lambda_wide = tau_wide / (1.0 / k2_ref)

fig, ax1 = plt.subplots(figsize=(7.6, 4.8))
ax1.plot(Lambda_wide, S_of_tau, color="#2E5597", linewidth=2.4, label="Sharpness S")
ax1.set_xlabel("\u039B = \u03C4_release / \u03C4_cure")
ax1.set_ylabel("Sharpness S = (t\u2089\u2080\u2212t\u2081\u2080)/t\u2081\u2080", color="#2E5597", fontsize=11.5)
ax1.tick_params(axis="y", labelcolor="#2E5597")
ax2 = ax1.twinx()
ax2.plot(Lambda_wide, t90_of_tau, color="#B5651D", linewidth=2.2, linestyle="--", label="t\u2089\u2080")
ax2.axhline(75, color="gray", linestyle=":", linewidth=1.4)
ax2.set_ylabel("Time to 90% conversion (min)", color="#B5651D", fontsize=11.5)
ax2.tick_params(axis="y", labelcolor="#B5651D")
ax1.set_title("Result 7: S and t\u2089\u2080 diverge \u2014 absolute timing, not S alone,\ndetermines when release becomes rate-limiting", fontsize=11.5)
fig.subplots_adjust(left=0.14, right=0.86, top=0.85, bottom=0.13)
plt.savefig(OUT + "t7_sharpness_vs_lambda.png", dpi=300, bbox_inches="tight"); plt.close()

# report the crossover point where S starts increasing (rate-limiting onset)
dS = np.diff(S_of_tau)
turn_idx = np.argmin(S_of_tau)  # sharpest point
print(f"\nResult 7: sharpness minimum (S={S_of_tau[turn_idx]:.2f}) occurs at Lambda={Lambda_wide[turn_idx]:.2f}; "
      f"S rises to {S_of_tau[-1]:.2f} by Lambda={Lambda_wide[-1]:.2f} (release becomes rate-limiting)")
crossing_75 = Lambda_wide[np.argmax(np.array(t90_of_tau) > 75)] if np.any(np.array(t90_of_tau) > 75) else None
print(f"Result 7: t90 first exceeds 75 min at Lambda \u2248 {crossing_75:.2f}")

# =====================================================================
# RESULT 8 (NEW): Realistic (non-ideal) thermal gate -- finite transition width
# Replaces the Heaviside step with a logistic gate of width w (degC), testing
# robustness to a non-ideal (real material) melting/activation range.
# =====================================================================
def simulate_smooth_gate(T_target_C, T_melt_C, width_C, tau_release, k2v, t_max=150.0, dt=0.02):
    n_steps = int(t_max / dt)
    alpha = 0.0
    T_K = T_target_C + 273.15
    k1 = k1_of_T(T_K)
    G_max = 1.0 / (1.0 + np.exp(-(T_target_C - T_melt_C) / max(width_C, 1e-6))) if width_C > 0 else float(T_target_C > T_melt_C)
    times = [0.0]; traj = [0.0]; t = 0.0
    for _ in range(n_steps):
        t += dt
        C = G_max * (1.0 - np.exp(-t / tau_release))
        dalpha = C * (k1 + k2v * max(alpha, 1e-6) ** M_EXP) * (1 - alpha) ** N_EXP
        alpha = min(max(alpha + dalpha * dt, 0.0), 0.999)
        times.append(t); traj.append(alpha)
    return np.array(times), np.array(traj)

widths = [0, 2, 5, 10, 15]
fig, ax = plt.subplots(figsize=(6.6, 4.2))
alpha25_widths = []
for w in widths:
    times, traj = simulate_smooth_gate(T_demo, T_melt, w, tau_release=6.0, k2v=k2_ref)
    ax.plot(times, traj, label=f"gate width = {w}\u00b0C", linewidth=2)
    idx25 = np.searchsorted(times, 25)
    alpha25_widths.append(traj[min(idx25, len(traj)-1)])
ax.axvspan(25, 75, color="gray", alpha=0.12)
ax.set_xlabel("Time (min)"); ax.set_ylabel("Conversion, \u03B1")
ax.set_title(f"Result 8: Effect of a realistic (non-ideal) thermal gate\n(T_target\u2212T_melt = {T_demo-T_melt:.0f}\u00b0C)", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT + "t8_gate_width.png", dpi=300, bbox_inches="tight"); plt.close()

print("\nResult 8: alpha at t=25 min vs. gate width:", {w: round(a, 3) for w, a in zip(widths, alpha25_widths)})

t90_widths = []
for w in widths:
    times, traj = simulate_smooth_gate(T_demo, T_melt, w, tau_release=6.0, k2v=k2_ref)
    t90_widths.append(crossing_time(times, traj, 0.90))
print("Result 8: t90 vs. gate width:", {w: (round(t, 1) if not np.isnan(t) else None) for w, t in zip(widths, t90_widths)})

# =====================================================================
# RESULT 9: Robustness under parameter uncertainty (Monte Carlo)
# Propagates uncertainty in kinetic/release parameters (literature-informed
# ranges, not fitted to any measured system) to test how sensitive the
# predicted transformation window is, and which parameter dominates.
# =====================================================================
N_MC = 20000
rng = np.random.default_rng(20260822)
k1_mc = rng.uniform(0.0005, 0.0025, N_MC)
k2max_mc = rng.uniform(0.04, 0.22, N_MC)
m_mc = rng.uniform(0.5, 1.5, N_MC)
n_mc = rng.uniform(1.0, 2.0, N_MC)
melt_offset_mc = rng.normal(12.0, 4.0, N_MC)
tau_release_mc = rng.uniform(2.0, 15.0, N_MC)
T_target_mc = T_demo
T_melt_mc = T_target_mc - melt_offset_mc

dt_mc = 0.05; t_max_mc = 150.0
alpha_mc = np.zeros(N_MC)
gate_open_mc = T_melt_mc < T_target_mc
t10_mc = np.full(N_MC, np.nan); t90_mc = np.full(N_MC, np.nan)
t_cur = 0.0
for _ in range(int(t_max_mc / dt_mc)):
    t_cur += dt_mc
    Cf = np.where(gate_open_mc, 1.0 - np.exp(-t_cur / tau_release_mc), 0.0)
    k2eff = k2max_mc * Cf
    da = (k1_mc + k2eff * np.power(np.clip(alpha_mc, 1e-6, None), m_mc)) * np.power(1 - alpha_mc, n_mc)
    alpha_mc = np.clip(alpha_mc + da * dt_mc, 0.0, 0.999)
    newly10 = np.isnan(t10_mc) & (alpha_mc >= 0.10); t10_mc[newly10] = t_cur
    newly90 = np.isnan(t90_mc) & (alpha_mc >= 0.90); t90_mc[newly90] = t_cur

valid_mc = ~np.isnan(t90_mc)
p_in_window = np.mean(valid_mc & (t90_mc <= 75) & (t10_mc >= 5))
print(f"\nResult 9 (Monte Carlo, N={N_MC}): P(t90<=75 and t10>=5) = {p_in_window:.1%}")
print(f"Result 9: median t90 = {np.nanmedian(t90_mc):.1f} min (IQR {np.nanpercentile(t90_mc[valid_mc],25):.1f}-{np.nanpercentile(t90_mc[valid_mc],75):.1f})")

params_mc = {"k1": k1_mc, "k2_max": k2max_mc, "m": m_mc, "n": n_mc,
             "T_melt offset": melt_offset_mc, "tau_release": tau_release_mc}
def spearman(x, y):
    mask = ~np.isnan(y)
    return np.corrcoef(rankdata(x[mask]), rankdata(y[mask]))[0, 1]
sens_mc = {k: spearman(v, t90_mc) for k, v in params_mc.items()}
print("Result 9 sensitivity (Spearman r with t90):", {k: round(v, 2) for k, v in sens_mc.items()})

fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
axes[0].hist(t90_mc[valid_mc], bins=60, color="#2E5597", alpha=0.8)
axes[0].axvspan(0, 75, color="gray", alpha=0.12, label="Meets 75 min target")
axes[0].set_xlabel("Time to 90% conversion (min)"); axes[0].set_ylabel("Monte Carlo draws")
axes[0].set_title(f"Distribution under parameter uncertainty\n(N={N_MC}, illustrative ranges)", fontsize=11.5)
axes[0].legend(fontsize=9)
names_sorted = sorted(sens_mc.keys(), key=lambda k: abs(sens_mc[k]))
vals_sorted = [sens_mc[k] for k in names_sorted]
colors = ["#B5651D" if v < 0 else "#2E5597" for v in vals_sorted]
axes[1].barh(names_sorted, vals_sorted, color=colors)
axes[1].axvline(0, color="black", linewidth=0.8)
axes[1].set_xlabel("Spearman rank correlation with t\u2089\u2080")
axes[1].set_title("Sensitivity ranking", fontsize=11.5)
plt.tight_layout(); plt.savefig(OUT + "t9_mc_robustness.png", dpi=300, bbox_inches="tight"); plt.close()

# =====================================================================
# RESULT 6 (micromechanics, illustrative): strength vs conversion & filler fraction
# Ryshkewitch-Duckworth-type porosity-strength scaling (illustrative constants)
# =====================================================================
alpha_range = np.linspace(0.5, 1.0, 60)
phi_range = np.linspace(0.10, 0.40, 60)
A_grid, P_grid = np.meshgrid(alpha_range, phi_range)
sigma0_fully_cured = 9000.0   # psi, illustrative fully-cured unfilled matrix strength
b_porosity = 5.0               # illustrative Ryshkewitch exponent
k_filler = 2.0                 # illustrative filler-reinforcement factor
porosity = 0.35 * (1 - A_grid)  # illustrative: porosity closes as conversion completes
strength = sigma0_fully_cured * np.exp(-b_porosity * porosity) * (1 + k_filler * P_grid)

fig, ax = plt.subplots(figsize=(6.8, 4.4))
im = ax.pcolormesh(A_grid, P_grid, strength, shading="auto", cmap="magma")
cbar = plt.colorbar(im, ax=ax); cbar.set_label("Predicted crush strength (psi, illustrative)")
cs = ax.contour(A_grid, P_grid, strength, levels=[4500, 6000], colors="white", linewidths=1.8)
ax.clabel(cs, fmt="%d psi")
ax.set_xlabel("Cure conversion, \u03B1"); ax.set_ylabel("Filler volume fraction, \u03C6")
ax.set_title("Result 6: Illustrative strength model \u2014 NOT calibrated to measured data", fontsize=9.5)
plt.tight_layout(); plt.savefig(OUT + "t6_micromechanics.png", dpi=300, bbox_inches="tight"); plt.close()

print("\nAll figures written to", OUT)
