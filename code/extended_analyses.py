"""
Extended analyses (v3) for the coupled thermal-activation / autocatalytic-cure model
=====================================================================================
Adds four analyses to the v2 study WITHOUT changing any v2 number (theoretical_study.py is untouched):

  E1  Capsule-population model: the gate is not postulated but emerges from a distribution of capsule melt
      temperatures (and, separately, release timescales). Gives a testable relation between the DSC melting-endotherm
      width and the effective gate width of the single-gate model.
  E2  Storage-stability axis and control-family robustness: the claim "transport stability vs transformation speed" needs a
      storage-stability objective that v2 did not measure. Three objectives (premature conversion at 25 min and t90 at the
      operating temperature; conversion after 150 min at a storage temperature below T_melt) are compared against several tuned
      control families. Reports whether ANY control dominates the coupled system and which family comes closest.
  E3  Design-viability map over (gate width, operating temperature): where the storage, induction and completion criteria hold together.
  E4  Variance-based (Sobol) attribution of t90, next to the Spearman ranking of v2.

All numbers are model outputs under the same illustrative parameters as v2; none is experimental.

    python extended_analyses.py   ->  ../results_v3.json  and  ../figures/t10..t13*.png
"""
from extended_core import *  # noqa: F401,F403
from extended_core import integrate, sharp, col, gate, np, plt, norm, qmc, brentq, OUT, HERE
import extended_core as ec
RES = {}
# ------------------------------------------------------------------ validation against the v2 adaptive-RK45 reference
t0 = time.time()
ref = integrate([T_REF], dt=0.02)
print("=== validation (reference case, T_inf = 140 C, dt = 0.02 min, w = 3 C) ===")
print("  batch RK4:  t10 = %.3f  t90 = %.3f  S = %.4f  alpha25 = %.4f  overshoot = %.1f" % (ref["t10"][0], ref["t90"][0], sharp(ref)[0], ref["a_q"][0], ref["overshoot"][0]))
print("  v2 reported (adaptive RK45): t10 = 26.75, t90 = 34.20, S = 0.279, alpha25 = 0.072, overshoot = 111.7")
RES["validation"] = dict(t10=float(ref["t10"][0]), t90=float(ref["t90"][0]), S=float(sharp(ref)[0]), alpha25=float(ref["a_q"][0]),
                         overshoot=float(ref["overshoot"][0]), v2=dict(t10=26.75, t90=34.20, S=0.279, alpha25=0.072, overshoot=111.7))
COUP = dict(a25=float(ref["a_q"][0]), t90=float(ref["t90"][0]))
st = integrate([T_STORE], dt=0.05)
COUP["a_store"] = float(st["a_end"][0])
print("  coupled system storage conversion (T_inf = %.0f C, 150 min): %.5f" % (T_STORE, COUP["a_store"]))
RES["coupled_point"] = COUP


# ungated autocatalytic baseline (same kinetics, no gate): what a formulation without encapsulation does at storage
ung_op = integrate([T_REF], gated=False, dt=0.02)
ung_st = integrate([T_STORE], gated=False, dt=0.05)
RES["ungated_baseline"] = dict(a25=float(ung_op["a_q"][0]), t90=float(ung_op["t90"][0]), S=float(sharp(ung_op)[0]), a_store=float(ung_st["a_end"][0]),
                               t50_store=float(ung_st["t50"][0]) if not np.isnan(ung_st["t50"][0]) else None)
print("  ungated autocatalytic baseline: alpha25 = %.4f  t90 = %.2f  storage conversion (150 min) = %.4f" % (RES["ungated_baseline"]["a25"], RES["ungated_baseline"]["t90"], RES["ungated_baseline"]["a_store"]))

# ================================================================== E1 capsule population
print("\n=== E1 capsule-population model ===")
K = 61
z = norm.ppf((np.arange(K) + 0.5) / K)                      # deterministic Gaussian quantiles of the melt temperature
W_CAP = 0.5                                                  # each capsule has a near-step gate
sigmas = [0.0, 1.0, 2.0, 4.0, 6.0]
pop = []
for s in sigmas:
    r_op = integrate([T_REF], Tmelt=T_MELT, width=W_CAP, off=s * z, dt=0.02)
    r_st = integrate([T_STORE], Tmelt=T_MELT, width=W_CAP, off=s * z, dt=0.05)
    pop.append(dict(sigma=s, a25=float(r_op["a_q"][0]), t90=float(r_op["t90"][0]), S=float(sharp(r_op)[0]), a_store=float(r_st["a_end"][0]),
                    overshoot=float(r_op["overshoot"][0])))
    print("  sigma_Tm = %.1f C: alpha25 = %.4f  t90 = %.2f  S = %.3f  storage = %.5f  overshoot = %.1f" % (s, pop[-1]["a25"], pop[-1]["t90"], pop[-1]["S"], pop[-1]["a_store"], pop[-1]["overshoot"]))


def single_gate_store(w):
    return float(integrate([T_STORE], width=w, dt=0.05)["a_end"][0])


def single_gate_a25(w):
    return float(integrate([T_REF], width=w, dt=0.05)["a_q"][0])


weff = []
for row in pop:
    s = row["sigma"]
    analytic = float(np.sqrt(W_CAP ** 2 + 3.0 * s ** 2 / np.pi ** 2))             # logistic-variance matching: pi^2 w^2 / 3 = sigma^2 + pi^2 W^2 / 3
    target = row["a_store"]
    lo, hi = 0.3, 30.0
    if single_gate_store(lo) > target:
        w_match = lo
    else:
        w_match = brentq(lambda w: single_gate_store(w) - target, lo, hi, xtol=1e-3)
    weff.append(dict(sigma=s, w_analytic=analytic, w_matched_storage=float(w_match)))
    print("  sigma = %.1f: analytic w_eff = %.2f C, w matched on storage leakage = %.2f C" % (s, analytic, w_match))
RES["population_sigma_Tm"] = dict(K=K, w_capsule=W_CAP, rows=pop, equivalent_width=weff)

# release-time dispersion (shell thickness): lognormal scatter of the release time constant, sharp gate
cvs = [0.0, 0.25, 0.5, 1.0]
tau_pop = []
zt = norm.ppf((np.arange(K) + 0.5) / K)
for cv in cvs:
    sg = np.sqrt(np.log(1 + cv ** 2))
    rs = np.exp(sg * zt - 0.5 * sg ** 2)[None, :]       # unit-mean lognormal multiplier of the release time
    r_op = integrate([T_REF], width=W_CAP, rel_scale=rs, off=np.zeros(K), dt=0.02)
    tau_pop.append(dict(cv=cv, a25=float(r_op["a_q"][0]), t90=float(r_op["t90"][0]), S=float(sharp(r_op)[0])))
    print("  release-time CV = %.2f: alpha25 = %.4f  t90 = %.2f  S = %.3f" % (cv, tau_pop[-1]["a25"], tau_pop[-1]["t90"], tau_pop[-1]["S"]))
RES["population_release_time_cv"] = tau_pop

fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
sg_arr = np.array([r["sigma"] for r in pop])
ax[0].plot(sg_arr, [r["a_store"] for r in pop], "o-", color="#2E5597", lw=2, label="capsule population (model)")
ax[0].set_xlabel("Spread of capsule melt temperature, σ (°C)"); ax[0].set_ylabel("Storage conversion at %.0f°C, 150 min" % T_STORE)
ax[0].set_title("(a) Storage leakage rises with melt-point spread"); ax[0].grid(alpha=0.25)
ax[1].plot(sg_arr, [r["w_analytic"] for r in weff], "s--", color="#B5651D", lw=2, label="analytic: sqrt(w₀²+3σ²/π²)")
ax[1].plot(sg_arr, [r["w_matched_storage"] for r in weff], "o-", color="#2E5597", lw=2, label="matched on storage leakage")
ax[1].set_xlabel("σ (°C)"); ax[1].set_ylabel("Equivalent single-gate width w_eff (°C)"); ax[1].legend(fontsize=9)
ax[1].set_title("(b) Emergent gate width vs. DSC-measurable σ"); ax[1].grid(alpha=0.25)
for r in pop:
    pass
ax[2].plot([r["sigma"] for r in pop], [r["t90"] for r in pop], "o-", color="#2E5597", lw=2, label="t₉₀ (min)")
ax[2].set_xlabel("σ (°C)"); ax[2].set_ylabel("t₉₀ (min)", color="#2E5597")
ax2 = ax[2].twinx(); ax2.plot([r["sigma"] for r in pop], [r["a25"] for r in pop], "s--", color="#B5651D", lw=2)
ax2.set_ylabel("α at 25 min", color="#B5651D"); ax[2].set_title("(c) Operating-temperature behaviour"); ax[2].grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT + "t10_population.png", bbox_inches="tight"); plt.close()

# ================================================================== E2 storage axis and control families
print("\n=== E2 storage-stability axis and control families ===")
FAMS = {}


def evaluate(label, Tsets, builder):
    """builder(T_inf, idx_arrays) not needed: pass parameter dict and gating info"""


def run_family(name, params, mode, gated, rel_scale=1.0):
    S = len(next(iter(params.values())))
    op = integrate(np.full(S, T_REF), mode=mode, gated=gated, p=params, rel_scale=rel_scale, dt=0.01)
    sto = integrate(np.full(S, T_STORE), mode=mode, gated=gated, p=params, rel_scale=rel_scale, dt=0.01)
    out = dict(a25=op["a_q"], t90=op["t90"], t10=op["t10"], a_store=sto["a_end"], params=params, mode=mode, gated=gated, rel_scale=np.broadcast_to(np.asarray(rel_scale, dtype=float), (S,)).copy(),
               physical=(op["overshoot"] <= DT_AD + 1.0) & (sto["overshoot"] <= DT_AD + 1.0))
    FAMS[name] = out
    return out


# F1 gated first-order (Arrhenius with Ea of the control = Ea2 or Ea1), k swept
kk = np.geomspace(0.05, 6.0, 50) * K2_REF
for nm, eac in (("F1a gated first-order, Ea_c = Ea2", EA2), ("F1b gated first-order, Ea_c = Ea1", EA1)):
    run_family(nm, dict(kfirst=kk, eac=np.full(len(kk), eac)), "nth", True)
# F2 gated n-th order, n in {0.5, 2, 3}
for nord in (0.5, 2.0, 3.0):
    run_family("F2 gated order n = %.1f" % nord, dict(kfirst=kk, nord=np.full(len(kk), nord)), "nth", True)
# F3 un-gated autocatalytic with the prefactor A2 tuned
a2s = np.geomspace(0.02, 4.0, 60) * A2
run_family("F3 ungated autocatalytic, A2 tuned", dict(a2=a2s), "auto", False)
# F4 gated first-order with tuned release time (slow release produces an induction by itself)
rsc = np.geomspace(0.3, 12.0, 8)
kk2 = np.geomspace(0.05, 6.0, 30) * K2_REF
KK, RS = np.meshgrid(kk2, rsc)
fam4 = run_family("F4 gated first-order, tuned release time", dict(kfirst=KK.ravel(), eac=np.full(KK.size, EA2)), "nth", True, rel_scale=RS.ravel())
# F5 gated autocatalytic with detuned kinetics (is the coupled point on the front of its own family?)
a2g = np.geomspace(0.2, 3.0, 25) * A2
rsg = np.geomspace(0.3, 6.0, 8)
AG, RG = np.meshgrid(a2g, rsg)
run_family("F5 gated autocatalytic, A2 and release time tuned", dict(a2=AG.ravel()), "auto", True, rel_scale=RG.ravel())

from scipy.integrate import solve_ivp


def verify_design(d, i, Tinf):
    """Independent check of one control design with the adaptive stiff solver LSODA (different integrator and step control)."""
    pr = {k: float(v[i]) for k, v in d["params"].items()}
    rs = float(d["rel_scale"][i])
    gated, mode = d["gated"], d["mode"]
    Tref_K = T_REF + 273.15

    def rhs(t, y):
        T, C, a = y
        TK = T + 273.15
        ac = min(max(a, 0.0), 0.999)
        krel = A_REL * np.exp(-EA_REL / (R_GAS * TK)) / rs
        dC = krel * (1.0 / (1.0 + np.exp(np.clip(-(T - T_MELT) / W_DEF, -60, 60))) - C) if gated else 0.0
        if mode == "auto":
            da = C * (A1 * np.exp(-EA1 / (R_GAS * TK)) + pr.get("a2", A2) * np.exp(-EA2 / (R_GAS * TK)) * max(ac, 1e-6) ** M_EXP) * (1 - ac) ** N_EXP
        else:
            da = C * pr["kfirst"] * np.exp(-pr.get("eac", EA2) / R_GAS * (1.0 / TK - 1.0 / Tref_K)) * (1 - ac) ** pr.get("nord", 1.0)
        if a >= 0.999:
            da = 0.0
        return [(Tinf - T) / TAU_TH + DT_AD * da, dC, da]
    sol = solve_ivp(rhs, [0, 150.0], [Tinf - 30.0, 0.0 if gated else 1.0, 0.0], method="LSODA", rtol=1e-8, atol=1e-10, dense_output=True)
    tt = np.linspace(0, 150, 3001)
    Y = sol.sol(tt)
    t90v = tt[np.argmax(Y[2] >= 0.9)] if np.any(Y[2] >= 0.9) else np.nan
    return dict(a25=float(sol.sol(25.0)[2]), a_end=float(Y[2][-1]), t90=float(t90v), Tmax=float(Y[0].max()))


fam_summary = {}
for nm, d in FAMS.items():
    ok = ~np.isnan(d["t90"]) & d["physical"]
    two = ok & (d["a25"] <= COUP["a25"]) & (d["t90"] <= COUP["t90"]) & ((d["a25"] < COUP["a25"]) | (d["t90"] < COUP["t90"]))
    three = two & (d["a_store"] <= COUP["a_store"])
    confirmed3 = 0
    for i in np.where(three)[0]:
        vo, vs = verify_design(d, i, T_REF), verify_design(d, i, T_STORE)
        if vo["Tmax"] <= T_REF + DT_AD + 1 and vo["a25"] <= COUP["a25"] and vo["t90"] <= COUP["t90"] and vs["a_end"] <= COUP["a_store"]:
            confirmed3 += 1
    gap = np.where(ok, np.maximum.reduce([(d["a25"] - COUP["a25"]) / max(COUP["a25"], 1e-9), (d["t90"] - COUP["t90"]) / COUP["t90"],
                                          (d["a_store"] - COUP["a_store"]) / max(COUP["a_store"], 1e-4)]), np.nan)
    best_gap = float(np.nanmin(gap)) if np.any(ok) else float("nan")
    fam_summary[nm] = dict(n=int(len(d["t90"])), n_physical=int((d["physical"]).sum()), n_converged=int(ok.sum()), dominates_on_2_axes=int(two.sum()),
                           dominates_on_3_axes=int(three.sum()), confirmed_by_LSODA=int(confirmed3), min_worst_relative_gap=best_gap)
    print("  %-52s n=%3d physical=%3d  dom(a25,t90)=%3d  dom(3 axes)=%3d (LSODA-confirmed %d)  best worst-axis gap=%+.2f" % (
        nm, len(d["t90"]), d["physical"].sum(), two.sum(), three.sum(), confirmed3, best_gap))
RES["control_families"] = fam_summary
RES["control_note"] = "Arrhenius first-order/nth-order controls use Ea_c as stated and the same gate and energy balance as the coupled system"

fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
cols = ["#999999", "#bbbbbb", "#B5651D", "#d98d3d", "#e8b25f", "#8a2be2", "#2a9d8f", "#555555"]
for (nm, d), c in zip(FAMS.items(), cols):
    ax[0].scatter(d["a25"], d["t90"], s=12, alpha=0.6, label=nm, color=c)
    ax[1].scatter(d["a_store"], d["t90"], s=12, alpha=0.6, color=c)
for a_ in ax:
    a_.axhline(75, color="gray", ls=":", lw=1.2); a_.grid(alpha=0.25); a_.set_ylim(0, 150)
ax[0].scatter([COUP["a25"]], [COUP["t90"]], marker="*", s=300, color="red", edgecolor="black", zorder=6, label="Coupled system")
ax[1].scatter([COUP["a_store"]], [COUP["t90"]], marker="*", s=300, color="red", edgecolor="black", zorder=6)
ax[0].set_xlabel("α at 25 min, T∞ = 140°C"); ax[0].set_ylabel("t₉₀ at 140°C (min)"); ax[0].set_title("(a) The two v2 objectives")
ax[1].set_xlabel("Storage conversion at %.0f°C, 150 min" % T_STORE); ax[1].set_ylabel("t₉₀ at 140°C (min)"); ax[1].set_title("(b) Storage-stability objective added")
ax[1].set_xscale("symlog", linthresh=1e-3)
ax[0].legend(fontsize=6.5, loc="upper right")
plt.tight_layout(); plt.savefig(OUT + "t11_storage_pareto.png", bbox_inches="tight"); plt.close()

# ================================================================== E3 design-viability map
print("\n=== E3 design-viability map over (gate width, operating temperature) ===")
widths = np.array([0.5, 1, 2, 3, 5, 8, 12, 16, 20])
Tops = np.arange(122.0, 171.0, 4.0)
WW, TT = np.meshgrid(widths, Tops)
op = integrate(TT.ravel(), width=WW.ravel(), dt=0.05, T0_off=30.0)
leak = np.array([float(integrate([T_STORE], width=w, dt=0.05)["a_end"][0]) for w in widths])
EPS_LEAK = 0.01
t90g = op["t90"].reshape(WW.shape); t10g = op["t10"].reshape(WW.shape); osg = op["overshoot"].reshape(WW.shape)
viable = (leak[None, :] <= EPS_LEAK) & (t90g <= 75.0) & (t10g >= 5.0)
print("  storage leakage vs width:", {float(w): round(float(l), 4) for w, l in zip(widths, leak)})
print("  viable fraction of the grid = %.1f%%" % (100 * viable.mean()))
Tmin = {}
for j, w in enumerate(widths):
    ok = np.where(viable[:, j])[0]
    Tmin[float(w)] = (float(Tops[ok].min()), float(Tops[ok].max())) if len(ok) else None
print("  viable operating-temperature window (min,max) by width:", Tmin)
w_crit = float(widths[np.max(np.where(leak <= EPS_LEAK))]) if np.any(leak <= EPS_LEAK) else None
RES["design_map"] = dict(widths=widths.tolist(), T_ops=Tops.tolist(), leakage=leak.tolist(), eps_leak=EPS_LEAK, viable_fraction=float(viable.mean()),
                         window_by_width=Tmin, widest_viable_gate=w_crit,
                         peak_overshoot_range=[float(np.nanmin(osg)), float(np.nanmax(osg))])

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
im = ax[0].pcolormesh(widths, Tops, np.where(viable, 1.0, 0.0) + 0.5 * ((leak[None, :] <= EPS_LEAK) & ~viable), shading="nearest", cmap="Greens", vmin=0, vmax=1.2)
cs = ax[0].contour(widths, Tops, t90g, levels=[75], colors="k", linewidths=1.5); ax[0].clabel(cs, fmt="t₉₀ = 75 min")
ax[0].axvline(widths[np.max(np.where(leak <= EPS_LEAK))] if w_crit else 0, color="red", ls="--", lw=1.5)
ax[0].set_xscale("log"); ax[0].set_xlabel("Gate width w (°C)"); ax[0].set_ylabel("Operating temperature T∞ (°C)")
ax[0].set_title("(a) Dark green: viable (leakage ≤ %.2f, t₁₀ ≥ 5, t₉₀ ≤ 75 min);
light green: leakage met, t₉₀ not" % EPS_LEAK, fontsize=9.5)
ax[1].plot(widths, np.maximum(leak, 1e-5), "o-", color="#2E5597", lw=2); ax[1].axhline(EPS_LEAK, color="red", ls="--")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_ylim(5e-5, 0.1); ax[1].set_xlabel("Gate width w (°C)"); ax[1].set_ylabel("Storage conversion at %.0f°C, 150 min" % T_STORE)
ax[1].set_title("(b) Storage criterion fixes the widest admissible gate"); ax[1].grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT + "t12_design_map.png", bbox_inches="tight"); plt.close()

# E3b design rule: widest admissible gate and largest admissible melt-point spread for a stated storage margin and leakage tolerance
print("\n--- E3b design rule: admissible spread of capsule melt temperature ---")
zq = z


def storage_pop(sigma, margin):
    return float(integrate([T_MELT - margin], Tmelt=T_MELT, width=W_CAP, off=sigma * zq, dt=0.05)["a_end"][0])


def storage_single(w, margin):
    return float(integrate([T_MELT - margin], width=w, dt=0.05)["a_end"][0])


table = []
for margin in (3.0, 6.0, 10.0):
    for eps in (0.001, 0.01, 0.05):
        w_star = brentq(lambda w: storage_single(w, margin) - eps, 0.3, 60.0, xtol=1e-2) if storage_single(60.0, margin) > eps else None
        s_star = brentq(lambda sg: storage_pop(sg, margin) - eps, 0.0, 40.0, xtol=1e-2) if (storage_pop(0.0, margin) < eps < storage_pop(40.0, margin)) else None
        table.append(dict(margin=margin, eps=eps, w_star=None if w_star is None else float(w_star), sigma_star=None if s_star is None else float(s_star)))
        print("  margin %4.1f C, tolerance %.3f: widest single gate w* = %s C ; largest melt-point spread sigma* = %s C" % (
            margin, eps, "n/a" if w_star is None else "%.2f" % w_star, "n/a" if s_star is None else "%.2f" % s_star))
RES["design_rule"] = table
T_lo = brentq(lambda T: float(np.nan_to_num(integrate([T], width=W_DEF, dt=0.05)["t90"][0], nan=1e3)) - 75.0, 122.0, 140.0, xtol=0.05)
print("  lowest operating temperature with t90 <= 75 min (w = 3 C): %.1f C" % T_lo)
RES["design_rule_T_min_t90_75"] = float(T_lo)

# ================================================================== E4 Sobol attribution
print("\n=== E4 Sobol attribution (Jansen estimators, Saltelli sampling) ===")
names = ["Ea1", "Ea2", "A2", "m", "n", "T_melt offset", "tau_release", "tau_thermal", "DT_ad"]
d = len(names)
NB = 2048
U = qmc.Sobol(d=2 * d, scramble=True, seed=20260906).random(NB)
UA, UB = U[:, :d], U[:, d:]


def to_params(Um):
    return dict(ea1=4.2e4 + Um[:, 0] * 1.6e4, ea2=5.8e4 + Um[:, 1] * 1.4e4, a2=1.0e7 + Um[:, 2] * 4.0e7, m=0.6 + Um[:, 3] * 0.8, n=1.1 + Um[:, 4] * 0.8,
                off=norm.ppf(np.clip(Um[:, 5], 1e-6, 1 - 1e-6)) * 5.0 + 22.0, trel=2.0 + Um[:, 6] * 13.0, tth=1.0 + Um[:, 7] * 7.0, dtad=150.0 + Um[:, 8] * 150.0)


def model_t90(Um):
    q = to_params(Um)
    r = integrate(np.full(len(Um), T_REF), tau_th=q["tth"], rel_scale=(q["trel"] / TAU_REL_REF).reshape(-1, 1), Tmelt=T_REF - q["off"], width=W_DEF,
                  p=dict(ea1=q["ea1"], ea2=q["ea2"], a2=q["a2"], m=q["m"], n=q["n"], dtad=q["dtad"]), dt=0.05)
    return np.where(np.isnan(r["t90"]), 150.0, r["t90"]), r["t90"]


t1 = time.time()
fA, rawA = model_t90(UA)
fB, _ = model_t90(UB)
fAB = []
for i in range(d):
    Uab = UA.copy(); Uab[:, i] = UB[:, i]
    fAB.append(model_t90(Uab)[0])
fAB = np.array(fAB)
V = np.var(np.concatenate([fA, fB]), ddof=1)
S1 = [float((V - 0.5 * np.mean((fB - fAB[i]) ** 2)) / V) for i in range(d)]
ST = [float(0.5 * np.mean((fA - fAB[i]) ** 2) / V) for i in range(d)]
cens = float(np.mean(np.isnan(rawA)))
print("  N_base = %d (%d model runs, %.0f s); censored (no t90 within 150 min) = %.1f%%" % (NB, NB * (d + 2), time.time() - t1, 100 * cens))
for nme, a_, b_ in sorted(zip(names, S1, ST), key=lambda x: -x[2]):
    print("   %-14s S1 = %+.3f   ST = %.3f" % (nme, a_, b_))
RES["sobol"] = dict(N_base=NB, runs=NB * (d + 2), censored_fraction=cens, names=names, S1=S1, ST=ST, sum_S1=float(sum(S1)), interaction_share=float(1 - sum(S1)))
print("  sum S1 = %.3f (interaction share ~ %.3f)" % (sum(S1), 1 - sum(S1)))

order = np.argsort(ST)
fig, ax = plt.subplots(figsize=(7.4, 4.6))
y = np.arange(d)
ax.barh(y + 0.2, np.array(ST)[order], height=0.38, color="#2E5597", label="total-order ST")
ax.barh(y - 0.2, np.array(S1)[order], height=0.38, color="#B5651D", label="first-order S1")
ax.set_yticks(y); ax.set_yticklabels(np.array(names)[order]); ax.set_xlabel("Variance share of t₉₀"); ax.legend(); ax.grid(alpha=0.25, axis="x")
ax.set_title("Sobol attribution of t₉₀ (N_base = %d)" % NB)
plt.tight_layout(); plt.savefig(OUT + "t13_sobol.png", bbox_inches="tight"); plt.close()

with open(os.path.join(HERE, "..", "results_v3.json"), "w", encoding="utf-8") as fh:
    json.dump(RES, fh, indent=2, default=str)
print("\ntotal runtime %.0f s; results_v3.json written" % (time.time() - t0))
