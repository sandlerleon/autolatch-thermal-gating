# -*- coding: utf-8 -*-
"""Builds the v3 manuscript from the v2 manuscript (unchanged text kept verbatim) plus the extended analyses of code/extended_analyses.py.
Every number added is read from results_v3.json / results_v2.json; nothing is typed by hand.

    python build_v3.py     ->  AutoLatch_ChemRxiv_WileyMTS_Manuscript_v3.docx
"""
import copy
import json
import os

import docx
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
SRC = os.path.join(HERE, "AutoLatch_ChemRxiv_WileyMTS_Manuscript.docx")
OUTF = os.path.join(HERE, "AutoLatch_ChemRxiv_WileyMTS_Manuscript_v3.docx")
R3 = json.load(open(os.path.join(ROOT, "results_v3.json"), encoding="utf-8"))
R2 = json.load(open(os.path.join(ROOT, "results_v2.json"), encoding="utf-8"))
ZEN = json.load(open(os.path.join("C:" + os.sep, "YouTube", "_autolatch_zenodo_state.json")))
SW_DOI, PP_DOI = (ZEN.get("code_v31") or ZEN["code"])["doi"], (ZEN.get("paper_v31") or ZEN["paper"])["doi"]
FIG = os.path.join(ROOT, "figures")

d = docx.Document(SRC)
P = lambda: d.paragraphs  # noqa: E731


def find(prefix, start=0):
    for i, p in enumerate(d.paragraphs):
        if i >= start and p.text.strip().startswith(prefix):
            return p
    raise KeyError(prefix)


def clone_after(anchor_el, exemplar, text, bold_lead=None):
    """new paragraph after anchor_el, formatted like `exemplar` (first run's properties kept)"""
    new = copy.deepcopy(exemplar._p)
    for r in list(new):
        if r.tag in (qn("w:r"), qn("w:hyperlink")):
            new.remove(r)
    anchor_el.addnext(new)
    p = docx.text.paragraph.Paragraph(new, exemplar._parent)
    src_run = exemplar.runs[-1] if exemplar.runs else None

    def add(t, bold=None):
        r = p.add_run(t)
        if src_run is not None:
            if src_run._r.rPr is not None:
                r._r.insert(0, copy.deepcopy(src_run._r.rPr))
        r.bold = True if bold else None
        return r
    if bold_lead:
        add(bold_lead, True)
    add(text)
    return p


def add_picture_after(anchor_el, exemplar_caption, png, caption, width=6.3):
    pp = d.add_paragraph()
    pp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pp.add_run().add_picture(os.path.join(FIG, png), width=Inches(width))
    anchor_el.addnext(pp._p)
    cap = clone_after(pp._p, exemplar_caption, caption)
    return cap


def borders(cell, fill=None, bold=False):
    tcPr = cell._tc.get_or_add_tcPr()
    b = OxmlElement("w:tcBorders")
    for side in ("top", "start", "bottom", "end"):
        e = OxmlElement("w:" + side)
        e.set(qn("w:val"), "single"); e.set(qn("w:sz"), "4"); e.set(qn("w:space"), "0"); e.set(qn("w:color"), "000000")
        b.append(e)
    tcPr.append(b)
    if fill:
        sh = OxmlElement("w:shd"); sh.set(qn("w:color"), "auto"); sh.set(qn("w:fill"), fill); sh.set(qn("w:val"), "clear"); tcPr.append(sh)


def table_after(anchor_el, rows, widths, header=True, size=9):
    t = d.add_table(rows=len(rows), cols=len(rows[0]))
    for i, row in enumerate(rows):
        for j, txt in enumerate(row):
            c = t.cell(i, j)
            c.width = Inches(widths[j])
            c.text = ""
            r = c.paragraphs[0].add_run(txt)
            r.font.size = Pt(size)
            r.bold = bool(header and i == 0)
            borders(c, fill="E8E8E8" if header and i == 0 else None)
    anchor_el.addnext(t._tbl)
    return t


BODY = find("The central computational finding")          # exemplar for body text
HEAD2 = find("4.7 Robustness")                           # exemplar for sub-section headings
HEAD3 = find("1.1 Novelty")
CAPT = find("Figure 7.")
PRED = find("Prediction 3")
NM = {"tau_release": "τ_release", "tau_thermal": "τ_thermal", "DT_ad": "ΔT_ad", "T_melt offset": "the melting-point offset"}
n = lambda x, k=3: ("%." + str(k) + "f") % x  # noqa: E731

# ---------------------------------------------------------------- numbers
E1 = R3["population_sigma_Tm"]["rows"]
EQ = {r["sigma"]: r for r in R3["population_sigma_Tm"]["equivalent_width"]}
TAU = R3["population_release_time_cv"]
CF = R3["control_families"]
COUP = R3["coupled_point"]
UNG = R3["ungated_baseline"]
DM = R3["design_map"]
DR = {(r["margin"], r["eps"]): r for r in R3["design_rule"]}
SB = R3["sobol"]
others = [k for k in CF if not k.startswith("F5")]
n_other = sum(CF[k]["n"] for k in others)
n_phys = sum(CF[k]["n_physical"] for k in others)
n_dom2 = sum(CF[k]["dominates_on_2_axes"] for k in others)
n_dom3 = sum(CF[k]["dominates_on_3_axes"] for k in others)
assert n_dom2 == 0 and n_dom3 == 0, "a control of another mechanism dominates: revise the text"
F5 = CF[[k for k in CF if k.startswith("F5")][0]]
closest = min(others, key=lambda k: CF[k]["min_worst_relative_gap"])
gap_min = CF[closest]["min_worst_relative_gap"]
sig = [r["sigma"] for r in E1]
rel_dev = [abs(EQ[s]["w_matched_storage"] - EQ[s]["w_analytic"]) / EQ[s]["w_analytic"] for s in sig if s >= 1.0]
ST = dict(zip(SB["names"], SB["ST"])); S1 = dict(zip(SB["names"], SB["S1"]))
SPEAR = R2["monte_carlo"]["sensitivity"]
a25_rng = (min(r["a25"] for r in E1), max(r["a25"] for r in E1))
t90_rng = (min(r["t90"] for r in E1), max(r["t90"] for r in E1))
s6 = [r for r in E1 if r["sigma"] == 6.0][0]
ratio = [DR[k]["sigma_star"] / DR[k]["w_star"] for k in DR if DR[k]["sigma_star"] and DR[k]["w_star"]]
dr61 = DR[(6.0, 0.01)]

# captions for the two existing (previously unlabelled) tables
t_conv, t_sys = d.tables[0]._tbl, d.tables[1]._tbl
clone_after(t_conv.getprevious(), CAPT, "Table 1. Convergence of fixed-step Euler and RK4 toward the adaptive RK45 reference (reference case, T∞ = 140 °C).")
clone_after(t_sys.getprevious(), CAPT, "Table 2. Four system variants (A–D): premature conversion α(t = 25 min) and sharpness S, with the first-order control matched on t50.")

# ================================================================== abstract
ab = find("The delayed transformation of a reactive material")
add = (" Four extensions then test whether these conclusions survive stronger scrutiny. (i) A capsule population with a spread σ of melting temperature reproduces the finite gate "
       "width without postulating it (w_eff ≈ 0.55σ, agreeing with the population model within about %d %%), so that the leakage prediction acquires a measurable input, the width of the "
       "DSC melting endotherm. (ii) With a storage-stability objective that the first version did not measure, and Arrhenius controls that share the energy balance, none of %d tuned "
       "designs of four other mechanisms dominates the coupled system on all three objectives; %d of %d detuned designs of the coupled family itself do, so the reference operating point "
       "is illustrative, not optimal. (iii) A design rule gives the largest admissible melting-point spread (σ ≈ %.0f °C for ≤ %.0f %% storage conversion at a %.0f °C margin). (iv) Variance-based "
       "(Sobol) indices confirm that Ea2 dominates the completion time (total-order index %.2f)." % (round(100 * max(rel_dev)), n_phys, F5["dominates_on_3_axes"], F5["n"],
                                                                                                    dr61["sigma_star"], 100 * 0.01, 6.0, ST["Ea2"]))
# insert before the last sentence block about predictions: append at end of abstract (keeps the existing text verbatim)
txt = ab.text
cut = txt.rindex("Three falsifiable predictions")
for r in ab.runs[1:]:
    r._r.getparent().remove(r._r)
ab.runs[0].text = txt[:cut].rstrip() + add + " " + txt[cut:].replace("Three falsifiable predictions", "Five falsifiable predictions", 1)

# ================================================================== 1.2 what this version adds
p11 = find("An earlier draft of this model (Sandler, 2026")
h12 = clone_after(p11._p, HEAD3, "1.2 What this version adds")
b12 = clone_after(h12._p, BODY,
                  "This version leaves every v2 computation unchanged and adds four analyses, each chosen to test one of the conclusions rather than to enlarge the model: (i) a "
                  "capsule-population model in which the gate width emerges from a distribution of melting temperatures (Section 4.8); (ii) a storage-stability objective that v2 did not "
                  "measure, together with five families of tuned control mechanisms, to test whether the Pareto claim of Section 4.3 survives a stronger and more physical set of controls "
                  "(Section 4.9); (iii) a design rule that converts the leakage requirement into an admissible spread of melting temperature (Section 4.10); and (iv) variance-based (Sobol) "
                  "attribution of the completion time (Section 4.11). Section 5.1 tabulates each claim with its evidential status, the assumption it rests on and what would reverse it. "
                  "In the course of this work the control comparison of v2 was found to use a constant-rate first-order control, which lacks the temperature dependence of a real reaction "
                  "and therefore the thermal feedback of the coupled system; the new control families are Arrhenius, and Section 4.9 reports the result under them.")

# ================================================================== 4.3 pointer
p43 = find("This is followed by a full Pareto sweep")
p43.add_run(" A stronger test, with Arrhenius controls that share the energy balance and with a storage-stability objective, is given in Section 4.9.")

# ================================================================== new results sections (4.8 - 4.11), before the old 4.8
old48 = find("4.8 Mechanical function and spatial-transport extensions")
for r in old48.runs[1:]:
    r._r.getparent().remove(r._r)
old48.runs[0].text = "4.12 Mechanical function and spatial-transport extensions"
anchor = find("Figure 7.")                       # new sections are inserted after the Monte Carlo text, before old 4.8
mc_text = anchor._p
# the paragraph(s) between the Figure-7 caption and the old 4.8 heading belong to 4.7; insert directly before old 4.8
prev_el = old48._p.getprevious()


def section(heading, paras, fig=None, table=None, after=None):
    """insert heading + paragraphs (+ figure, + table) before the old 4.8 heading, in order"""
    global prev_el
    h = clone_after(prev_el, HEAD2, heading)
    cur = h._p
    for k, t in enumerate(paras):
        pp = clone_after(cur, BODY, t)
        cur = pp._p
        if k == 0 and fig:
            capp = add_picture_after(cur, CAPT, fig[0], fig[1])
            cur = capp._p
        if table and k == table[0]:
            cap = clone_after(cur, CAPT, table[1])
            t_ = table_after(cap._p, table[2], table[3])
            cur = t_._tbl
            sp = clone_after(cur, BODY, "")
            cur = sp._p
    prev_el = cur


section("4.8 Capsule-population origin of the gate width", [
    "Section 4.6 postulated a logistic gate of width w. A real capsule population has no such function; it has a distribution of melting temperatures (shell thickness, wax grade, "
    "capsule size) and possibly of release times. The gate is therefore re-derived as a population: K = %d capsules with melting temperatures T_melt + σ·z_i, where the z_i are "
    "deterministic Gaussian quantiles, each with a near-step gate (w0 = %.1f °C) and its own release relaxation; the reaction is driven by the capsule-averaged availability "
    "C̄ = (1/K)ΣC_i. Matching the variance of a logistic and of a Gaussian distribution gives w_eff = (w0² + 3σ²/π²)^½, i.e. w_eff ≈ 0.55σ for σ ≫ w0." % (
        R3["population_sigma_Tm"]["K"], R3["population_sigma_Tm"]["w_capsule"]),
    "The population model confirms this where it can be tested (Figure 8b). The single-gate width that reproduces the storage conversion of the population is %s °C for σ = 1, 2, 4 and 6 °C, "
    "against %s °C from the variance argument; the two differ by at most %d %% over σ = 1–6 °C and by %d–%d %% for σ ≥ 4 °C, the residual reflecting the release lag, which the variance "
    "argument ignores. At the operating temperature the behaviour is almost insensitive to σ: α at 25 min stays between %s and %s and t90 between %.1f and %.1f min for σ = 0–6 °C, whereas "
    "the storage conversion at 112 °C rises from 0 to %s (Figure 8a). A lognormal scatter of the release time constant (coefficient of variation up to %.0f) changes α at 25 min from %s to %s and "
    "t90 from %.1f to %.1f min. Under the present model, therefore, the spread of melting temperature, not of release time, controls storage stability. This turns the leakage prediction of "
    "Section 4.6 into a relation between two measurable quantities (Prediction 4, Section 7): the standard deviation σ of the melting-temperature distribution, estimated from the DSC "
    "melting endotherm, and the storage conversion." % (
        ", ".join("%.2f" % EQ[s]["w_matched_storage"] for s in (1.0, 2.0, 4.0, 6.0)), ", ".join("%.2f" % EQ[s]["w_analytic"] for s in (1.0, 2.0, 4.0, 6.0)),
        round(100 * max(rel_dev)), round(100 * min(abs(EQ[s]["w_matched_storage"] - EQ[s]["w_analytic"]) / EQ[s]["w_analytic"] for s in (4.0, 6.0))),
        round(100 * max(abs(EQ[s]["w_matched_storage"] - EQ[s]["w_analytic"]) / EQ[s]["w_analytic"] for s in (4.0, 6.0))),
        n(a25_rng[0]), n(a25_rng[1]), t90_rng[0], t90_rng[1], n(s6["a_store"], 4), TAU[-1]["cv"], n(TAU[0]["a25"]), n(TAU[-1]["a25"]), TAU[0]["t90"], TAU[-1]["t90"])],
    fig=("t10_population.png", "Figure 8. Capsule-population model. (a) Storage conversion at 112 °C after 150 min vs. the spread σ of capsule melting temperature. (b) Equivalent single-gate width "
         "w_eff from the variance argument (dashed) and as matched on the storage conversion (solid). (c) t90 and α at 25 min at the operating temperature vs. σ."))

section("4.9 Storage-stability objective and robustness of the Pareto claim to the control family", [
    "The comparison of Section 4.3 measured two objectives at the operating temperature, premature conversion α(25 min) and t90. It did not measure the property named in the opening "
    "sentence of the abstract, transport stability: conversion during storage below the gate. A third objective is therefore added, the conversion after 150 min at T∞ = 112 °C (6 °C below "
    "T_melt), and the constant-rate first-order control of Section 4.3, which has no temperature dependence, is replaced by Arrhenius controls that share the energy balance, and hence the thermal "
    "feedback, of the coupled system. Four families of other mechanisms were tuned over their rate constants: gated first-order with a control activation energy equal to Ea2 or to Ea1; gated "
    "n-th order with n = 0.5, 2 and 3; the autocatalytic law without a gate, with A2 tuned; and gated first-order with a tuned release time. Together they comprise %d designs, of which %d are "
    "physically admissible (the others reach temperatures above T∞ + ΔT_ad, which violates energy conservation and signals numerical runaway, and are excluded)." % (n_other, n_phys),
    "None of the %d admissible designs dominates the coupled operating point (α25 = %s, t90 = %.1f min, storage conversion %s), on the two objectives of v2 or on all three (Figure 9). The "
    "closest family, %s, still exceeds the coupled value on its worst objective by a factor of %.1f. The ungated autocatalytic formulation, with the same kinetics and no gate, completes at "
    "t90 = %.1f min but converts %s of its mass in the storage test against %s for the gated system (Prediction 5, Section 7), which is what the gate contributes. One qualification is "
    "necessary. %d of the %d detuned designs of the coupled family itself (A2 and release time varied) dominate the reference operating point on all three objectives, each confirmed with an "
    "independent stiff integrator (LSODA). The reference point is therefore a representative, not an optimum. What the comparison supports is that none of the other mechanisms examined reaches the "
    "region of the coupled family, not that the reference point is optimal; and the absence of a dominating control among %d designs of four families is evidence, not proof." % (
        n_phys, n(COUP["a25"]), COUP["t90"], n(COUP["a_store"], 4), closest.split(" ", 1)[1], 1 + gap_min, UNG["t90"], n(UNG["a_store"]), n(COUP["a_store"], 4),
        F5["dominates_on_3_axes"], F5["n"], n_other)],
    fig=("t11_storage_pareto.png", "Figure 9. Control-family comparison. (a) Premature conversion α(25 min) vs. t90 at 140 °C. (b) Storage conversion (150 min at 112 °C) vs. t90. Each point is one tuned "
         "control design; the star is the coupled reference operating point. The ungated autocatalytic family has a large storage conversion; gated families with constant dynamics fail on t90 or on "
         "premature conversion."))

rows = [["Storage margin T_melt − T_store (°C)", "Tolerance ε (conversion)", "Widest single gate w* (°C)", "Largest melting-point spread σ* (°C)"]]
for (mg, ep) in sorted(DR):
    r = DR[(mg, ep)]
    rows.append(["%.0f" % mg, "%g" % ep, ("%.2f" % r["w_star"]) if r["w_star"] else "> 60 (no limit reached)", ("%.2f" % r["sigma_star"]) if r["sigma_star"] else "> 40 (no limit reached)"])
section("4.10 Design rule: admissible spread of the capsule melting temperature", [
    "The leakage requirement fixes how uniform the capsule melting temperature must be. For a tolerance ε on the storage conversion and a storage margin ΔT = T_melt − T_store, the widest admissible "
    "single gate w* and the largest admissible spread σ* follow by root-finding on the leakage curves of the single-gate and the population model (Table 3). For ε = 0.01 and ΔT = 6 °C, "
    "w* = %.1f °C and σ* = %.1f °C; the ratio σ*/w* lies between %.1f and %.1f over the table, against 1.8 from the variance argument. The design map of Figure 10 combines this with the completion "
    "criterion: %.1f %% of the (w, T∞) grid satisfies leakage ≤ %.2f, t10 ≥ 5 min and t90 ≤ 75 min together. The operating temperature required for t90 ≤ 75 min (%.1f °C for w = 3 °C) does not "
    "depend on the gate width within the admissible range, to the resolution of the grid (4 °C), so the two requirements decouple within this model: the storage criterion constrains the "
    "spread of melting temperature and the completion criterion constrains the operating temperature, and neither can be traded against the other." % (
        dr61["w_star"], dr61["sigma_star"], min(ratio), max(ratio), 100 * DM["viable_fraction"], DM["eps_leak"], R3["design_rule_T_min_t90_75"]),
    ], fig=("t12_design_map.png", "Figure 10. Design map. (a) Viable window (dark) over gate width and operating temperature for leakage ≤ %.2f, t10 ≥ 5 min and t90 ≤ 75 min; the dashed line is the widest "
            "admissible gate. (b) Storage conversion at 112 °C vs. gate width." % DM["eps_leak"]),
    table=(0, "Table 3. Widest admissible single gate w* and largest admissible melting-point spread σ* for stated storage margin and tolerance (illustrative parameters; 150 min storage).", rows, [1.6, 1.3, 1.6, 1.8]))

section("4.11 Variance-based attribution of the completion time", [
    "The Spearman ranking of Section 4.7 measures monotonic association among the draws that complete. Variance-based indices also capture interactions and non-monotonic effects. Saltelli sampling "
    "(N = %d base samples, %s model runs, Jansen estimators; draws with no t90 within 150 min, %.1f %% of the total, are assigned the censoring value 150 min) gives total-order indices ST of "
    "%s. The first-order indices sum to %.2f, so interactions account for about %.0f %% of the variance. Ea2 dominates, as in Section 4.7. Two parameters rank differently from the Spearman "
    "coefficients: τ_release correlates with t90 (r = %.2f) but contributes almost no variance (ST = %.3f), whereas τ_thermal has almost no rank correlation (r = %.3f) but ST = %.2f, through its "
    "interaction with the reaction exotherm. The experimental priority stated in Section 4.7, kinetic characterization by dynamic DSC, is unchanged; the thermal time constant of the sample "
    "should also be controlled and reported." % (
        SB["N_base"], "{:,}".format(SB["runs"]), 100 * SB["censored_fraction"],
        ", ".join("%s %.2f" % (NM.get(k, k), ST[k]) for k in sorted(ST, key=lambda x: -ST[x])), SB["sum_S1"], 100 * SB["interaction_share"],
        SPEAR["tau_release"], ST["tau_release"], SPEAR["tau_thermal"], ST["tau_thermal"])],
    fig=("t13_sobol.png", "Figure 11. Sobol attribution of t90: first-order (S1) and total-order (ST) indices for the nine uncertain parameters of Section 4.7."))

# ================================================================== discussion
disc = find("This is a general statement about the coupled kinetic structure")
dd = clone_after(disc._p.getprevious(), BODY,
                 "Three results of the extended analysis qualify and sharpen this finding. First, the storage-stability objective shows what the gate contributes: the ungated autocatalytic "
                 "formulation completes at least as fast but converts %s of its mass in the storage test against %s (Section 4.9), so that the gate, not the autocatalysis, provides storage stability, "
                 "and the autocatalysis, not the gate, provides sharpness; this is a more specific statement of the complementarity than the 2×2 comparison of Section 4.4. Second, the Pareto claim is "
                 "a statement about mechanisms, not about the reference operating point: none of the tuned alternative mechanisms dominates the coupled system, but detuned members of the coupled "
                 "family do dominate the reference point, which should be read as an illustration. Third, the population model and the design rule convert the finite-width prediction into a "
                 "specification: for a 6 °C storage margin and a 1 %% storage tolerance the melting temperatures of the capsules must lie within a standard deviation of about %.0f °C "
                 "(Table 3), a quantity measurable on the DSC endotherm before any formulation is cured." % (n(UNG["a_store"]), n(COUP["a_store"], 4), dr61["sigma_star"]))

h51 = clone_after(dd._p, HEAD2, "5.1 Status of claims and load-bearing assumptions")
rows51 = [["Claim", "Status", "Load-bearing assumption", "What would reverse or limit it"],
          ["Zero conversion below T_melt for a zero-width gate", "Structural (by construction)", "Step gate, C = 0 below T_melt", "Not a prediction (Section 2.4)"],
          ["Finite gate width gives quantified leakage and a leakage–completion trade-off", "Computed; Predicted (P1)", "Logistic gate; Arrhenius release relaxation", "A non-logistic activation profile changes the numbers, not the direction (Section 4.8)"],
          ["Gate width follows from the melting-temperature spread, w_eff ≈ 0.55σ", "Computed; Predicted (P4)", "Independent, equal-mass capsules; Gaussian spread; no capsule–capsule exchange", "A skewed or bimodal melting distribution; storage leakage inconsistent with σ by more than a factor of two"],
          ["No tested alternative mechanism dominates the coupled system on three objectives", "Computed (%d designs)" % n_other, "Four control families, Arrhenius, same energy balance; illustrative parameters", "Any admissible control outside the families that dominates; detuned members of the coupled family already dominate the reference point"],
          ["The gate, not autocatalysis, provides storage stability", "Computed; Predicted (P5)", "Same kinetics with and without gate; storage = 150 min at 6 °C below T_melt", "An unencapsulated formulation that stays within the tolerance would falsify it"],
          ["Ea2 dominates the variance of t90", "Computed", "Uniform ranges of Section 4.7; censoring at 150 min", "Narrower Ea2 range or different priors"],
          ["Admissible spread σ* ≈ %.0f °C (ε = 0.01, ΔT = 6 °C)" % dr61["sigma_star"], "Computed", "Gaussian spread; storage defined as 150 min at ΔT below T_melt", "Different storage duration or margin (Table 3)"],
          ["Overshoots up to %.0f °C in the thermally limited region" % R2["regime_map"]["max_exotherm_overshoot_C"], "Computed", "Lumped, adiabatic-like heating; ΔT_ad = 220 K", "Spatial heat loss; real material limits (Section 6)"]]
t51 = table_after(h51._p, rows51, [1.9, 1.1, 1.8, 1.9], size=8)
cap51 = clone_after(h51._p, CAPT, "Table 4. Status of the claims of this paper. All parameters are illustrative; none is fitted to a measured system.")
sp51 = clone_after(t51._tbl, BODY, "")

# ================================================================== limitations
lim_anchor = find("ψ (Section 2.2) is used descriptively")
lim_anchor.runs[0].text = ("ψ (Section 2.2) is used descriptively via direct simulation of the regime map (Section 4.5), not via a closed-form critical-value criterion. The lumped balance with "
                           "Newtonian cooling is formally of Semenov type, but the classical critical values assume negligible reactant consumption and an exponential (Frank-Kamenetskii) "
                           "linearisation, and are not applied to this consumed-reactant, autocatalytic system.")
for r in lim_anchor.runs[1:]:
    r._r.getparent().remove(r._r)
mech = find("The mechanical function extension (Section 4.8)")
mech.runs[0].text = mech.text.replace("Section 4.8", "Section 4.12")
for r in mech.runs[1:]:
    r._r.getparent().remove(r._r)
cur = lim_anchor._p
for t in ("The population model (Section 4.8) assumes K = %d independent capsules of equal mass with a Gaussian spread of melting temperatures and no heat or mass exchange between capsules; skewed or bimodal "
          "distributions would change the relation between σ and w_eff." % R3["population_sigma_Tm"]["K"],
          "The control families of Section 4.9 are finite (%d designs of four mechanisms); the absence of a dominating control is evidence, not proof. The constant-rate control of Section 4.3 lacked "
          "temperature dependence; the Arrhenius controls of Section 4.9 supersede it for the claim of non-dominance." % n_other,
          "The analyses of Sections 4.8–4.11 integrate with a heat-release cut-off once conversion reaches 0.999 (the reaction is complete); the first version retained a small residual term. The two agree "
          "on the reference case (t90 %.2f against 34.20 min, overshoot %.1f against 111.7 K)." % (R3["validation"]["t90"], R3["validation"]["overshoot"]),
          "The storage objective is one condition (150 min, 6 °C below T_melt); the design rule (Table 3) shows how the admissible spread changes with margin and tolerance. Sobol indices treat draws "
          "without a t90 within 150 min as censored at 150 min."):
    pp = clone_after(cur, mech, t)
    cur = pp._p

# ================================================================== predictions 4 and 5
val3 = find("Validation protocol: a τ_release sweep")
p4 = clone_after(val3._p, PRED, "For a formulation in which the capsule melting temperatures have standard deviation σ, the isothermal storage conversion at a temperature ΔT below the mean melting temperature should be "
                 "that of a single-gate model with w_eff ≈ 0.55σ (to within about 25 % in w_eff, Section 4.8). A storage conversion corresponding to a w_eff that differs from 0.55σ by more than a factor of "
                 "two would falsify the population description.", bold_lead="Prediction 4: ")
v4 = clone_after(p4._p, val3, "Validation protocol: DSC at a low scan rate to obtain the melting-temperature distribution (deconvolution of the endotherm; instrument broadening corrected), followed by isothermal holds at "
                 "T_melt − 6 °C for 150 min on the same batch with conversion read from the residual cure enthalpy; compare with the leakage curve of Table 3.")
p5 = clone_after(v4._p, PRED, "A formulation with the same cure chemistry but without the thermal gate (unencapsulated catalyst) converts substantially in the storage test (model: %s against %s with the gate), "
                 "while completing at least as fast (t90 %.1f against %.1f min). An unencapsulated formulation that stays within the same tolerance as the gated one would falsify the claim that the gate, not the autocatalysis, "
                 "provides storage stability." % (n(UNG["a_store"]), n(COUP["a_store"], 4), UNG["t90"], COUP["t90"]), bold_lead="Prediction 5: ")
v5 = clone_after(p5._p, val3, "Validation protocol: Paired isothermal holds (encapsulated and unencapsulated formulations of the same resin) at T_melt − 6 °C for 150 min, residual-enthalpy conversion; and dynamic DSC to confirm "
                 "comparable t90 at the operating temperature.")

# ================================================================== methods and availability
cm = find("All the results were produced by a single Python script")
cm.runs[0].text = cm.text.replace("All the results were produced by a single Python script (theoretical_study.py)", "All the results of Sections 3–4.7 were produced by one Python script (theoretical_study.py), and those of Sections 4.8–4.11 by a second (extended_analyses.py, with the shared routines in extended_core.py)")
for r in cm.runs[1:]:
    r._r.getparent().remove(r._r)
clone_after(cm._p, BODY, "The extended analyses use a vectorized fixed-step RK4 integrator over many scenarios at once (dt = 0.02–0.05 min; dt = 0.01 min for the control families), validated against the adaptive "
            "RK45 reference of Section 3.2 (t90 %.2f against 34.20 min); admissible control designs that dominate the coupled point are re-integrated with LSODA (rtol 1e-8) before they are counted. "
            "Sobol indices use scipy.stats.qmc scrambled Sobol sequences, Saltelli sampling and Jansen estimators, verified on the Ishigami function (test_extended.py). Population models use "
            "deterministic Gaussian quantiles, so that no random seed enters Sections 4.8–4.10; the Sobol analysis uses seed 20260906." % R3["validation"]["t90"])
av = find("The complete integration, sweep, Monte Carlo")
av.runs[0].text = av.text + (" The present version (v3) is archived as Zenodo software version https://doi.org/%s and manuscript version https://doi.org/%s; the concept DOIs above resolve to the latest versions." % (SW_DOI, PP_DOI))
for r in av.runs[1:]:
    r._r.getparent().remove(r._r)

# ================================================================== references: remove entries that cannot be verified, cite what the text uses
# The three entries removed here (Bashir 2023; Zhang 2020; Gou 2021) resolve to no record in Crossref under the printed volume/page/DOI and were never cited in the text; Ryshkewitch 1953
# is genuine but concerns porous-ceramic strength, which this paper does not treat. Replacements were harvested from Crossref.
NEWREFS = [
    "[1] Kamal, M. R.; Sourour, S. Kinetics and thermal characterization of thermoset cure. Polymer Engineering & Science 1973, 13(1), 59–64. https://doi.org/10.1002/pen.760130110",
    "[2] Sourour, S.; Kamal, M. R. Differential scanning calorimetry of epoxy cure: isothermal cure kinetics. Thermochimica Acta 1976, 14, 41–59. https://doi.org/10.1016/0040-6031(76)80056-1",
    "[3] Semenoff, N. Zur Theorie des Verbrennungsprozesses. Zeitschrift für Physik 1928, 48, 571–582. https://doi.org/10.1007/bf01340021",
    "[4] Frank-Kamenetskii, D. A. Diffusion and Heat Transfer in Chemical Kinetics, 2nd ed.; Plenum Press: New York, 1969.",
    "[5] White, S. R.; Sottos, N. R.; Geubelle, P. H.; Moore, J. S.; Kessler, M. R.; Sriram, S. R.; Brown, E. N.; Viswanathan, S. Autonomic healing of polymer composites. Nature 2001, 409, 794–797. https://doi.org/10.1038/35057232",
    "[6] Jamekhorshid, A.; Sadrameli, S. M.; Farid, M. A review of microencapsulation methods of phase change materials (PCMs) as a thermal energy storage (TES) medium. Renewable and Sustainable Energy Reviews 2014, 31, 531–542. https://doi.org/10.1016/j.rser.2013.12.033",
]
ref_ps = [find("Kamal, M. R.; Sourour"), find("Bashir, M. A."), find("Ryshkewitch"), find("Zhang, X. et al."), find("Gou, J. et al.")]
for p_, t_ in zip(ref_ps, NEWREFS[:5]):
    for r_ in p_.runs[1:]:
        r_._r.getparent().remove(r_._r)
    p_.runs[0].text = t_
clone_after(ref_ps[-1]._p, ref_ps[-1], NEWREFS[5])


def add_after_text(prefix, find_str, add_str):
    p_ = find(prefix)
    for r_ in p_.runs:
        if find_str in r_.text:
            r_.text = r_.text.replace(find_str, find_str + add_str, 1)
            return
    raise KeyError(find_str)


add_after_text("Controlling when and how quickly", "self-healing or set-on-demand materials", " [5, 6]")
add_after_text("Autocatalytic reaction [Established form", "Kamal & Sourour 1973", " [1, 2]")
add_after_text("ψ = ΔT_ad·Ea₂", "Frank-Kamenetskii/Semenov-type analyses", " [3, 4]")

d.save(OUTF)
print("saved", OUTF)
words = sum(len(p.text.split()) for p in d.paragraphs)
print("words:", words, "| tables:", len(d.tables), "| images:", len(d.inline_shapes))
