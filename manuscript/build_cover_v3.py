# -*- coding: utf-8 -*-
"""Cover letter v3: the v2 letter with the date updated, the Pareto sentence qualified, and a paragraph on the four extended analyses."""
import copy
import json
import os

import docx

HERE = os.path.dirname(os.path.abspath(__file__))
R3 = json.load(open(os.path.join(HERE, "..", "results_v3.json"), encoding="utf-8"))
ZEN = json.load(open(os.path.join("C:" + os.sep, "YouTube", "_autolatch_zenodo_state.json")))
d = docx.Document(os.path.join(HERE, "AutoLatch_WileyMTS_Cover_Letter_Sandler.docx"))
CF = R3["control_families"]
others = [k for k in CF if not k.startswith("F5")]
n_phys = sum(CF[k]["n_physical"] for k in others)
dr = [r for r in R3["design_rule"] if r["margin"] == 6.0 and r["eps"] == 0.01][0]


def setp(p, text):
    for h in p._p.findall("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}hyperlink"):
        p._p.remove(h)
    for r in p.runs[1:]:
        r._r.getparent().remove(r._r)
    p.runs[0].text = text


for p in d.paragraphs:
    t = p.text
    if t.startswith("September 8, 2026"):
        setp(p, "October 6, 2026")
    elif t.startswith("The manuscript formulates and computationally analyzes"):
        setp(p, t.replace("a full 60-point Pareto sweep confirms the coupled system is non-dominated — no sampled first-order control matches or beats it simultaneously on both premature conversion and completion time.",
                          "a full 60-point Pareto sweep shows the coupled system is non-dominated by first-order controls on premature conversion and completion time, and the extended analysis described below tests this claim against stronger controls and a storage-stability objective."))
        new = copy.deepcopy(p._p)
        p._p.addnext(new)
        q = docx.text.paragraph.Paragraph(new, p._parent)
        setp(q, "This version adds four analyses, each designed to test a conclusion rather than to enlarge the model. A capsule-population model derives the gate width from a distribution of "
                "melting temperatures (w_eff ≈ 0.55σ), so that the leakage prediction acquires a measurable input, the width of the DSC melting endotherm. A storage-stability objective, which the "
                "earlier draft did not measure, is added and the constant-rate first-order control is replaced by Arrhenius controls that share the energy balance: none of %d admissible designs of four other "
                "mechanisms dominates the coupled system, whereas some detuned designs of the coupled family itself dominate its illustrative reference point, which is stated explicitly. A design rule "
                "gives the largest admissible melting-point spread (about %.0f °C for 1 %% storage conversion at a 6 °C margin), and Sobol indices confirm the dominance of the kinetic activation energy. "
                "A table lists every claim with its evidential status, the assumption on which it rests and what would reverse it, and two further falsifiable predictions are added. All numbers are model outputs under illustrative parameters; "
                "no experimental data are presented." % (n_phys, dr["sigma_star"]))
    elif t.startswith("The manuscript is explicit throughout"):
        setp(p, t.replace("three falsifiable predictions", "five falsifiable predictions"))
    elif t.startswith("The complete simulation code"):
        setp(p, t + " The present version is archived as https://doi.org/%s (code) and https://doi.org/%s (manuscript)." % (ZEN["code"]["doi"], ZEN["paper"]["doi"]))
out = os.path.join(HERE, "AutoLatch_WileyMTS_Cover_Letter_Sandler_v3.docx")
d.save(out)
print("saved", out)
