# -*- coding: utf-8 -*-
"""Publish the new-version drafts reserved by zenodo_reserve_v3.py (code and paper). Existing metadata of the previous version is kept; version,
publication date and a version note are updated. Token from ZENODO_TOKEN only.

    python zenodo_publish_v3.py code|paper [--dry]
"""
import json
import os
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request

TOKEN = os.environ["ZENODO_TOKEN"]
API = "https://zenodo.org/api"
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
STATE = os.path.join("C:" + os.sep, "YouTube", "_autolatch_zenodo_state.json")
TAG, DRY = "v3.0.0", "--dry" in sys.argv
NOTE = ("<p><strong>Version 3.</strong> Leaves every v2 number unchanged and adds four analyses that test the v2 conclusions: a capsule-population model in which the gate width "
        "emerges from a spread of melting temperatures (w_eff ~ 0.55 sigma); a storage-stability objective with five tuned control families (Arrhenius, same energy balance) to test the "
        "Pareto claim; a design rule for the admissible melting-point spread; and Sobol variance attribution. The v2 constant-rate first-order control is superseded for the non-dominance "
        "claim. A status-of-claims table and two further falsifiable predictions are added. All results are model outputs under illustrative parameters.</p>")


def req(method, url, data=None, raw=None, headers=None):
    h = {"Authorization": "Bearer " + TOKEN}
    h.update(headers or {})
    body = raw if raw is not None else (json.dumps(data).encode() if data is not None else None)
    if data is not None:
        h["Content-Type"] = "application/json"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, data=body, headers=h, method=method), timeout=600) as r:
            t = r.read()
            return json.loads(t) if t else {}
    except urllib.error.HTTPError as e:
        raise SystemExit("%s %s -> %s %s" % (method, url, e.code, e.read().decode()[:600]))


def run(kind):
    d = json.load(open(STATE))[kind]
    dep = req("GET", "%s/deposit/depositions/%s" % (API, d["id"]))
    for f in dep.get("files", []):
        try:
            req("DELETE", "%s/deposit/depositions/%s/files/%s" % (API, d["id"], f["id"]))
        except SystemExit:
            pass
    if kind == "code":
        tmp = os.path.join(os.environ.get("TEMP", "."), "autolatch-thermal-gating-%s.zip" % TAG)
        subprocess.check_call(["git", "-C", REPO, "archive", "--format=zip", "--prefix=autolatch-thermal-gating-%s/" % TAG, "-o", tmp, TAG])
        files = [tmp]
    else:
        files = [os.path.join(REPO, "manuscript", "AutoLatch_ChemRxiv_WileyMTS_Manuscript_v3.docx"), os.path.join(REPO, "manuscript", "AutoLatch_ChemRxiv_WileyMTS_Manuscript_v3.pdf")]
    for f in files:
        with open(f, "rb") as fh:
            req("PUT", "%s/%s" % (d["bucket"], urllib.parse.quote(os.path.basename(f))), raw=fh.read(), headers={"Content-Type": "application/octet-stream"})
        print("   uploaded", os.path.basename(f), "%.0f kB" % (os.path.getsize(f) / 1024))
    meta = dep["metadata"]
    meta["version"] = "3.0.0"
    meta["publication_date"] = "2026-10-06"
    if "Version 3." not in meta.get("description", ""):
        meta["description"] = NOTE + meta.get("description", "")
    meta["prereserve_doi"] = {"doi": d["doi"]}
    req("PUT", "%s/deposit/depositions/%s" % (API, d["id"]), data={"metadata": meta})
    print("   metadata written")
    if DRY:
        print("   DRY RUN - left unpublished")
        return
    pub = req("POST", "%s/deposit/depositions/%s/actions/publish" % (API, d["id"]))
    print("   PUBLISHED", pub.get("doi"), "concept", pub.get("conceptdoi"))


if __name__ == "__main__":
    run([a for a in sys.argv[1:] if not a.startswith("--")][0])
