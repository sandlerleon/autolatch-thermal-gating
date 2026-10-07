# -*- coding: utf-8 -*-
"""Open new-version drafts of the AutoLatch code and paper records and pre-reserve their DOIs, so the DOIs can be written into the
manuscript before anything is published. State (no token) goes to C:\YouTube\_autolatch_zenodo_state.json. Token from ZENODO_TOKEN only.

    python zenodo_reserve_v3.py
"""
import json
import os
import urllib.request

TOKEN = os.environ["ZENODO_TOKEN"]
STATE = os.path.join("C:" + os.sep, "YouTube", "_autolatch_zenodo_state.json")
API = "https://zenodo.org/api"
CONCEPTS = {"code": 22073392, "paper": 22073390}
SUFFIX = os.environ.get("KEY_SUFFIX", "")
st = json.load(open(STATE)) if os.path.exists(STATE) else {}


def req(method, url):
    r = urllib.request.Request(url, method=method, headers={"Authorization": "Bearer " + TOKEN})
    with urllib.request.urlopen(r, timeout=60) as resp:
        return json.load(resp)


for kind0, concept in CONCEPTS.items():
    kind = kind0 + SUFFIX
    if kind in st:
        print("already reserved:", kind, st[kind]["doi"])
        continue
    latest = req("GET", "%s/records/%d/versions/latest" % (API, concept))
    dep = req("POST", "%s/deposit/depositions/%s/actions/newversion" % (API, latest["id"]))
    draft = req("GET", dep["links"]["latest_draft"])
    st[kind] = {"id": draft["id"], "doi": draft["metadata"]["prereserve_doi"]["doi"], "bucket": draft["links"]["bucket"],
                "inherited_files": [f["id"] for f in draft.get("files", [])], "parent": latest["id"], "concept_doi": "10.5281/zenodo.%d" % concept}
    json.dump(st, open(STATE, "w"), indent=1)
    print(kind, st[kind]["id"], st[kind]["doi"], "inherited files:", len(st[kind]["inherited_files"]))
