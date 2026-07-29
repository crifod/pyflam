"""Robust fire-name matching: portal slugs vs Zenodo sonde filenames.

Substring containment fails on real pairs -- 'StaColomadeQueralt' against
'santa-coloma-de-queralt' -- which silently dropped the one fire carrying per-sonde
pyroconvection types. Similarity matching with an expansion table handles it.
"""
import difflib, re
EXPAND = {"sta": "santa", "sant": "santa", "d": "de"}
def norm(s):
    s = re.sub(r"[^a-z]", " ", s.lower())
    return "".join(EXPAND.get(w, w) for w in s.split())
def best(name, candidates, cutoff=0.82):
    n = norm(name)
    scored = [(difflib.SequenceMatcher(None, n, norm(c)).ratio(), c) for c in candidates]
    scored += [(1.0, c) for c in candidates if n in norm(c) or norm(c) in n]
    scored.sort(reverse=True)
    return scored[0][1] if scored and scored[0][0] >= cutoff else None
