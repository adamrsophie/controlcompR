"""Compute reference values with python-control and slycot.

Run from the repository root:
    python benchmarks/reference/reference.py
It writes benchmarks/reference/reference-values.json, which
benchmarks/reference/compare.R checks against controlcompR.
"""
import json
import pathlib

import control
import numpy as np

here = pathlib.Path(__file__).parent
problems = json.loads((here / "problems.json").read_text())
out = {"versions": {"control": control.__version__, "numpy": np.__version__}}

def arr(x):
    return np.array(x, dtype=float)

def tolist(x):
    x = np.asarray(x)
    if np.iscomplexobj(x):
        return {"re": x.real.tolist(), "im": x.imag.tolist()}
    return x.tolist()

out["care"] = [tolist(control.care(arr(p["A"]), arr(p["B"]), arr(p["Q"]), arr(p["R"]))[0])
               for p in problems["care"]]
out["dare"] = [tolist(control.dare(arr(p["A"]), arr(p["B"]), arr(p["Q"]), arr(p["R"]))[0])
               for p in problems["dare"]]
# python-control solves A X + X A' + Q = 0, the convention of ctrl_lyap().
out["lyap"] = [tolist(control.lyap(arr(p["A"]), arr(p["Q"]))) for p in problems["lyap"]]
out["dlyap"] = [tolist(control.dlyap(arr(p["A"]), arr(p["Q"]))) for p in problems["dlyap"]]

def ss(p):
    return control.ss(arr(p["A"]), arr(p["B"]), arr(p["C"]), arr(p["D"]))

out["hinfnorm"] = [float(control.norm(ss(p), "inf")) for p in problems["hinfnorm"]]
out["h2norm"] = [float(control.norm(ss(p), 2)) for p in problems["h2norm"]]
out["hsv"] = [tolist(control.hsvd(ss(p))) for p in problems["hsv"]]
c2d = []
for p in problems["c2d"]:
    d = control.c2d(ss(p), p["dt"], method="zoh")
    c2d.append({"A": tolist(d.A), "B": tolist(d.B)})
out["c2d"] = c2d
place = []
for p in problems["place"]:
    K = control.place(arr(p["A"]), arr(p["B"]), p["poles"])
    place.append(sorted(np.linalg.eigvals(arr(p["A"]) - arr(p["B"]) @ K).real.tolist()))
out["place"] = place
hinf = []
for p in problems["hinfsyn"]:
    try:
        K, CL, gam, rcond = control.hinfsyn(ss(p), p["nmeas"], p["ncon"])
        hinf.append({"gamma": float(gam),
                     "closed_loop_norm": float(control.norm(CL, "inf"))})
    except Exception as e:  # slycot missing or synthesis failure
        hinf.append({"error": str(e)})
out["hinfsyn"] = hinf

(here / "reference-values.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out["versions"]))
