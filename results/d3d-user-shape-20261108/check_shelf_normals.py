"""Cross-check shelf incidence with dense, finite-difference tangents on fixed geometries."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from d3d_optimization.tangent_geometry import Curve
from d3d_optimization.user_input import load_user_input

POINTS = 262145


def scan(field, endpoints, coefficients):
    u = np.linspace(0., 1., POINTS)
    rz = Curve(*endpoints, coefficients).point(u)
    tangent = np.gradient(rz, u, axis=0, edge_order=2)
    normal = np.column_stack([-tangent[:, 1], tangent[:, 0]])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    angles = field.incidence(rz, normal)
    i = int(np.argmin(angles))
    return {"minimum_deg": float(angles[i]), "u": float(u[i]), "RZ_m": rz[i].tolist()}


def check(study):
    report = {}
    for case in ("e806f47778", "1a3155370c"):
        source = study/case
        full = json.loads((source/"full"/"result.json").read_text())
        provenance = full["inputs"]
        user = load_user_input(provenance["equilibrium_path"], provenance["time_ms"])
        if user.equilibrium_sha256 != provenance["equilibrium_sha256"]:
            raise ValueError("equilibrium hash mismatch")
        endpoints = provenance["omp_endpoints_RZ_m"]["shelf"]
        cases = [("original", None), ("full", full)]
        other = [("floor-only", source/"floor-only")]
        if case == "e806f47778":
            other.append(("entry-demo", study/"entry-demo"))
        cases.extend((label, json.loads((run/"result.json").read_text())) for label, run in other)
        report[case] = {}
        for method in ("cubic", "quintic"):
            field, _, _ = user.field(method)
            out = {}
            for label, result in cases:
                vector = result["candidate"]["vector"][4:8] if result else [0.]*4
                out[label] = scan(field, endpoints, vector)
                if result:
                    out[label]["candidate_id"] = result["candidate"]["id"]
                if out[label]["minimum_deg"] < full["fixed_shelf"]["floor_deg"]:
                    raise ValueError(f"dense shelf check failed: {case}, {method}, {label}")
            report[case][method] = out
    return dict(points_per_curve=POINTS, normal="finite-difference derivative of R,Z; second-order endpoint formula",
                results=report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.study), indent=2, allow_nan=False))
