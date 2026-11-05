"""First-hit sensitivity of the negative OMP label tail on two downward private-flux legs.

Run as a module from the repository root. Uses archived positive-label power, without
renormalization. The bounds describe only this finite-ring, prescribed-release model.
"""
import argparse
import csv
import json
import os
from itertools import combinations
from pathlib import Path
from time import perf_counter

import numpy as np
import shapely
from scipy.optimize import brentq, linprog, minimize
from shapely.geometry import LineString, Polygon

from d3d_optimization.broadening_study import digest, load_json
from d3d_optimization.case import load_case, separatrix
from d3d_optimization.frozen_mapping import first_hits, trace_paths
from verification.nf_preimpact_comparison import Flux, ROOT, ranges
from verification.preimpact_routing import GaussianMass, Wall, merged, ring_grid, scores

DATA = Path("D:/DivertorOptimization")
STUDY = DATA / "studies/preimpact-firstwall-20261105"
SOURCE = DATA / "studies/firefly-baseline-20261105/OMP.csv"
EVIDENCE = DATA / "analyses/submission-adversarial-20261105/audit.json"
SHOT = ROOT / "input/D3D/d3d/d3d_shot_1a3155370c.parquet"
GEOMETRY = ROOT / "input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json"
CUTS = (-1.285, -1.305)


class Leg:
    def __init__(self, z, lo, hi, flux, inner, limit_is_sep):
        self.z, self.lo, self.hi, self.flux = z, lo, hi, flux
        self.inner = inner
        sample = np.linspace(lo, hi, 257)
        vals = flux.values(np.c_[sample, np.full(257, z)])
        direction = -1 if inner else 1
        assert np.all(direction * np.diff(vals) > 0)
        # At fixed Z the cubic spline is piecewise cubic in R. Check the extrema
        # of its quadratic radial derivative between every neighboring knot pair.
        knots = flux.field.spline.get_knots()[1]
        edges = np.unique(np.r_[lo, knots[(knots > lo) & (knots < hi)], hi])
        for a, b in zip(edges, edges[1:]):
            radii = np.array([a, (a + b) / 2, b])
            slopes = direction * flux.gradient(np.c_[radii, [z]*3])[:, 0]
            quad, linear, c = np.polyfit([0., .5, 1.], slopes, 2)
            minimum = min(slopes)
            if quad > 0 and 0 < -linear/(2*quad) < 1:
                minimum = min(minimum, c - linear**2/(4*quad))
            assert minimum > -1e-11, (z, a, b, minimum)
        self.umin = float(min(vals[0], vals[-1]))
        self.umax = 0. if limit_is_sep else float(max(vals[0], vals[-1]))
        assert self.umin < self.umax <= 0

    def starts(self, labels):
        roots = [brentq(lambda r: self.flux.values([[r, self.z]])[0] - u,
                        self.lo, self.hi, xtol=6e-15) for u in labels]
        return np.c_[roots, np.full(len(roots), self.z)]


def make_legs(wall, xp, z):
    flux = wall.flux
    assert z < xp[1]
    crit = brentq(lambda r: flux.gradient([[r, z]])[0, 0], xp[0]-.12, xp[0]+.1)
    val = lambda r: flux.values([[r, z]])[0]
    left = brentq(val, max(flux.field.r[0], xp[0]-.25), crit, xtol=6e-15)
    right = brentq(val, crit, xp[0]+.19, xtol=6e-15)
    line = LineString([(flux.field.r[0], z), (flux.field.r[-1], z)])
    poly = Polygon(wall.p)
    assert poly.is_valid
    domain = poly.intersection(line)
    parts = domain.geoms if hasattr(domain, "geoms") else [domain]
    interval, = [part.bounds for part in parts if part.bounds[0] < crit < part.bounds[2]]
    inner = Leg(z, max(left, interval[0]), crit, flux, True, left >= interval[0])
    outer = Leg(z, crit, min(right, interval[2]), flux, False, right <= interval[2])
    info = dict(critical_R_m=crit, magnetic_separatrix_R_m=[left,right],
                interior_R_m=[interval[0],interval[2]], inner_label_m=[inner.umin,inner.umax],
                outer_label_m=[outer.umin,outer.umax])
    return (inner, outer), info


def routed_pieces(grid, good):
    pieces = []
    for a, b, arc in merged(good):
        lo, hi = np.maximum(a, grid["lo"]), np.minimum(b, grid["hi"])
        valid = (grid["arc"] == arc) & (hi > lo)
        pieces.extend(zip(lo[valid], hi[valid], grid["ring"][valid]))
    pieces.sort()
    for left, right in zip(pieces, pieces[1:]):
        assert left[1] <= right[0]+2e-14
    return pieces


def bucket_table(routes, ends):
    edges = np.unique(np.r_[ends, [u for parts in routes for a,b,_ in parts for u in (a,b)]])
    mid = (edges[:-1]+edges[1:])/2
    bins = []
    for pieces in routes:
        dst = np.full(len(mid), -1, int)
        for lo, hi, ring in pieces:
            ok = (mid >= lo) & (mid < hi)
            assert np.all(dst[ok] == -1)
            dst[ok] = ring
        bins.append(dst)
    return edges, bins[0], bins[1]


def bin_mass(rings, mass, n):
    return np.bincount(rings, weights=mass, minlength=n)


def write_csv(file, rows):
    with file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def extremize(fixed, area, inner, outer, mass):
    """Sharp peak and A2 on common atoms, with a convex optimum certificate for max A2."""
    n = len(area)
    in_add, out_add = (bin_mass(idx, mass, n) for idx in (inner, outer))
    # Here the two destinations illuminate disjoint material rings. Their separate
    # depositions give a sharp lower and upper power for every individual ring.
    assert not np.intersect1d(inner, outer).size
    all_out = fixed + out_add
    amplitude = in_add - out_add
    total = all_out.sum()
    touched = np.unique(np.r_[inner, outer])
    q0, slope = all_out/area, amplitude/area
    # Peak of affine functions. All minimum candidates are endpoints or intersections.
    hot = np.argmax(q0)
    lines = np.unique(np.r_[touched, hot, np.argmax(fixed/area)])
    crossings = [0., 1.]
    for i,j in combinations(lines, 2):
        if slope[i] != slope[j]:
            t = (q0[j]-q0[i])/(slope[i]-slope[j])
            if 0 < t < 1:
                crossings.append(t)
    uniform_min_peak = min(float(np.max(q0 + t*slope)) for t in crossings)
    d_coef = np.sum(amplitude**2 / area)
    t_opt = float(np.clip(-np.dot(q0, amplitude)/d_coef, 0, 1))
    samples = {str(t):scores(all_out+t*amplitude, area) for t in (0., .25, .5, .75, 1., t_opt)}
    uniform = dict(samples=samples, inner_fraction_at_max_A2=t_opt,
                   peak=[uniform_min_peak, max(samples["0.0"]["peak"],samples["1.0"]["peak"])],
                   A2=[min(samples["0.0"]["A2"],samples["1.0"]["A2"]),samples[str(t_opt)]["A2"]])

    def profile(x):
        return fixed+bin_mass(inner, mass*x, n)+bin_mass(outer, mass*(1-x), n)

    def objective(x):
        power = profile(x)
        q = power/area
        return np.dot(power, q), 2*mass*(q[inner]-q[outer])

    # Use the variable-dependent portion to avoid stopping from constant-background
    # rounding. Verify quality from the global convex tangent bound on the box.
    def variable_objective(x):
        p = profile(x)[touched]
        q = np.zeros_like(area)
        q[touched] = p/area[touched]
        return np.dot(p, q[touched])*1e4, 2e4*mass*(q[inner]-q[outer])

    sol = minimize(variable_objective, np.full(len(mass),.5), jac=True, method="L-BFGS-B",
                   bounds=[(0,1)]*len(mass), options=dict(ftol=1e-15, gtol=1e-11, maxiter=4000))
    best, grad = objective(sol.x)
    gap = float(np.sum(np.maximum(grad*sol.x, -grad*(1-sol.x))))
    d_min_bound = best-gap
    assert gap < 1e-8 and d_min_bound > 0

    # Convex D attains a maximum at a vertex of the box. These fixed shapes have
    # few common atoms, so enumerate every binary split in memory-bounded batches.
    # Restrict calculations to touched rings; the other terms are constant.
    transfer = ((touched[:,None] == inner).astype(float) - (touched[:,None] == outer))*mass
    assert len(mass) <= 20, "exhaustive vertex search is restricted to the studied shapes"
    states = 1 << len(mass)
    best_variable_D = -np.inf
    for first in range(0,states,1024):
        ids = np.arange(first,min(first+1024,states),dtype=np.uint64)
        flags = ((ids[:,None] >> np.arange(len(mass),dtype=np.uint64)) & 1).astype(float)
        local_power = all_out[touched] + flags @ transfer.T
        D = np.sum(local_power**2/area[touched],axis=1)
        index = np.argmax(D)
        if D[index] > best_variable_D:
            best_variable_D = D[index]
            binary = flags[index].copy()
    d_max = objective(binary)[0]

    ring_upper = fixed+in_add+out_add
    upper_peak = float(np.max(ring_upper/area))  # Each such ring bound is attainable.
    fixed_peak = float(np.max(fixed/area))
    lower_witness = min(samples["0.0"]["peak"],samples["1.0"]["peak"])
    if abs(lower_witness-fixed_peak) < 1e-12:
        lower_peak = fixed_peak
    else:
        # All untouched constraints reduce to one constant. Use fraction variables.
        mat = transfer / area[touched,None]
        lp = linprog(np.r_[np.zeros(len(mass)),1],A_ub=np.c_[mat,-np.ones(len(touched))],
                     b_ub=-q0[touched],bounds=[(0.,1.)]*len(mass)+[(fixed_peak,None)],method="highs")
        assert lp.success
        lower_peak = float(lp.fun)
        np.testing.assert_allclose(scores(profile(lp.x[:-1]),area)["peak"],lower_peak,atol=1e-8)
    return (uniform,dict(peak=[lower_peak,upper_peak],A2=[total**2/d_max,total**2/d_min_bound],
                         attained_max_A2=total**2/best, exact_binary_vertices_checked=states, min_A2_inner_MW=float(np.dot(mass,binary)),
                         max_A2_inner_MW=float(np.dot(mass,sol.x)), max_A2_D_tangent_gap=gap),
            ring_upper, d_max, d_min_bound, binary, sol.x)


def residual_bounds(fixed, area, ring_upper, missing, Dmax, Dmin, sharp):
    """Any fraction of residual may reach any ring (arbitrary further transport)."""
    qhi=ring_upper/area
    return dict(peak=[sharp["peak"][0],float(np.max((ring_upper+missing)/area))],
                A2=[float(fixed.sum()**2/(Dmax+missing*np.max(2*qhi+missing/area))),float(1/Dmin)])


def surface_ranges(wall):
    original = np.array(load_json(GEOMETRY)["limiter"])
    def vertex(point):
        ids = np.flatnonzero(np.linalg.norm(wall.p-point,axis=1) < 1e-11)
        assert len(ids)==1
        return int(ids[0])
    return {j: (vertex(original[j]),vertex(original[j+1])) for j in (68,69)}


def ray_control(section, a, b, mass, wall, grid, surface, good, count, step):
    """Midpoint quadrature of exact masses, with no use of wall labels to find hits."""
    # Cosine auxiliary mesh resolves the fold; quantile edges lower histogram noise.
    t = np.linspace(0,1,1025)
    mesh = a+(b-a)*(.5-.5*np.cos(np.pi*t))
    cum = np.r_[0,np.cumsum(mass(mesh[:-1],mesh[1:]))]
    edges=np.interp(np.linspace(0,cum[-1],count+1),cum,mesh)
    edges[[0,-1]] = a,b
    mids=(edges[:-1]+edges[1:])/2
    starts=section.starts(mids)
    powers=mass(edges[:-1],edges[1:])
    poly=Polygon(wall.p)
    assert shapely.contains_xy(poly,starts[:,0],starts[:,1]).all()
    bsign=np.where(section.flux.field.magnetic(starts)[:,1] < 0,1,-1)
    assert np.all(bsign == (1 if section.inner else -1))
    points=[]; segments=[]; distances=[]
    for chunk in np.array_split(starts,(len(starts)+95)//96):
        rays=trace_paths(section.flux.field,chunk,step=step,max_length=1.)
        assert all(path[1,1] < path[0,1] for path in rays)
        hit=first_hits(rays,wall.p)
        assert np.all(hit.segment>=0)
        points.extend(hit.points);segments.extend(hit.segment);distances.extend(hit.distance)
    points,segments=np.asarray(points),np.asarray(segments)
    which = 69 if section.inner else 68
    start,end = surface[which]
    assert ((segments>=start)&(segments<end)).all()
    s=wall.hit_s(points,segments);arcs=wall.arc_id(s)
    for low,high,j in good:
        assert (arcs[(mids>=low)&(mids<high)]==j).all()
    actual = bin_mass((np.floor(s*len(grid["area"])/wall.s[-1]).astype(int))%len(grid["area"]),powers,len(grid["area"]))
    return actual,dict(rays=count,step_m=step,along_B_sign=int(bsign[0]),surface_original_segment=which,
                       hit_RZ_min_m=points.min(axis=0).tolist(),hit_RZ_max_m=points.max(axis=0).tolist(),
                       weighted_hit_RZ_m=(powers@points/powers.sum()).tolist(),
                       distance_to_hit_range_m=[min(distances),max(distances)],
                       hit_segment_range=[int(segments.min()),int(segments.max())],
                       max_hit_flux_label_error_m=float(np.max(abs(section.flux.values(points)-mids))))


def run_case(case, model, table, flux, xp, out, source_hashes):
    t=perf_counter()
    assert digest(model["path"]) == model["sha256"]
    with np.load(model["path"]) as raw:
        wall=Wall(raw["wall"],flux)
        np.testing.assert_allclose(raw["donor_region_power"].sum(axis=1),table["power_fraction"],atol=1e-16)
        np.testing.assert_allclose(raw["source_x_m"],table["x_mid_m"],atol=1e-14,rtol=0)
        archived_hit_s=wall.hit_s(raw["hit_points"],raw["hit_segment"])
        assert np.all(raw["hit_segment"]>=0)
        centers=flux.values(np.c_[table["R_m"],table["Z_m"]])
        mass=GaussianMass(centers,raw["thermal_sigma_m"],table["power_fraction"])
    grid=ring_grid(wall,.002,0);area=grid["area"]
    per_case=STUDY/(case+"_fine")
    def inherited(rel):
        path=per_case/rel
        assert digest(path)==source_hashes[str(path)]
        return path
    with np.load(inherited("omp_2mm_0.npz")) as old:
        omp=old["resolved_MW"]
        np.testing.assert_array_equal(area,old["area"])
    previous=load_json(inherited("result.json"))
    previous_ranges=load_json(inherited("routes.json"))
    assert previous["archived_hit_arc_matches"]==len(table) and not previous["mismatch_intervals"]
    positive_upper=previous_ranges["omp"]["assigned"][-1][1]
    positive_unresolved=float(mass(positive_upper,np.inf))
    input_total=table["power_W"].sum()/1e6
    assert len(table)==1280 and np.all(table["sign_along_B"]==-1)
    assert np.all(centers>0) and np.all(centers<positive_upper)
    zero=bin_mass((np.floor(archived_hit_s*len(area)/wall.s[-1]).astype(int))%len(area),
                  table["power_fraction"],len(area))
    assert abs(zero.sum()-input_total)<2e-13
    surface=surface_ranges(wall)
    info=dict(input_sha256=model["sha256"],input_MW=input_total,OMP_resolved_MW=float(omp.sum()),
              negative_MW=float(mass(-np.inf,0)),positive_high_unresolved_MW=positive_unresolved,
              zero_width=dict(negative_MW=0,unresolved_MW=0,source_MW=float(table["power_fraction"].sum()),
                              direct_wall_MW=float(zero.sum()),archived_first_hit_2mm=scores(zero,area)),
              cuts={})
    rows=[]
    for z in CUTS:
        legs,record=make_legs(wall,xp,z)
        checked=[]; all_pieces=[]
        for section in legs:
            good,bad=ranges(wall,section,section.umin,section.umax,.0005,True)
            refined,refined_bad=ranges(wall,section,section.umin,section.umax,.00025,False)
            assert merged(good)==merged(refined) and bad==refined_bad
            # No positive-power arc ambiguity is accepted by this run.
            assert not bad, (case,z,bad)
            checked.append(merged(good))
            all_pieces.append(routed_pieces(grid,good))
        span=[min(s.umin for s in legs),0.]
        edges,idx_in,idx_out=bucket_table(all_pieces,span)
        weights=mass(edges[:-1],edges[1:])
        common=(idx_in>=0)&(idx_out>=0)
        alone=(idx_in<0)&(idx_out>=0)
        assert not ((idx_in>=0)&(idx_out<0)).any()
        insufficient=(idx_in<0)&(idx_out<0)
        unresolved=float(mass(-np.inf,span[0])+positive_unresolved+weights[insufficient].sum())
        record.update(routes={"inner":checked[0],"outer":checked[1]}, common_MW=float(weights[common].sum()),
                      outer_only_MW=float(weights[alone].sum()),uncovered_internal_MW=float(weights[insufficient].sum()),
                      deeper_unresolved_MW=float(mass(-np.inf,span[0])),unresolved_MW=unresolved)
        assert np.any(common) and weights[common].sum()>0
        m=weights[common]; ii=idx_in[common]; oo=idx_out[common]
        only=bin_mass(idx_out[alone],weights[alone],len(area))
        inner=bin_mass(ii,m,len(area));outer=bin_mass(oo,m,len(area))
        inputs={"fixed_single_leg":(omp+only,unresolved)}
        if np.any(alone):
            inputs["single_leg_unassigned"]=(omp,unresolved+weights[alone].sum())
        arrays=dict(s_m=grid["s_mid"],area_m2=area,rz_m=wall.points_at(grid["s_mid"]),OMP_MW=omp,zero_width_OMP_MW=zero,
                    outer_only_MW=only, inner_common_MW=inner,outer_common_MW=outer,
                    atom_edges_m=edges,atom_mass_MW=weights,atom_inner_ring=idx_in,atom_outer_ring=idx_out)
        record["scenarios"]={}
        for key,(fixed,missing) in inputs.items():
            uniform,free,ring_hi,dmax,dmin,binary,x=extremize(fixed,area,ii,oo,m)
            power=fixed+outer
            balance=float(power.sum()+missing-input_total)
            assert abs(balance)<3e-12
            row=dict(uniform=uniform,free_label=free,arbitrary_residual_bounds=residual_bounds(
                power,area,ring_hi,missing,dmax,dmin,free),resolved_MW=float(power.sum()),
                unresolved_MW=missing,balance_error_MW=balance)
            record["scenarios"][key]=row
            arrays.update({key+"_min_A2_inner_fraction":binary,key+"_max_A2_inner_fraction":x,
                           key+"_ring_lower_MW":fixed,key+"_ring_upper_MW":ring_hi})
            for label,s in uniform["samples"].items():
                a=float(label)
                placed=fixed+(1-a)*outer+a*inner
                imax=np.argmax(placed/area)
                rows.append(dict(case=case,cut_Z_m=z,single_leg_policy=key,inner_common_fraction=a,
                                 inner_common_MW=a*m.sum(),outer_common_MW=(1-a)*m.sum(),resolved_MW=power.sum(),
                                 unresolved_MW=missing,peak_MW_m2=s["peak"],A2_m2=s["A2"],
                                 peak_R_m=arrays["rz_m"][imax,0],peak_Z_m=arrays["rz_m"][imax,1]))
        # Reproduce the full archived outboard endpoint exactly before replacing its routing.
        with np.load(inherited(f"leg{z:g}_2mm_0.npz")) as old:
            record["archived_outboard_max_ring_delta_MW"]=float(np.max(abs(old["resolved_MW"]-omp-only-outer)))
            assert record["archived_outboard_max_ring_delta_MW"] < 5e-12
        # Zero width leaves original OMP donors only and no power on either hypothetical leg.
        a,b=max(s.umin for s in legs),min(s.umax for s in legs)
        assert abs(mass(a,b)-m.sum())<2e-13
        record["direct"]={}
        for side,section,j,ref in zip(("inner","outer"),legs,(0,1),(inner,outer)):
            controls={}
            for n,step in ((256,.0005),(1024,.00025)):
                direct,c=ray_control(section,a,b,mass,wall,grid,surface,checked[j],n,step)
                fixed=omp+only
                stats=scores(fixed+direct,area);target=scores(fixed+ref,area)
                c.update(TV_MW=float(abs(direct-ref).sum()/2),ledger_error_MW=float(direct.sum()-m.sum()),
                         peak_relative_delta=stats["peak"]/target["peak"]-1,A2_relative_delta=stats["A2"]/target["A2"]-1,
                         max_ring_delta_MW=float(np.max(abs(direct-ref))),private_only_peak_relative_delta=float(
                             np.max(direct/area)/np.max(ref/area)-1))
                assert abs(c["ledger_error_MW"]) < 2e-13
                controls[str(n)]=c
                arrays[f"direct_{side}_{n}_MW"]=direct
            record["direct"][side]=controls
        # All reported ranges for a geometry are generated from its own wall.
        np.savez_compressed(out/f"leg-{case}-z{z:g}.npz",**arrays)
        info["cuts"][str(z)]=record
        print(case,z,"common MW",float(m.sum()),"free",record["scenarios"]["fixed_single_leg"]["free_label"],flush=True)
    info["elapsed_s"]=perf_counter()-t
    return info,rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/nf-finalist-robustness-20261106")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    field,lcfs,axis=load_case(SHOT,1400);r,xp,_=separatrix(field,lcfs,axis)
    table=np.genfromtxt(SOURCE,delimiter=",",names=True)
    assert np.isclose(table["power_W"].sum(),1e6,rtol=1e-13) and np.all(table["Z_m"]==axis[1])
    np.testing.assert_allclose(table["power_fraction"],table["power_W"]/1e6,rtol=2e-14)
    prior=load_json(STUDY/"verified.json")["inputs"]
    for p,k in ((SOURCE,"source_sha256"),(EVIDENCE,"audit_sha256"),(SHOT,"shot_sha256")):
        assert digest(p)==prior[k]
    source_hashes=load_json(STUDY/"sha256.json")
    report=dict(inputs=prior,sha256_script=digest(__file__),xpoint_RZ_m=xp.tolist(),geometry_sha256=digest(GEOMETRY),
                grid="2mm_0",cases={},cpu_threads={k:os.getenv(k) for k in
                ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")})
    rows=[]
    for case,entry in load_json(EVIDENCE)["inputs"].items():
        report["cases"][case],rws=run_case(case,entry,table,Flux(field,(r,axis[1])),xp,args.output,source_hashes)
        rows.extend(rws)
    write_csv(args.output/"leg-mixtures.csv",rows)
    (args.output/"leg-results.json").write_text(json.dumps(report,indent=2,allow_nan=False),encoding="utf-8")


if __name__=="__main__":
    main()
