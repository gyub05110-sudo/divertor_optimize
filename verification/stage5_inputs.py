"""Register the fixed Stage 4 material walls, equilibria and OMP sources for FIREFLY."""
import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import LineString, Polygon

from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.case import separatrix
from d3d_optimization.frozen_field import exponential_source
from d3d_optimization.user_input import GEOMETRY, POLICY, REFERENCE, _parquet
from d3d_optimization.user_output import file_hash
from verification.nf_prospective_inputs import frame_fingerprint, publish, require


STUDIES = Path("D:/DivertorOptimization/studies")
STAGE4 = STUDIES / "nf-stage4-20261109"
DEFAULT = STUDIES / "nf-stage5-firefly-20261110"
CASES = ("3264fdd71f", "7ce879f1b3", "70481becc0", "5d7f95a730", "1255ab8051",
         "d880609264", "f317e77cf7")
FORCED = ("floor-only", "floor-entry", "floor-shelf", "floor-connector")
SOURCE_HASHES = {"inputs/registration.json": "629906a4c416b43c80fc77a97d56ee166c66a439a8f63c61bfab99817e271e36",
                 "expected.json": "a93beab78292dae9ac282272290cf8533d629bcb0fbe7a888957ddfe72354115",
                 "stage4_report.json": "ea51d468e649f5e1cdef8015167c47c923564b0d958da37b1aca5df3e3520b29"}
REFERENCE_SI = dict(power_W=1e6, density_m3=1e19, chi_perp_m2_s=.3, lambda_q_m=.003, gamma=7., kappa0=2000.)
SOURCE_HEADER = "donor_id,x_low_m,x_high_m,x_mid_m,R_m,Z_m,sign_along_B,power_fraction,power_W"
MEASUREMENT_FILES = {"geometry.csv", "boundary.npz", "donors.csv.gz", "profile.csv.gz"}


def requested_roles(case):
    require(case in CASES, "unregistered Stage 5 case")
    roles = [(30713, role) for role in ("baseline", *FORCED, "hotspot", "full")]
    if case == "70481becc0":
        roles.extend([(30757, "hotspot"), (30757, "full")])
    return roles


def array_hash(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def numeric_text(values, **kwargs):
    stream = io.StringIO()
    np.savetxt(stream, values, fmt="%.17g", **kwargs)
    return stream.getvalue().encode("ascii")


def geqdsk(field, axis, psi_sep, lcfs, limiter, label):
    """Historical native FLARE convention: psi[Z,R] Wb/rad and signed F=R Bphi."""
    r, z = field.r, field.z
    require(max(np.max(np.abs(v-np.linspace(v[0],v[-1],len(v)))) for v in (r,z)) <= 2e-7,
            "GEQDSK requires a uniform rectangular grid")
    r0, z0 = axis
    psi_axis = field.values([axis])[0]
    f = r0 * field.magnetic([axis])[0,2]
    out = io.StringIO()
    out.write(f"{label[:48]:48s} 0 {len(r)} {len(z)}\n")

    def rows(numbers):
        flat = np.asarray(numbers).ravel(order="C")
        for i in range(0,len(flat),5):
            out.write("".join(f"{v:16.9E}" for v in flat[i:i+5])+"\n")

    rows([np.ptp(r), np.ptp(z), r0, r[0], (z[-1]+z[0])/2])
    rows([r0, z0, psi_axis, psi_sep, f/r0])
    rows([0., psi_axis, 0., r0, 0.])
    rows([z0, 0., psi_sep, 0., 0.])
    zeros = np.zeros_like(r)
    for a in (np.full_like(r,f), zeros, zeros, zeros, field.psi, zeros):
        rows(a)
    out.write(f"{len(lcfs)} {len(limiter)}\n")
    rows(lcfs)
    rows(limiter)
    return out.getvalue().encode("ascii")


def donor_table(field, rsep, zaxis, wall, width, power, offsets, watts):
    """Check the Stage 4 donor quadrature on the original wall; export the archived watts."""
    ray = LineString([(rsep,zaxis), (field.r[-1],zaxis)])
    intersections = LineString(wall).intersection(ray)
    points = getattr(intersections, "geoms", [intersections])
    outward = [p.x-rsep for p in points if p.geom_type=="Point" and not p.is_empty and p.x>rsep]
    require(outward, "missing outboard midplane wall crossing")
    edges = min(.1, min(outward))*np.linspace(0,1,len(offsets)+1)**2
    source = exponential_source(field,rsep,zaxis,edges,width)
    require(np.allclose(source.x,offsets,rtol=0,atol=1e-15), "registered OMP offset mismatch")
    require(np.allclose(power*source.weights,watts,rtol=2e-12,atol=0) and np.all(watts>0)
            and abs(watts.sum()/power-1)<5e-14, "registered OMP assigned power mismatch")
    rz = np.column_stack([rsep+offsets,np.full(len(offsets),zaxis)])
    require(np.all(shapely.contains_xy(Polygon(wall),*rz.T)), "OMP donor outside material wall")
    bz = field.magnetic(rz)[:,1]
    require(np.all(bz!=0), "undetermined downward source branch")
    return np.column_stack([np.arange(len(offsets)),edges[:-1],edges[1:],offsets,rz,
                            np.where(bz<0,1,-1),watts/power,watts])


def probe_table(field, lcfs, wall, axis, xp, donors):
    r, z = np.meshgrid(np.linspace(field.r[0],field.r[-1],91),np.linspace(field.z[0],field.z[-1],113))
    inside = shapely.contains_xy(Polygon(wall),r,z)
    points = np.vstack([donors, axis, xp, lcfs, np.column_stack([r[inside],z[inside]])])
    a, b = field.values([axis, xp])
    return np.c_[points,field.magnetic(points),(field.values(points)-a)/(b-a)]


def native_cdf(weights):
    cdf = np.cumsum(weights)
    require(np.all(weights>0) and abs(cdf[-1]-1)<5e-14,"invalid prescribed source weights")
    # Preserve valid registered representations. A rounded penultimate sum above
    # one needs endpoint scaling before applying the native exact endpoint.
    if cdf[-2]>1.:
        cdf /= cdf[-1]
    cdf[-1]=1.
    require(np.all(np.diff(cdf)>=0),"nonmonotonic native source CDF")
    return cdf


def metadata_fields(field, lcfs, axis, raw, registration):
    topo = raw["equilibrium_topology"]
    rsep, xp, psi_sep = separatrix(field,lcfs,axis)
    require(np.array_equal(topo["axis_RZ_m"],axis) and np.allclose(topo["xpoint_RZ_m"],xp,atol=1e-12,rtol=0)
            and raw["fixed_wall"]["LCFS_RZ_little_endian_float64_sha256"]==array_hash(lcfs)
            and registration["fingerprint"]==frame_fingerprint(field.psi), "original topology or frame differs")
    endpoints = raw["fixed_shelf"]["endpoints"]
    require(all(np.allclose(field.magnetic([pt["RZ_m"]])[0],pt["B_RZphi_T"],rtol=5e-13,atol=1e-14)
                for pt in endpoints), "signed field differs from archived shelf samples")
    data = dict(axis_RZ_m=axis.tolist(),xpoint_RZ_m=xp.tolist(),omp_separatrix_RZ_m=[rsep,float(axis[1])],
                psi_axis_Wb_per_rad=float(field.values([axis])[0]),psi_sep_Wb_per_rad=psi_sep,
                F_T_m=float(field.magnetic([axis])[0,2]*axis[0]),tf_current_kA=field.tf_current_ka,
                lcfs_vertices=len(lcfs),lcfs_sha256=array_hash(lcfs),grid_shape_ZR=list(field.psi.shape),
                grid_uniform_approx_max_m=[float(np.max(np.abs(v-np.linspace(v[0],v[-1],len(v)))))
                                           for v in (field.r,field.z)])
    return data


def checked_measurement(directory, measurement):
    files = measurement["sha256"]
    require(directory.name==measurement["directory"] and set(files)==MEASUREMENT_FILES and
            all(file_hash(directory/name)==sha for name,sha in files.items()), "original measurement digest differs")
    geometry = np.genfromtxt(directory/"geometry.csv", delimiter=",", names=True, usecols=(0,1,2))
    donors = np.genfromtxt(directory/"donors.csv.gz", delimiter=",", names=True,usecols=(0,1,3,4,7))
    with np.load(directory/"boundary.npz", allow_pickle=False) as saved:
        wall = saved["wall"]
        offsets = donors["source_x_m"]
        hits = np.c_[donors["first_R_m"],donors["first_Z_m"],donors["parallel_length_m"]]
        require(np.array_equal(np.column_stack([geometry["R_m"],geometry["Z_m"]]),wall) and
                np.array_equal(geometry["vertex"],np.arange(len(wall))) and
                np.array_equal(wall[0],wall[-1]) and np.all(np.linalg.norm(np.diff(wall,axis=0),axis=1)>0)
                and Polygon(wall).is_valid, "original material geometry differs")
        require(np.array_equal(offsets,saved["source_x_m"]) and len(donors)==1280
                and np.array_equal(hits[:,:2],saved["hit_points"]) and
                np.array_equal(hits[:,2],saved["thermal_length_m"]) and np.all(np.isfinite(hits))
                and np.all(saved["hit_segment"]>=0) and float(saved["assigned_power_W"])==REFERENCE_SI["power_W"],
                "original donor or first hit reference differs")
    return wall,offsets,donors["assigned_W"],hits


def check_run(study, registration, run, arm, seed):
    case = registration["id"]
    folder = study/"runs"/case/str(seed)/arm
    require(run["relative_output"]==folder.relative_to(study).as_posix() and run["status"]=="qualified"
            and run["receipt_status"]=="completed", "registered winner unavailable")
    for name, key in (("result.json","result_sha256"),("receipt.json","receipt_sha256"),("events.jsonl","events_sha256")):
        require(file_hash(folder/name)==run[key], "timed run digest differs: "+str(folder/name))
    raw = json.loads((folder/"result.json").read_bytes())
    args, source = raw["arguments"], raw["inputs"]
    candidate = raw["candidate"]
    checkpoint = run["checkpoints"]["600"]
    require(args["family"]==arm and args["seed"]==seed and run["case_id"]==case and
            args["time_ms"]==source["time_ms"]==registration["time_ms"]==1400 and
            source["equilibrium_sha256"]==registration["sha256"] and source["equilibrium_path"]==registration["path"]
            and raw["reference_SI"]==REFERENCE_SI and raw["levels"]["fine"]==["cubic",1280,.000125,641,.00025],
            "Stage 4 operating point or source differs")
    require(raw["status"]==run["status"] and raw["candidate"]==run["candidate"]
            and Candidate.from_record(candidate).identity==candidate["id"]==checkpoint["id"]
            and raw["selected"]["candidate_eligible"] and raw["selected"]["physics_measurable"]
            and raw["baseline"]["physics_measurable"] and run["baseline_peak_MW_m2"]==raw["baseline"]["peak_MW_m2"]
            and checkpoint["peak_MW_m2"]==raw["selected"]["peak_MW_m2"], "timed winner identity or peak differs")
    return folder,raw


def register_case(study, out, reg, source, emit):
    case = reg["id"]
    runs = {(r["seed"],r["arm"]):r for r in source["runs"] if r["case_id"]==case}
    common = runs[30713,"hotspot"]
    _, first = check_run(study,reg,common,"hotspot",30713)
    require(file_hash(reg["path"])==reg["sha256"] and reg["status"]=="ready" and reg["exact_frame_count"]==1,
            "raw input registration changed")
    field,lcfs,axis,metadata = _parquet(Path(reg["path"]).read_bytes(),1400)
    provenance = first["inputs"]
    require(reg["frame_index"]==metadata["frame_index"] and all(provenance[k]==v for k,v in metadata.items()),
            "parquet magnetic or frame settings differ from original")
    for path,key in ((GEOMETRY,"geometry"),(POLICY,"design_policy"),(REFERENCE,"design_reference")):
        require(file_hash(path)==provenance[key+"_sha256"],"original fixed input differs: "+key)
    limiter = np.asarray(json.loads(REFERENCE.read_bytes())["wall"])
    material = np.asarray(json.loads(GEOMETRY.read_bytes())["limiter"])
    require(np.array_equal(limiter,np.vstack([material,material[0]])) and array_hash(limiter)==provenance["wall_sha256"],
            "supplied material limiter differs")
    topology = metadata_fields(field,lcfs,axis,first,reg)
    xp = topology["xpoint_RZ_m"]
    rsep,zaxis = topology["omp_separatrix_RZ_m"]
    roles, hashes, donors = [], {}, None

    def put(name, data):
        emit(out/case/name,data)
        hashes[name]=hashlib.sha256(data).hexdigest()

    # Keep the supplied physical limiter in the model. The transport mesh has a separate casing.
    put("model/fixed.geqdsk",geqdsk(field,axis,topology["psi_sep_Wb_per_rad"],lcfs,limiter,f"DIII-D {case} frozen 1400 ms"))
    put("model/.boundary",b"[axisurf]\nfilename: wall.txt\nunits: m\n")
    put("model/.bfield",b"[equi2d_geqdsk]\nfilename: fixed.geqdsk\nspline_order: 3\n")
    put("model/.equi2d",b"1\n"+numeric_text([xp]))
    put("model/wall.txt",numeric_text(limiter))
    for seed,role in requested_roles(case):
        arm = "hotspot" if role=="baseline" else role
        run = runs[seed,arm]
        folder,raw = check_run(study,reg,run,arm,seed)
        which = "baseline" if role=="baseline" else "selected"
        measurement = raw[which]
        directory = folder/which
        require(raw["inputs"]==provenance and raw["equilibrium_topology"]==first["equilibrium_topology"],
                "physical input differs between roles")
        wall,offsets,watts,hits = checked_measurement(directory,measurement)
        wall_id = array_hash(wall)
        wall_file = "walls/"+wall_id+".wall"
        hit_file = "walls/"+wall_id+".hits"
        put(wall_file,numeric_text(wall))
        put(hit_file,numeric_text(hits))
        if donors is None:
            require(LineString(wall).hausdorff_distance(LineString(limiter))<1e-12 and
                    abs(LineString(wall).length-LineString(limiter).length)<1e-10,"baseline material wall differs")
            donors = donor_table(field,rsep,zaxis,limiter,REFERENCE_SI["lambda_q_m"],REFERENCE_SI["power_W"],offsets,watts)
            put("OMP.csv",numeric_text(donors,delimiter=",",header=SOURCE_HEADER,comments=""))
            cdf=native_cdf(donors[:,7])
            put("OMP.source",f"{len(donors)}\n".encode()+numeric_text(np.c_[donors[:,4:6],cdf]))
            put("field_probes.txt",numeric_text(probe_table(field,lcfs,limiter,axis,xp,donors[:,4:6])))
        require(np.array_equal(donors[:,3],offsets) and np.array_equal(donors[:,8],watts),
                "per-geometry assigned donor powers or offsets differ")
        require(np.all(shapely.contains_xy(Polygon(wall),donors[:,4],donors[:,5])),
                "one or more source starts outside registered candidate")
        roles.append(dict(role=role,seed=seed,wall_id=wall_id,wall_file=wall_file,hit_file=hit_file,
                          candidate_id=raw["candidate"]["id"] if which=="selected" else None,
                          stage4_directory=str(directory.resolve()),stage4_result_sha256=run["result_sha256"],
                          stage4_receipt_sha256=run["receipt_sha256"],stage4_events_sha256=run["events_sha256"],
                          stage4_peak_MW_m2=measurement["peak_MW_m2"],stage4_A2_m2=measurement["A2_m2"],
                          stage4_peak_region=measurement["peak_region"],files_sha256=measurement["sha256"]))
    result = dict(case_id=case,registration_rank=reg["rank"],time_ms=1400,frame_index=reg["frame_index"],
                  raw_parquet=reg["path"],raw_sha256=reg["sha256"],frame_fingerprint=reg["fingerprint"],
                  predicted_family=common["predicted_family"],predicted_region=common["predicted_region"],
                  predicted_floor_sufficient=common["predicted_floor_sufficient"],abstained=common["abstained"],
                  baseline_peak_region=first["baseline"]["peak_region"],reference_SI=REFERENCE_SI,topology=topology,
                  source_input_provenance=provenance,donor_count=len(donors),roles=roles,files_sha256=hashes)
    emit(out/case/"case.json",json_bytes(result))
    return result


def json_bytes(value):
    return (json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n").encode()


def register(study, stage4=STAGE4, verify=False):
    out=study/"inputs"
    for name,sha in SOURCE_HASHES.items():
        require(file_hash(stage4/name)==sha,"sealed Stage 4 source differs: "+name)
    report=json.loads((stage4/"stage4_report.json").read_bytes())
    reg=json.loads((stage4/"inputs/registration.json").read_bytes())
    require(report["complete"] and report["registration_sha256"]==SOURCE_HASHES["inputs/registration.json"] and
            report["plan_sha256"]==SOURCE_HASHES["expected.json"],"Stage 4 report provenance differs")
    eligible = [r["case_id"] for r in report["runs"] if r["arm"]=="hotspot" and r["seed"]==30713 and r["status"]=="qualified"]
    require(len(eligible)==len(CASES) and set(eligible)==set(CASES),"Stage 4 eligible cohort differs")

    def emit(path, data):
        if verify:
            require(path.read_bytes()==data, "registered Stage 5 artifact differs: "+str(path))
        else:
            publish(path,data)

    cases=[]
    for row in reg["cases"]:
        if row["id"] in CASES:
            record=register_case(stage4,out,row,report,emit)
            cases.append(record)
            print(row["id"],len(record["roles"]),"roles",len({r["wall_id"] for r in record["roles"]}),"walls",flush=True)
    require(tuple(c["case_id"] for c in cases)==CASES and sum(len(c["roles"]) for c in cases)==51,
            "case or role inventory incomplete")
    result=dict(schema="nf-stage5-inputs.v1",stage4_source_sha256=SOURCE_HASHES,
                driver_sha256=file_hash(__file__),source_columns=SOURCE_HEADER.split(","),
                native_source_columns=["R_m","Z_m","CDF"],field_probe_columns=["R_m","Z_m","BR_T","BZ_T","Bphi_T","psiN"],
                hit_columns=["R_m","Z_m","first_parallel_length_m"],field_units="SI, poloidal flux Wb/rad",
                eqdsk_convention="psi[Z,R] in C order, signed RBphi constant, archived ideal 144-turn coil",
                case_record_sha256={c["case_id"]:file_hash(out/c["case_id"]/"case.json") for c in cases},cases=cases)
    emit(out/"panel.json",json_bytes(result))
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,default=DEFAULT)
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    register(args.study.resolve(),verify=args.verify)
