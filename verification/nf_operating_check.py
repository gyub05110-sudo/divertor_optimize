"""Compare the reused and independently rebuilt magnetic workspace at a new operating point."""
import argparse
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import numpy as np
from d3d_optimization.broadening_study import save,load_json
from d3d_optimization.boundary_evaluator import export_boundary
from d3d_optimization.campaign_geometry import Candidate,evaluate_candidate,evaluator
from d3d_optimization.operating_sweep import reuse_boundary
from d3d_optimization.physical_profile import read_physical_profile
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_search_benchmark import ANCHORS,ROOT,prep


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    reference=ReferencePoint(power_W=2e6,density_m3=1.3e19,chi_perp_m2_s=.6)
    evidence={}
    for source in ('OMP','IMP'):
        nominal,initial=prep(ROOT,source,ReferencePoint(),0)
        before=nominal.baseline['peak_q']
        t=perf_counter();shared=reuse_boundary(nominal,reference);shared_t=perf_counter()-t
        fresh,fresh_t=prep(ROOT,source,reference,0)
        c=Candidate.from_record(load_json(ANCHORS)[source][0]['candidate'])
        results={}
        candidate_times={}
        for label in ('baseline','candidate'):
            raw=[]
            for setup,obj in [('shared',shared),('fresh',fresh)]:
                t=perf_counter()
                result=obj.baseline if label=='baseline' else evaluate_candidate(evaluator(obj,c.radius_mm,'gpu'),c)
                if label=='candidate':candidate_times[setup+'_gpu_s']=perf_counter()-t
                if not result['feasible']:raise ValueError(f'{source} {label}: {result["errors"]}')
                raw.append(result)
            left,right=raw
            if label=='candidate':
                t=perf_counter()
                cpu=evaluate_candidate(evaluator(fresh,c.radius_mm,'cpu'),c)
                candidate_times['fresh_cpu_s']=perf_counter()-t
                if not cpu['feasible']:raise ValueError('CPU parity candidate infeasible')
                qgpu=right['deposition']['q'];qcpu=cpu['deposition']['q']
                np.testing.assert_allclose(qgpu,qcpu,rtol=2e-11,atol=1e-12)
                candidate_times['max_q_relative_error']=float(np.max(abs(qgpu-qcpu))/np.max(qcpu))
            diff={}
            for name,arr1,arr2 in [('hit_RZ',left['hits'].points,right['hits'].points),
                                  ('q',left['deposition']['q'],right['deposition']['q']),
                                  ('cell_power',left['deposition']['cell_power'],right['deposition']['cell_power']),
                                  ('Tu',left['thermal']['Tu_eV'],right['thermal']['Tu_eV']),
                                  ('sigma',left['thermal']['sigma_m'],right['thermal']['sigma_m'])]:
                np.testing.assert_allclose(arr1,arr2,rtol=1e-12,atol=1e-14)
                diff[name]=float(np.nanmax(abs(arr1-arr2)))
            output=a.output/f'{source}_{label}.npz'
            export_boundary(output,shared,left)
            profile=read_physical_profile(output,reference.power_W)
            np.testing.assert_allclose(profile['continuous_q_MW_m2'],right['deposition']['q']*reference.power_W/1e6)
            np.testing.assert_allclose(profile['cell_power_MW'].sum(),2*right['deposition']['ledger']['deposited'])
            results[label]=dict(peak_MW_m2=float(profile['continuous_q_MW_m2'].max()),power_MW=float(profile['cell_power_MW'].sum()),
                                max_absolute_differences=diff)
        assert before==nominal.baseline['peak_q'] and fresh.thermal_paths is not shared.thermal_paths
        evidence[source]=dict(nominal_setup_s=initial,reused_setup_s=shared_t,independent_fresh_setup_s=fresh_t,
                              candidate_parity=candidate_times,results=results)
        save(a.output/f'{source}.json',evidence[source],exclusive=True)
    save(a.output/'summary.json',dict(reference=asdict(reference),cases=evidence),exclusive=True)

if __name__=='__main__':main()

