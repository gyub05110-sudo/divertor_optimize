"""Independent quadrature control: trace displaced rays, histogram first contacts directly.
Uses no wall-flux extrema, branch partition, or inversion of a wall flux coordinate.
"""
import argparse
from pathlib import Path
from time import perf_counter
import numpy as np
from d3d_optimization.case import load_case,separatrix
from d3d_optimization.frozen_mapping import trace_paths,first_hits
from d3d_optimization.broadening_study import load_json,save,digest
from verification.nf_preimpact_comparison import Flux,Section,private_section,ROOT
from verification.preimpact_routing import GaussianMass,Wall,scores


class Weights:
    def __init__(self,mixtures,lower,upper):
        self.mixes=mixtures;self.cp=None
        self.timing={'using_gpu':False}
        try:
            from d3d_optimization.gpu_footprint import cuda
            cp=cuda()
        except (ImportError,RuntimeError,OSError) as e:
            self.timing['fallback']=str(e)
            return
        self.cp=cp
        self.kernel=cp.RawKernel(r'''
        extern "C" __global__ void flux_mass(const double* lo,const double* hi,const double* c,
            const double* s,const double* w,double* out,const int n) {
            __shared__ double block[128];
            int b=blockIdx.x, t=threadIdx.x;
            double sum=0;
            for(int j=t;j<n;j+=blockDim.x){
                double a=(lo[b]-c[j])/(s[j]*1.4142135623730951), h=(hi[b]-c[j])/(s[j]*1.4142135623730951);
                double p=(a>=0.0)?erfc(a)-erfc(h):erfc(-h)-erfc(-a);
                sum+=w[j]*0.5*p;
            }
            block[t]=sum; __syncthreads();
            for(int k=64;k>0;k/=2){if(t<k) block[t]+=block[t+k]; __syncthreads();}
            if(t==0) out[b]=block[0];
        }''','flux_mass')
        self.params={k:tuple(cp.asarray(a) for a in (m.c,m.s,m.w)) for k,m in mixtures.items()}
        x=np.linspace(lower,upper,2049);k=next(iter(mixtures))
        self.gpu(k,x[:-1],x[1:])
        t=perf_counter();ref=mixtures[k](x[:-1],x[1:]);cpu_t=perf_counter()-t
        t=perf_counter();y=self.gpu(k,x[:-1],x[1:]);gpu_t=perf_counter()-t
        np.testing.assert_allclose(ref,y,rtol=1e-10,atol=2e-15)
        self.timing.update(cpu_s=cpu_t,gpu_s=gpu_t,max_mass_difference=float(np.max(abs(ref-y))),using_gpu=gpu_t<cpu_t)

    def gpu(self,k,a,b):
        cp=self.cp
        aa=cp.asarray(np.ascontiguousarray(a));bb=cp.asarray(np.ascontiguousarray(b));out=cp.empty_like(aa)
        self.kernel((len(a),),(128,),(aa,bb,*self.params[k],out,np.int32(len(self.mixes[k].c))))
        return cp.asnumpy(out)

    def mass(self,k,a,b):
        return self.gpu(k,a,b) if self.timing['using_gpu'] else self.mixes[k](a,b)


def ray_quadrature(section,n,weights,walls,arrays,step):
    # A composite mixture chooses only the integration mesh. Each geometry keeps its original weights.
    t=np.linspace(0,1,2049)
    fine=section.umin+(section.umax-section.umin)*(.5-.5*np.cos(np.pi*t))
    blend=np.mean([weights.mass(k,fine[:-1],fine[1:]) for k in walls],axis=0)
    cum=np.r_[0.,np.cumsum(blend)]
    edges=np.interp(np.linspace(0,cum[-1],n+1),cum,fine);edges[[0,-1]]=[section.umin,section.umax]
    middle=.5*(edges[1:]+edges[:-1]);starts=section.starts(middle)
    masses={k:weights.mass(k,edges[:-1],edges[1:]) for k in walls}
    for i in range(0,n,96):
        cut=slice(i,i+96);paths=trace_paths(section.flux.field,starts[cut],step=step,max_length=8.)
        for k,wall in walls.items():
            h=first_hits(paths,wall.p)
            if not np.all(h.segment>=0):raise ValueError('missed wall in direct ray check')
            s=wall.hit_s(h.points,h.segment)
            bins=(np.floor(s*len(arrays[k])/wall.s[-1]).astype(int))%len(arrays[k])
            arrays[k]+=np.bincount(bins,weights=masses[k][cut],minlength=len(arrays[k]))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True);p.add_argument('--evidence',type=Path,required=True)
    p.add_argument('--source-csv',type=Path,required=True);a=p.parse_args()
    out=a.study/'raycheck';out.mkdir(exist_ok=False)
    source=np.genfromtxt(a.source_csv,names=True,delimiter=',')
    field,lcfs,ax=load_case(ROOT/'input/D3D/d3d/d3d_shot_1a3155370c.parquet',1400)
    rsep,xp,_=separatrix(field,lcfs,ax);z=ax[1];flux=Flux(field,(rsep,z))
    centers=flux.values(np.c_[source['R_m'],source['Z_m']]);mixes={};walls={};expected={};areas={}
    for k,v in load_json(a.evidence)['inputs'].items():
        assert digest(v['path'])==v['sha256']
        with np.load(v['path']) as f:
            mixes[k]=GaussianMass(centers,f['thermal_sigma_m'],source['power_fraction'])
            walls[k]=Wall(f['wall'],flux)
        with np.load(a.study/(k+'_fine')/'leg-1.305_2mm_0.npz') as f:
            expected[k]=f['resolved_MW'];areas[k]=f['area']
    sections=[private_section(-1.305,w,xp) for w in walls.values()]
    assert np.ptp([s.lo for s in sections])<1e-11 and np.ptp([s.hi for s in sections])<1e-11
    omp=Section(z,rsep,rsep+source['x_high_m'][-1],flux)
    weights=Weights(mixes,sections[0].umin,omp.umax)
    report=dict(weight_integration=weights.timing,levels={})
    for n in (2048,8192):
        t=perf_counter();arrays={k:np.zeros_like(p) for k,p in expected.items()}
        ray_quadrature(omp,n,weights,walls,arrays,.00075)
        ray_quadrature(sections[0],n//8,weights,walls,arrays,.00075)
        rows={}
        for k,x in arrays.items():
            a0=areas[k];control=scores(x,a0);ref=scores(expected[k],a0)
            rows[k]=dict(direct=control,peak_relative_delta=control['peak']/ref['peak']-1,
                         A2_relative_delta=control['A2']/ref['A2']-1,mass_delta=float(x.sum()-expected[k].sum()),
                         TV_MW=float(np.abs(expected[k]-x).sum()/2))
            np.savez_compressed(out/f'{k}_{n}.npz',first_hit_MW=x,area=a0)
        report['levels'][n]=dict(seconds=perf_counter()-t,rows=rows)
        print(n,rows,flush=True)
    save(out/'result.json',report,exclusive=True)

if __name__=='__main__':main()

