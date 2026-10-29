"""Optional FP64 CUDA sum; only density arithmetic differs from the CPU model."""
import os,sys
from pathlib import Path
import numpy as np
from .boundary_transport import wrapped_density

_cp=None
_kernel=None
_dll_handles=[]


def cuda():
    global _cp,_kernel
    if _cp is not None:return _cp
    runtime=Path(os.environ.get('DIVERTOR_GPU_RUNTIME','D:/DivertorOptimization/runtime/cupy'))
    if runtime.is_dir():
        sys.path.append(str(runtime))
        for folder in runtime.glob('nvidia/*/bin'):
            if hasattr(os,'add_dll_directory'):_dll_handles.append(os.add_dll_directory(str(folder)))
    cache=Path(os.environ.get('CUPY_CACHE_DIR','D:/DivertorOptimization/runtime/cupy-cache'))
    cache.mkdir(parents=True,exist_ok=True);os.environ['CUPY_CACHE_DIR']=str(cache)
    import cupy as cp
    _kernel=cp.RawKernel(r'''
    extern "C" __global__ void footprint(const double* a,const double* mu,const double* sigma,
        const double* weight,double period,int na,int nd,double* q) {
        int j=blockDim.x*blockIdx.x+threadIdx.x;
        if(j>=na)return;
        double value=0.;
        for(int i=0;i<nd;i++){
            double d=a[j]-mu[i];
            if(d>period/2.)d-=period;else if(d< -period/2.)d+=period;
            double v=0.;
            for(int k=-1;k<=1;k++){
                double z=(d+k*period)/sigma[i];v+=exp(-0.5*z*z);
            }
            value+=weight[i]*v/(2.5066282746310005024*sigma[i]);
        }
        q[j]=value;
    }''','footprint',options=('--fmad=false',))
    _cp=cp
    return cp


def density_sum(area,centers,widths,weights,period,backend='gpu'):
    area=np.asarray(area,float)%period;centers=np.asarray(centers,float)%period
    widths=np.asarray(widths,float);weights=np.asarray(weights,float)
    if backend=='cpu':return sum((w*wrapped_density(area,c,s,period) for c,s,w in zip(centers,widths,weights)),np.zeros_like(area))
    if backend!='gpu':raise ValueError('backend must be cpu or gpu')
    narrow=widths/period<.12;q=np.zeros_like(area)
    if np.any(narrow):
        cp=cuda();a=cp.asarray(area);mu=cp.asarray(centers[narrow]);s=cp.asarray(widths[narrow]);w=cp.asarray(weights[narrow]);result=cp.empty_like(a)
        _kernel(((len(area)+127)//128,),(128,),(a,mu,s,w,np.float64(period),np.int32(len(area)),np.int32(narrow.sum()),result))
        q=cp.asnumpy(result)
    for c,s,w in zip(centers[~narrow],widths[~narrow],weights[~narrow]):q+=w*wrapped_density(area,c,s,period)
    if not np.all(np.isfinite(q)):raise FloatingPointError('nonfinite GPU density')
    return q
