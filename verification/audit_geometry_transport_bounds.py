"""Independent numerical falsification checks, NOT an interval proof oracle."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline,bisplev
from scipy.optimize import brentq
from d3d_optimization.reference_transport import solve_two_point,ReferencePoint,E_CHARGE


def numeric(model,x):
    surfaces=[];areas=[];gn,gw=np.polynomial.legendre.leggauss(48)
    for s in model['surfaces']:
        a,b,n=map(np.asarray,(s['a'],s['b'],s['normal']));v=b-a;cs=[]
        for c in s['controls']:
            if not isinstance(c,dict):cs.append(float(c))
            elif 'variable' in c:cs.append(x[c['variable']])
            else:
                q=c['junction'];theta=float(q['angle_low'])+x[q['variable']]*(float(q['angle_high'])-float(q['angle_low']))
                t=np.array([np.cos(theta),np.sin(theta)])
                cs.append(q['sign']*(v@v)*(t@n)/(t@v)/15)
        sp=BSpline(s['knots'],cs,s['degree']);dsp=sp.derivative()
        def area(end,a=a,v=v,n=n,sp=sp,dsp=dsp):
            cuts=sorted(set([0.,end]+[float(t) for t in sp.t if 0<float(t)<end]));total=0.
            for lo,hi in zip(cuts,cuts[1:]):
                u=lo+(hi-lo)*(gn+1)/2
                p=a+u[:,None]*v+sp(u)[:,None]*n
                tang=v+dsp(u)[:,None]*n
                total+=(hi-lo)/2*np.dot(gw,2*np.pi*p[:,0]*np.linalg.norm(tang,axis=1))
            return total
        surfaces.append((a,v,n,sp,area));areas.append(area(1))
    prefix=np.r_[0.,np.cumsum(areas)]
    field=model['field'];tck=(np.asarray(field['z_knots']),np.asarray(field['r_knots']),np.asarray(field['coefficients']).ravel(),field['degree_z'],field['degree_r'])
    rows=[];uquad,wquad=np.polynomial.legendre.leggauss(192);uquad=(uquad+1)/2
    for donor in model['donors']:
        p=np.asarray(donor['path']);delta=np.diff(p,axis=0);h=np.linalg.norm(delta,axis=1)
        hits=[]
        for si,(a,v,n,sp,area) in enumerate(surfaces):
            inverse=np.linalg.inv(np.column_stack([v,n]));coordinates=(p-a)@inverse.T
            uu,dd=coordinates.T
            # Independent scalar root solve on every path edge that can cross.
            fs=dd-sp(np.clip(uu,0,1))
            candidates=np.flatnonzero((fs[:-1]*fs[1:]<=0)&(np.maximum(uu[:-1],uu[1:])>=0)&(np.minimum(uu[:-1],uu[1:])<=1))
            for j in candidates:
                lo,hi=0.,1.;du=uu[j+1]-uu[j]
                if du:
                    roots=sorted([-uu[j]/du,(1-uu[j])/du]);lo=max(lo,roots[0]);hi=min(hi,roots[1])
                if lo>hi:continue
                def f(t):return dd[j]+t*(dd[j+1]-dd[j])-sp(uu[j]+t*du)
                if f(lo)*f(hi)>0:continue
                t=brentq(f,lo,hi,xtol=1e-14);u=uu[j]+t*du
                hits.append((j,t,si,u))
        if not hits:raise AssertionError('independent numerical no-hit')
        j,t,si,u=min(hits);point=p[j]+t*delta[j]
        metric=np.asarray(donor['metric'][:j+2],float);f=np.asarray(donor['integrand'][:j+2],float)
        mapping=np.asarray(donor['mapping'][:j+2],float)
        def cumulative(values):return np.r_[0,np.cumsum(h[:j+1]*(values[:-1]+values[1:])/2)]
        cm=cumulative(metric);cf=cumulative(f)
        L=cm[j]+h[j]*(t*(1-t/2)*metric[j]+t*t/2*metric[j+1])
        J=cf[j]+h[j]*(t*(1-t/2)*f[j]+t*t/2*f[j+1])
        g=np.array([bisplev(point[1],point[0],tck,dy=1),bisplev(point[1],point[0],tck,dx=1)],float)
        gnrm=np.linalg.norm(g);B=np.hypot(gnrm,field['toroidal_K'])/point[0]
        ref=ReferencePoint(**model['reference']);state=solve_two_point(donor['qu'],B/donor['b0'],J,ref)
        y=state['Tt_eV']**1.75;yu=state['Tu_eV']**1.75;S=y+yu
        partial=J*(1-uquad)*(1+uquad*(yu-y)/S)
        idx=np.clip(np.searchsorted(cf,partial,side='right')-1,0,j);d=partial-cf[idx]
        slope=(f[idx+1]-f[idx])/h[idx]
        dx=2*d/(f[idx]+np.sqrt(f[idx]**2+2*slope*d))
        values=mapping[idx]+(dx/h[idx])*(mapping[idx+1]-mapping[idx])
        sigma=np.sqrt(3*E_CHARGE*ref.density_m3*state['Tu_eV']*ref.chi_perp_m2_s/ref.kappa0)*J/S*np.dot(wquad/2,values)
        center=(prefix[si]+surfaces[si][4](u))/prefix[-1]
        physical=sigma*model['grad_sep']/gnrm
        rows.append({'length_m':L,'integral_Bratio_m':J,'sigma_flux_m':sigma,
            'sigma_physical_m':physical,'sigma_area_m2':2*np.pi*point[0]*physical,
            'center_fraction':center,'arrival_RZ_m':point,'edge':j,'fraction':t})
    return prefix[-1],rows


def audit(certificate,samples=3):
    if not certificate['transport_parameter_enclosure_proven']:raise ValueError('resolved certificate required')
    model=certificate['model'];bounds=np.asarray(certificate['box'],float);rng=np.random.default_rng(4157)
    points=[bounds.mean(axis=1)]+[rng.uniform(bounds[:,0],bounds[:,1]) for _ in range(samples-1)]
    checks=0
    def contained(value,interval):
        nonlocal checks
        lo,hi=map(float,interval);pad=2e-8*max(1.,abs(value))
        if not lo-pad<=value<=hi+pad:raise AssertionError((value,lo,hi))
        checks+=1
    for x in points:
        area,rows=numeric(model,x);contained(area,certificate['total_area_m2'])
        for numeric_row,bound in zip(rows,certificate['donors']):
            for key in ('length_m','integral_Bratio_m','sigma_flux_m','sigma_physical_m','sigma_area_m2','center_fraction'):
                contained(numeric_row[key],bound[key])
            for k in range(2):contained(numeric_row['arrival_RZ_m'][k],bound['arrival_RZ_m'][k])
            branches=bound['first_hit']['possible_hits']
            if not any(h['edge']==numeric_row['edge'] and float(h['fraction'][0])-1e-8<=numeric_row['fraction']<=float(h['fraction'][1])+1e-8 for h in branches):
                raise AssertionError('independent first hit outside retained union')
            checks+=1
    return {'samples':len(points),'checks':checks,'passed':True,
            'scope':'independent floating-point falsification only; not the inclusion proof',
            'comparison_tolerance':'2e-8 * max(1,abs(value)); hit fraction 1e-8'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('certificate');parser.add_argument('--samples',type=int,default=3)
    args=parser.parse_args();r=json.loads(Path(args.certificate).read_text());print(json.dumps(audit(r,args.samples),indent=2))
