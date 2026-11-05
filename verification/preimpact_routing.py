"""Axisymmetric first-hit transfer on identified monotone wall arcs."""
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr


class Wall:
    def __init__(self,points,flux):
        self.p=np.asarray(points,float)
        self.flux=flux
        dz=np.diff(self.p,axis=0)
        ds=np.linalg.norm(dz,axis=1)
        if not np.all(ds>0) or not np.allclose(self.p[0],self.p[-1],atol=1e-12,rtol=0):
            raise ValueError('a closed, ordered wall with nonzero edges is required')
        self.t=dz/ds[:,None]
        self.s=np.r_[0.,np.cumsum(ds)]
        self.area=np.r_[0.,np.cumsum(np.pi*(self.p[:-1,0]+self.p[1:,0])*ds)]
        # Include all geometric corners and isolate stationary values along long edges.
        ns=np.maximum(1,np.ceil(ds/.0005).astype(int))
        seg=np.repeat(np.arange(len(ds)),ns)
        at=np.arange(len(seg))-np.repeat(np.r_[0,np.cumsum(ns)[:-1]],ns)
        lo=self.s[seg]+at*ds[seg]/ns[seg]
        hi=self.s[seg]+(at+1)*ds[seg]/ns[seg]
        g0=np.sum(flux.gradient(self.points_at(lo))*self.t[seg],axis=1)
        g1=np.sum(flux.gradient(self.points_at(hi))*self.t[seg],axis=1)
        critical=[]
        for i in np.flatnonzero(g0*g1<0):
            tangent=self.t[seg[i]]
            derivative=lambda x:float(flux.gradient(self.points_at([x]))[0]@tangent)
            critical.append(brentq(derivative,lo[i],hi[i],xtol=2e-14))
        probe=np.unique(np.r_[self.s,lo[g0==0],critical])
        u=flux.values(self.points_at(probe))
        # Explicitly split flat arcs as well. They carry zero continuous label measure.
        directions=np.sign(np.diff(u))
        split=np.flatnonzero(directions[1:]!=directions[:-1])+1
        self.turn=probe[np.r_[0,split,len(probe)-1].astype(int)]
        self.level=flux.values(self.points_at(self.turn))

    def points_at(self,s):
        return np.column_stack([np.interp(s,self.s,col) for col in self.p.T])

    def area_at(self,s):
        j=np.clip(np.searchsorted(self.s,s,side='right')-1,0,len(self.s)-2)
        ds=s-self.s[j]
        return self.area[j]+np.pi*ds*(2*self.p[j,0]+ds*self.t[j,0])

    def arc_id(self,s):
        return np.clip(np.searchsorted(self.turn,s,side='right')-1,0,len(self.turn)-2)

    def hit_s(self,points,seg):
        return self.s[seg]+np.sum((points-self.p[seg])*self.t[seg],axis=1)


def ring_grid(wall,spacing,phase):
    length=wall.s[-1];n=round(length/spacing)
    if n<2 or not 0<=phase<1:raise ValueError('invalid wall binning')
    edges=(np.arange(n)+phase)*length/n
    cuts=np.unique(np.r_[0.,edges,wall.turn,length])
    mid=(cuts[:-1]+cuts[1:])/2
    ring=np.floor(mid*n/length-phase).astype(int)%n
    u=wall.flux.values(wall.points_at(cuts))
    lo=np.minimum(u[:-1],u[1:]);hi=np.maximum(u[:-1],u[1:])
    lmin=np.full(n,np.inf);lmax=np.full(n,-np.inf)
    np.minimum.at(lmin,ring,lo);np.maximum.at(lmax,ring,hi)
    area=np.bincount(ring,weights=np.diff(wall.area_at(cuts)),minlength=n)
    if not np.all(area>0):raise ValueError('zero ring area')
    return dict(area=area,s_mid=((np.arange(n)+phase+.5)*length/n)%length,
                ring=ring,lo=lo,hi=hi,arc=wall.arc_id(mid),u_min=lmin,u_max=lmax)


class GaussianMass:
    def __init__(self,centers,widths,weights):
        self.c,self.s,self.w=map(np.asarray,(centers,widths,weights))
        if (not np.all(np.isfinite([self.c,self.s,self.w])) or np.any(self.s<=0)
                or np.any(self.w<0) or abs(self.w.sum()-1)>1e-12):
            raise ValueError('finite unit-power mixture with positive widths required')

    def __call__(self,lo,hi):
        lo,hi=np.broadcast_arrays(np.asarray(lo,float),np.asarray(hi,float))
        if np.any(hi<lo):raise ValueError('unordered interval')
        flat=np.empty(lo.size)
        for i in range(0,lo.size,256):
            a=(lo.ravel()[i:i+256,None]-self.c)/self.s
            b=(hi.ravel()[i:i+256,None]-self.c)/self.s
            flat[i:i+256]=np.where(a>=0,ndtr(-a)-ndtr(-b),ndtr(b)-ndtr(a))@self.w
        return flat.reshape(lo.shape)


def split_routes(cuts,trace,depth=6,fractions=(.002,.5,.998)):
    lo=np.asarray(cuts[:-1]);hi=np.asarray(cuts[1:])
    good=[];bad=[]
    for level in range(depth+1):
        if not len(lo):break
        x=lo[:,None]+(hi-lo)[:,None]*fractions
        index=np.asarray(trace(x.ravel())).reshape(x.shape)
        same=np.all(index==index[:,:1],axis=1)
        for l,h,k in zip(lo[same],hi[same],index[same,0]):
            if k<0:bad.append((l,h))
            else:good.append((l,h,int(k)))
        a=lo[~same];b=hi[~same];m=(a+b)/2
        if level==depth:
            bad.extend(zip(a,b));break
        lo=np.r_[a,m];hi=np.r_[m,b]
    return sorted(good),sorted(bad)


def merged(good):
    result=[]
    for lo,hi,j in sorted(good):
        if result and result[-1][1]==lo and result[-1][2]==j:
            result[-1]=(result[-1][0],hi,j)
        else:result.append((lo,hi,j))
    return result


def route_bins(wall,grid,good,unknown,mass):
    del wall  # grid includes exact label extrema and area on this material wall
    p=np.zeros_like(grid['area']);cap=np.zeros_like(p)
    for l,h,arc in merged(good):
        a=np.maximum(grid['lo'],l);b=np.minimum(grid['hi'],h)
        ok=(grid['arc']==arc)&(b>a)
        np.add.at(p,grid['ring'][ok],mass(a[ok],b[ok]))
    for l,h in unknown:
        a=np.maximum(grid['u_min'],l);b=np.minimum(grid['u_max'],h)
        ok=b>a
        cap[ok]+=mass(a[ok],b[ok])
    return p,cap


def scores(p,a):
    return dict(peak=float(np.max(p/a)),A2=float(p.sum()**2/np.sum(p*p/a)))


def complete_bounds(p,a,unknown,cap=None):
    """Conservative bounds allowing any fraction of the unassigned power to reach the wall."""
    if unknown<0 or not (np.all(a>0) and np.all(p>=0) and p.sum()>0):raise ValueError('invalid ledger')
    c=np.full_like(p,unknown) if cap is None else np.minimum(unknown,cap)
    if np.any(c<0):raise ValueError('negative bin cap')
    q=p/a
    d=float(np.sum(p*q))
    extra=min(float(np.sum((2*p+c)*c/a)),unknown*float(np.max(2*q+c/a)))
    return dict(peak=[float(q.max()),float(np.max((p+c)/a))],
                A2=[float(p.sum()**2/(d+extra)),float((p.sum()+unknown)**2/d)])

