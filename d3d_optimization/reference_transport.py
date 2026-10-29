"""Fixed-input, one-way effective two-point conduction closure.

Temperatures are pre-redistribution tube states. No particle/momentum equation
or feedback from transverse deposition is solved. Electron chi uses C=1.5 e n.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.optimize import brentq
from .broadening import PathLengths

E_CHARGE = 1.602176634e-19
M_DEUTERIUM = 3.3435837724e-27
M_ELECTRON = 9.1093837139e-31


@dataclass(frozen=True)
class ReferencePoint:
    power_W: float = 1e6
    density_m3: float = 1e19
    chi_perp_m2_s: float = .3
    lambda_q_m: float = .003
    gamma: float = 7.
    kappa0: float = 2000.

    def __post_init__(self):
        if not all(np.isfinite(x) and x > 0 for x in asdict(self).values()):
            raise ValueError('reference inputs must be positive and finite')


def solve_two_point(qu, field_ratio, integral_Bratio_m, reference):
    """Unique positive log-temperature root of energy and sheath conditions."""
    if not all(np.isfinite(x) and x > 0 for x in (qu, field_ratio, integral_Bratio_m)):
        raise ValueError('positive finite tube power, field ratio and integral required')
    r=reference
    qt=qu*field_ratio
    energy=3.5*qu*integral_Bratio_m/r.kappa0
    coefficient=r.gamma*E_CHARGE*r.density_m3/2*np.sqrt(2*E_CHARGE/M_DEUTERIUM)
    log_rhs=np.log(qt/coefficient)
    def residual(log_t):
        return (2/7)*np.logaddexp(3.5*log_t,np.log(energy))+.5*log_t-log_rhs
    upper=(2/3)*log_rhs
    lower=min(upper-10,2*(log_rhs-(2/7)*np.log(energy))-10)
    log_t=brentq(residual,lower,upper,xtol=1e-13)
    tt=float(np.exp(log_t))
    tu=float(np.exp((2/7)*np.logaddexp(3.5*log_t,np.log(energy))))
    nt=r.density_m3*tu/(2*tt)
    sheath=r.gamma*E_CHARGE*nt*tt*np.sqrt(2*E_CHARGE*tt/M_DEUTERIUM)
    # Compare complete positive terms, avoiding cancellation at small L.
    energy_error=abs(tu**3.5-(tt**3.5+energy))/(tt**3.5+energy)
    return dict(Tu_eV=tu,Tt_eV=tt,nt_m3=nt,qu_W_m2=float(qu),qt_W_m2=float(qt),
                conduction_relative_residual=float(energy_error),
                sheath_relative_residual=float(abs(sheath-qt)/qt))


def cumulative_trapezoid(values, arc):
    return np.r_[0.,np.cumsum(np.diff(arc)*(values[:-1]+values[1:])/2)]


class ThermalPaths:
    """Cached magnetic integrals; thermal states recomputed at every first hit.

    Field interpolation on the traced polyline is a numerical approximation.
    Quadrature order and path refinement are independently exposed for audits.
    The static bulk pressure applies through the last quadrature node; n_t is
    a separate sonic sheath quantity, never inserted into the bulk integral.
    """
    def __init__(self, field, paths, source_rz, source, reference, lengths=None,
                 quadrature_order=128):
        if not isinstance(quadrature_order,int) or quadrature_order<8:
            raise ValueError('quadrature order must be an integer >=8')
        if not np.isclose(source.width,reference.lambda_q_m,rtol=0,atol=1e-14):
            raise ValueError('reference and imposed source widths differ')
        self.field,self.paths,self.source,self.reference=field,paths,source,reference
        self.lengths=lengths if lengths is not None else PathLengths(field,paths)
        self.order=quadrature_order
        grad_sep=abs(float(field.gradient([source_rz])[0,0]))
        if not np.isfinite(grad_sep) or grad_sep<=0:
            raise ValueError('nonzero finite source flux gradient required')
        self.grad_sep=grad_sep
        self.data=[]
        for p,arc,metric in zip(self.lengths.paths,self.lengths.s,self.lengths.metric):
            b=np.linalg.norm(field.magnetic(p),axis=1)
            mapping=np.linalg.norm(field.gradient(p),axis=1)/grad_sep
            # g ds_parallel is regular at an X-point in a consistent field.
            spatial=mapping*metric
            integral=cumulative_trapezoid(b/b[0]*metric,arc)
            self.data.append((b,spatial,integral))
        if len(source.x)!=len(paths):
            raise ValueError('source/path count mismatch')
        self.qu=reference.power_W*np.exp(-np.asarray(source.x)/source.width)/source.normalization

    def evaluate(self,hits,quadrature_order=None):
        n=len(self.paths)
        if np.shape(hits.distance)!=(n,) or np.shape(hits.segment)!=(n,):
            raise ValueError('one hit per donor required')
        keys=['Tu_eV','Tt_eV','nt_m3','qu_W_m2','qt_W_m2',
              'conduction_relative_residual','sheath_relative_residual',
              'sigma_m','length_m','integral_Bratio_m','max_free_streaming_ratio',
              'conductive_temperature_drop']
        out={key:np.full(n,np.nan) for key in keys}
        good=np.flatnonzero(np.asarray(hits.segment)>=0)
        distances=np.where(np.asarray(hits.segment)>=0,hits.distance,np.inf)
        out['length_m']=self.lengths.at(distances)
        if not len(good):return out
        bhit=np.linalg.norm(self.field.magnetic(np.asarray(hits.points)[good]),axis=1)
        nodes,weights=np.polynomial.legendre.leggauss(quadrature_order or self.order)
        unit=(nodes+1)/2
        r=self.reference
        for k,i in enumerate(good):
            distance=float(hits.distance[i]);arc=self.lengths.s[i]
            b,spatial,cumulative=self.data[i]
            if distance<=0 or distance>arc[-1]+1e-10:
                raise ValueError('positive hit inside traced path required')
            # Integral within the last segment uses a linearly varying metric.
            j=max(0,min(int(np.searchsorted(arc,distance,side='left')-1),len(arc)-2))
            h=distance-arc[j];fraction=h/(arc[j+1]-arc[j])
            integrand=b/b[0]*self.lengths.metric[i]
            total=cumulative[j]+h*(integrand[j]+.5*fraction*(integrand[j+1]-integrand[j]))
            state=solve_two_point(float(self.qu[i]),bhit[k]/b[0],total,r)
            # y=T^(7/4) removes the integrable cold-target singularity.
            # Compute yu-yt through the energy identity to avoid cancellation.
            yt=state['Tt_eV']**1.75;yu=state['Tu_eV']**1.75
            a=3.5*self.qu[i]/r.kappa0
            y=yt+unit*(a*total/(yu+yt))
            partial=total*(1-unit*(y+yt)/(yu+yt))
            index=np.clip(np.searchsorted(cumulative,partial,side='right')-1,0,len(arc)-2)
            delta=partial-cumulative[index]
            slope=(integrand[index+1]-integrand[index])/(arc[index+1]-arc[index])
            dx=2*delta/(integrand[index]+np.sqrt(np.maximum(integrand[index]**2+2*slope*delta,0)))
            s=arc[index]+dx
            temp=y**(4/7)
            density=r.density_m3*state['Tu_eV']/temp
            mapping_over_b=spatial/self.lengths.metric[i]*b[0]/b
            prefactor=np.sqrt(1.5*E_CHARGE*r.density_m3*state['Tu_eV']*r.chi_perp_m2_s/r.kappa0)
            spreading=2*prefactor*total/(yu+yt)*np.dot(weights/2,np.interp(s,arc,mapping_over_b))
            local_q=self.qu[i]*np.interp(s,arc,b)/b[0]
            q_free=density*E_CHARGE*temp*np.sqrt(2*E_CHARGE*temp/M_ELECTRON)
            upstream_ratio=self.qu[i]/(r.density_m3*E_CHARGE*state['Tu_eV']*np.sqrt(2*E_CHARGE*state['Tu_eV']/M_ELECTRON))
            # Exact bulk target limit uses 2*n_t, not sheath n_t.
            target_ratio=r.gamma/2*np.sqrt(M_ELECTRON/M_DEUTERIUM)
            state.update(sigma_m=float(spreading/np.sqrt(2)),
                         integral_Bratio_m=float(total),
                         max_free_streaming_ratio=float(max(np.max(local_q/q_free),upstream_ratio,target_ratio)),
                         conductive_temperature_drop=float(1-(state['Tt_eV']/state['Tu_eV'])**3.5))
            if not all(np.isfinite(v) and v>=0 for v in state.values()) or state['sigma_m']<=0:
                raise ValueError('nonfinite/nonpositive thermal closure; no clipping allowed')
            for key,value in state.items():out[key][i]=value
        return out


def thermal_summary(state,weights,reference):
    """Power-weighted diagnostic fractions; warnings are not regime proofs."""
    weights=np.asarray(weights,float)
    good=np.isfinite(state['Tt_eV'])
    total=float(weights.sum())
    hit=float(weights[good].sum())
    if total<=0 or hit<=0:return {'resolved_power_fraction':0.}
    result=dict(resolved_power_fraction=hit/total,
        cold_target_power_fraction=float(weights[good&(state['Tt_eV']<5)].sum()/total),
        free_streaming_ratio_above_0p1_power_fraction=float(weights[good&(state['max_free_streaming_ratio']>.1)].sum()/total),
        weak_conduction_power_fraction=float(weights[good&(state['conductive_temperature_drop']<.1)].sum()/total),
        sigma_above_0p3_lambda_power_fraction=float(weights[good&(state['sigma_m']>.3*reference.lambda_q_m)].sum()/total))
    for key in ['Tu_eV','Tt_eV','nt_m3','sigma_m','length_m']:
        result[key]={'min':float(np.min(state[key][good])), 'max':float(np.max(state[key][good])),
                     'power_weighted_mean':float(np.dot(weights[good],state[key][good])/hit)}
    for key in ['conduction_relative_residual','sheath_relative_residual']:
        result[key]=float(np.max(state[key][good]))
    return result
