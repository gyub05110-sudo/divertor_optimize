"""Budgeted spatial peak bounds; optional finalist layer, not a search gate.

Old certificate algorithms are unchanged. Every operation is replayable; wall
time stops between probes/splits (one split evaluates two cells and one point).
"""
import copy
import math
from time import perf_counter
from .transport_peak_bounds import I,D,DOWN,UP,KernelBounds,rational_interval
from .optimality_bounds import pi_bounds
from .geometry_transport_certificate import universal_peak_bounds,digest


class SpatialBounds:
    """Split narrow/broad width families, retaining both at their boundary."""
    def __init__(self,model):
        self.model=model;self.parts=[]
        area=I(*model['period_m2']);pi=rational_interval(*pi_bounds())
        for donor in model['donors']:
            weight=I(*donor['weight_MW']);sigma=I(*donor['sigma_m2'])
            # Validate all centers, including the Fourier-only branch.
            center=donor.get('center_fraction');absolute=donor.get('center_m2')
            if (center is None)==(absolute is None):raise ValueError('one center representation required')
            c=I(*(center if center is not None else absolute))
            if c.lo<0 or c.hi>(D(1) if center is not None else area.lo):raise ValueError('center outside period')
            if weight.hi==0:continue
            threshold=UP.multiply(D('.2'),area.hi);branches=[]
            if sigma.lo<threshold:
                narrow=copy.deepcopy(donor);narrow['sigma_m2']=[str(sigma.lo),str(min(sigma.hi,threshold))]
                try:branches.append(KernelBounds({'period_m2':model['period_m2'],'donors':[narrow]}))
                except ValueError as exc:
                    if 'image-series budget' not in str(exc):raise
                    # Extremely uncertain period: safe fallback, never drop donor.
                    branches.append(I(0,I(*universal_peak_bounds({'period_m2':model['period_m2'],'donors':[narrow]})).hi))
            if sigma.hi>=threshold:
                broad=I(max(sigma.lo,threshold),sigma.hi)
                exponent=-2*pi.square()*(broad/area).square()
                mass=I(exponent.exp().hi)/(1-I((3*exponent).exp().hi))
                bound=weight/area*(1+I(-1,1)*2*mass)
                branches.append(I(max(D(0),bound.lo),bound.hi))
            self.parts.append(branches)

    def density(self,t):
        total=I(0)
        for branches in self.parts:
            values=[b if isinstance(b,I) else b.density(t) for b in branches]
            total+=I(min(v.lo for v in values),max(v.hi for v in values))
        return total


def _initial(model):
    lo,hi=map(D,universal_peak_bounds(model))
    return {'lower':lo,'leaves':[{'id':1,'a':D(0),'b':D(1),'upper':hi}]}


def _operation(state,evaluator,op):
    if op['kind']=='probe':
        state['lower']=max(state['lower'],evaluator.density(I(op['point'])).lo)
        return
    active=next(x for x in state['leaves'] if x['id']==op['id'])
    a,b=active['a'],active['b'];mid=DOWN.divide(DOWN.add(a,b),D(2))
    if not a<mid<b:raise ValueError('coordinate precision exhausted')
    children=[]
    for identity,x,y in [(2*active['id'],a,mid),(2*active['id']+1,mid,b)]:
        children.append({'id':identity,'a':x,'b':y,'upper':min(active['upper'],evaluator.density(I(x,y)).hi)})
    state['lower']=max(state['lower'],evaluator.density(I(mid)).lo)
    state['leaves'].remove(active);state['leaves'].extend(children)


def _evidence(state):
    upper=max(x['upper'] for x in state['leaves']);lower=state['lower']
    if lower>upper:raise ArithmeticError('inconsistent peak enclosure')
    return {'peak_bracket_MW_m2':[str(lower),str(upper)],
            'gap_upper_MW_m2':str(UP.subtract(upper,lower)),
            'leaves':[{'id':x['id'],'interval':[str(x['a']),str(x['b'])],'upper':str(x['upper'])}
                      for x in sorted(state['leaves'],key=lambda x:x['id'])]}


def certify(model,tolerance='.01',max_splits=128,seconds=10,problem_id=None):
    if not math.isfinite(seconds) or seconds<0:raise ValueError('finite nonnegative time budget required')
    tol=D(tolerance)
    if not tol.is_finite() or tol<=0 or not isinstance(max_splits,int) or not 0<=max_splits<=10000:
        raise ValueError('positive tolerance and bounded split count required')
    if problem_id is not None and (not isinstance(problem_id,str) or not problem_id):raise ValueError('nonempty problem identity required')
    start=perf_counter();state=_initial(model);ops=[];reason='time_budget';largest=perf_counter()-start;splits=0
    evaluator=None
    if perf_counter()-start<seconds:
        t=perf_counter();evaluator=SpatialBounds(model)
        ordered=sorted(model['donors'],key=lambda r:UP.divide(I(*r['weight_MW']).hi,I(*r['sigma_m2']).lo),reverse=True)
        area=I(*model['period_m2']);points=[D(0),D('.5')]
        for row in ordered[:8]:
            c=I(*row['center_fraction']) if 'center_fraction' in row else I(*row['center_m2'])/area
            points.append(max(D(0),min(D(1),DOWN.divide(DOWN.add(c.lo,c.hi),D(2)))))
        largest=max(largest,perf_counter()-t)
        for point in dict.fromkeys(points):
            if perf_counter()-start>=seconds:break
            op={'kind':'probe','point':str(point)};t=perf_counter();_operation(state,evaluator,op)
            largest=max(largest,perf_counter()-t);ops.append(op)
        while True:
            gap=D(_evidence(state)['gap_upper_MW_m2'])
            if gap<=tol:reason='tolerance';break
            if perf_counter()-start>=seconds:reason='time_budget';break
            if splits>=max_splits:reason='split_budget';break
            active=max(state['leaves'],key=lambda x:(x['upper'],-x['id']))
            mid=DOWN.divide(DOWN.add(active['a'],active['b']),D(2))
            if not active['a']<mid<active['b']:reason='coordinate_precision';break
            op={'kind':'split','id':active['id']};t=perf_counter();_operation(state,evaluator,op)
            largest=max(largest,perf_counter()-t);ops.append(op);splits+=1
    t=perf_counter();evidence=_evidence(state)
    result=dict(schema='selective-spatial-peak.v1',model=copy.deepcopy(model),model_sha256=digest(model),
                problem_id=problem_id,identity_sha256=digest({'model':model,'problem_id':problem_id}),
                tolerance=str(tol),max_splits=max_splits,seconds_budget=seconds,
                stop_reason=reason,completed_splits=splits,operations=ops,**evidence,
                tolerance_met=D(evidence['gap_upper_MW_m2'])<=tol,
                scope='declared finite periodic-kernel family only',geometry_global_optimum_proven=False)
    largest=max(largest,perf_counter()-t);elapsed=perf_counter()-start
    result.update(elapsed_seconds=elapsed,budget_overrun_seconds=max(0.,elapsed-seconds),largest_atomic_operation_seconds=largest)
    return result


def replay(record):
    if (record.get('schema')!='selective-spatial-peak.v1' or digest(record['model'])!=record['model_sha256']
        or digest({'model':record['model'],'problem_id':record.get('problem_id')})!=record.get('identity_sha256')):
        raise ValueError('peak replay model mismatch')
    state=_initial(record['model']);evaluator=SpatialBounds(record['model']) if record['operations'] else None
    for op in record['operations']:
        if op['kind'] not in ('probe','split'):raise ValueError('peak replay operation mismatch')
        _operation(state,evaluator,op)
    expected=_evidence(state)
    if (any(record[k]!=v for k,v in expected.items())
        or record['completed_splits']!=sum(op['kind']=='split' for op in record['operations'])
        or record['tolerance_met']!=(D(expected['gap_upper_MW_m2'])<=D(record['tolerance']))):
        raise ValueError('peak replay mismatch')
    return True


def compare(a,b,*,problem_id_a,problem_id_b):
    if (not problem_id_a or problem_id_a!=problem_id_b
        or a.get('problem_id')!=problem_id_a or b.get('problem_id')!=problem_id_b):
        raise ValueError('same problem identity required')
    started=perf_counter();replay(a);replay(b)
    la,ua=map(D,a['peak_bracket_MW_m2']);lb,ub=map(D,b['peak_bracket_MW_m2'])
    return {'status':'a_lower' if ua<lb else 'b_lower' if ub<la else 'unresolved',
            'problem_id':problem_id_a,'replay_seconds':perf_counter()-started,'scope':'peak ranking only; feasibility is separate',
            'geometry_global_optimum_proven':False}
