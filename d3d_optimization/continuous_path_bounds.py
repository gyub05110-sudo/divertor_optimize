"""Validated graph-chart continuation for an exact axisymmetric spline field.

The guide chooses coordinate planes only. Picard inclusion proves the path tube;
the conserved flux contracts endpoints without resetting to the guide.
"""
import copy
from collections import Counter
from .enclosure_math import I,D,DOWN,UP,interval,record,hull,intersect,norm,SplineField,dot
from .geometry_transport_certificate import digest


def absolute(x):
    if x.lo>=0:return x
    if x.hi<=0:return I(x.hi.copy_negate(),x.lo.copy_negate())
    return I(0,max(x.lo.copy_negate(),x.hi))


class ContinuousField(SplineField):
    def __init__(self,model):
        for key in ('r','z'):
            k=model['degree_'+key];t=[interval(v) for v in model[key+'_knots']]
            if not isinstance(k,int) or k<1 or any(x.lo!=x.hi for x in t):
                raise ValueError('exact positive-degree C1 field required')
            if any(a.lo>b.lo for a,b in zip(t,t[1:])):raise ValueError('ordered field knots required')
            counts=Counter(x.lo for x in t if t[k].lo<x.lo<t[-k-1].lo)
            if any(v>k-1 for v in counts.values()):raise ValueError('C1 flux spline required')
        for row in model['coefficients']:
            if any(interval(v).lo!=interval(v).hi for v in row):raise ValueError('exact field coefficients required')
        if interval(model['toroidal_K']).lo!=interval(model['toroidal_K']).hi:raise ValueError('exact toroidal constant required')
        super().__init__(model)

    def value(self,p):
        p=list(map(interval,p));self.gradient(p) # includes domain checks
        return self._tensor(self.c,self.r,self.z,self.kr,self.kz,p)

    def chart(self,p,axis,sign,direction,b0,grad_sep):
        g=self.gradient(p);cross_axis=1-axis
        denom=g[cross_axis]
        oriented=direction*sign*(g[0] if axis==1 else -g[1])
        if p[0].lo<=0 or oriented.lo<=0:raise ValueError('chart direction or nonzero denominator not proven')
        f=-sign*g[axis]/denom
        gn=norm(g);bn=(gn.square()+self.K.square()).sqrt();B=bn/p[0]
        metric=bn/absolute(denom)
        return f,metric,B/b0*metric,gn/grad_sep*b0/B


def chart_step(field,start,target,axis,direction,psi0,b0,grad_sep):
    """Enclose a variable-time passage to one exact coordinate plane."""
    start=list(map(interval,start));target=interval(target)
    if axis not in (0,1) or direction not in (-1,1) or target.lo!=target.hi:
        raise ValueError('exact target plane, axis and branch required')
    other=1-axis;delta=target-start[axis]
    if delta.lo>0:sign=1
    elif delta.hi<0:sign=-1
    else:raise ValueError('target plane not strictly ahead of all incoming states')
    H=sign*delta;time=I(0,H.hi);axis_range=start[axis]+sign*time
    f0=field.chart(start,axis,sign,direction,b0,grad_sep)[0]
    seed=hull([start[other],start[other]+time*f0])
    # A vanishing initial slope can still curve immediately. This is only a
    # trial radius; the strict Picard test, never this heuristic, proves it.
    radius=max(UP.subtract(start[other].hi,start[other].lo),UP.divide(H.hi,D(4)),D('1e-12'))
    accepted=None
    for _ in range(7):
        q=seed+I(radius.copy_negate(),radius)
        tube=[None,None];tube[axis]=axis_range;tube[other]=q
        f,metric,integrand,mapping=field.chart(tube,axis,sign,direction,b0,grad_sep)
        image=start[other]+time*f
        if q.lo<image.lo and image.hi<q.hi:
            accepted=(tube,f,metric,integrand,mapping);break
        radius=UP.multiply(radius,D(2))
    if accepted is None:raise ValueError('Picard tube inclusion not established')
    tube,f,metric,integrand,mapping=accepted
    if getattr(field,'tight_tubes',False):
        # Picard already proved the entire trajectory lies in image. Intersect
        # that enclosure, then reevaluate; no guide/nominal state substitution.
        for _ in range(3):
            tube[other]=intersect(tube[other],start[other]+time*f)
            if tube[other] is None:raise ArithmeticError('proved tube disappeared')
            f,metric,integrand,mapping=field.chart(tube,axis,sign,direction,b0,grad_sep)
    endpoint=[None,None];endpoint[axis]=target
    endpoint[other]=intersect(start[other]+H*f,tube[other])
    if endpoint[other] is None:raise ArithmeticError('Picard endpoint disappeared')
    for _ in range(4):
        q=endpoint[other];mid=DOWN.divide(DOWN.add(q.lo,q.hi),D(2))
        center=list(endpoint);center[other]=I(mid)
        deriv=field.gradient(endpoint)[other]
        if deriv.lo<=0<=deriv.hi:break
        refined=intersect(q,I(mid)-(field.value(center)-psi0)/deriv)
        if refined is None:raise ArithmeticError('flux invariant lost the proved solution')
        endpoint[other]=refined
    return {'start':list(map(record,start)),'end':list(map(record,endpoint)),
            'tube':list(map(record,tube)),'axis':axis,'sign':sign,'time':record(H),
            'length_increment':record(H*metric),'integral_increment':record(H*integrand),
            'mapping':record(mapping),'picard_inclusion_proven':True}


def trace(field_model,guide,grad_sep,max_cells=2048,max_depth=10,field_enclosure='legacy'):
    if not isinstance(max_cells,int) or not 1<=max_cells<=100000:raise ValueError('positive bounded cell budget required')
    if not isinstance(max_depth,int) or not 0<=max_depth<=20:raise ValueError('subdivision depth 0..20 required')
    p=[list(map(interval,point)) for point in guide]
    if len(p)<2 or any(len(q)!=2 or any(x.lo!=x.hi for x in q) for q in p):raise ValueError('exact 2D guide vertices required')
    if field_enclosure not in ('legacy','centered-v1'):raise ValueError('unknown field enclosure')
    if field_enclosure=='centered-v1':
        from .selective_transport import TightContinuousField
        field=TightContinuousField(field_model)
    else:field=ContinuousField(field_model)
    grad_sep=interval(grad_sep)
    if grad_sep.lo<=0:raise ValueError('positive source gradient required')
    start=p[0];g=field.gradient(start);bpol=[-g[1],g[0]]
    orientation=dot(bpol,[b-a for a,b in zip(p[0],p[1])])
    if orientation.lo>0:direction=1
    elif orientation.hi<0:direction=-1
    else:raise ValueError('initial branch direction not determined')
    _,b0=field.quantities(start);psi0=field.value(start)
    inputs={'field':copy.deepcopy(field_model),'guide':copy.deepcopy(guide),'grad_sep':record(grad_sep),
            'max_cells':max_cells,'max_depth':max_depth}
    if field_enclosure!='legacy':inputs['field_enclosure']=field_enclosure
    out={'schema':'continuous-fieldline-path.v1','inputs':inputs,'input_sha256':digest(inputs),
         'cells':[],'complete':False,'continuous_prefix_proven':False,'unresolved':[],
         'direction':direction,'b0':record(b0),'psi0':record(psi0)}
    pending=[(a,b,0) for a,b in zip(p,p[1:])][::-1]
    L=I(0);J=I(0);attempts=0
    while pending:
        if len(out['cells'])>=max_cells or attempts>=max_cells*16:
            out['unresolved'].append('cell or attempt budget exhausted');break
        a,b,depth=pending.pop();axis=max(range(2),key=lambda j:absolute(b[j]-a[j]).hi)
        attempts+=1
        try:cell=chart_step(field,start,b[axis],axis,direction,psi0,b0,grad_sep)
        except ValueError as exc:
            if depth>=max_depth:
                out['unresolved'].append(str(exc));break
            mid=[I(DOWN.divide(DOWN.add(x.lo,y.lo),D(2))) for x,y in zip(a,b)]
            if any(x.lo!=y.lo for x,y in zip(mid,a)) and any(x.lo!=y.lo for x,y in zip(mid,b)):
                pending.extend([(mid,b,depth+1),(a,mid,depth+1)]);continue
            out['unresolved'].append(str(exc));break
        cell['length_start']=record(L);cell['integral_start']=record(J)
        L+=interval(cell['length_increment']);J+=interval(cell['integral_increment'])
        cell['length_end']=record(L);cell['integral_end']=record(J)
        out['cells'].append(cell);start=list(map(interval,cell['end']))
    out.update(complete=not pending and not out['unresolved'],continuous_prefix_proven=bool(out['cells']),attempts=attempts)
    return out


def replay(recorded):
    i=recorded['inputs']
    expected=trace(i['field'],i['guide'],i['grad_sep'],i['max_cells'],i['max_depth'],i.get('field_enclosure','legacy'))
    if expected!=recorded:raise ValueError('continuous path replay mismatch')
    return True
