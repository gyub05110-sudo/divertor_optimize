"""Tighter field ranges and selective continuation with unchanged exact inputs."""
from .enclosure_math import I,D,UP,interval,hull,intersect,record
from .continuous_path_bounds import ContinuousField,trace,replay


def blend(a,b,alpha):
    # Both algebraic forms include the same correlated affine combination.
    # Intersecting their ranges cannot remove its actual value.
    return intersect((1-alpha)*a+alpha*b,a+alpha*(b-a))


def tight_spline_range(knots,controls,degree,u):
    t=[interval(x) for x in knots];u=interval(u);c=list(map(interval,controls))
    if (any(x.lo!=x.hi for x in t) or len(t)!=len(c)+degree+1 or degree<0
        or any(a.lo>b.lo for a,b in zip(t,t[1:])) or u.lo<t[degree].lo or u.hi>t[len(c)].hi):
        raise ValueError('invalid exact spline domain')
    t=[x.lo for x in t];values=[]
    for span in range(degree,len(c)):
        if t[span]>=t[span+1] or u.hi<t[span] or u.lo>t[span+1]:continue
        v=I(max(u.lo,t[span]),min(u.hi,t[span+1]));work=c[span-degree:span+1]
        for r in range(1,degree+1):
            for j in range(degree,r-1,-1):
                i=span-degree+j;alpha=(v-t[i])/(I(t[i+degree-r+1])-t[i])
                work[j]=blend(work[j-1],work[j],alpha)
        values.append(intersect(work[degree],hull(c[span-degree:span+1])))
    return hull(values)


class TightContinuousField(ContinuousField):
    tight_tubes=True

    def _tensor(self,c,r,z,kr,kz,p):
        tz=[interval(x).lo for x in z];zz=p[1];values=[]
        for span in range(kz,len(c)):
            if tz[span]>=tz[span+1] or zz.hi<tz[span] or zz.lo>tz[span+1]:continue
            local=[tight_spline_range(r,c[j],kr,p[0]) for j in range(span-kz,span+1)]
            work=local[:];v=I(max(zz.lo,tz[span]),min(zz.hi,tz[span+1]))
            for q in range(1,kz+1):
                for j in range(kz,q-1,-1):
                    i=span-kz+j;alpha=(v-tz[i])/(I(tz[i+kz-q+1])-tz[i])
                    work[j]=blend(work[j-1],work[j],alpha)
            values.append(intersect(work[kz],hull(local)))
        return hull(values)


def trace_tight(field_model,guide,grad_sep,max_cells=4096,max_depth=10):
    return trace(field_model,guide,grad_sep,max_cells,max_depth,field_enclosure='centered-v1')


def replay_tight(path):
    if path['inputs'].get('field_enclosure')!='centered-v1':raise ValueError('tight replay mode missing')
    return replay(path)


def select_cells(path,hit_cells,limit=24):
    """Heuristic allocation only; each newly proposed guide is proved afresh."""
    if not isinstance(limit,int) or limit<1:raise ValueError('positive cell budget required')
    if any(not 0<=j<len(path['cells']) for j in hit_cells):raise ValueError('invalid hit cell')
    last=max(hit_cells) if hit_cells else len(path['cells'])-1
    selected={j for h in hit_cells for j in (h-1,h,h+1) if 0<=j<=last}
    def width(j):
        x=interval(path['cells'][j]['integral_increment']);return UP.subtract(x.hi,x.lo)
    ordered=sorted(range(last+1),key=lambda j:(j in selected,width(j)),reverse=True)
    return set(ordered[:limit])


def refine_donor(path,selected,factor=4):
    from verification.run_continuous_path_study import refined_guide
    if not path['complete']:raise ValueError('refinement requires complete original horizon')
    guide=refined_guide(path,selected,factor)
    # Guide midpoints may round long exact decimals. Preserve every original
    # chart's terminal plane verbatim, in particular the final trace horizon.
    position=0
    for index,cell in enumerate(path['cells']):
        position+=factor if index in selected else 1
        axis=cell['axis'];guide[position][axis]=cell['end'][axis][0]
    i=path['inputs']
    result=trace_tight(i['field'],guide,i['grad_sep'],max_cells=max(4096,len(guide)*2),max_depth=10)
    if result['complete']:
        old=path['cells'][-1];new=result['cells'][-1];axis=old['axis']
        if new['axis']!=axis or new['end'][axis]!=old['end'][axis]:
            raise ValueError('refinement must preserve exact terminal plane')
    return result


def next_donor(rows,visits,max_visits=2):
    """Schedule by conservative amplitude contribution, not just power weight."""
    choices=[]
    for row in rows:
        index=row['index'];count=visits.get(index,0)
        if count>=max_visits or row['status']=='zero_input_power':continue
        if row['status']=='unresolved':score=float('inf')
        else:
            w=interval(row['weight_MW']);s=interval(row['sigma_area_m2'])
            score=float(w.hi/s.lo)/(1+count)
        choices.append((score,-index))
    return -max(choices)[1] if choices else None
