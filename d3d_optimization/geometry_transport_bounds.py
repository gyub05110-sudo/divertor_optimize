"""All-root inclusion and first-hit existence for exact normal-graph walls.

The paths are declared polylines, not certified magnetic ODE solutions.
"""
from .enclosure_math import (I,D,DOWN,UP,interval,record,hull,intersect,
    spline_range,derivative,sin_cos,cross,dot,norm,rational_interval)
from .optimality_bounds import pi_bounds


class Surface:
    def __init__(self,model,box,panels=16):
        if not isinstance(panels,int) or not 1<=panels<=256:raise ValueError('panels 1..256 required')
        self.a=list(map(interval,model['a']));self.b=list(map(interval,model['b']))
        self.v=[b-a for a,b in zip(self.a,self.b)];self.n=list(map(interval,model['normal']))
        if len(self.a)!=2 or len(self.b)!=2 or len(self.n)!=2:raise ValueError('2D geometry required')
        self.det=cross(self.v,self.n)
        if self.det.lo<=0<=self.det.hi:raise ValueError('singular normal graph')
        self.t=model['knots'];self.degree=model['degree'];self.c=[]
        for value in model['controls']:
            if not isinstance(value,dict):self.c.append(interval(value));continue
            if 'variable' in value:self.c.append(box[value['variable']]);continue
            spec=value['junction'];fraction=box[spec['variable']]
            if fraction.lo<0 or fraction.hi>1:raise ValueError('junction fraction outside [0,1]')
            theta=interval(spec['angle_low'])+fraction*(interval(spec['angle_high'])-interval(spec['angle_low']))
            s,c=sin_cos(theta);tangent=[c,s];den=dot(tangent,self.v)
            if den.lo<=0:raise ValueError('junction denominator not proven positive; subdivide box')
            self.c.append(interval(spec['sign'])*dot(self.v,self.v)*dot(tangent,self.n)/den/15)
        # Validate domain and exact endpoints before any coverage assertion.
        spline_range(self.t,self.c,self.degree,I(0,1))
        self.dt,self.dc,self.dd=derivative(self.t,self.c,self.degree)
        self.exact_straight=bool(model.get('exact_straight_area',False)) and all(c.lo==c.hi==0 for c in self.c)
        self.panels=panels
        cuts=sorted(set([D(0),D(1)]+[interval(x).lo for x in self.t if 0<interval(x).lo<1]))
        self.cells=[]
        for a,b in zip(cuts,cuts[1:]):
            inner=[a]+[DOWN.add(a,DOWN.divide(DOWN.multiply(UP.subtract(b,a),D(j)),D(panels))) for j in range(1,panels)]+[b]
            self.cells.extend((l,h) for l,h in zip(inner,inner[1:]) if h>l)
        self.bounds=[hull(self.point(I(a,b))[k] for a,b in self.cells) for k in range(2)]
        self.pi=rational_interval(*pi_bounds())
        self.area_cells=[self._integral(a,b) for a,b in self.cells]
        self.total_area=sum(self.area_cells,I(0))

    def point(self,u):
        d=spline_range(self.t,self.c,self.degree,u)
        return [a+u*v+n*d for a,v,n in zip(self.a,self.v,self.n)]

    def _integral(self,a,b):
        u=I(a,b);p=self.point(u)
        if p[0].lo<=0:raise ValueError('positive radius not proven; subdivide box')
        if self.exact_straight:
            ra=self.a[0]+I(a)*self.v[0];rb=self.a[0]+I(b)*self.v[0]
            return self.pi*(ra+rb)*norm(self.v)*(I(b)-a)
        dp=spline_range(self.dt,self.dc,self.dd,u)
        speed=norm([v+n*dp for v,n in zip(self.v,self.n)])
        return 2*self.pi*p[0]*speed*(I(b)-a)

    def area(self,u):
        u=interval(u)
        if u.lo<0 or u.hi>1:raise ValueError('area parameter outside [0,1]')
        def upto(x):
            value=I(0)
            for (a,b),full in zip(self.cells,self.area_cells):
                if x<=a:break
                value+=full if x>=b else self._integral(a,x)
            return value
        # Positive radius/speed makes cumulative area monotone.
        return I(upto(u.lo).lo,upto(u.hi).hi)

    def coordinate(self,p):
        delta=[q-a for q,a in zip(p,self.a)]
        return cross(delta,self.n)/self.det,cross(self.v,delta)/self.det


def first_hits(path,surfaces,depth=6):
    if not isinstance(depth,int) or not 0<=depth<=12:raise ValueError('root depth 0..12 required')
    p=[list(map(interval,point)) for point in path]
    if len(p)<2 or any(len(x)!=2 for x in p):raise ValueError('2D path required')
    if any(any(x.lo!=x.hi for x in q) for q in p):raise ValueError('exact path vertices required')
    if any(all(x.lo==y.lo for x,y in zip(a,b)) for a,b in zip(p,p[1:])):raise ValueError('duplicate path vertex')
    possible=[];guarantees=[]
    for index,surface in enumerate(surfaces):
        previous_sign=None;previous_end=None
        for edge,(a,b) in enumerate(zip(p,p[1:])):
            # Exact scalar comparisons: no approximate broadphase exclusions.
            if any(max(a[k].hi,b[k].hi)<surface.bounds[k].lo or
                   min(a[k].lo,b[k].lo)>surface.bounds[k].hi for k in range(2)):
                previous_sign=None;previous_end=None;continue
            u0,d0=surface.coordinate(a);u1,d1=surface.coordinate(b)
            du=u1-u0;dd=d1-d0
            domain=I(0,1)
            if du.lo>0 or du.hi<0:
                endpoints=[(I(q)-u0)/du for q in (0,1)]
                domain=intersect(domain,hull(endpoints))
            elif hull([u0,u1]).hi<0 or hull([u0,u1]).lo>1:domain=None
            if domain is None:
                previous_sign=None;previous_end=None;continue
            def residual(t):
                u=intersect(u0+t*du,I(0,1))
                return None if u is None else d0+t*dd-spline_range(surface.t,surface.c,surface.degree,u)
            # A sign bracket across vertices is valid only along a connected
            # portion wholly within this SAME material graph's domain.
            inside=all(0<=u.lo and u.hi<=1 for u in (u0,u1))
            if not inside or previous_end!=edge:
                previous_sign=None
            if inside:
                for t in (D(0),D(1)):
                    r=residual(I(t))
                    sign=1 if r.lo>0 else -1 if r.hi<0 else 0
                    if r.lo==r.hi==0:
                        guarantees.append((edge,t))
                    if sign:
                        if previous_sign is not None and sign!=previous_sign:
                            guarantees.append((edge,t))
                        previous_sign=sign
                previous_end=edge+1
            else:previous_end=None
            pending=[(domain,0)]
            while pending:
                t,level=pending.pop();r=residual(t)
                if r is None or r.lo>0 or r.hi<0:continue
                if level<depth and t.lo<t.hi:
                    mid=DOWN.divide(DOWN.add(t.lo,t.hi),D(2))
                    if t.lo<mid<t.hi:
                        pending.extend([(I(t.lo,mid),level+1),(I(mid,t.hi),level+1)]);continue
                u=intersect(u0+t*du,I(0,1))
                full_u=u0+t*du
                if full_u.lo>=0 and full_u.hi<=1:
                    left=residual(I(t.lo));right=residual(I(t.hi))
                    if left.lo==left.hi==0:guarantees.append((edge,t.lo))
                    if right.lo==right.hi==0:guarantees.append((edge,t.hi))
                    if (left.hi<0 and right.lo>0) or (left.lo>0 and right.hi<0):
                        guarantees.append((edge,t.hi))
                point=[a[k]+t*(b[k]-a[k]) for k in range(2)]
                curve_point=surface.point(u)
                point=[intersect(x,y) for x,y in zip(point,curve_point)]
                if any(x is None for x in point):continue
                possible.append({'edge':edge,'fraction':record(t),'surface':index,
                                 'u':record(u),'point':list(map(record,point))})
    cutoff=min(guarantees) if guarantees else None
    if cutoff:
        kept=[]
        for h in possible:
            if (h['edge'],D(h['fraction'][0]))>cutoff:continue
            if h['edge']==cutoff[0]:h['fraction'][1]=str(min(D(h['fraction'][1]),cutoff[1]))
            kept.append(h)
        possible=kept
    if cutoff and not possible:raise ArithmeticError('existence proof without any retained root')
    return {'guaranteed_hit':cutoff is not None,'possible_hits':possible,
            'guaranteed_before':None if cutoff is None else [cutoff[0],str(cutoff[1])],
            'status':'enclosed' if cutoff is not None else 'unresolved_finite_horizon',
            'unique_surface_proven':bool(cutoff and len({h['surface'] for h in possible})==1)}
