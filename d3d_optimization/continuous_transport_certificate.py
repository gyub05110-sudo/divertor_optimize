"""Geometry verification with replayed continuous-fieldline caches.

Caches are independent of wall variables and reusable across coefficient boxes.
The fast optimizer and its stored incumbents are not modified here.
"""
import copy
from .enclosure_math import I,D,DOWN,UP,interval,record,hull,intersect,dot
from .geometry_transport_bounds import Surface
from .geometry_transport_certificate import digest,_validate,universal_peak_bounds
from .continuous_path_bounds import ContinuousField,replay as replay_path
from .thermal_box_bounds import thermal_width


def branches(model):
    field=ContinuousField(model['field']);result=[]
    for donor in model['donors']:
        if len(donor['path'])<2:raise ValueError('source branch requires two guide points')
        a,b=[list(map(interval,p)) for p in donor['path'][:2]]
        g=field.gradient(a);orientation=dot([-g[1],g[0]],[y-x for x,y in zip(a,b)])
        if orientation.lo>0:result.append(1)
        elif orientation.hi<0:result.append(-1)
        else:raise ValueError('model source branch is ambiguous')
    return result


def cache_key(model):
    return digest({'field':model['field'],'grad_sep':record(interval(model['grad_sep'])),
                   'starts':[d['path'][0] for d in model['donors']],'branches':branches(model)})


class ValidatedPaths:
    """Replay once, own a private copy, then reuse for many geometry boxes."""
    def __init__(self,model,records):
        if len(records)!=len(model['donors']):raise ValueError('one continuous path per donor required')
        owned=copy.deepcopy(records)
        for donor,path,direction in zip(model['donors'],owned,branches(model)):
            i=path['inputs']
            if (i['field']!=model['field'] or list(map(record,map(interval,i['guide'][0])))!=list(map(record,map(interval,donor['path'][0])))
                    or record(interval(i['grad_sep']))!=record(interval(model['grad_sep']))):
                raise ValueError('continuous cache/model mismatch')
            if path['direction']!=direction:raise ValueError('continuous cache branch mismatch')
            replay_path(path)
        self._paths=owned;self._key=cache_key(model)
        self._hashes=[digest(p) for p in owned]

    def check(self,model):
        if cache_key(model)!=self._key or [digest(p) for p in self._paths]!=self._hashes:
            raise ValueError('continuous cache/model mismatch')


def first_hits(path,surfaces):
    possible=[];guarantees=[]
    for si,s in enumerate(surfaces):
        previous_sign=None;previous_cell=None
        for ci,cell in enumerate(path['cells']):
            tube=list(map(interval,cell['tube']))
            if any(tube[k].hi<s.bounds[k].lo or tube[k].lo>s.bounds[k].hi for k in range(2)):
                previous_sign=None;previous_cell=None;continue
            u,d=s.coordinate(tube);clipped=intersect(u,I(0,1))
            if clipped is None:
                previous_sign=None;previous_cell=None;continue
            from .enclosure_math import spline_range
            residual=d-spline_range(s.t,s.c,s.degree,clipped)
            if residual.lo<=0<=residual.hi:
                curve=s.point(clipped);point=[intersect(x,y) for x,y in zip(tube,curve)]
                if all(x is not None for x in point):
                    possible.append({'cell':ci,'surface':si,'u':record(clipped),'point':list(map(record,point))})
            # Whole un-clipped tube must stay in the same material domain.
            inside=u.lo>=0 and u.hi<=1
            if not inside or previous_cell!=ci-1:previous_sign=None
            if inside:
                for endpoint in ('start','end'):
                    eu,ed=s.coordinate(list(map(interval,cell[endpoint])))
                    rr=ed-spline_range(s.t,s.c,s.degree,eu)
                    sign=1 if rr.lo>0 else -1 if rr.hi<0 else 0
                    if rr.lo==rr.hi==0:guarantees.append(ci)
                    if sign:
                        if previous_sign is not None and sign!=previous_sign:guarantees.append(ci)
                        previous_sign=sign
                previous_cell=ci
            else:previous_cell=None
    cutoff=min(guarantees) if guarantees else None
    if cutoff is not None:possible=[p for p in possible if p['cell']<=cutoff]
    if cutoff is not None and not possible:raise ArithmeticError('continuous existence without retained root')
    return {'guaranteed_hit':cutoff is not None,'guaranteed_cell':cutoff,'possible_hits':possible}


def mapping_average(cells,J,state,panels):
    S=(state['target_y'].square()+state['energy']).sqrt()+state['target_y']
    ratio=state['energy']/S.square();total=I(0)
    cuts=[D(0)]+[DOWN.divide(D(j),D(panels)) for j in range(1,panels)]+[D(1)]
    for a,b in zip(cuts,cuts[1:]):
        u=I(a,b);p=intersect(J*(1-u)*(1+u*ratio),I(0,J.hi))
        values=[interval(c['mapping']) for c in cells
                if interval(c['integral_start']).lo<=p.hi and interval(c['integral_end']).hi>=p.lo]
        if not values:raise ArithmeticError('continuous inverse-integral branches disappeared')
        total+=hull(values)*(I(b)-a)
    return intersect(total,hull(interval(c['mapping']) for c in cells))


def certify_box(model,domain,cache,area_panels=8,thermal_panels=64):
    if not isinstance(cache,ValidatedPaths):raise ValueError('replayed continuous cache required')
    if not isinstance(thermal_panels,int) or not 1<=thermal_panels<=512:raise ValueError('thermal panels 1..512 required')
    cache.check(model);box=list(map(interval,domain));_validate(model,box)
    out={'schema':'continuous-geometry-transport.v1','model':copy.deepcopy(model),'box':list(map(record,box)),
         'path_sha256':list(cache._hashes),'area_panels':area_panels,'thermal_panels':thermal_panels,
         'transport_parameter_enclosure_proven':False,'continuous_fieldline_enclosure_proven':False,
         'geometry_feasibility_proven':False,'geometry_global_optimum_proven':False,
         'continuous_source_enclosure_proven':False,'field_reconstruction_error_enclosed':False,
         'certified_geometry_optimality_gap':None,'donors':[],'kernel':None,'unresolved':[]}
    try:surfaces=[Surface(s,box,area_panels) for s in model['surfaces']]
    except ValueError as exc:
        out['unresolved'].append({'stage':'geometry','reason':str(exc)});return out
    areas=[I(0)]
    for s in surfaces:areas.append(areas[-1]+s.total_area)
    total=areas[-1];out['total_area_m2']=record(total)
    field=ContinuousField(model['field']);grad_sep=interval(model['grad_sep']);kernel=[]
    for index,(donor,path) in enumerate(zip(model['donors'],cache._paths)):
        w=interval(donor['weight_MW']);row={'index':index,'weight_MW':record(w)};out['donors'].append(row)
        if w.hi==0:row['status']='zero_input_power';continue
        try:
            hit=first_hits(path,surfaces);row['first_hit']=hit
            if not hit['guaranteed_hit']:raise ValueError('no guaranteed first hit in proved continuous prefix')
            possible=hit['possible_hits'];indices=sorted({h['cell'] for h in possible})
            L=hull(I(interval(path['cells'][j]['length_start']).lo,interval(path['cells'][j]['length_end']).hi) for j in indices)
            J=hull(I(interval(path['cells'][j]['integral_start']).lo,interval(path['cells'][j]['integral_end']).hi) for j in indices)
            cells=path['cells'][:max(indices)+1];avg=hull(interval(c['mapping']) for c in cells)
            points=[hull(interval(h['point'][k]) for h in possible) for k in range(2)]
            gradient,B=field.quantities(points)
            state=thermal_width(donor['qu'],B/interval(path['b0']),J,avg,model['reference'])
            avg=mapping_average(cells,J,state,thermal_panels)
            sigma=intersect(state['sigma_m'],state['sigma_factor']*avg)
            physical=sigma*grad_sep/gradient;sigma_area=2*surfaces[0].pi*points[0]*physical
            centers=[]
            for si in sorted({h['surface'] for h in possible}):
                u=hull(interval(h['u']) for h in possible if h['surface']==si)
                centers.append(areas[si]+surfaces[si].area(u))
            center=intersect(hull(centers)/total,I(0,1))
            row.update(status='enclosed',arrival_RZ_m=list(map(record,points)),length_m=record(L),
                       integral_Bratio_m=record(J),sigma_flux_m=record(sigma),sigma_physical_m=record(physical),
                       sigma_area_m2=record(sigma_area),center_fraction=record(center),
                       path_horizon_complete=path['complete'],proved_prefix_cells=len(path['cells']))
            kernel.append({'weight_MW':record(w),'sigma_m2':record(sigma_area),'center_fraction':record(center)})
        except ValueError as exc:
            row.update(status='unresolved',reason=str(exc));out['unresolved'].append({'stage':'donor','index':index,'reason':str(exc)})
    out['prescribed_power_MW']=record(sum((interval(d['weight_MW']) for d in model['donors']),I(0)))
    if not out['unresolved']:
        out.update(transport_parameter_enclosure_proven=True,continuous_fieldline_enclosure_proven=True,
                   kernel={'period_m2':record(total),'donors':kernel})
        out['peak_bracket_MW_m2']=universal_peak_bounds(out['kernel'])
    return out


def replay(recorded,cache):
    expected=certify_box(recorded['model'],recorded['box'],cache,recorded['area_panels'],recorded['thermal_panels'])
    if expected!=recorded:raise ValueError('continuous transport replay mismatch')
    return True


def cover(model,cache,splits=2,area_panels=4,thermal_panels=16):
    if not isinstance(splits,int) or not 0<=splits<=64:raise ValueError('split budget 0..64 required')
    domain=list(map(interval,model['domains']));history=[]
    def leaf(identity,box):return {'id':identity,'certificate':certify_box(model,list(map(record,box)),cache,area_panels,thermal_panels)}
    leaves=[leaf(1,domain)]
    for _ in range(splits):
        choices=[]
        for i,entry in enumerate(leaves):
            c=entry['certificate'];b=list(map(interval,c['box']))
            widths=[UP.divide(UP.subtract(x.hi,x.lo),UP.subtract(d.hi,d.lo)) if d.lo<d.hi else D(0) for x,d in zip(b,domain)]
            dimension=max(range(len(b)),key=lambda j:widths[j])
            if (model.get('source')=='IMP' and len(b)==9 and
                    any('junction denominator' in q['reason'] for q in c['unresolved']) and widths[8]>0):dimension=8
            choices.append((bool(c['unresolved']),widths[dimension],-i,dimension))
        _,width,negative_index,dimension=max(choices);index=-negative_index
        if width<=0:break
        entry=leaves[index];b=list(map(interval,entry['certificate']['box']));x=b[dimension]
        mid=DOWN.divide(DOWN.add(x.lo,x.hi),D(2))
        if not x.lo<mid<x.hi:break
        left=list(b);right=list(b);left[dimension]=I(x.lo,mid);right[dimension]=I(mid,x.hi)
        history.append({'parent':entry['id'],'dimension':dimension,'boundary':str(mid)})
        leaves[index:index+1]=[leaf(2*entry['id'],left),leaf(2*entry['id']+1,right)]
    return {'schema':'continuous-geometry-cover.v1','model':copy.deepcopy(model),'splits':splits,
            'area_panels':area_panels,'thermal_panels':thermal_panels,'history':history,'leaves':leaves,
            'all_boxes_resolved':all(e['certificate']['transport_parameter_enclosure_proven'] for e in leaves),
            'geometry_global_optimum_proven':False,'certified_geometry_optimality_gap':None}


def replay_cover(recorded,cache):
    expected=cover(recorded['model'],cache,recorded['splits'],recorded['area_panels'],recorded['thermal_panels'])
    if expected!=recorded:raise ValueError('continuous cover replay mismatch')
    return True
