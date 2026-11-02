"""Replayable inclusion certificates for the declared frozen-path family.

Feasibility, continuous magnetic trajectories, source integration and bitwise
production roundoff are separate obligations, never inferred from inclusion.
"""
import copy
import hashlib
import json
from .enclosure_math import (I,D,DOWN,UP,interval,record,hull,intersect,SplineField,rational_interval)
from .optimality_bounds import pi_bounds
from .geometry_transport_bounds import Surface,first_hits
from .thermal_box_bounds import PathThermal,thermal_width

SCHEMA='geometry-transport-enclosure.v1'


def digest(model):
    return hashlib.sha256(json.dumps(model,sort_keys=True,allow_nan=False).encode()).hexdigest()


def universal_peak_bounds(kernel):
    """P/A <= peak <= sum_i w_i*(1/A + 1/(sqrt(2*pi)*sigma_i)).

    A single wrapped Gaussian is maximized at its center (positive Fourier
    coefficients). Its centered image sum is bounded by g(0) plus the integral
    of the monotone tails / period. No image truncation budget is needed.
    """
    area=interval(kernel['period_m2']);pi=rational_interval(*pi_bounds())
    if area.lo<=0:raise ValueError('positive kernel area required')
    total=I(0);upper=I(0)
    for row in kernel['donors']:
        weight=interval(row['weight_MW']);sigma=interval(row['sigma_m2'])
        if weight.lo<0 or sigma.lo<=0:raise ValueError('positive width and nonnegative weight required')
        total+=weight;upper+=weight*(I(1)/area+I(1)/((2*pi).sqrt()*sigma))
    return [str((total/area).lo),str(upper.hi)]


def _validate(model,box):
    if model['schema']!='geometry-transport-input.v1':raise ValueError('unsupported model schema')
    domains=[interval(x) for x in model['domains']]
    if not domains or len(box)!=len(domains):raise ValueError('one interval per variable required')
    for x,d in zip(box,domains):
        if x.lo<d.lo or x.hi>d.hi:raise ValueError('box outside declared domain')
    if not model['surfaces'] or not model['donors']:raise ValueError('wall and donors required')
    for s,next_s in zip(model['surfaces'],model['surfaces'][1:]+model['surfaces'][:1]):
        k=s['degree'];knots=[interval(t) for t in s['knots']]
        if (not isinstance(k,int) or k<1 or len(knots)!=len(s['controls'])+k+1
                or any(t.lo!=t.hi for t in knots)
                or any(t.lo!=0 for t in knots[:k+1])
                or any(t.hi!=1 for t in knots[-k-1:])
                or any(t.lo<0 or t.hi>1 for t in knots)):
            raise ValueError('clamped spline on [0,1] required for closed wall')
        if s['b']!=next_s['a']:raise ValueError('ordered closed wall required')
        if any(isinstance(c,dict) or interval(c).lo!=0 or interval(c).hi!=0 for c in (s['controls'][0],s['controls'][-1])):
            raise ValueError('closed wall requires fixed zero displacement endpoints')
    if interval(model['grad_sep']).lo<=0:raise ValueError('positive source gradient required')
    power=[interval(t['weight_MW']) for t in model['donors']]
    if any(w.lo<0 for w in power) or sum((w.hi for w in power),D(0))<=0:
        raise ValueError('nonnegative donor powers with positive total required')


def certify_box(model,domain,root_depth=6,area_panels=16,thermal_panels=0):
    if not isinstance(thermal_panels,int) or not 0<=thermal_panels<=512:raise ValueError('thermal panels 0..512 required')
    box=[interval(x) for x in domain];_validate(model,box)
    out={'schema':SCHEMA,'model':copy.deepcopy(model),'model_sha256':digest(model),
         'box':list(map(record,box)),'root_depth':root_depth,'area_panels':area_panels,
         'scope':'exact continuous spline wall; frozen poloidal paths and nodal thermal model',
         'transport_parameter_enclosure_proven':False,'continuous_fieldline_enclosure_proven':False,
         'production_roundoff_enclosed':False,'geometry_feasibility_proven':False,
         'geometry_global_optimum_proven':False,'certified_geometry_optimality_gap':None,
         'donors':[],'kernel':None,'unresolved':[]}
    if thermal_panels:out['thermal_panels']=thermal_panels
    try:
        surfaces=[Surface(s,box,area_panels) for s in model['surfaces']]
        field=SplineField(model['field'])
    except ValueError as exc:
        out['unresolved'].append({'stage':'geometry_or_field','reason':str(exc)})
        return out
    area=[I(0)]
    for s in surfaces:area.append(area[-1]+s.total_area)
    total=area[-1];out['total_area_m2']=record(total)
    if total.lo<=0:
        out['unresolved'].append({'stage':'area','reason':'nonpositive total area lower bound'})
        return out
    grad_sep=interval(model['grad_sep']);kernel=[]
    for index,tube in enumerate(model['donors']):
        w=interval(tube['weight_MW'])
        if w.hi==0:
            out['donors'].append({'index':index,'status':'zero_input_power'});continue
        row={'index':index,'weight_MW':record(w)};out['donors'].append(row)
        try:
            hit=first_hits(tube['path'],surfaces,root_depth);row['first_hit']=hit
            if not hit['guaranteed_hit']:raise ValueError('not every geometry proven to hit before finite path end')
            possible=hit['possible_hits'];last=max(h['edge'] for h in possible)
            thermal=PathThermal(tube,last)
            # Merge local fractions per edge; all retained first-hit branches
            # remain included. Correlations may be relaxed, never guessed.
            groups={}
            for h in possible:groups.setdefault(h['edge'],[]).append(interval(h['fraction']))
            states=[thermal.at(j,hull(ts)) for j,ts in sorted(groups.items())]
            L=hull(s[0] for s in states);J=hull(s[1] for s in states);avg=hull(s[2] for s in states)
            points=[hull(interval(h['point'][k]) for h in possible) for k in range(2)]
            gradient,B=field.quantities(points)
            state=thermal_width(interval(tube['qu']),B/interval(tube['b0']),J,avg,model['reference'])
            sigma=state['sigma_m']
            if thermal_panels:
                avg=thermal.average(J,state['target_y'],state['energy'],thermal_panels)
                sigma=intersect(sigma,state['sigma_factor']*avg)
            physical=sigma*grad_sep/gradient
            sigma_area=2*surfaces[0].pi*points[0]*physical
            centers=[]
            by_surface={}
            for h in possible:by_surface.setdefault(h['surface'],[]).append(interval(h['u']))
            for j,us in sorted(by_surface.items()):centers.append(area[j]+surfaces[j].area(hull(us)))
            center=intersect(hull(centers)/total,I(0,1))
            if center is None:raise ArithmeticError('area coordinate outside period')
            row.update(status='enclosed',length_m=record(L),integral_Bratio_m=record(J),
                       arrival_RZ_m=list(map(record,points)),sigma_flux_m=record(sigma),
                       sigma_physical_m=record(physical),sigma_area_m2=record(sigma_area),
                       center_fraction=record(center),target_gradient=record(gradient),
                       Tu_eV=record(state['Tu_eV']),Tt_eV=record(state['Tt_eV']))
            kernel.append({'weight_MW':record(w),'sigma_m2':record(sigma_area),'center_fraction':record(center)})
        except ValueError as exc:
            row.update(status='unresolved',reason=str(exc))
            out['unresolved'].append({'stage':'donor','index':index,'reason':str(exc)})
    if not out['unresolved']:
        out['kernel']={'period_m2':record(total),'donors':kernel}
        out['transport_parameter_enclosure_proven']=True
        out['peak_bracket_MW_m2']=universal_peak_bounds(out['kernel'])
    out['prescribed_power_MW']=record(sum((interval(t['weight_MW']) for t in model['donors']),I(0)))
    return out


def replay(recorded):
    expected=certify_box(recorded['model'],recorded['box'],recorded['root_depth'],recorded['area_panels'],recorded.get('thermal_panels',0))
    if expected!=recorded:raise ValueError('geometry transport certificate replay mismatch')
    return True


def cover(model,splits=0,root_depth=6,area_panels=16):
    if not isinstance(splits,int) or not 0<=splits<=64:raise ValueError('split budget 0..64 required')
    domains=[interval(x) for x in model['domains']]
    leaves=[(1,domains)];history=[]
    for _ in range(splits):
        index=max(range(len(leaves)),key=lambda j:max(UP.subtract(x.hi,x.lo) for x in leaves[j][1]))
        identity,box=leaves[index]
        dimension=max(range(len(box)),key=lambda j:UP.subtract(box[j].hi,box[j].lo))
        x=box[dimension];mid=DOWN.divide(DOWN.add(x.lo,x.hi),D(2))
        if not x.lo<mid<x.hi:break
        left=list(box);right=list(box);left[dimension]=I(x.lo,mid);right[dimension]=I(mid,x.hi)
        leaves[index:index+1]=[(2*identity,left),(2*identity+1,right)]
        history.append({'parent':identity,'dimension':dimension,'boundary':str(mid)})
    records=[{'id':identity,'certificate':certify_box(model,list(map(record,box)),root_depth,area_panels)} for identity,box in leaves]
    return {'schema':'geometry-transport-cover.v1','model':copy.deepcopy(model),'splits':splits,
            'root_depth':root_depth,'area_panels':area_panels,'history':history,'leaves':records,
            'all_boxes_resolved':all(x['certificate']['transport_parameter_enclosure_proven'] for x in records),
            'continuous_fieldline_enclosure_proven':False,'geometry_global_optimum_proven':False}


def replay_cover(recorded):
    expected=cover(recorded['model'],recorded['splits'],recorded['root_depth'],recorded['area_panels'])
    if expected!=recorded:raise ValueError('geometry cover replay mismatch')
    return True
