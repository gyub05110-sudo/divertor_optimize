"""Numerical coverage helpers; missing connectivity never implies zero heat."""
import numpy as np

def refinement_cells(edges, arclength, components, q, gap=.002):
    edges=np.asarray(edges);s=np.asarray(arclength);q=np.asarray(q)
    selected=set()
    for i in range(len(s)-1):
        if (components[i]!=components[i+1] or
            (components[i] in ('floor','shelf') and abs(s[i+1]-s[i])>gap)):
            selected.update((i,i+1))
    for row in q:
        if np.max(row)>0:
            for i in np.flatnonzero(row>=.8*np.max(row)):
                selected.update((max(0,i-1),i,min(len(s)-1,i+1)))
    return sorted(i for i in selected if edges[i+1]-edges[i]>1e-7)

def connection_status(in_scope, reached, is_first):
    if not in_scope:return 'outside_source_scope'
    if not reached:return 'trace_unresolved'
    return 'direct_hit' if is_first else 'blocked_first_hit'

def surface_integral(vertices, q):
    vertices=np.asarray(vertices);q=np.asarray(q)
    areas=2*np.pi*(vertices[:-1,0]+vertices[1:,0])/2*np.linalg.norm(np.diff(vertices,axis=0),axis=1)
    return np.sum(q*areas,axis=-1)

def compare_levels(previous, latest, peak_tol=.01, location_tol=.002, ratio_tol=.01):
    peak=np.abs(np.asarray(latest['peaks'])/np.asarray(previous['peaks'])-1)
    position=np.linalg.norm(np.asarray(latest['locations'])-np.asarray(previous['locations']),axis=1)
    ratio=np.abs(np.asarray(latest['ratios'])/np.asarray(previous['ratios'])-1)
    return {'peak_relative_change':peak.tolist(),'peak_location_change_m':position.tolist(),
            'ratio_relative_change':ratio.tolist(),
            'passed':bool(np.max(peak)<=peak_tol and np.max(position)<=location_tol and np.max(ratio)<=ratio_tol)}
